//	BLMPI.cc
/*
This source file is part of G4beamline, http://g4beamline.muonsinc.com
Copyright (C) 2002-2013 by Tom Roberts, all rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
*/

//	This is the only G4beamline file that uses MPI routines.
//	It needs G4BL_MPI defined, or it is all a big no-op.
/*
	While it sounds simple to have rank 0 generate events, parcel
	them out to the workers, and collect output data from the workers,
	in practice it is rather complicated to build a program that is
	robust.

	Because processes can hang, rank 0 keeps a an alarm set to
	fire if ever 2*eventTimeLimit passes without any MPI message.
	Normally the logic in nPackTracks() ensures the GetBeamTracks
	messages prevent that, but at the end of events it is possible
	for just a few workers to be busy, and fail to send enough
	messages. So KeepAlive messages are enabled in each worker
	when it receives a KeepAlive message; each worker then sends
	a KeepAlive message every event, ensuring sufficient messages
	are received by rank 0.
*/

#ifdef G4BL_MPI
#include <mpi.h>
#include <unistd.h>
#endif

#include <stdio.h>
#include <stdarg.h>
#include <signal.h>
#include <vector>
#include <map>

#include "G4ApplicationState.hh"
#include "G4StateManager.hh"
#include "G4PrimaryTransformer.hh"

#include "BLAssert.hh"
#include "BLAlarm.hh"
#include "BLSignal.hh"
#include "BLTime.hh"
#include "BLMPI.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLRunManager.hh"
#include "BLBeam.hh"
#include "BLElement.hh"
#include "BLNTuple.hh"
#include "BLRealTimeMonitor.hh"
#include "BLWriteAsciiFile.hh"

// MPI parameter values
int MPI_MaxTracksPerMsg=100;
int MPI_GoalMsgPerSec=200;
int MPI_PackNTuples=1000;
int MPI_IdleSleep = 0;
int MPI_WorkerOutput = 1;
int MPI_Debug = 0;

// class BLMPI static values
int BLMPI::rank=1;	// makes isRank1() return true if not MPI
int BLMPI::nRanks=0;
int BLMPI::alarmTime=0;
bool BLMPI::processingFatalException=false;
bool BLMPI::readyForBeam=false;
int BLMPI::nExiting=0;
bool BLMPI::inScan=false;
int BLMPI::inMPI=0;
double BLMPI::startupTime=0.0;
double BLMPI::startupWaitTime=0.0;
double BLMPI::computeTime=0.0;
double BLMPI::computeWaitTime=0.0;

#ifdef G4BL_MPI // extends to the end of this file, except for dummy routines

// local static variables
static MPI_Comm allWorkers;
static bool sendKeepAlive = false;
static long startTime=0;
static BLNTuple *allTracksNTuple=0;

/**	enum MPITags defines the tags for MPI messages.
 *	Collective routines use no tag.
 *	NTuple handles are also tags >= FIRST_NTUPLE_HANDLE.
 **/
enum MPITags { CREATE_NTUPLE=1274, NTUPLE_HANDLE, 
		GET_GLOBAL_PARAMETERS,  GLOBAL_PARAMETERS,
		GET_BEAM_TRACKS, BEAM_TRACKS, TRACK, 
		WRITEFILE_MSG, EXCEPTIONMSG, NEXT_TIME_STEP,
		KEEP_ALIVE_MSG, EXITING_MSG, ABORT_MSG };
const int FIRST_SUM_TAG=2000;
const int FIRST_NTUPLE_HANDLE=10000;

const int MAX_NTUPLE_DOUBLES=4096;	// Max # doubles in an NTuple message
const int MAX_FILENAME_CHARS=128;	// Max # chars in a filename
const int MAX_DATA_CHARS=32768;		// Max # data chars in WRITEFILE_MSG.
const int MAX_PRINTF_CHARS=1024;	// Max # chars in rank0BufferPrintf();


static void log(const char *fmt, ...)
{
	if(MPI_Debug <= 0) return;
	
	va_list ap;
	va_start(ap,fmt);

	printf("%08ld(%d): ",BLTime::timems()-startTime,BLMPI::getRank());
	vprintf(fmt,ap);
	fflush(stdout);
}

/**	struct Msg is a base class for MPI messages.
 *	There must be no virtual functions in any Msg struct.
 *
 *	Each derived struct has datatype(), sendTo0() or send(), and recv()
 *	functions.
 **/
struct Msg {
	/// Constructor (no data).
	Msg() { }

	/// Constructs an MPI_Datatype for a struct consisting of chars,
	/// ints, and doubles, in that order. Be sure the # chars is a
	/// multiple of 4.
	static MPI_Datatype charIntDouble(int nchar, int nint, int ndouble);

	static void error(const char *type, int dest, int tag) {
		char tmp1[64], tmp2[64];
		sprintf(tmp1,"MPI Error in %s",type);
		sprintf(tmp2,"dest=%d src=%d tag=%d",dest,BLMPI::getRank(),tag);
		G4Exception("BLMPI",tmp1,JustWarning,tmp2);
	}

	/// Sends 1 struct of type datatype to rank 0.
	static void sendTo0(void *data, MPI_Datatype datatype, int tag)
		{ if(MPI_Send(data,1,datatype,0,tag,MPI_COMM_WORLD) !=
		     						MPI_SUCCESS) 
		     	error("MPI_Send",0,tag);
		  log("sendTo0(tag=%d)\n",tag);
		}

	/// Receives a message with specified tag from any source.
	/// The source can be obtained from status.
	static void recv(void *data, MPI_Datatype datatype,int tag,
							MPI_Status &status)
		{ MPI_Recv(data,1,datatype,MPI_ANY_SOURCE,tag, MPI_COMM_WORLD,
								&status); 
		  log("recv(tag=%d from=%d)\n",status.MPI_TAG,status.MPI_SOURCE);
		}

	/// sends an NTuple row to rank 0 (nData should be an integer multiple
	/// of the # doubles in each row)
	static void sendNTupleRows(double *data, int nData, int handle)
		{ MPI_Send(data,nData,MPI_DOUBLE,0,handle,MPI_COMM_WORLD); 
		  log("sendNTupleRows(%d)\n",handle);
		}

	/// Sends a message to any rank, without waiting.
	static void sendNoWait(int dest, void *data,
						MPI_Datatype datatype, int tag)
		{ MPI_Request retval;
		  MPI_Isend(data,1,datatype,dest,tag,MPI_COMM_WORLD,&retval);
		  log("sendNoWait(dest=%d tag=%d)\n",dest,tag);
		}
};
struct CreateNTuple : public Msg {
	char type[32];
	char category[128];
	char name[128];
	char filename[128];
	char fields[4000-32-3*128];
	MPI_Datatype datatype() { return charIntDouble(4000,0,0); }
	void sendTo0() { Msg::sendTo0(this,datatype(),CREATE_NTUPLE); }
	void recv(MPI_Status &status) 
		{ Msg::recv(this,datatype(),CREATE_NTUPLE,status); }
};
struct NTupleHandle : public Msg {
	int handle;
	MPI_Datatype datatype() { return charIntDouble(0,1,0); }
	void send(int dest) 
	      { MPI_Send(this,1,datatype(),dest,NTUPLE_HANDLE,MPI_COMM_WORLD); 
		  log("send(dest=%d tag=%d)\n",dest,NTUPLE_HANDLE);
	      }
	void recv(MPI_Status &status) 
		{ Msg::recv(this,datatype(),NTUPLE_HANDLE,status); }
};
struct GlobalParameters : public Msg {
	int maxTracks;	// max # tracks in a BEAM_TRACKS message
	int nodeTracks;	// max # tracks for each node (collective mode)
	MPI_Datatype datatype() { return charIntDouble(0,2,0); }
	void send(int dest) 
	      { MPI_Send(this,1,datatype(),dest,GLOBAL_PARAMETERS,
	      						MPI_COMM_WORLD); 
		  log("send(dest=%d tag=%d\n",dest,GLOBAL_PARAMETERS);
	      }
	void recv(MPI_Status &status) 
		{ Msg::recv(this,datatype(),GLOBAL_PARAMETERS,status); }
};
struct Track : public Msg {
	int PDGid;
	int eventID;
	int trackID;		// external value
	int parentID;		// external value
	int secondaryTrackID;
	int filler;
	double position[3];	// mm
	double direction[3];
	double polarization[3];
	double kineticEnergy;	// MeV
	double time;		// ns
	double weight;
	MPI_Datatype datatype() { return charIntDouble(0,6,12); }
	int tag() { return TRACK; }
	void send(int dest) 
	      { MPI_Send(this,1,datatype(),dest,TRACK,MPI_COMM_WORLD); 
		log("send(dest=%d tag=%d)\n",dest,TRACK);
	      }
	void recv(MPI_Status &status) 
		{ Msg::recv(this,datatype(),TRACK,status); }
};
struct GetBeamTracks : public Msg {
	int sequence;
	GetBeamTracks() { sequence = 0; }
	MPI_Datatype datatype() { return charIntDouble(0,1,0); }
	void sendTo0() { 
		++sequence;
		log("sending GetBeamTracks(%d) to 0\n",sequence);
		Msg::sendTo0(this,datatype(),GET_BEAM_TRACKS); }
	void recv(MPI_Status &status) 
		{ Msg::recv(this,datatype(),GET_BEAM_TRACKS,status); 
		  log("received GetBeamTracks(%d) from %d\n",sequence,
		  					status.MPI_SOURCE);
		}
};
struct ExceptionMsg : public Msg {
	char origin[64];
	char code[64];
	char severity[32];
	char description[256];
	char particleName[64];
	int eventID;
	int trackID;
	int abortProgram;
	int padding;
	double kineticEnergy;
	MPI_Datatype datatype() { return charIntDouble(64+64+32+256+64,4,1); }
	void sendTo0() 
		{ log("ExceptionMsg code=%s\n",code);
		  Msg::sendTo0(this,datatype(),EXCEPTIONMSG); 
		}
	void recv(MPI_Status &status) 
		{ Msg::recv(this,datatype(),EXCEPTIONMSG,status); 
		  log("ExceptionMsg code=%s\n",code);
		}
};
struct BeamTracksMsg : public Msg {
	Track *tracks;
	int nTracks;
	bool busy;
	MPI_Request request;
	BeamTracksMsg() : busy(false), request()
		{ tracks=new Track[MPI_MaxTracksPerMsg];
		  nTracks=MPI_MaxTracksPerMsg; }
	~BeamTracksMsg() { if(tracks) delete[] tracks; tracks=0; nTracks=0; }
	MPI_Datatype datatype() { return tracks->datatype(); }
	void send(int dest, int n) 
		{ BLAssert(isIdle());
		  busy = true;
		  MPI_Isend(tracks,n,tracks->datatype(),dest,BEAM_TRACKS,
						MPI_COMM_WORLD,&request); 
		  log("send BeamTracksMsg(dest=%d tag=%d n=%d ev=%d)\n",dest,BEAM_TRACKS,n,(n>0?tracks[0].eventID:-1));
		}
	void recv(MPI_Status &status)	// get # Tracks from status
		{ MPI_Recv(tracks,nTracks,tracks->datatype(),MPI_ANY_SOURCE,
					BEAM_TRACKS, MPI_COMM_WORLD, &status); 
		  log("recvBeamTracks()\n");
		}
	bool isIdle() { // rank 0 MUST call isIdle() before using tracks
		if(!busy) return true;
		int flag=0;
		MPI_Test(&request,&flag,MPI_STATUS_IGNORE);
		busy = (flag == 0);
		return !busy;
	}
};
struct NextTimeStep : public Msg {
	double suggestedDeltaT;	// ns
	MPI_Datatype datatype() { return charIntDouble(0,0,1); }
	void send(int dest) 
	     { MPI_Send(this,1,datatype(),dest,NEXT_TIME_STEP,MPI_COMM_WORLD); 
	       log("send(dest=%d tag=%d\n",dest,NEXT_TIME_STEP);
	     }
	void recv(MPI_Status &status) 
		{ Msg::recv(this,datatype(),NEXT_TIME_STEP,status); }
};
struct WriteFileMsg : public Msg {
	char filename[MAX_FILENAME_CHARS];
	char data[MAX_DATA_CHARS];
	int ndata;
	WriteFileMsg() { filename[0]='\0'; data[0]='\0'; ndata=0; }
	void sendTo0() {
		MPI_Send(this,MAX_FILENAME_CHARS+ndata+1,MPI_CHAR,0,
						WRITEFILE_MSG,MPI_COMM_WORLD);
		data[0] = '\0';
		ndata = 0;
		log("sendTo0(tag=%d)\n",WRITEFILE_MSG);
	}
	void recv(MPI_Status &status) { 
		MPI_Recv(this,MAX_FILENAME_CHARS+MAX_DATA_CHARS,MPI_CHAR,
			MPI_ANY_SOURCE,WRITEFILE_MSG, MPI_COMM_WORLD, &status); 
		data[MAX_DATA_CHARS-1] = '\0';
		ndata = strlen(data);
	}
};
struct KeepAliveMsg : public Msg {
	int sequence;
	KeepAliveMsg() { sequence = 0; }
	MPI_Datatype datatype() { return charIntDouble(0,1,0); }
	void sendTo0() { 
		  ++sequence;
		  Msg::sendTo0(this,datatype(),KEEP_ALIVE_MSG); 
		  log("send KEEP_ALIVE_MSG(%d) to 0\n",sequence);
		}
	void recv(MPI_Status &status) 
		{ Msg::recv(this,datatype(),KEEP_ALIVE_MSG,status); 
		  log("received KEEP_ALIVE_MSG(%d) from %d\n",sequence,
		  					status.MPI_SOURCE);
		}
	void sendNoWait(int dest) 
		{
		  ++sequence;
		  Msg::sendNoWait(dest,this,datatype(),KEEP_ALIVE_MSG); 
		  log("send KEEP_ALIVE_MSG(%d) to %d\n",sequence,dest);
		}
};
struct ExitingMsg : public Msg {
	double startupTime;
	double startupWaitTime;
	double computeTime;
	double computeWaitTime;
	MPI_Datatype datatype() { return charIntDouble(0,0,4); }
	void sendTo0() { Msg::sendTo0(this,datatype(),EXITING_MSG); }
	void recv(MPI_Status &status) 
		{ Msg::recv(this,datatype(),EXITING_MSG,status); }
};
struct AbortProgram : public Msg {
	int dummy;
	MPI_Datatype datatype() { return charIntDouble(0,1,0); }
	void sendTo0() { Msg::sendTo0(this,datatype(),ABORT_MSG); }
	void recv(MPI_Status &status) 
		{ Msg::recv(this,datatype(),ABORT_MSG,status); }
	void sendNoWait(int dest)
		{ Msg::sendNoWait(dest,this,datatype(),ABORT_MSG); }
};

class SumData {
	double *data;
	int nData;
	int nSum;
public:
	SumData() { data=0; nData=0; nSum=0; }
	void add(double d[], int nd);
	void getSum(double d[], int nd);
};

static std::map<int,SumData> sumData;
static std::map<G4String,WriteFileMsg> writeFiles;

/**	class MPIBeam implements a beam for worker ranks.
 **/
class MPIBeam : public BLBeam {
	G4ParticleGun *particleGun;
public:
	MPIBeam() : BLBeam() { particleGun = new G4ParticleGun(1); }

	virtual void init() { }

	virtual bool generateReferenceParticle(G4Event *event) { return false; }

	virtual bool nextBeamEvent(G4Event *event);
};

/**	class MPINTuple implements a generic MPI NTuple.
 **/
class MPINTuple : public BLNTuple {
	int handle;
	int nfields;
	G4String fields;
	double *buffer;
	int iBuffer;	// current position in buffer[]
	int nBuffer;	// # rows in buffer[]
public:
	MPINTuple(int _handle, G4String _name, G4String _fields);
	~MPINTuple() { close(); }
	virtual void appendRow(double data[], int n);
	virtual bool readRow(double data[], int ndata) { return false; }
	virtual int getNData() { return nfields; }
	virtual G4String getFieldNames() { return fields; }
	virtual void annotate(G4String line) { }
	void close() { flush(); }
	void flush() { 
		if(iBuffer > 0) 
			Msg::sendNTupleRows(buffer,nfields*iBuffer,handle); 
		iBuffer = 0;
	}
	void doSummary() { }
};
/**	class MPINTupleHandler implements the handler for MPI NTuples.
 *	NOTE: this class forces itself as the handler for all NTuple types.
 *	So it must be explicitly constructed in MPI mode (i.e. no default 
 *	instance exists).
 **/
class MPINTupleHandler : public BLNTupleHandler {
	static std::map<G4String,int> name2handle;
	static std::map<int,BLNTuple*> handle2NTuple;
	static int nextHandle;
public:
	MPINTupleHandler() { 
		BLNTuple::registerHandler("MPI",this);
		BLNTuple::setForceHandler(this);
	}
	BLNTuple *create(G4String type, G4String category,
			G4String name, G4String fields, G4String filename);
	BLNTuple *read(G4String type, G4String category, G4String name,
					G4String fields,G4String filename) {
		return 0;
	}
	virtual bool ignoreMultipleCreates() { return true; }
	static int getHandle(const char *name) { return name2handle[name]; }
	static BLNTuple *getNTuple(int handle) { return handle2NTuple[handle]; }
};
std::map<G4String,int> MPINTupleHandler::name2handle;
std::map<int,BLNTuple*> MPINTupleHandler::handle2NTuple;
int MPINTupleHandler::nextHandle=FIRST_NTUPLE_HANDLE;


// utility routine delarations

/// handleNTupleMessage() will receive and handle any message related to NTuples
/// in rank 0. Returns true if the message was handled, false otherwise.
static bool handleNTupleMessage(MPI_Status &status);

/// handleTraceNtuple() keeps a buffer for each worker, adding data to it
/// until the current track ends, at which time it flushes the buffer to
/// the ntuple. This keeps all rows from each track together, but the order
/// of tracks is arbitrary. Applies only to the single ntuple whose name
/// starts "AllTracks"; called from handleNTupleMessage().
/// Called with data=0 and nData=0 to flush all internal buffers.
static void handleTraceNtuple(BLNTuple *ntuple,double *data,int nData,int src);

/// handleSumMessage() will receive and handle any Sum message 
/// in rank 0. Returns true if the message was handled, false otherwise.
static bool handleSumMessage(MPI_Status &status);

/// sendBeamTracks() is used by rank 0 to send a BEAM_TRACKS message to dest.
static void sendBeamTracks(int dest);

/// processOneTrackAndAllSecondaries() is used in rank nonzero processes to 
/// track a Track and any secondaries it might create (recursively)
static void processOneTrackAndAllSecondaries(Track &track);

/// plotMsg() generates an NTuple of messages/sec received by rank 0.
static void plotMsg(bool force=false);

//	monitor the real-time usage. States:
//		0 startup
//		1 startup-wait
//		2 compute
//		3 compute-wait
//		4 closeup
BLRealTimeMonitor monitor;

void BLMPI::scan(bool wait)
{
	log("BLMPI::scan() entered inScan=%d\n",(inScan?1:0));
	if(inScan || !isMPI()) return;
	inScan = true;

	MPI_Status status;

	if(rank != 0) { // never waits, but might call closeupAndExit()
		static KeepAliveMsg keepAlive;
		if(sendKeepAlive) {
			keepAlive.sendNoWait(0);
		}
		int flag=0;
		MPI_Iprobe(0,KEEP_ALIVE_MSG,MPI_COMM_WORLD,&flag,&status);
		if(flag != 0) {
			log("received KEEP_ALIVE_MSG\n");
			keepAlive.sendNoWait(0);
			sendKeepAlive = true;
			KeepAliveMsg m;
			m.recv(status);
			flag = 0;
		}
		MPI_Iprobe(0,ABORT_MSG,MPI_COMM_WORLD,&flag,&status);
		if(flag == 0) {
			inScan = false;
			return;
		}
		// We just received an AbortProgram message
		AbortProgram abt;
		abt.recv(status);
		fflush(stdout);
		// cannot use G4Exception, as rank 0 is not receiving messages
		fprintf(stderr,"AbortProgram Message Received -- Exiting\n");
		closeupAndExit(0);
	}

	// rank 0
	for(;;) {
		if(BLSignal::sigusr1()) {
			log("sigusr1 detected - Ntuples flushed\n");
			G4Exception("BLMPI","Ntuples Flushed",JustWarning,"");
			BLNTuple::flushAll();
		}
		BLAlarm::set(alarmTime);
		if(!wait) {
			int flag=0;
			MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,
								&flag, &status);
			log("scan(wait=false) returns no message\n");
			if(flag == 0) break;
		} else {
			monitor.incrState(1);
			log("scan(wait=true) waits for message\n");
			// blocking probe for any message from anybody
			if(MPI_IdleSleep <= 0) {
			   MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
			} else {
			   for(;;) {
				int flag;
				MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&flag,&status);
				if(flag != 0) break;
				BLTime::sleepms(MPI_IdleSleep);
			   }
			}
			monitor.incrState(-1);
		}

		log("scan() received message tag=%d source=%d\n",
					status.MPI_TAG,status.MPI_SOURCE);
		// Avoid infinite loop if !readyForBeam
		if(!readyForBeam && status.MPI_TAG == GET_BEAM_TRACKS) {
			BLAssert(!wait);
			log("scan() not readyForBeam, returns\n");
			break;
		}

		plotMsg();

		if(handleNTupleMessage(status)) continue;
		if(handleSumMessage(status)) continue;

		switch(status.MPI_TAG) {
		case GET_BEAM_TRACKS:
			BLAssert(readyForBeam);
			{ struct GetBeamTracks get;
			  get.recv(status);
			  if(processingFatalException) {
				BeamTracksMsg msg;
				msg.send(status.MPI_SOURCE,0);
			  	break;
			  }
			  sendBeamTracks(status.MPI_SOURCE);
			}
			break;
		case EXITING_MSG:
			{ struct ExitingMsg exiting;
			  exiting.recv(status);
			  startupTime += exiting.startupTime;
			  startupWaitTime += exiting.startupWaitTime;
			  computeTime += exiting.computeTime;
			  computeWaitTime += exiting.computeWaitTime;
			  log("ExitingMsg received, have %d/%d\n",
			  				nExiting+1,nRanks-1);
			  if(++nExiting == nRanks-1) {
			  	inScan = false;
				return;
			  }
			}
			break;
		case EXCEPTIONMSG:
			{ struct ExceptionMsg e;
			  e.recv(status);
			  BLManager::getObject()->printException(e.origin,
			  		e.code,e.severity,e.description,
					e.eventID,e.trackID,e.particleName,
					e.kineticEnergy,e.abortProgram,true);
			  if(e.abortProgram != 0) {
				log("Fatal ExceptionMsg received\n");
				if(!processingFatalException) {
				    processingFatalException = true;
				    monitor.setState(4);
				    log("Sending AbortProgram to Workers\n");
				    AbortProgram abt;
				    for(int i=1; i<nRanks; ++i)
					abt.sendNoWait(i);
				}
			  }
			}
			break;
		case WRITEFILE_MSG:
			{ WriteFileMsg m;
			  m.recv(status);
			  FILE *f = BLWriteAsciiFile::fopen(m.filename);
			  fputs(m.data,f);
			}
			break;
		case KEEP_ALIVE_MSG:
			{ KeepAliveMsg m;
			  m.recv(status);
			}
			break;	// nothing to do
		default:
			BLAssert(false && "Illegal incoming message tag in rank 0");
		}
	}
	inScan = false;
}

/**	rank0AlarmHandler() will abort the MPI job, protecting against infinite
 *	hangs. Note that scan() sets an alarm after receiving every message.
 *	used in rank 0 only.
 **/
void BLMPI::rank0AlarmHandler()
{
	fflush(stdout);
	BLManager::getObject()->printException("BLMPI","MPI is Hung",
			  	"Fatal Exception","",0,0,"",0,true);
	fprintf(stderr,"\n*******************************************\n");
	fprintf(stderr,  "***  MPI is hung! Aborting entire job.  ***\n");
	fprintf(stderr,  "*******************************************\n");
	MPI_Abort(MPI_COMM_WORLD,99);
	_exit(99);
}

/**	mainRankZero() is the main program for the rank 0 process
 *	Note this is called after the input.file has been read and the
 *	reference particle has been tracked.
 **/
void BLMPI::mainRankZero()
{
	log("mainRankZero entered\n");
	if(!BLElement::allOK()) {
		G4Exception("BLMPI::main","Element Error",FatalException,
				"Not all BLElement-s are OK");
	}

	BLManager::getObject()->setState(BEAM);
	BLRunManager::getObject()->beginRun();

	BLAlarm::clear();
	BLAlarm::setCustomHandler(rank0AlarmHandler);
	BLAlarm::set(alarmTime);

	printf("================== Begin Tracking Beam ===============\n");
	readyForBeam = true;
	monitor.setState(2);	// compute 
	scan(true);		// handles messages until all EXITING_MSG rcvd.

	closeupAndExit(0);
}

bool MPIBeam::nextBeamEvent(G4Event *event)
{
	static bool requestSent=false;
	static GetBeamTracks getMsg;
	static BeamTracksMsg tracksMsg;
	static int count=0;
	static int next=0;
	static int prevEventID = -999;

	log("nextBeamEvent() entered\n");
	// return next track in tracksMsg
nextEv:	while(next < count) {
		BLMPI::scan();
		Track &track=tracksMsg.tracks[next++];
		G4ParticleDefinition *particleDefinition =
					G4ParticleTable::GetParticleTable()->
		    				FindParticle(track.PDGid);
		if(!particleDefinition) {
			char tmp[64];
			sprintf(tmp,"PDGid=%d",track.PDGid);
			G4Exception("MPIBeam","Invalid PDGid - track abandoned",
				JustWarning,tmp);
			continue;
		}
		particleGun->SetParticleDefinition(particleDefinition);
		particleGun->SetParticleTime(track.time);
		G4ThreeVector pos(track.position[0],track.position[1],track.position[2]);
		particleGun->SetParticlePosition(pos);
		particleGun->SetParticleEnergy(track.kineticEnergy);
		G4ThreeVector dir(track.direction[0],track.direction[1],track.direction[2]);
		particleGun->SetParticleMomentumDirection(dir);
		G4ThreeVector pol(track.polarization[0],track.polarization[1],track.polarization[2]);
		particleGun->SetParticlePolarization(pol);
		particleGun->GeneratePrimaryVertex(event);
		event->SetEventID(track.eventID);
		event->GetPrimaryVertex()->SetWeight(track.weight);
		if(track.eventID != prevEventID) {
			setRandomSeedToTrack(track.eventID);
			BLManager::getObject()->clearTrackIDMap();
			BLManager::getObject()->setNextSecondaryTrackID(track.secondaryTrackID);
			prevEventID = track.eventID;
		}
		BLManager::getObject()->setPrimaryTrackID(track.trackID,track.parentID);
		if(track.trackID >= track.secondaryTrackID) {
			G4Exception("beam command","Large Primary TrackID",
				JustWarning,
				"Confusion with secondary tracks is likely");
		}
		log("nextBeamEvent() returns event %d\n",track.eventID);
		return true;
	}

	// request another BeamTracksMsg, if necessary
	if(!requestSent) {
		monitor.incrState(1);
		log("nextBeamEvent sending initial getMsg\n");
		getMsg.sendTo0();
		requestSent = true;
		monitor.incrState(-1);
	}

	// wait for and receive a tracksMsg
	BLAssert(requestSent);
	MPI_Status status;
	monitor.incrState(1);
	BLMPI::scan();
	requestSent = false; // because we are receiving the request we sent
	log("nextBeamEvent receiving tracksMsg\n");
	tracksMsg.recv(status);
	log("nextBeamEvent sending getMsg\n");
	getMsg.sendTo0(); // request next block of events while we work
	requestSent = true; // because we just sent it
	monitor.incrState(-1);
	next = 0;
	MPI_Get_count(&status,tracksMsg.datatype(),&count);
	if(count == 0) {
		log("nextBeamEvent returns, end of events received\n");
		fflush(stdout);
		fprintf(stderr,"End of Events Message Received -- Exiting\n");
		return false;
	}
	log("nextBeamEvent count=%d firstEv=%d\n",count,
						tracksMsg.tracks[0].eventID);

	// go back and process this tracksMsg
	goto nextEv;
}

/**	mainRankNonZero is the main program for ranks other than 0
 *	Note this is called after the input.file has been read, and the
 *	reference particle has been tracked.
 **/
void BLMPI::mainRankNonZero()
{
	log("mainRankNonZero entered\n");
	BLManager *manager = BLManager::getObject();

	manager->clearBeamVector();
	manager->registerBeam(new MPIBeam());

	monitor.setState(2);	// compute

	// return to have the normal main() process events.
	return;
}

void BLMPI::init(std::vector<G4String> &argV)
{
	if(!isMPI()) {
		Param.setParam("MPI_rank",-1);
		printf("G4beamline built for MPI, but not in MPI mode\n");
		return;
	}
	startTime = BLTime::timems();
	monitor.setState(0);	// startup

	// MPI_Init -- convert G4String array to char* array
	int argc = argV.size();
	char **argv = new char *[argc+1];
	for(int i=0; i<argc; ++i)
		argv[i] = strdup(argV[i].c_str());
	argv[argc] = 0;
	MPI_Init(&argc,&argv);
	for(int i=0; i<argc; ++i)
		delete argv[i];
	delete argv;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nRanks);
	MPI_Comm_split(MPI_COMM_WORLD,(rank==0?MPI_UNDEFINED:0),rank,
								&allWorkers);
	BLAssert(nRanks > 1);

	if(rank != 0) {
		if(MPI_WorkerOutput != 0) {
			char tmp[32];
			sprintf(tmp,"rank%d.out",rank);
			stdout = freopen(tmp,"w",stdout);
			fclose(stderr);
			stderr = fdopen(1,"w");
			setvbuf(stderr,0,_IONBF,0);
		} else {
			fclose(stdout);
			fclose(stderr);
		}
	}

	// This causes rank 0 to collapse multiple NTUple definitions.
	// It also causes other ranks to create NTuples via MPI messages.
	new MPINTupleHandler();

	// set default values of MPI parameters
	Param.setParam("MPI_rank",rank);
	if(!Param.isDefined("MPI_MaxTracksPerMsg"))
		Param.setParam("MPI_MaxTracksPerMsg",MPI_MaxTracksPerMsg);
	if(!Param.isDefined("MPI_GoalMsgPerSec"))
		Param.setParam("MPI_GoalMsgPerSec",MPI_GoalMsgPerSec);
	if(!Param.isDefined("MPI_PackNTuples"))
		Param.setParam("MPI_PackNTuples",MPI_PackNTuples);
	if(!Param.isDefined("MPI_IdleSleep"))
		Param.setParam("MPI_IdleSleep",MPI_IdleSleep);
	if(!Param.isDefined("MPI_WorkerOutput"))
		Param.setParam("MPI_WorkerOutput",MPI_WorkerOutput);
	if(!Param.isDefined("MPI_Debug"))
		Param.setParam("MPI_Debug",MPI_Debug);

	printf("MPI MODE rank=%d #ranks=%d\n",rank,nRanks);
}

bool BLMPI::isMPI()
{
	if(inMPI != 0) return inMPI > 0;

	// NOTE: isMPI() must determine whether we are in MPI mode _BEFORE_
	// MPI_Init() is called (because we must not call it if not in MPI
	// mode). So we look for variables in the environment.

	// documented variable for OpenMPI
	if(getenv("OMPI_UNIVERSE_SIZE") != 0)
		inMPI = 1;
	// documented variable for aprun on a Cray at NERSC
	else if(getenv("ALPS_APP_DEPTH") != 0) 
		inMPI = 1;
	// GUESSED variable for srun on a Cray at NERSC
	else if(getenv("SLURM_SRUN_COMM_HOST") != 0)
		inMPI = 1;
	// GUESSED variable for mpirun on heimdall
	else if(getenv("MPICH_INTERFACE_HOSTNAME") != 0)
		inMPI = 1;
	else
		inMPI = -1;

	return inMPI > 0;
}

void BLMPI::rank0BufferPrintf(G4String filename, const char *format, ...)
{
	BLAssert(filename.size() < MAX_FILENAME_CHARS);
	va_list ap;
	va_start(ap,format);
	char tmp[MAX_PRINTF_CHARS];
	int n=vsnprintf(tmp,MAX_PRINTF_CHARS,format,ap);
	if(n < 0) return;
	if(n > MAX_PRINTF_CHARS-1) n = MAX_PRINTF_CHARS-1;
	WriteFileMsg &m = writeFiles[filename];
	if(m.filename[0] == '\0') strcpy(m.filename,filename.c_str());
	if(n+m.ndata >= MAX_DATA_CHARS) m.sendTo0();
	strcpy(m.data+m.ndata,tmp);
	m.ndata += n;
}

void BLMPI::rank0BufferFlush()
{
	std::map<G4String,WriteFileMsg>::iterator i;
	for(i=writeFiles.begin(); i!=writeFiles.end(); ++i) {
		if(i->second.ndata > 0) {
			if(!isMPI() || rank == 0) {
				FILE *f=BLWriteAsciiFile::fopen(i->second.filename);
				fputs(i->second.data,f);
			} else {
				i->second.sendTo0();
			}
			i->second.ndata = 0;
		}
	}
}

bool BLMPI::sumAllWorkers(double data[], int nData)
{
	static int tag=FIRST_SUM_TAG;

	if(!isMPI()) return true;
	if(rank == 0) {
		sumData[tag++].getSum(data,nData); // create if not found
		return true;
	} else {
		MPI_Send(data,nData,MPI_DOUBLE,0,tag++,MPI_COMM_WORLD);
		log("sendTo0(SumData tag=%d nData=%d)\n", tag-1,nData);
	}
	return false;
}

void closeupAndExit99()
{
	BLMPI::closeupAndExit(99);
}

void BLMPI::main()
{
	// Note this is called after input.file has been read and the 
	// reference particle has been tracked.

	if(!isMPI()) return;

	// This should be the last call to atexit(), hence the first callback;
	// it short-circuits all others (from libraries) and exits.
	atexit(closeupAndExit99);

	// get current (final) values of parameters
	MPI_MaxTracksPerMsg = Param.getInt("MPI_MaxTracksPerMsg");
	MPI_GoalMsgPerSec = Param.getInt("MPI_GoalMsgPerSec");
	MPI_PackNTuples = Param.getInt("MPI_PackNTuples");
	MPI_IdleSleep = Param.getInt("MPI_IdleSleep");
	MPI_WorkerOutput = Param.getInt("MPI_WorkerOutput");
	MPI_Debug = Param.getInt("MPI_Debug");
	alarmTime = Param.getInt("eventTimeLimit") * 2;
	if(alarmTime < 60) alarmTime = 60;

	if(rank == 0)
		mainRankZero();
	else
		mainRankNonZero();
	// return means the normal main() will process events
}

void closeupAlarmhandler()
{
	G4Exception("closeupAndExit","Alarm Timeout",JustWarning,"");
	BLMPI::closeupAndExit(99);
}

void BLMPI::closeupAndExit(int value)
{
	BLAlarm::clear();
	BLAlarm::setCustomHandler(closeupAlarmhandler); // comes right back here
	if(alarmTime < 600) alarmTime = 600;
	BLAlarm::set(alarmTime);

	BLRunManager *runManager = BLRunManager::getObject();
	monitor.setState(4);
//MPI_Debug = 1;
//log("closeupAndExit forcing MPI_Debug=1\n");

	static int state=0;
	if(state != 0) {
		//@ char tmp[32];
		//@ sprintf(tmp,"state=%d",state);
		//@ G4Exception("BLMPI::closeupAndExit","closeupAndExit re-entered",
		//@ 					JustWarning,tmp);
		log("closeupAndExit re-entered, state=%d\n",state);
	}

	switch(state) { // all states flow into the next
	case 0:	++state;
		log("closeupAndExit state=0, discard incoming messages\n");
		if(rank != 0) {
			// discard any incoming messages
			MPI_Status status;
			for(;;) {
			    int flag=0;
			    MPI_Iprobe(0,ABORT_MSG,MPI_COMM_WORLD,&flag,&status);
			    if(flag == 0) break;
			    char data[32000];
			    MPI_Recv(data,32000,MPI_CHAR,MPI_ANY_SOURCE,
			    		MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			}
		}
	case 1: ++state;
		log("closeupAndExit state=%d, endRun\n",state-1);
		BLAlarm::set(alarmTime);
		// In rank 0, these will do the real thing, writing to files,
		// etc. By the time rank 0 is here, all workers are waiting
		// in state 6 inside MPI_Finalize().
		// On worker ranks, some of these might call 
		// MPI_Reduce(...AllWorkers...) or other collective MPI calls,
		// but they will be in the same order for all ranks. Note that
		// rank 0 is waiting to receive an ExitingMsg from all workers
		// before calling closeupAndExit, so rank 0 must NOT 
		// participate in any collective MPI calls.
		runManager->endRun();
	case 2: ++state;
		log("closeupAndExit state=%d, BLNTuple::summary\n",state-1);
		BLAlarm::set(alarmTime);
		plotMsg(true);
		handleTraceNtuple(allTracksNTuple,0,0,0); // flush all buffers
		BLNTuple::summary();
	case 3: ++state;
		log("closeupAndExit state=%d, BLNTuple::closeAll\n",state-1);
		BLAlarm::set(alarmTime);
		BLNTuple::closeAll();
	case 4: ++state;
		log("closeupAndExit state=%d, handleCallbacks\n",state-1);
		BLAlarm::set(alarmTime);
		BLManager::getObject()->handleCallbacks(2);
		fflush(stdout);
	case 5: ++state;
		BLAlarm::set(alarmTime);
		if(rank != 0) {
			log("closeupAndExit state=%d, send ExitingMsg\n",
								state-1);
			// Send ExitingMsg to rank 0
			ExitingMsg exiting;
			exiting.startupTime = monitor.getTime(0);
			exiting.startupWaitTime = monitor.getTime(1);
			exiting.computeTime = monitor.getTime(2);
			exiting.computeWaitTime = monitor.getTime(3);
			exiting.sendTo0();
		} else {
			log("closeupAndExit state=%d, no-op\n",state-1);
		}
	case 6:	++state;
		log("closeupAndExit state=%d, MPI_Finalize\n",state-1);
		BLAlarm::clear(); // MPI_Finalize() can take a long time
		MPI_Finalize();
	case 7:	++state;
		log("closeupAndExit state=%d, BLWriteAsciiFile::closeAll\n",state-1);
		BLAlarm::set(alarmTime);
		BLWriteAsciiFile::closeAll();
	case 8: ++state;
		log("closeupAndExit state=%d, MPI summary\n",state-1);
		BLAlarm::set(alarmTime);
		monitor.setState(-1);
		if(rank == 0) {
		  double t=((double)BLTime::timems() - (double)startTime)/1.0E3;
		  printf("MPI COMPUTATION COMPLETE  elapsed time=%.3f sec.  %d workers\n",
			t,nRanks-1);
		  printf(" Rank0 Setup:%-9.3f Wait:%-9.3f    Run:%-9.3f Wait:%-9.3f\n",
			monitor.getTime(0),monitor.getTime(1),monitor.getTime(2),
			monitor.getTime(3));
		  printf(" 1-N   Setup:%-9.3f Wait:%-9.3f    Run:%-9.3f Wait:%-9.3f\n",
			startupTime,startupWaitTime,computeTime,computeWaitTime);
		  printf(" Agv.  Setup:%-9.3f Wait:%-9.3f    Run:%-9.3f Wait:%-9.3f\n",
			startupTime/(nRanks-1),startupWaitTime/(nRanks-1),
			computeTime/(nRanks-1),computeWaitTime/(nRanks-1));
		  startupTime += monitor.getTime(0);
		  startupWaitTime += monitor.getTime(1);
		  computeTime += monitor.getTime(2);
		  computeWaitTime += monitor.getTime(3);
		  printf(" Total Setup:%-9.3f Wait:%-9.3f    Run:%-9.3f Wait:%-9.3f\n",
			startupTime,startupWaitTime,computeTime,computeWaitTime);
		  printf(" Closeup took %.3f seconds.\n",monitor.getTime(4));
		}
		fflush(stdout);
	case 9:	++state;
		log("closeupAndExit state=%d\n",state-1);
		BLAlarm::set(alarmTime);
		if(value == 0)
			printf("g4beamline: simulation complete\n");
		else
			printf("g4beamline: simulation aborted\n");
		fflush(stdout);
	default: ++state;
		log("closeupAndExit state=%d _exit\n",state-1);
		BLAlarm::set(alarmTime);
		fflush(stdout);
	}

	_exit(value);
}

void BLMPI::sendExceptionMessage(const char *origin, const char *code,
		const char *severity, const char *description,
		int eventID, int trackID, const char *particleName,
		G4double kineticEnergy, bool abortProgram)
{
	BLAssert(rank > 0);

	ExceptionMsg e;
	strncpy(e.origin,origin,sizeof(e.origin));
	strncpy(e.code,code,sizeof(e.code));
	strncpy(e.severity,severity,sizeof(e.severity));
	strncpy(e.description,description,sizeof(e.description));
	e.eventID = eventID;
	e.trackID = trackID;
	strncpy(e.particleName,particleName,sizeof(e.particleName));
	e.kineticEnergy = kineticEnergy;
	e.abortProgram = (abortProgram ? 1 : 0);

	e.sendTo0();
}

MPI_Datatype Msg::charIntDouble(int nchar, int nint, int ndouble)
{
	static std::map<unsigned long,MPI_Datatype> known;

	BLAssert(nchar<4096 && nint<1024 && ndouble<1024 && (nchar&3)==0);
	unsigned long key = ((unsigned long)nchar<<20) | (nint<<10) | ndouble;
	MPI_Datatype v = known[key];
	if(v != 0) return v;

	int len[3]; MPI_Aint offset[3]; MPI_Datatype type[3];
	int count=0, addr=0;
	if(nchar > 0) {
		len[count] = nchar;
		offset[count] = addr;
		type[count++] = MPI_CHAR;
		addr += nchar * sizeof(char);
	}
	if(nint > 0) {
		BLAssert(addr%sizeof(int) == 0);
		len[count] = nint;
		offset[count] = addr;
		type[count++] = MPI_INT;
		addr += nint * sizeof(int);
	}
	if(ndouble > 0) {
		BLAssert(addr%sizeof(double) == 0);
		len[count] = ndouble;
		offset[count] = addr;
		type[count++] = MPI_DOUBLE;
		addr += ndouble * sizeof(double);
	}
	MPI_Type_struct(count,len,offset,type,&v);
	MPI_Type_commit(&v);
	known[key] = v;
	return v;
}

bool handleNTupleMessage(MPI_Status &status)
{
	if(status.MPI_TAG >= FIRST_NTUPLE_HANDLE) {
		static double *data=0;
		static int ndata=0;
		int count;
		MPI_Get_count(&status,MPI_DOUBLE,&count);
		if(count > ndata) {
			if(data) delete[] data;
			data = new double[count];
			ndata = count;
		}
		MPI_Recv(data,count,MPI_DOUBLE,status.MPI_SOURCE,status.MPI_TAG,
							MPI_COMM_WORLD,&status);
		BLNTuple *ntuple = MPINTupleHandler::getNTuple(status.MPI_TAG);
		if(ntuple) {
			if(ntuple == allTracksNTuple) {
				handleTraceNtuple(ntuple,data,count,
							status.MPI_SOURCE);
			} else {
				int k=ntuple->getNData();
				for(int j=0; j<count; j+=k)
					ntuple->appendRow(&data[j],k);
			}
		} else {
			char tmp[64];
			sprintf(tmp,"handle=%d",status.MPI_TAG);
			G4Exception("MPImsgRcv", "Unknown NTuple handle",
							JustWarning,tmp);
		}
		return true;
	} else if(status.MPI_TAG == CREATE_NTUPLE) {
		CreateNTuple c;
		c.recv(status);
		G4String name = c.category;
		if(name.size() > 0) name += "/";
		name += c.name;
		NTupleHandle h;
		h.handle = MPINTupleHandler::getHandle(name);
		if(h.handle == 0) {
			BLNTuple *ntuple = BLNTuple::create(c.type,c.category,
						c.name,c.fields,c.filename);
			BLAssert(ntuple != 0);
			h.handle = MPINTupleHandler::getHandle(name);
			BLAssert(h.handle != 0);
			if(strncmp(c.name,"AllTracks",9) == 0 || 
						name.find("AllTracks") == 0)
				allTracksNTuple = ntuple;
		}
		monitor.incrState(1);
		h.send(status.MPI_SOURCE);
		monitor.incrState(-1);
		return true;
	}

	return false;
}



const int NROW=18, NEROW=28;         // trace NTuple can have these two sizes
struct Row { double data[NROW]; };
struct ERow { double data[NEROW]; };
static int fileEvent = -999; // previous event # written to NTuple file

static void dumpRowBuf(BLNTuple *ntuple, std::vector<Row> &buf)
{
	char tmp[64];

	if(ntuple == 0 || buf.size() == 0) return;

	// buf contains rows all with the same event and track
	int event=buf[0].data[8];
	int track=buf[0].data[9];

	if(fileEvent != event && fileEvent != -999) {
		sprintf(tmp,"# End of Event %d",fileEvent);
		ntuple->annotate(tmp);
		ntuple->annotate("");
	}
	fileEvent = event;

	sprintf(tmp,"# Event %d Track %d",event,track);
	ntuple->annotate(tmp);
	for(int j=0; j<buf.size(); ++j)
		ntuple->appendRow(buf[j].data,NROW);
	buf.clear();
	ntuple->annotate("");
}

static void dumpERowBuf(BLNTuple *ntuple, std::vector<ERow> &buf)
{
	char tmp[64];

	if(ntuple == 0 || buf.size() == 0) return;

	// buf contains rows all with the same event and track
	int event=buf[0].data[8];
	int track=buf[0].data[9];

	if(fileEvent != event && fileEvent != -999) {
		sprintf(tmp,"# End of Event %d",fileEvent);
		ntuple->annotate(tmp);
		ntuple->annotate("");
	}
	fileEvent = event;

	sprintf(tmp,"# Event %d Track %d",event,track);
	ntuple->annotate(tmp);
	for(int j=0; j<buf.size(); ++j)
		ntuple->appendRow(buf[j].data,NEROW);
	buf.clear();
	ntuple->annotate("");
}

static void handleTraceNtuple(BLNTuple *ntuple,double *data,int nData,int src)
{
	static std::vector< std::vector<Row> > rowBuf;
	static std::vector< std::vector<ERow> > eRowBuf;
	static int ranks=0;
	static int ntupleSize=0;

	if(!BLMPI::isRank0() || !ntuple) return;

	if(ranks < 2) {
		ranks = BLMPI::getNranks();
		if(ranks < 2) return;
		ntupleSize = ntuple->getNData();
		if(ntupleSize == sizeof(Row)/sizeof(double)) {
			// entry for rank 0 is never used
			for(int i=0; i<ranks; ++i)
			    rowBuf.push_back(std::vector<Row>());
		} else if(ntupleSize == sizeof(ERow)/sizeof(double)) {
			// entry for rank 0 is never used
			for(int i=0; i<ranks; ++i)
			    eRowBuf.push_back(std::vector<ERow>());
		}
	}

	// flush all buffers
	if(data == 0 && nData == 0) {
		if(rowBuf.size() > 0) {
			for(int i=0; i<ranks; ++i) {
				std::vector<Row> &buf = rowBuf[i];
				dumpRowBuf(ntuple,buf);
			}
		}
		if(eRowBuf.size() > 0) {
			for(int i=0; i<ranks; ++i) {
				std::vector<ERow> &buf = eRowBuf[i];
				dumpERowBuf(ntuple,buf);
			}
		}
		return;
	}

	// Ensure valid src and at least one row in data.
	if(src < 1 || src >= ranks || nData < ntupleSize) return;

	if(ntupleSize == NROW) {
		// All Row-s in buf always have the same event and track values.
		std::vector<Row> &buf = rowBuf[src];
		// If buf is empty, set event and track to values in first 
		// data row, otherwise set them to values in buf.
		int event, track;
		if(buf.size() > 0) {
			event = (int)buf[0].data[8];
			track = (int)buf[0].data[9];
		} else {
			event = (int)data[8];
			track = (int)data[9];
		}
		// loop over data, flushing buf whenever event or track changes
		for(int i=0; i<nData; i+=ntupleSize) {
			double *r = &data[i];
			if((int)r[8] != event || (int)r[9] != track) {
				// never come here when buf is empty, so
				// now buf holds all rows of this event & track
				dumpRowBuf(ntuple,buf);
				event = r[8]; // values from what will become
				track = r[9]; // the first Row in buf
			}
			buf.push_back(*(Row*)r);
		}
	} else if(ntupleSize == NEROW) {
		// All Row-s in buf always have the same event and track values.
		std::vector<ERow> &buf = eRowBuf[src];
		// If buf is empty, set event and track to values in first 
		// data row, otherwise set them to values in buf.
		int event, track;
		if(buf.size() > 0) {
			event = (int)buf[0].data[8];
			track = (int)buf[0].data[9];
		} else {
			event = (int)data[8];
			track = (int)data[9];
		}
		// loop over data, flushing buf whenever event or track changes
		for(int i=0; i<nData; i+=ntupleSize) {
			double *r = &data[i];
			if((int)r[8] != event || (int)r[9] != track) {
				// never come here when buf is empty, so
				// now buf holds all rows of this event & track
				dumpERowBuf(ntuple,buf);
				event = r[8]; // values from what will become
				track = r[9]; // the first ERow in buf
			}
			buf.push_back(*(ERow*)r);
		}
	}
}



int BeamTracksMsgCount = 0;

/// determine the number of tracks to pack into a BeamTrackMsg
int nPackTracks()
{
	static int current=1; // start with 1 track/msg
	static long last=0L;
	static int increase=0, decrease=0;
	if(last == 0L) {
		last = BLTime::timems();
		increase = MPI_GoalMsgPerSec + MPI_GoalMsgPerSec/1;
		decrease = MPI_GoalMsgPerSec - MPI_GoalMsgPerSec/2;
		BeamTracksMsgCount = 0;
	}

	// Update current value of packing once per second, trying to keep
	// close to the goal in msgs/sec.
	if(BLTime::timems()-last >= 1000 || BeamTracksMsgCount > increase) {
		int prev=current;
		if(BeamTracksMsgCount > increase) 
			current *= 2;
		else if(BeamTracksMsgCount < decrease)
			current -= 1;
		if(current < 1) current = 1;
		if(current > MPI_MaxTracksPerMsg) current = MPI_MaxTracksPerMsg;
		if(current != prev)
		    printf("Now packing %d Tracks in BeamTracksMsg\n",current);
		BeamTracksMsgCount = 0;
		last = BLTime::timems();
	}
	
	return current;
}

void sendBeamTracks(int dest)
{
	static int nev=0;
	static std::vector<BeamTracksMsg*> msgVector;
	if(msgVector.size() == 0) {
		for(int i=0; i<BLMPI::getNranks(); ++i)
			msgVector.push_back(new BeamTracksMsg());
	}
	BLAssert(dest > 0 && dest < BLMPI::getNranks());
	BeamTracksMsg &msg = *msgVector[dest];

	// make sure the previous message send is complete
	int idle = msg.isIdle();
	BLAssert(idle);

	static BLManager *manager=BLManager::getObject();
	static BLRunManager *runManager = BLRunManager::getObject();

	int n, nPack=nPackTracks();
	for(n=0; n<nPack; ++n) {
		G4Event *event;
		G4Track *track;
		if(!runManager->getNextBeamEventAndTrack(&event,&track))
			break;
		manager->incrEventsProcessed(event->GetEventID());
		msg.tracks[n].PDGid = track->GetDefinition()->GetPDGEncoding();
		msg.tracks[n].eventID = event->GetEventID();
		msg.tracks[n].trackID = manager->getPrimaryTrackID();
		msg.tracks[n].parentID = manager->getPrimaryParentID();
		msg.tracks[n].secondaryTrackID = 
					manager->getNextSecondaryTrackID();
		G4ThreeVector pos = track->GetPosition();
		msg.tracks[n].position[0] = pos[0];
		msg.tracks[n].position[1] = pos[1];
		msg.tracks[n].position[2] = pos[2];
		G4ThreeVector dir = track->GetMomentumDirection();
		msg.tracks[n].direction[0] = dir[0];
		msg.tracks[n].direction[1] = dir[1];
		msg.tracks[n].direction[2] = dir[2];
		G4ThreeVector pol = track->GetPolarization();
		msg.tracks[n].polarization[0] = pol[0];
		msg.tracks[n].polarization[1] = pol[1];
		msg.tracks[n].polarization[2] = pol[2];
		msg.tracks[n].kineticEnergy = track->GetKineticEnergy();
		msg.tracks[n].time = track->GetGlobalTime();
		msg.tracks[n].weight = track->GetWeight();
		delete event;
		delete track;
	}

	log("Send %d Beam Tracks to %d\n",n,dest);
	monitor.incrState(1);
	msg.send(dest,n);	// n=0 is OK, indicates no more beam tracks
	++BeamTracksMsgCount;
	monitor.incrState(-1);
	if(n == 0) {
		monitor.setState(4);
		if(!sendKeepAlive) {
			// tell all workers to send KeepAliveMsg-s to rank 0
			for(int i=1; i<BLMPI::getNranks(); ++i) {
				KeepAliveMsg m;
				m.sendNoWait(i);
			}
			sendKeepAlive = true;
		}
	}
}

void processOneTrackAndAllSecondaries(Track &track)
{
	static BLManager *manager=0;
	static BLRunManager *runManager = 0;
	if(manager == 0) {
		manager = BLManager::getObject();
		runManager = BLRunManager::getObject();
	}

	G4ThreeVector pos(track.position[0],track.position[1],
							track.position[2]);
	G4ThreeVector dir(track.direction[0],track.direction[1],
							track.direction[2]);
	G4ParticleDefinition *particle = 
		G4ParticleTable::GetParticleTable()->FindParticle(track.PDGid);
	G4DynamicParticle *dyn = new G4DynamicParticle(particle,dir,
							track.kineticEnergy);
	G4Track *g4track = new G4Track(dyn,track.time,pos);
	g4track->SetTrackID(track.trackID);
	g4track->SetParentID(track.parentID);
	BLManager::getObject()->setExternalTrackID(g4track,track.trackID,
							track.parentID);
	g4track->SetWeight(track.weight);
	manager->setEventID(track.eventID);
	runManager->processOneTrack(g4track);
 	runManager->processAllSecondariesAndDeferredTracks(10000,track.trackID);
}

void plotMsg(bool force)
{
	if(BLMPI::getRank() != 0) return;

	static long start = 0L;
	static long current = 0L;
	static long nMsg = 0L;
	static BLNTuple *ntuple=0;
	const char fields[] = "t:MsgPerSec";

	if(ntuple == 0) {
		start = current = BLTime::time();
		ntuple = BLNTuple::create("","","MPI_MsgPerSec",fields,"");
		nMsg = 0;
	}

	++nMsg;
	long now = BLTime::time();
	if(force && now == current) ++now;	// force=true is final call
	if(now == current) return;

	while(current < now) {
		double data[2];
		data[0] = current - start;
		data[1] = nMsg - 1;
		ntuple->appendRow(data,2);
		++current;
		nMsg = 1L;
	}
}


void SumData::add(double d[], int nd)
{
	if(data == 0) {
		data = new double[nd];
		nData = nd;
		for(int i=0; i<nd; ++i)
			data[i] = 0.0;
	}
	BLAssert(nData == nd);
	for(int i=0; i<nData; ++i)
		data[i] += d[i];
	++nSum;
}

void SumData::getSum(double d[], int nd)
{
	if(data == 0) {
		G4Exception("SumData","No Workers Contributed Data",
								JustWarning,"");
		for(int i=0; i<nd; ++i)
			d[i] = 0;
		return;
	}
	BLAssert(nData == nd);
	if(BLMPI::getNranks()-1 != nSum) {
		char tmp[64];
		sprintf(tmp,"%d workers out of %d",nSum,BLMPI::getNranks()-1);
		G4Exception("SumData","Not All Workers Contributed Data",
							JustWarning,tmp);
	}
	for(int i=0; i<nd; ++i)
		d[i] = data[i];
}

bool handleSumMessage(MPI_Status &status)
{
	int tag=status.MPI_TAG;
	if(tag >= FIRST_SUM_TAG && tag < FIRST_NTUPLE_HANDLE) {
		int nData;
		MPI_Get_count(&status,MPI_DOUBLE,&nData);
		double *data = new double[nData];
		MPI_Recv(data,nData,MPI_DOUBLE,status.MPI_SOURCE,tag,
						MPI_COMM_WORLD, &status); 
		sumData[tag].add(data,nData); // create if not found
		delete data;
		return true;
	}

	return false;
}


MPINTuple::MPINTuple(int _handle, G4String _name, G4String _fields) 
				: BLNTuple(_name,_fields)
{
	handle = _handle;
	fields = _fields;
	nfields = 1;
	for(unsigned i=0; i<fields.size(); ++i)
		if(fields(i) == ':') ++nfields;
	nBuffer = MPI_PackNTuples;
	if(nBuffer > MAX_NTUPLE_DOUBLES/nfields)
		nBuffer = MAX_NTUPLE_DOUBLES/nfields;
	if(nBuffer < 1) nBuffer = 1;
	buffer = new double[nBuffer*nfields];
	iBuffer = 0;
}

void MPINTuple::appendRow(double data[], int n)
{
	static BLManager *manager=BLManager::getObject();

	if(handle == 0 || BLMPI::isRank0() || manager->getState() == TUNE || 
					manager->getState() == REFERENCE) 
		return;
	BLAssert(iBuffer < nBuffer);
	int j = iBuffer*nfields;
	for(int i=0; i<n; ++i)
		buffer[j++] = data[i];
	if(++iBuffer >= nBuffer) {
		Msg::sendNTupleRows(buffer,nBuffer*nfields,handle);
		iBuffer = 0;
	}
}

BLNTuple *MPINTupleHandler::create(G4String type, G4String category,
			G4String name, G4String fields, G4String filename)
{
	if(BLMPI::isRank0()) {
		G4String n = category;
		if(n.size() > 0) n += "/";
		n += name;
		int h = getHandle(n);
		if(h != 0) {
			return getNTuple(h);
		}
		BLNTupleHandler *prev = BLNTuple::getForceHandler();
		BLNTuple::setForceHandler(0);
		BLNTuple *ntuple = BLNTuple::create(type,category,
							name,fields,filename);
		BLAssert(ntuple != 0);
		BLNTuple::setForceHandler(prev);
		name2handle[n] = nextHandle;
		handle2NTuple[nextHandle] = ntuple;
		++nextHandle;
		return ntuple;
	}

	// send a CreateNTuple message to rank 0
	CreateNTuple c;
	char tmp[32];
	if(strcpy(tmp,"type")==0 || type.size() >= sizeof(c.type) ||
	   strcpy(tmp,"category")==0 || category.size() >= sizeof(c.category) ||
	   strcpy(tmp,"name")==0 || name.size() >= sizeof(c.name) ||
	   strcpy(tmp,"filename")==0 || filename.size() >= sizeof(c.filename) ||
	   strcpy(tmp,"fields")==0 || fields.size() >= sizeof(c.fields)) {
		strcat(tmp," too large");
		G4Exception("MPINTupleHandler::create",
				"Invalid NTuple Definition",FatalException,tmp);
	}
	strcpy(c.type,type);
	strcpy(c.category,category);
	strcpy(c.name,name);
	strcpy(c.filename,filename);
	strcpy(c.fields,fields);
	c.sendTo0();

	// wait for and receive the NTupleHandle response
	MPI_Status status;
	NTupleHandle h;
	h.recv(status);
	return new MPINTuple(h.handle,name,fields);
}

#else // G4BL_MPI

void BLMPI::init(std::vector<G4String> &argV)
{
}

bool BLMPI::isMPI()
{
	return false;
}

void BLMPI::scan(bool wait) {
}

void BLMPI::rank0BufferPrintf(G4String filename, const char *format, ...)
{
	va_list ap;
	va_start(ap,format);
	FILE *f=BLWriteAsciiFile::fopen(filename); // only calls fopen once
	vfprintf(f,format,ap);
}

void BLMPI::rank0BufferFlush()
{
}

bool BLMPI::sumAllWorkers(double data[], int nData)
{
	return true;
}

void BLMPI::main()
{
}

void BLMPI::closeupAndExit(int value)
{
}

void BLMPI::sendExceptionMessage(const char *origin, const char *code,
		const char *severity, const char *description,
		int eventID, int trackID, const char *particleName,
		G4double kineticEnergy, bool abortProgram)
{
}

#endif // G4BL_MPI
