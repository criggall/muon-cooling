//	BLManager.hh
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

#ifndef BLMANAGER_HH
#define BLMANAGER_HH

#include <vector>
#include <map>
#include <set>
#include "G4VUserDetectorConstruction.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VExceptionHandler.hh"
#include "G4Material.hh"

#include "BLCallback.hh"
#include "BLRunManager.hh"
#include "BLPhysics.hh"
#include "BLZStep.hh"
#include "BLUserCode.hh"
#include "BLSourceRun.hh"

class BLBeam;

#undef IDLE
#undef VISUAL
#undef TUNE
#undef REFERENCE
#undef BEAM
#undef SPECIAL
#undef SOURCE

// State SPECIAL is used for pre-tracking in collective mode -- all
// entries into NTuples should be omitted in SPECIAL mode.
enum BLManagerState { IDLE, VISUAL, TUNE, REFERENCE, BEAM, SPECIAL, SOURCE };
enum VerboseFormat { TAG,NSTEP,GLOBAL,CL,CLX,KE,STEP,VOL,PROCESS,B,E,MAT,
		     P, ID, PART, SEG, WT, POLAR, NEWLINE, EXT };
enum PRNGSeedMethod { EVENT_NUMBER, NO_SEED, TIME_US };

/**	BLManager is the overall manager for g4beamline, managing all aspects
 *	of the run.
 *
 * 	Note it is uses all the G4 user action classes, implemented as classes
 *	under BLManager.
 *
 *	TrackID-s: It is non-trivial to preserve TrackID-s from input files,
 *	because Geant4 considers TrackID to be an internal variable, and it
 *	assigns them in a non-customizable manner. So BLManager keeps a
 *	trackIDMap[] that converts from internal (Geant4) to external (User)
 *	TrackID-s. It also uses BLTrackInfo to associate the external IDs
 *	with the track -- the external trackID is determined either in BLBeam
 *	(when the track is created), or in BLManager::PreUserTrackingAction
 *	for secondaries.
 **/
class BLManager : public G4VUserDetectorConstruction, 
			G4VExceptionHandler {
public:	// user Action classes
	class RunAction {
	public:	virtual void BeginOfRunAction(const G4Run *run) = 0;
	public:	virtual void EndOfRunAction(const G4Run *run) = 0;
	};
	class EventAction {
	public:	virtual void BeginOfEventAction(const G4Event* event) = 0;
	public:	virtual void EndOfEventAction(const G4Event* event) = 0;
	};
	class TrackingAction {
	protected: G4TrackingManager *fpTrackingManager;
	public: void SetTrackingManagerPointer(G4TrackingManager *p) 
		{fpTrackingManager=p;}
	public:	virtual void PreUserTrackingAction(const G4Track *track) = 0;
	public:	virtual void PostUserTrackingAction(const G4Track *track) = 0;
	};
	class SteppingAction {
	public:	virtual void UserSteppingAction(const G4Step *step) = 0;
	};
	class ZSteppingAction {
	public:	virtual void UserZSteppingAction(const G4Track *track) = 0;
	};
	class PrimaryGeneratorAction {
	public:	virtual void GeneratePrimaries(G4Event *event) = 0;
	};
	class StackingAction {
	// returns only fUrgent or fKill.
	public: virtual G4ClassificationOfNewTrack 
					ClassifyNewTrack(const G4Track*) = 0;
	public:	virtual void NewStage() { }
	public:	virtual void PrepareNewEvent() { }
	};
private:
	struct ZStep {
		G4double z;
		ZSteppingAction *action;
		ZStep(G4double _z, ZSteppingAction *a) { z=_z; action=a; }
	};
	static BLManager *blManager;
	static bool initialized;
	BLRunManager *runManager;
	BLManagerState state;
	G4int steppingVerbose;
	unsigned int  beamIndex;
	G4int histoUpdate;
	time_t startProgram;
	time_t startRun;
	time_t startEvent;
	G4int wallClockLimit;
	G4int eventTimeLimit;
	BLPhysics *physics;
	G4VPhysicalVolume *worldPhysicalVolume;
	G4int eventID;
	G4int trackID;
	const G4Track *currentTrack;
	G4bool allExceptions;
	G4int fatalExceptions;
	G4int eventsAborted;
	G4int stuckTracks;
	G4int warnings;
	G4int prevEventID;
	G4int eventsProcessed;
	bool endRun;
	PRNGSeedMethod seedMethod;
	std::set<int> keepEventList;
	std::set<int> skipEventList;
	G4SteppingManager *fpSteppingManager;
	G4TrackingManager *fpTrackingManager;
	std::vector<BLBeam*> beamVector;
	std::vector<BLBeam*> referenceVector;
	std::vector<RunAction*> runActionVector;
	std::vector<RunAction*> beamRunActionVector;
	std::vector<EventAction*> eventActionVector;
	std::vector<EventAction*> beamEventActionVector;
	std::vector<TrackingAction*> trackingActionVector;
	std::vector<BLCallback*> preReferenceCallbackVector;
	std::vector<BLCallback*> postReferenceCallbackVector;
	std::vector<BLCallback*> postTrackingCallbackVector;
	std::vector<BLCallback*> replaceMainLoopCallbackVector;
	std::vector<BLCallback*> visualizationCallbackVector;
	std::vector<BLCallback*> physicsCallbackVector;
	std::vector<SteppingAction*> allStepVector;
	std::map<G4VPhysicalVolume*,SteppingAction*> allStepMap;
	std::map<G4VPhysicalVolume*,SteppingAction*> tpStepMap;
	std::map<G4VPhysicalVolume*,SteppingAction*> rpStepMap;
	std::vector<SteppingAction*> tpStepVector;
	std::vector<SteppingAction*> rpStepVector;
	std::map<G4VPhysicalVolume*,SteppingAction*> beamStepMap;
	std::vector<SteppingAction*> beamStepVector;
	std::vector<BLSourceRun*> sourceRunVector;
	std::vector<int> verboseFormat;
	G4double zTolerance;
	std::vector<ZStep> tuneZStep;
	std::vector<ZStep> referenceZStep;
	std::vector<ZStep> beamZStep;
	std::vector<ZStep> *currentZStep;
	unsigned indexZStep;
	G4double prevZ;
	G4int nStuckSteps;
	std::vector<StackingAction*> stackingActionVector;
	int nextSecondaryTrackID;
	std::map<G4int,G4int> trackIDMap;
	int primaryTrackID;
	int primaryParentID;
	std::vector<BLUserCode*> userCodeVector;
	std::map<G4String,int> exceptionCount;
	BLSourceRun *sourceRun;
	/// private constructor -- immediate construction only.
	BLManager();
	void insertZStep(std::vector<ZStep>& vector, G4double z, ZSteppingAction *action);
public:
	/// getObject() will return a pointer to the single BLManager object,
	/// creating it if necessary. Does only the immediate constructor, not
	/// delayedConstruction(). that means it is Ok to register capabilities,
	/// but not much else.
	static BLManager *getObject();

	/// delayedConstruction() performs things which must wait until all
	/// static initializers have executed (e.g. in Geant4 routines).
	void delayedConstruction();

	/// Destructor.
	~BLManager();

	/// getState() returns the current state.
	BLManagerState getState() { return state; }

	/// setState sets the state
	void setState(BLManagerState _state) { state = _state; }

	/// getSteppingVerbose() returns steppingVerbose. NOTE: use this
	/// during tracking instead of Param.getInt("steppingVerbose");
	int getSteppingVerbose() { return steppingVerbose; }

	/// setSteppingVerbose() updates steppingVerbose. -- NOTE: many
	/// other classes relay on the Parameter, not the valud in this class.
	void setSteppingVerbose(int v) { steppingVerbose = v; }

	/// getEventTimeLimit() returns the CPU time limit for events (seconds).
	/// -1 mean infinite.
	int getEventTimeLimit() { return eventTimeLimit; }

	/// setEventTimeLimit() sets the CPU time limit for events (seconds).
	/// -1 mean infinite.
	void setEventTimeLimit(int sec) { eventTimeLimit = sec; }

	/// getEventID() gets the current eventID;
	G4int getEventID() const { return eventID; }

	/// setEventID() sets the current eventID;
	void setEventID(int evId) { eventID=evId; prevEventID=evId-1; }

	/// skipEvent() determines if this EventID should be skipped. It
	/// checks keepEventList and skipEventList.
	bool skipEvent(int EventID) {
	  if(skipEventList.count(EventID) > 0) return true;
	  return keepEventList.size() > 0 && keepEventList.count(EventID) == 0;
	}

	/// setKeepEvent() puts EventID into the keepEventList.
	/// (if keepEventList is empty, all events are kept.)
	void setKeepEvent(int EventID) { keepEventList.insert(EventID); }

	/// setSkipEvent() puts EventID into the skipEventList.
	void setSkipEvent(int EventID) { skipEventList.insert(EventID); }

	/// incrEventsProcessed() will increment eventsProcessed.
	/// For special uses only (e.g.MPI).
	void incrEventsProcessed(int eventID);

	/// showAllExceptions() sets the flag that prevents thinning out
	/// multiple exception printouts. Returns the previous value.
	bool showAllExceptions(bool value=true)
		{ bool tmp=allExceptions;  allExceptions = value; return tmp; }

	/// registerSteppingAction() registers a SteppingAction to
	/// be called for each step (regardless of state).
	void registerSteppingAction(SteppingAction *sa)
		{ allStepVector.push_back(sa); }

	/// registerSteppingAction() registers a SteppingAction to
	/// be called for each step (regardless of state), for every step
	/// involving the physicalVol.
	/// The callback will be called if physicalVol is either
	/// the pre- or post-step physical volume (once if both).
	/// If physicalVol==0 it is called every reference particle step.
	/// LIMITATION: only one callback can be registered for a given 
	/// physicalVol (except 0).
	void registerSteppingAction(G4VPhysicalVolume *physicalVol,
						SteppingAction *sa)
		{ if(physicalVol != 0) allStepMap[physicalVol] = sa;
		  else allStepVector.push_back(sa); }

	/// registerTuneParticleStep() registers a SteppingAction to
	/// be called for every tune particle step involving the physicalVol.
	/// The callback will be called if physicalVol is either
	/// the pre- or post-step physical volume (once if both).
	/// If physicalVol==0 it is called every tune particle step.
	/// LIMITATION: only one callback can be registered for a given 
	/// physicalVol (except 0).
	void registerTuneParticleStep(G4VPhysicalVolume *physicalVol,
					SteppingAction *sa)
		{ if(physicalVol != 0) tpStepMap[physicalVol] = sa; 
		  else tpStepVector.push_back(sa); }

	/// registerReferenceParticleStep() registers a SteppingAction to
	/// be called for every reference particle step involving the
	/// physicalVol.
	/// The callback will be called if physicalVol is either
	/// the pre- or post-step physical volume (once if both).
	/// If physicalVol==0 it is called every reference particle step.
	/// LIMITATION: only one callback can be registered for a given 
	/// physicalVol (except 0).
	void registerReferenceParticleStep(G4VPhysicalVolume *physicalVol,
					SteppingAction *sa)
		{ if(physicalVol != 0) rpStepMap[physicalVol] = sa; 
		  else rpStepVector.push_back(sa); }

	/// registerBeamStep() registers a SteppingAction to
	/// be called for every beam step involving the physicalVol.
	/// (beam particles are everything except reference and tune.)
	/// The callback will be called if physicalVol is either
	/// the pre- or post-step physical volume (once if both). 
	/// if physicalVol==0 it is called every beam step.
	/// LIMITATION: only one callback can be registered for a given 
	/// physicalVol (except 0).
	void registerBeamStep(G4VPhysicalVolume *physicalVol,
					SteppingAction *sa)
		{ if(physicalVol != 0) beamStepMap[physicalVol] = sa; 
		  else beamStepVector.push_back(sa); }

	/// registerZStep() will force a step to occur near the given z
	/// position, and will call the ZSteppingAction for it, interpolating
	/// to the desired z value (Centerline coords).
	/// when is a bitwise OR of 1=tune, 2=reference, 4=beam.
	void registerZStep(G4double z, ZSteppingAction *sa, G4int when=7);

	/// registerStackingAction() registers a StackingAction to be called
	/// by the Geant4 stacking manager.
	void registerStackingAction(StackingAction *sa)
		{ stackingActionVector.push_back(sa); }

	/// registerSourceRun() will arrange to perform a source run.
	void registerSourceRun(BLSourceRun *sr)
		{ sourceRunVector.push_back(sr); }

	/// setSteppingFormat() sets the verbose printing format according
	/// to parameter "steppingFormat".
	void setSteppingFormat();

	/// getFormatHelp() returns a string with help text about valid format 
	/// items.
	G4String getFormatHelp();

	/// appendVerboseFormat() appends fmt to the format for printing when
	/// param steppingVerbose is nonzero
	void appendVerboseFormat(G4String fmt);

	/// steppingVerbosePrint() will print this step according to the current
	/// verboseFormat.
	/// Prints header if header != 0).
	void steppingVerbosePrint(const G4Step *step, const G4Track *track, int
								header=0);

	/// setSteppingManager() sets the pointer to the current
	/// G4SteppingManager.
	void setSteppingManager(G4SteppingManager *p) { fpSteppingManager = p; }

	/// getSteppingManager() returns a pointer to the current
	/// G4SteppingManager.
	G4SteppingManager *getSteppingManager() { return fpSteppingManager; }

	/// setTrackingManager() sets the pointer to the current
	/// G4TrackingManager.
	void setTrackingManager(G4TrackingManager *p) { fpTrackingManager = p; }

	/// initialize() will initialize the BLManager object, and the 
	/// geant4 kernel, thus constructing the geometry in the
	/// world group. Note that registerPhysics() must be called
	/// before initialize() (normally done by a "physics" command
	/// in the input file).
	void initialize();

	/// isInitialized() returns true if the BLManager has been initialized.
	static bool isInitialized() { return initialized; }

	/// trackTuneAndReferenceParticles() will generate and track the tune
	/// particle and then the reference particle.
	void trackTuneAndReferenceParticles();

	/// handleSourceRun() will perform any source runs that have been
	/// registered;
	void handleSourceRun();

	/// trackBeam() will generate the defined beam and track each event.
	void trackBeam();

	/// displayVisual() will display the detector visually.
	void displayVisual();

	/// displayGeometry() will display the geant4 geometry.
	/// This is a hierarchical ASCII listing of all volumes.
	/// if phys==0 then use the worldPhysicalVolume.
	void displayGeometry(G4VPhysicalVolume *phys=0, int level=0);

	/// registerPhysics() registers the BLPhysics object.
	/// It also sets the physics list to the BLRunManager, so following
	/// commands can find particle by name.
	void registerPhysics(BLPhysics *_physics)
		{ physics = _physics;
		  runManager->SetUserInitialization(physics->getPhysicsList());
		  handleCallbacks(-1);
		}
	
	/// getPhysics returns the registered BLPhysics object.
	BLPhysics *getPhysics() { return physics; }

	/// registerBeam() registers a BLBeam object for beam generation.
	/// Multiple BLBeam objects can be registered, used in order.
	void registerBeam(BLBeam *_beam) { beamVector.push_back(_beam); }

	/// clearBeamVector() will clear the beamVector. used in BLMPI.
	void clearBeamVector() { beamVector.clear(); }

	std::vector<BLBeam*> *getBeamVector() { return &beamVector; }

	/// registerReference() registers a BLBeam object for reference 
	/// particle generation.
	/// Multiple BLBeam objects can be registered, used in order.
	void registerReference(BLBeam *_beam)
		{ referenceVector.push_back(_beam); }

	/// nReference() returns the number of reference particles registered.
	int nReference() { return referenceVector.size(); }

	/// registerRunAction() registers a UserRunAction.
	/// If beamOnly is true (the default), the callback is made only
	/// if state==BEAM.
	void registerRunAction(RunAction *a, G4bool beamOnly=true)
		{ if(beamOnly)
			beamRunActionVector.push_back(a);
		  else
			runActionVector.push_back(a);
		}

	/// registerEventAction() registers a UserEventAction.
	/// If beamOnly is true (the default), the callback is made only
	/// if state==BEAM.
	void registerEventAction(EventAction *a, G4bool beamOnly=true)
		{ if(beamOnly)
			beamEventActionVector.push_back(a);
		  else
			eventActionVector.push_back(a);
		}

	/// registerTrackingAction() registers a UserTackingAction.
	/// By default puts new entry at the back end of the vector;
	/// if front is true, puts new entry at the front.
	void registerTrackingAction(TrackingAction *a, bool front=false)
	    { if(front)
		trackingActionVector.insert(trackingActionVector.begin(),a); 
	      else
		trackingActionVector.push_back(a); 
	    }

	/// registerUserCode() registers a BLUserCode instance.
	void registerUserCode(BLUserCode *instance)
		{ userCodeVector.push_back(instance); }

	/// getUserCodeInstances() returns a vector of all registered
	/// instances of BLUserCode with the specified type.
	std::vector<BLUserCode*> getUserCodeInstances(G4String type)
		{ std::vector<BLUserCode*> ret;
		  for(unsigned i=0; i<userCodeVector.size(); ++i) {
			if(type == userCodeVector[i]->getType())
				ret.push_back(userCodeVector[i]);
		  }
		  return ret;
		}


	/// registerCallback() registers a BLCallback.
	/// type=0 for pre-Tune particle,
	/// type=1 for post-Reference (pre-beam tracking),
	/// type=2 for post-beam tracking.
	/// type=3 for replacing the main program loop.
	/// type=4 for visualization.
	/// type=-1 for inside physics definition
	/// NOTE: if there are type=3 callbacks, when the last one returns
	/// the closes up by summarizing NTuples and callng handleCallbacks(2),
	/// and then the program exits. This prevents the main program loop from
	/// executing. type=3 callbacks are called just after the type=1
	/// callbacks (i.e. after the post-Reference callbacks).
	void registerCallback(BLCallback *cb, int type) {
		if(type==0) preReferenceCallbackVector.push_back(cb);
		else if(type==1) postReferenceCallbackVector.push_back(cb);
		else if(type==2) postTrackingCallbackVector.push_back(cb);
		else if(type==3) replaceMainLoopCallbackVector.push_back(cb);
		else if(type==4) visualizationCallbackVector.push_back(cb);
		else if(type==-1) physicsCallbackVector.push_back(cb);
	}

	/// handleCallbacks() calls all applicable registered callbacks.
	/// type=0 for pre-reference particle,
	/// type=1 for post-center (pre-beam tracking),
	/// type=2 for post-beam tracking (just before program exit).
	/// type=3 for replacing the main program loop.
	/// type=4 for visualization.
	/// type=-1 for inside physics definition
	void handleCallbacks(int type);

	/// getWorldPhysicalVolume() returns a pointer to the world PV.
	/// Note that it must already have been consructed, so this function
	/// returns NULL before construct() is called (by main).
	G4VPhysicalVolume *getWorldPhysicalVolume()
		{ return worldPhysicalVolume; }

	/// setPRNGSeedMethod() will set the method used to seed the
	/// pseudo random number generator at the start of each event.
	void setPRNGSeedMethod(PRNGSeedMethod method)
		{ seedMethod = method; }

	/// getPRNGSeedMethod() will return the method used to seed the
	/// pseudo random number generator at the start of each event.
	PRNGSeedMethod getPRNGSeedMethod() { return seedMethod; }

	// virtual functions from the geant4 base classes. All of them
	// call the registered actions, as appropriate for state and
	// the current physical volume (if applicable). Due to the design
	// of the Geant4 callback classes, these are mediated by classes of
	// the form BLManager_name.

	/// UserSteppingAction() from G4UserSteppingAction.
	void UserSteppingAction(const G4Step *step);

	/// Construct() from G4VUserDetectorConstruction.
	G4VPhysicalVolume *Construct();

	/// PreUserTrackingAction() from G4UserTrackingAction.
	void PreUserTrackingAction(const G4Track *track);

	/// PostUserTrackingAction() from G4UserTrackingAction.
	void PostUserTrackingAction(const G4Track *track);

	/// BeginOfRunAction() from G4UserRunAction.
	void BeginOfRunAction(const G4Run *run);

	/// EndOfRunAction() from G4UserRunAction.
	void EndOfRunAction(const G4Run *run);

	/// BeginOfEventAction() from G4UserEventAction.
	void BeginOfEventAction(const G4Event* event);

	/// EndOfEventAction() from G4UserEventAction.
	void EndOfEventAction(const G4Event* event);

	/// GeneratePrimaries() from G4VUserPrimaryGeneratorAction.
	void GeneratePrimaries(G4Event *event);

	/// ClassifyNewTrack() from G4UserStackingAction.
	G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);

	/// NewStage() from G4userStackingAction.
	void NewStage();

	/// PrepareNewEvent() from G4StackingAction.
	void PrepareNewEvent();

	/// clearTrackIDMap() clears the TrackID map.
	void clearTrackIDMap() { trackIDMap.clear(); }

	/// setNextSecondaryTrackID() sets the external TrackID for the
	/// next secondary track. Automatically incremented for subsequent
	/// secondaries.
	void setNextSecondaryTrackID(int next) { nextSecondaryTrackID = next; }
	int getNextSecondaryTrackID() { return nextSecondaryTrackID; }

	/// getExternalTrackID() returns the external TrackID for the given
	/// track.
	/// In collective mode, internal and external trackID-s are the same.
	int getExternalTrackID(const G4Track *track);

	/// getExternalParentID() returns the external ParentID for the given
	/// track.
	/// In collective mode, internal and external trackID-s are the same.
	int getExternalParentID(const G4Track *track);

	/// setExternalTrackID() will set the external trackID and parentID.
	/// if trackID<0 uses nextSecondarTrackID++.
	void setExternalTrackID(G4Track *track, int trackID, int parentID);

	/// getPrimaryTrackID() returns the primaryTrackID set by the beam command
	int getPrimaryTrackID() { return primaryTrackID; }

	/// getPrimaryParentID() returns the primaryParentID set by the beam command
	int getPrimaryParentID() { return primaryParentID; }

	/// setPrimaryTrackID() sets track and parent IDs for a primary track.
	void setPrimaryTrackID(int t, int p)
		{ primaryTrackID = t; primaryParentID = p; 
		  trackIDMap[1] = t;
		}

	/// Notify() from G4VExceptionHandler.
	G4bool Notify(const char* originOfException, const char* exceptionCode,
                        G4ExceptionSeverity severity, const char* description);

	/// printException() does the actual printing of a G4Exception.
	/// calls tallyException().
	void printException(const char *origin, const char *code,
		const char *severity, const char *description,
		int evid, int trkid, const char *particleName,
		G4double kineticEnergy, bool toBeAborted, bool printOnce=false);

	/// tallyException() does the actual counting of a G4Exception, for
	/// printing at program exit.
	int tallyException(const char *code);

	/// exceptionSummary() prints a summary of all exceptions
	void exceptionSummary();

#ifdef G4BL_GSL
	/// GSL error handler
	static void gsl_error_handler(const char *reason, const char *file,
						int lineno, int gsl_errno);
#endif
};

#endif // BLMANAGER_HH
