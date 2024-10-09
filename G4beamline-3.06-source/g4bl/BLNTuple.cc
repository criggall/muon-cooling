//	BLNTuple.cc
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

#include <vector>
#include <map>
#include <set>
#include <time.h>

#include "BLAssert.hh"
#include "BLNTuple.hh"
#include "BLParam.hh"
#include "BLTrackFile.hh"
#include "BLFOR009.hh"
#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLWriteAsciiFile.hh"
#ifdef USE_SHARED_OBJECTS
#include "BLLoad.hh"
#endif // USE_SHARED_OBJECTS

bool BLNTuple::enableOutput = true;
#ifdef G4BL_ROOT
G4String BLNTuple::defaultType = "root";
#else
G4String BLNTuple::defaultType = "ascii";
#endif
std::map<G4String,BLNTupleHandler*> *BLNTuple::handler = 0;
std::vector<BLNTuple *> BLNTuple::list;
std::set<G4String> BLNTuple::nameSet;
BLNTupleHandler *BLNTuple::forceHandler = 0;

const char TrackFields[] =
    "x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight";
const int NTrackFields=12;
const char TraceFields[] =
    "x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight:Bx:By:Bz:Ex:Ey:Ez";
const int NTraceFields=18;


void BLNTuple::registerHandler(G4String typeName, BLNTupleHandler *h)
{
	// permit this to be called during initializations, without
	// worrying about ordering.
	if(!handler)
		handler = new std::map<G4String,BLNTupleHandler*>();

	for(unsigned i=0; i<typeName.size(); ++i)
		typeName[i] = tolower(typeName[i]);

	(*handler)[typeName] = h;
}

BLNTupleHandler *BLNTuple::getHandler(G4String type)
{
	if(forceHandler) return forceHandler;

	for(unsigned i=0; i<type.size(); ++i)
		type[i] = tolower(type[i]);

	if(!handler) return 0;
	return (*handler)[type]; // returns 0 for unknown type
}

BLNTuple *BLNTuple::create(G4String type, G4String category, G4String name,
					G4String fields, G4String filename)
{
	if(type == "") type = defaultType;
	if(!enableOutput) type = "Dummy";

	BLNTupleHandler *h = getHandler(type);
	if(h) {
		G4String fullName=category+"/"+name;
		if(category == "") fullName=name;
		if(!h->ignoreMultipleCreates() && !isUniqueName(fullName))
			G4Exception("BLNTuple","Duplicate NTuple",JustWarning,
								fullName);

		BLNTuple *p = h->create(type,category,name,fields,filename);
		if(p) list.push_back(p);
		return p;
	}

	char tmp[128];
	sprintf(tmp,"unknown format '%s'",type.c_str());
	G4Exception("BLNTuple",tmp,FatalException, "");

	return 0;
}

BLNTuple *BLNTuple::read(G4String type, G4String category, G4String name,
					G4String fields, G4String filename)
{
	if(type == "") type = defaultType;

	BLNTupleHandler *h = getHandler(type);
	if(h) {
		BLNTuple *p = h->read(type,category,name,fields,filename);
		return p;
	}

	char tmp[128];
	sprintf(tmp,"unknown format '%s'",type.c_str());
	G4Exception("BLNTuple",tmp,FatalException, "");

	return 0;
}

void BLNTuple::flushAll()
{
	for(unsigned i=0; i<list.size(); ++i)
		list[i]->flush();
}

void BLNTuple::closeAll()
{
	for(unsigned i=0; i<list.size(); ++i)
		list[i]->close();
}

void BLNTuple::summary()
{
	for(unsigned i=0; i<list.size(); ++i)
		list[i]->doSummary();
}

void BLNTuple::update()
{
	// no-op for now....
}

G4String BLNTuple::getFormatList()
{
	G4String ret;
	if(handler == 0) return ret;
	std::map<G4String,BLNTupleHandler*>::iterator i;
	for(i=(*handler).begin(); i!=(*handler).end(); ++i) {
		if(i->first == "Dummy") continue;
		if(i != (*handler).begin()) ret += " ";
		ret += i->first;
	}
	return ret;
}

bool BLNTuple::isUniqueName(G4String name)
{
	if(nameSet.count(name) != 0) return false;
	nameSet.insert(name);
	return true;
}

void BLNTuple::doCallbacks(double data[], int ndata)
{
	for(unsigned i=0; i<callback.size(); ++i)
		callback[i]->ntupleCallback(this,data,ndata);
}

#ifdef USE_SHARED_OBJECTS

// auto-loading handler for HistoScope NTuples
class HistoHandlerLoader : public BLNTupleHandler {
	bool loadHistoSharedObject(G4String type) {
		if(!BLLoad::load("g4bl-histoscope")) {
		    G4Exception("BLNTuple","HistoScope not supported",
		    					FatalException, "");
		    return false;
		}
		// verify the loaded .so registered this type
		if(BLNTuple::getHandler(type) == this) {
		    G4Exception("BLNTuple","HistoScope INTERNAL ERROR",
		    					FatalException, "");
		    return false;
		}
		return true;
	}
public:
	HistoHandlerLoader() {
		if(BLNTuple::getHandler("histo") == 0) {
			BLNTuple::registerHandler("histo",this); 
			BLNTuple::registerHandler("HistoScope",this); 
			BLNTuple::registerHandler("histoscope",this); 
		}
	}
	BLNTuple *create(G4String type, G4String category,
			G4String name, G4String fields, G4String filename) {
		if(!loadHistoSharedObject(type)) return 0;
		return BLNTuple::getHandler(type)->create(type,category,name,
						fields,filename);
	}
	BLNTuple *read(G4String type, G4String category, G4String name,
					G4String fields, G4String filename) {
		if(!loadHistoSharedObject(type)) return 0;
		return BLNTuple::getHandler(type)->read(type,category,name,
						fields, filename);
	}
};
HistoHandlerLoader defaultHistoHandlerLoader;


// auto-loading handler for Root NTuples
class RootHandlerLoader : public BLNTupleHandler {
	bool loadRootSharedObject(G4String type) {
		if(!BLLoad::load("g4bl-root")) {
		    G4Exception("BLNTuple","Root not supported",
		    					FatalException, "");
		    return false;
		}
		// verify the loaded .so registered this type
		if(BLNTuple::getHandler(type) == this) {
		    G4Exception("BLNTuple","Root INTERNAL ERROR",
		    					FatalException, "");
		    return false;
		}
		return true;
	}
public:
	RootHandlerLoader() {
		if(BLNTuple::getHandler("root") == 0) {
			BLNTuple::registerHandler("root",this); 
			BLNTuple::registerHandler("Root",this); 
		}
	}
	BLNTuple *create(G4String type, G4String category,
			G4String name, G4String fields, G4String filename) {
		if(!loadRootSharedObject(type)) return 0;
		return BLNTuple::getHandler(type)->create(type,category,name,
						fields,filename);
	}
	BLNTuple *read(G4String type, G4String category, G4String name,
					G4String fields, G4String filename) {
		if(!loadRootSharedObject(type)) return 0;
		return BLNTuple::getHandler(type)->read(type,category,name,
						fields,filename);
	}
};
RootHandlerLoader defaultRootHandlerLoader;

#endif // USE_SHARED_OBJECTS



/**	class DummyNTuple implements a Dummy NTuple.
 **/
class DummyNTuple : public BLNTuple {
public:
	DummyNTuple() : BLNTuple("","") { }
	~DummyNTuple() { }
	virtual void appendRow(double data[], int n) { }
	virtual bool readRow(double data[], int ndata) { return false; }
	virtual int getNData() { return 1; }
	void flush() { }
	void close() { }
	void doSummary() { }
};

class DummyNTupleHandler : public BLNTupleHandler {
public:
	DummyNTupleHandler() { 
		BLNTuple::registerHandler("Dummy",this);
	}
	BLNTuple *create(G4String type, G4String category,
			G4String name, G4String fields, G4String filename) {
		return new DummyNTuple();
	}
	BLNTuple *read(G4String type, G4String category, G4String name,
					G4String fields,G4String filename) {
		return 0;
	}
};
DummyNTupleHandler defaultDummyNTupleHandler;



/**	class TrackFileNTuple implements a BLTrackFile NTuple.
 **/
class TrackFileNTuple : public BLNTuple {
	BLTrackFile *tf;
	G4String title;
	int entries;
	static std::vector<TrackFileNTuple*> list;
public:
	TrackFileNTuple(G4String category, G4String name, G4String fields,
							G4String filename);
	~TrackFileNTuple() { close(); }
	virtual void appendRow(double data[], int n);
	virtual bool readRow(double data[], int ndata) { return false; }
	void flush() { if(tf) tf->flush(); }
	void close() { if(tf) delete tf; tf = 0; }
	void doSummary() {
		printf("NTuple %-15s %8d entries\n",title.c_str(), entries);
	}
};
std::vector<TrackFileNTuple*> TrackFileNTuple::list;

class TrackFileNTupleHandler : public BLNTupleHandler {
public:
	TrackFileNTupleHandler() { 
		BLNTuple::registerHandler("TrackFile",this);
		BLNTuple::registerHandler("BLTrackFile",this);
		BLNTuple::registerHandler("trackfile",this);
	}
	BLNTuple *create(G4String type, G4String category,
			G4String name, G4String fields, G4String filename) {
		return new TrackFileNTuple(category,name,fields,filename);
	}
	BLNTuple *read(G4String type, G4String category, G4String name,
					G4String fields,G4String filename) {
		return 0;
	}
};
TrackFileNTupleHandler defaultTrackFileNTupleHandler;

TrackFileNTuple::TrackFileNTuple(G4String category, G4String name,
		G4String _fields, G4String filename) : BLNTuple(name,_fields)
{
	nData = NTrackFields;
	fields = TrackFields;
	if(filename == "") filename = name+".txt";
	if(filename.find(".txt") == filename.npos)
		filename += ".txt";
	tf = new BLTrackFile(filename,category+"/"+name,"w");
	title = name;
	entries = 0;
	list.push_back(this);
}

void TrackFileNTuple::appendRow(double data[], int n)
{
	BLAssert(n == NTrackFields);

	// make values too small for a float be exactly zero
	for(int i=0; i<n; ++i) {
		if(fabs(data[i]) < 1.0E-44) data[i] = 0.0;
	}

	G4ThreeVector pos(data[0]*mm,data[1]*mm,data[2]*mm);
	G4ThreeVector momentum(data[3]*MeV,data[4]*MeV,data[5]*MeV);
	G4double time = data[6]*ns;
	int PDGid = (int)data[7];
	int eventID = (int)data[8];
	int trackID = (int)data[9];
	int parentID = (int)data[10];
	double weight = data[11];

	tf->write(pos,time,momentum,PDGid,eventID,trackID,parentID,weight);

	++entries;

	doCallbacks(data,n);
}





/**	class FOR009NTuple implements a FOR009.DAT NTuple.
 *
 *	If multiple instances specify different filename-s, the first one wins.
 **/
class FOR009NTuple : public BLNTuple {
	G4String title;
	int entries;
	static BLFOR009 *for009file;
	static G4String for009filename;
	static std::vector<FOR009NTuple*> for009list;
public:
	FOR009NTuple(G4String category, G4String name, G4String fields,
							G4String filename);
	~FOR009NTuple() { }
	virtual void appendRow(double data[], int n);
	virtual bool readRow(double data[], int ndata) { return false; }
	virtual bool needsReference() { return true; }
	void flush() { }
	void close();
	void doSummary() { /* close provides summary per region */ }
};
BLFOR009 *FOR009NTuple::for009file = 0;
G4String FOR009NTuple::for009filename = "";
std::vector<FOR009NTuple*> FOR009NTuple::for009list;

class FOR009NTupleHandler : public BLNTupleHandler {
public:
	FOR009NTupleHandler() { 
		BLNTuple::registerHandler("FOR009",this);
		BLNTuple::registerHandler("FOR009.DAT",this);
		BLNTuple::registerHandler("for009",this);
		BLNTuple::registerHandler("for009.dat",this);
	}
	BLNTuple *create(G4String type, G4String category,
			G4String name, G4String fields, G4String filename) {
		return new FOR009NTuple(category,name,fields,filename);
	}
	BLNTuple *read(G4String type, G4String category, G4String name,
					G4String fields,G4String filename) {
		return 0;
	}
};
FOR009NTupleHandler defaultFOR009NTupleHandler;


FOR009NTuple::FOR009NTuple(G4String category, G4String name,
		G4String _fields, G4String filename) : BLNTuple(name,_fields)
{
	nData = NTraceFields;
	fields = TraceFields;
	G4String file=filename;
	if(file == "") file = name;
	if(file.find(".txt") == file.npos)
		file += ".txt";
	if(for009file == 0) {
		for009filename = file;
		for009file = new BLFOR009(file,category+"/"+name,"w");
	} else if(filename != "" && file != for009filename) {
		G4Exception("BLNTuple-FOR009","Multiple Filenames",JustWarning,
								"");
	}
	title = name;
	entries = 0;
	for009list.push_back(this);
}

void FOR009NTuple::appendRow(double data[], int n)
{
	BLAssert(n == NTraceFields);

	G4ThreeVector pos(data[0]*mm,data[1]*mm,data[2]*mm);
	G4ThreeVector momentum(data[3]*MeV,data[4]*MeV,data[5]*MeV);
	G4ThreeVector Bfield(data[12]*tesla,data[13]*tesla,data[14]*tesla);
	G4ThreeVector Efield(data[15]*(megavolt/meter),
			data[16]*(megavolt/meter),data[17]*(megavolt/meter));
	G4double time = data[6]*ns;
	int PDGid = (int)data[7];
	int eventID = (int)data[8];
	int trackID = (int)data[9];
	int parentID = (int)data[10];
	double weight = data[11];

	for009file->write(pos,time,momentum,Bfield,Efield,PDGid,eventID,trackID,
							parentID,weight);

	++entries;

	doCallbacks(data,n);
}

void FOR009NTuple::close()
{
	if(for009file) for009file->close();
	for009file = 0;
}




/**	class AsciiNTuple implements a generic ASCII NTuple.
 **/
const char BLTRACKFILE[]="x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight";
const char BLTRACKFILE2[]="x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight:Bx:By:Bz:Ex:Ey:Ez:ProperTime:PathLength:PolX:PolY:PolZ:InitX:InitY:InitZ:InitT:InitKE";
class AsciiNTuple : public BLNTuple {
	FILE *fd;
	G4String title;
	int entries;
	bool *intFields;	// array indicating which fields are integers
public:
	AsciiNTuple(G4String category, G4String name, G4String fields,
							G4String filename);
	~AsciiNTuple() { close(); }
	virtual void appendRow(double data[], int n);
	virtual bool readRow(double data[], int ndata) { return false; }
	virtual void annotate(G4String line) 
		{ fprintf(fd,"%s\n",line.c_str()); }
	void flush() { fflush(fd); }
	void close() { if(fd != 0) BLWriteAsciiFile::fclose(fd); fd = 0; }
	void doSummary() {
		printf("NTuple %-15s %8d entries\n",title.c_str(), entries);
	}
};

class AsciiNTupleHandler : public BLNTupleHandler {
public:
	AsciiNTupleHandler() { 
		BLNTuple::registerHandler("Ascii",this);
		BLNTuple::registerHandler("ascii",this);
		BLNTuple::registerHandler("ASCII",this);
	}
	BLNTuple *create(G4String type, G4String category,
			G4String name, G4String fields, G4String filename) {
		return new AsciiNTuple(category,name,fields,filename);
	}
	BLNTuple *read(G4String type, G4String category, G4String name,
					G4String fields,G4String filename) {
		return 0;
	}
};
AsciiNTupleHandler defaultAsciiNTupleHandler;

static bool endsWith(G4String s, G4String v)
{
	return s.size() >= v.size() && s.substr(s.size()-v.size()) == v;
}

AsciiNTuple::AsciiNTuple(G4String category, G4String name,
		G4String _fields, G4String filename) : BLNTuple(name,_fields)
{
	if(filename == "") filename = name+".txt";
	if(filename.find(".txt") == filename.npos)
		filename += ".txt";
	fd = BLWriteAsciiFile::fopen(filename);
	if(!fd) G4Exception("AsciiNTuple","Cannot create file",FatalException,
					filename);
	title = name;
	entries = 0;
	fields = _fields;

	if(fields == BLTRACKFILE) {
		fprintf(fd,"#BLTrackFile %s\n",name.c_str());
	} else if(fields == BLTRACKFILE2) {
		fprintf(fd,"#BLTrackFile2 %s\n",name.c_str());
	} else {
		fprintf(fd,"# %s\n",name.c_str());
	}

	G4String s(fields);
	nData = 1;
	for(;;) {
		G4String::size_type p = s.find(":");
		if(p == s.npos) break;
		s.replace(p,1," ");
		++nData;
	}
	fprintf(fd,"#%s\n",s.c_str());

	// initialize intFields[] (compare via endsWith() for ntuple command)
	std::vector<G4String> v=BLCommand::splitString(fields,":");
	assert(v.size() == (unsigned)nData);
	intFields = new bool[nData];
	for(int i=0; i<nData; ++i)
		intFields[i] = (endsWith(v[i],"EventID") ||
				endsWith(v[i],"TrackID") ||
				endsWith(v[i],"ParentID") || 
				endsWith(v[i],"PDGid"));
}

void AsciiNTuple::appendRow(double data[], int n)
{
	assert(n == nData);
	for(int i=0; i<n; ++i) {
		if(intFields[i]) {
			fprintf(fd,"%d ",(int)data[i]);
		} else {
			// make values too small for a float be exactly zero
			if(fabs(data[i]) < 1.0E-44) data[i] = 0.0;
			fprintf(fd,"%.6g ",data[i]);
		}
	}
	fprintf(fd,"\n");
	++entries;

	doCallbacks(data,n);
}

