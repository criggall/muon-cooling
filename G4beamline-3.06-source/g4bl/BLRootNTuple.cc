//	RootNTuple.cc
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

#ifdef G4BL_ROOT

#include <stdio.h>
#include <vector>
#include <time.h>

#include "BLAssert.hh"
#include "BLNTuple.hh"
#include "BLParam.hh"
#include "BLTrackFile.hh"
#include "BLFOR009.hh"
#include "BLCommand.hh"

#include "TROOT.h"
#include "TPluginManager.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TDirectory.h"


/**	class RootNTuple implements a Root NTuple.
 **/
class RootNTuple : public BLNTuple {
	G4String name;
	int nEntries;
	G4String title;
	TNtuple *ntuple;
	TFile *readFile;
	int readIndex;
	float *readData;
	static TFile *file;
	static G4String filename;
	static RootNTuple *first;
	friend class RootNTupleHandler;
public:
	RootNTuple(G4String category, G4String name, G4String fields);
	RootNTuple(TFile *file, TNtuple *nt, G4String fields);
	~RootNTuple() { 
		if(readFile) {
			readFile->Close();
			delete readFile;
			readFile=0; 
		}
		if(ntuple) {
			delete ntuple;
			ntuple = 0;
		}
	}
	virtual void appendRow(double data[], int n);
	virtual bool readRow(double data[], int ndata);
	void flush();
	void close();
	void doSummary();
};
TFile *RootNTuple::file=0;
G4String RootNTuple::filename;
RootNTuple *RootNTuple::first=0;

class RootNTupleHandler : public BLNTupleHandler {
	static void init() {
	    static bool firstCall=true;
	    if(firstCall) {
		/* magic line from Rene - for future reference! */
		gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo",
			"*",
			"TStreamerInfo",
			"RIO",
			"TStreamerInfo()");
		firstCall = false;
	    }
	}
public:
	RootNTupleHandler() {
		// cannot init() here -- gROOT not setup
		BLNTuple::registerHandler("root",this); 
		BLNTuple::registerHandler("Root",this); 
	}
	BLNTuple *create(G4String type, G4String category,
			G4String name, G4String fields, G4String filename) {
		init();
		if(Param.getString("histoFile") == "") return 0;
		return new RootNTuple(category,name,fields);
	}
	BLNTuple *read(G4String type, G4String category, G4String name,
					G4String fields, G4String filename);
};
RootNTupleHandler defaultRootNTupleHandler;


RootNTuple::RootNTuple(G4String category, G4String _name, G4String fields)
			: BLNTuple(_name,fields)
{
	name = _name;
	nEntries = 0;
	title = "";
	ntuple = 0;
	readFile = 0;
	readIndex = 0;
	readData = 0;

	if(file == 0) {
		filename = Param.getString("histoFile");
		if(filename == "") return;
		if(filename.find(".root") == filename.npos)
			filename.append(".root");
		file = new TFile(filename.c_str(),"recreate");
		if(file->IsZombie()) {
			G4Exception("RootNTuple","Cannot create file",
				FatalException,filename);
		}
	}
	if(file->IsZombie()) {
		return;
	}

	// create the directory(s) for category, and cd() to it
	TDirectory *dir = file;
	if(!dir->cd()) {
		G4Exception("RootNTuple","Cannot cd into file",
				FatalException,filename);
	}
	for(G4String::size_type place=0,next=0; next!=category.npos; place=next+1) {
		next = category.find("/",place);
		G4String tmp = category.substr(place,next-place);
		if(tmp.size() == 0) continue;
		if(dir->GetDirectory(tmp) && dir->cd(tmp)) {
			continue; // exists
		}
		// created in current directory:
		TDirectory *d = dir->mkdir(tmp,tmp);
		if(!d || !d->cd()) {
			G4Exception("RootNTuple","Cannot create directory",
				FatalException,tmp);
		}
		dir = d;
	}
	
	title = name;
	if(category.size() != 0) {
		title = category;
		title += "/";
		title += name;
	}
	// created in current directory of the Root file:
	ntuple = new TNtuple(name,title,fields);
	if(!ntuple) return;

	nData = ntuple->GetNvar();

	if(first == 0) first = this;
}

void RootNTuple::appendRow(double data[], int n)
{
	if(ntuple == 0) return; // Zombie output file, be silent here but warn 
				// in flush(), close(), and doSummary().
	BLAssert(n==nData && n>0);

	// convert data[] to float-s, and Fill the TNtuple.
	static float *floatData=0;
	static int nFloatData=0;
	if(nFloatData < n) {
		if(floatData) delete floatData;
		floatData = new float[n];
		nFloatData = n;
	}
	for(int i=0; i<n; ++i)
		floatData[i] = data[i];
	ntuple->Fill(floatData);

	++nEntries;
	doCallbacks(data,n);
}

BLNTuple *RootNTupleHandler::read(G4String type, G4String category,
			G4String name, G4String fields, G4String filename)
{
	init();

	TFile *readFile = new TFile(filename,"read");
	if(readFile->IsZombie()) {
		G4Exception("RootNTuple","Cannot read file",
				JustWarning,filename);
		return 0;
	}

	if(category != "") name = category+"/"+name;
	TNtuple *ntuple = (TNtuple *)readFile->Get(name);
	if(!ntuple) {
		name = filename + ":" + name;
		G4Exception("RootNTuple","Cannot find NTuple",
				JustWarning,name);
		return 0;
	}

	return new RootNTuple(readFile,ntuple,fields);
}

RootNTuple::RootNTuple(TFile *file, TNtuple *nt, G4String fields)
					: BLNTuple(nt->GetName(),fields)
{
	ntuple = nt;
	name = ntuple->GetName();
	nData = ntuple->GetNvar();
	nEntries = ntuple->GetEntries();
	title = ntuple->GetTitle();
	readFile = file;
	readIndex = 0;
	readData = new float[nData];
	for(int i=0; i<nData; ++i) readData[i] = 0.0;

	// set Branch addresses to our readData[]
	std::vector<G4String> list = BLCommand::splitString(fields,":");
	int ifield=0;
	for(int i=0; i<list.size(); ++i) {
		if(ifield >= nData) break;
	 	if(list[i].size() == 0) continue;
		if(ntuple->SetBranchAddress(list[i],&readData[ifield]) >= 0) {
			++ifield;
		}
	}
	nData = ifield;
}

bool RootNTuple::readRow(double data[], int ndata)
{
	BLAssert(ntuple!=0 && ndata==nData);
	if(readIndex >= nEntries) return false;
	ntuple->GetEntry(readIndex++);
	for(int i=0; i<ndata; ++i)
		data[i] = readData[i];
	return true;
}

void RootNTuple::flush()
{
	if(ntuple == 0) {
		return;
	}

	if(file == 0) return;

	ntuple->AutoSave("SaveSelf,Overwrite");
}

void RootNTuple::close()
{
	// close the file only once, not evey NTuple, only if some exist
	if(file == 0 || first != this) return;
	if(ntuple == 0) {
		return;
	}

	file->Write();
	file->Close();
	printf("NTuple wrote Root File '%s'\n",filename.c_str());
	file = 0;
}

void RootNTuple::doSummary()
{
	if(ntuple == 0) 
		G4Exception("Root output","No ntuple",JustWarning, 
						"Unwritable file?");
	printf("NTuple %-15s %8d entries\n",name.c_str(),nEntries);
}

#else
int noRoot = 0;
#endif // G4BL_ROOT
