//	BLCMDsourceonly.cc
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

#include <stdio.h>
#include <vector>
#include <map>
#include <cmath>	// std::isnan()

#include "BLCommand.hh"
#include "BLCallback.hh"
#include "BLManager.hh"
#include "BLTrackFile.hh"
#include "BLNTuple.hh"
#include "BLMPI.hh"

/**	class BLCMDsourceonly - Convert the simulation to source only --
 *	generate events and write them to a file withut simulating them.
 *
 *	format=BLTrackFile or format=G4beamline
 *		Output file is a BLTrakckFile, version 1.
 *
 *	format=MCNP
 *! File Format: x y z u v w erg tme wgt ipt
 *!         (x,y,z in cm; erg in MeV; tme in shakes; ipt is particle id)
 *!         (u,v,w) are direction cosines (px/ptot,py/ptot,pz/ptot)
 *!         wgt is the particle's weight
 *!         lines beginning with # are comments and are ignored.
 *
 **/
class BLCMDsourceonly : public BLCommand, public BLCallback, 
					public BLManager::StackingAction {
	G4String filename;
	G4String format;
	BLTrackFile *trackFile;
	FILE *out;
public:
	BLCMDsourceonly();

	G4String commandName() { return "sourceonly"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	virtual void callback(int type);

	// from BLManager::StackingAction -- output track to out and kill it.
	virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*);

	int iptFromPDGid(int PDGid);
};

BLCMDsourceonly defaultSourceonlyCommand;

BLCMDsourceonly::BLCMDsourceonly() : BLCommand(), BLCallback()
{
	registerCommand(BLCMDTYPE_OTHER);
	setSynopsis("Run a source-only simulation.");
	setDescription("More to come...");

	// initialize class variables here
	filename = "source.txt";
	format = "BLTrackFile";
	trackFile = 0;
	out = 0;
}

int BLCMDsourceonly::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(BLMPI::isMPI())
	    G4Exception("sourceonly","Incompatible with MPI",FatalException,"");

	handleNamedArgs(namedArgs);

	print("");

	// replace the main run loop, register our StackingAction
	BLManager::getObject()->registerCallback(this,3);
	BLManager::getObject()->registerStackingAction(this);

	if(format == "MCNP") {
		out = fopen(filename.c_str(),"w");
		if(!out) {
			printError("sourceonly: cannot write '%s'\n",
							filename.c_str());
			return -1;
		}
	} else if(format == "BLTrackFile" || format == "G4beamline") {
		trackFile = new BLTrackFile(filename,"source","w",1);
	} else {
		printError("sourceonly: invalid format '%s'\n",format.c_str());
		return -1;
	}

	return 0;
}

void BLCMDsourceonly::defineNamedArgs()
{
	argString(filename,"filename","Filename of output file.");
	argString(format,"format","Format of output file (BLTrackFile).");
	argString(filename,"file","Synonym for filename.");
}	

void BLCMDsourceonly::callback(int type)
{
	G4Exception("sourceonly","Source-Only run",JustWarning,"");
	BLManager::getObject()->trackBeam();

	printf("\nEnd of sourceonly run.\n");
	if(trackFile) delete trackFile;
	trackFile = 0;
	if(out) fclose(out);
	out = 0;

	BLNTuple::summary();
	BLNTuple::closeAll();

	// handle post-tracking callbacks
	//BLManager::getObject()->handleCallbacks(2);
}

G4ClassificationOfNewTrack BLCMDsourceonly::ClassifyNewTrack(
							const G4Track *track)
{
	G4ThreeVector pos = track->GetPosition();
	G4ThreeVector mom = track->GetMomentum();
	G4ThreeVector dir = track->GetMomentumDirection();
	double t = track->GetGlobalTime();
	double ke = track->GetKineticEnergy();
	double w = track->GetWeight();
	int pdgid = track->GetDefinition()->GetPDGEncoding();
	int eventID = BLManager::getObject()->getEventID();
	int trackID = BLManager::getObject()->getExternalTrackID(track);
	int parentID = track->GetParentID();

	if(trackFile) {
		if(std::isnan(mom[0]) || std::isnan(mom[1]) || 
							std::isnan(mom[2]))
			mom = G4ThreeVector(0.0,0.0,0.0);
		trackFile->write(pos,t,mom,pdgid,eventID,trackID,parentID,w);
	} else if(out) {
		// File Format: x y z u v w erg tme wgt ipt
		int ipt = iptFromPDGid(pdgid);
		fprintf(out,"%.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %.6g %d\n",
			pos[0]/cm,pos[1]/cm,pos[2]/cm,dir[0],dir[1],dir[2],
			ke/MeV,t/(ns*10),w,ipt);
	} else {
		G4Exception("sourceonly","No Output",FatalException,"");
	}

	return fKill;
}

int BLCMDsourceonly::iptFromPDGid(int PDGid)
{
	static std::map<int,int> hash;
	if(hash.size() == 0) {
		hash[2112] = 1; hash[-2112] = 5;	// neutron
		hash[22] = 2;				// photon / gamma
		hash[11] = 3; hash[-11] = 8;		// electron
		hash[13] = 4; hash[-13] = 16;		// muon
		hash[12] = 6; hash[-12] = 17;		// nu_e
		hash[14] = 7; hash[-14] = 18;		// nu_mu
		hash[2212] = 9; hash[-2212] = 19;	// proton
		hash[3122] = 10; hash[-3122] = 25;	// lambda
		hash[3222] = 11; hash[-3222] = 26;	// sigma+
		hash[3112] = 12; hash[-3112] = 27;	// sigma-
		hash[3322] = 13; hash[-3322] = 28;	// xi0
		hash[3312] = 14; hash[-3312] = 29;	// xi-
		hash[3334] = 15; hash[-3334] = 30;	// omega-
		hash[211] = 20; hash[-211] = 35;	// pi+
		hash[111] = 21;				// pi0
		hash[321] = 22; hash[-321] = 36;	// kaon+
		hash[310] = 23; hash[130] = 24;		// kaon0S, kaon0L
		hash[1000010020] = 31;			// deuteron
		hash[1000010030] = 32;			// triton
		hash[1000020030] = 33;			// helion (He3)
		hash[1000020040] = 34;			// alpha
	}

	int v = hash[PDGid];

	if(v == 0 && PDGid > 1000000000) {
		int Z = (PDGid/10000) % 1000;
		int A = (PDGid/10) % 1000;
		v = Z*1000 + A;
		if(PDGid%10 != 0)
			G4Exception("sourceonly",
				"Excited ion, excitation omitted",
				JustWarning,"");
	}

	return v;
}

