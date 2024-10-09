//	BLCMDtotalenergy.cc
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
#include <map>
#include <vector>
#include "fnmatch.h"

#include "BLManager.hh"
#include "BLCommand.hh"
#include "BLWriteAsciiFile.hh"
#include "BLMPI.hh"

/**	class BLCMDtotalenergy prints the total energy deposited into each
 *	selected volume, at the end of the run.
 **/
class BLCMDtotalenergy : public BLCommand, public BLCallback, 
			public BLManager::SteppingAction,
			public BLManager::TrackingAction {
	G4String volumes;
	G4int ancestors;
	G4String filename;
	G4int verbose;
	bool alreadyRegistered;
	struct Totalizer {
		double total;
		Totalizer() { total=0.0; }
	};
	BLManager *manager;
	std::map<G4VPhysicalVolume*,Totalizer> total;
	std::vector<G4VPhysicalVolume*> printOrder;
	std::vector<G4String> patterns;
	G4VPhysicalVolume *world;
public:
	/// Constructor.
	BLCMDtotalenergy();

	/// commandName() returns "totalenergy".
	virtual G4String commandName() { return "totalenergy"; }
	
	/// command() implements the totalenergy command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// callback() from BLCallback.
	void callback(int type);

	/// walk the tree of physical volumes
	void walkPVTree(G4VPhysicalVolume *pv, G4String indent="");

	/// isMatch() returns true if the name of this PV matches some pattern
	bool isMatch(G4VPhysicalVolume *pv);

	/// sum() will sum the energy into the relevant volumes.
	void sum(const G4TouchableHandle &th, double energy);

	/// UserSteppingAction() from BLManager::SteppingAction.
	void UserSteppingAction(const G4Step *step);

	/// PreUserTrackingAction() from BLManager::TrackingAction.
	void PreUserTrackingAction(const G4Track *track);

	/// PostUserTrackingAction() from BLManager::TrackingAction.
	void PostUserTrackingAction(const G4Track *track);
};
BLCMDtotalenergy defaultTotalenergy;


BLCMDtotalenergy::BLCMDtotalenergy() : volumes("*"), total(), printOrder(), 
							patterns()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Print total energy deposited in selected volumes.");
	setDescription("At end of run, prints the total energy deposited in "
		"the selected volumes.\n\n"
		"Volume-name patterns are like UNIX filename patterns: "
		"e.g. '*[AB]*' matches any name containing an A or a B.\n\n"
		"Tracks that are killed have their kinetic energy summed into "
		"the volume where they were killed, unless they are killed "
		"because they leave the World.\n\n"
		"With ancestors=1, energy deposited in matching volumes is "
		"added into their ancestors; energy deposited directly into "
		"those ancestors is not summed into them unless their names "
		"also match. "
		"That is, if A is a child of B, but only A matches the list "
		"of volume-names, energy "
		"deposited into A will be reported in both A and B, but energy "
		"deposited directly into B is ignored.\n\n"
		"This command is not placed into the geometry.");

	ancestors = 0;
	filename = "";
	verbose = 0;
	alreadyRegistered = false;
	world = 0; // initialized in callback()
}

int BLCMDtotalenergy::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	manager = BLManager::getObject();
	int retval = handleNamedArgs(namedArgs);

	if(alreadyRegistered) {
		printError("Duplicate totalenergy command");
		return -1;
	}

	manager->registerCallback(this,0);
	manager->registerCallback(this,2);
	manager->registerSteppingAction(this);
	manager->registerTrackingAction(this);
	alreadyRegistered = true;

	print("");

	patterns = splitString(volumes, ",", true);

	return retval;
}

void BLCMDtotalenergy::defineNamedArgs()
{
	argString(volumes,"volumes","Comma-separated list of Volume-Name patterns (*)");
	argInt(ancestors,"ancestors","Set nonzero to sum energy into all ancestor (enclosing) volumess (0).");
	argInt(ancestors,"enclosing","Synonym for ancestors.");
	argString(filename,"filename","Filename to write summary (stdout).");
	argString(filename,"file","Synonym for filename.");
	argInt(verbose,"verbose","Verbosity level (0).");
}

void BLCMDtotalenergy::callback(int type)
{
	if(type == 0) { 
		// Walk the tree of physical volumes
		printf("\nPhysical volumes for totalenergy (* = totaled):\n");
		world = manager->getWorldPhysicalVolume();
		walkPVTree(world);
	} else if(type == 2) {
		if(manager->getSteppingVerbose() > 0 && verbose == 0)
			verbose = manager->getSteppingVerbose();
		FILE *out=stdout;
		bool printing = (!BLMPI::isMPI() || BLMPI::isRank0());
		if(printing && filename.size() > 0) {
			out = BLWriteAsciiFile::fopen(filename);
			if(!out) {
				G4Exception("totalenergy",
					"Cannot write to file - stdout used",
					JustWarning,filename);
				out = stdout;
			}
		}
		if(printing)
			fprintf(out,"# Total energy deposited in selected volumes (MeV) (Joules):\n");
		if(BLMPI::isMPI()) {
			unsigned nd = printOrder.size();
			double *data = new double[nd];
			for(unsigned i=0; i<nd; ++i) {
				G4VPhysicalVolume *pv = printOrder[i];
				data[i] = total[pv].total;
			}
			BLMPI::sumAllWorkers(data,nd);
			for(unsigned i=0; i<nd; ++i) {
				G4VPhysicalVolume *pv = printOrder[i];
				total[pv].total = data[i];
			}
		}
		if(!printing) return;
		for(unsigned i=0; i<printOrder.size(); ++i) {
			G4VPhysicalVolume *pv = printOrder[i];
			double t=total[pv].total/MeV;
			fprintf(out,"%s\t%.4g\t%.4g\n",pv->GetName().c_str(),
				total[pv].total/MeV,total[pv].total/joule);
		}
		if(out != stdout) {
			BLWriteAsciiFile::fclose(out);
			printf("totalenergy: wrote '%s'\n",filename.c_str());
		}
	}
}

void BLCMDtotalenergy::walkPVTree(G4VPhysicalVolume *pv, G4String indent)
{
	printf("%s%s",indent.c_str(),pv->GetName().c_str());
	if(isMatch(pv)) {
		printf(" *");
		if(total.count(pv) != 0)
			G4Exception("totalenergy command",
				"Duplicate PhysicalVolume cut",JustWarning, "");
		total[pv] = Totalizer();
		printOrder.push_back(pv);
	}
	printf("\n");

	G4LogicalVolume *lv = pv->GetLogicalVolume();
	int n=lv->GetNoDaughters();
	for(int i=0; i<n; ++i)
		walkPVTree(lv->GetDaughter(i),indent+"    ");
}

bool BLCMDtotalenergy::isMatch(G4VPhysicalVolume *pv)
{
	const char *name = pv->GetName().c_str();
	for(unsigned i=0; i<patterns.size(); ++i) {
		if(patterns[i].size() == 0) continue;
		if(fnmatch(patterns[i].c_str(),name,0) == 0)
			return true;
	}
	return false;
}

void BLCMDtotalenergy::sum(const G4TouchableHandle &th, double energy)
{
	if(!th) return;

	for(int depth=0; depth<=th->GetHistoryDepth(); ++depth) {
		G4VPhysicalVolume *pv = th->GetVolume(depth);
		if(!pv) break;
		if(depth == 0 && total.count(pv) == 0) break;
		total[pv].total += energy; // create if needed
		if(verbose > 0)
		    printf("totalenergy: adding %.4f MeV to '%s'\n",
		    		energy/MeV,pv->GetName().c_str());
		if(pv == world) break;
		if(!ancestors) break;
	}
}

void BLCMDtotalenergy::UserSteppingAction(const G4Step *step)
{
	if(BLManager::getObject()->getState() != BEAM) return;
	G4double energy = step->GetTotalEnergyDeposit();
	G4Track *track = step->GetTrack();
	sum(track->GetTouchableHandle(),energy*track->GetWeight());
}

void BLCMDtotalenergy::PreUserTrackingAction(const G4Track *track)
{
}

void BLCMDtotalenergy::PostUserTrackingAction(const G4Track *track)
{
	//if(track->GetTrackStatus() == fSuspend) return;
	if(track->GetNextVolume() == 0) return; // don't sum if leaving World
	sum(track->GetTouchableHandle(),
				track->GetKineticEnergy()*track->GetWeight());
}

