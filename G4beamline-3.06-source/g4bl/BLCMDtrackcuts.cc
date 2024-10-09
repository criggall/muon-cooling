//	BLCMDtrackcuts.cc
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

#include <set>

#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "BLCommand.hh"
#include "BLParam.hh"
#include "BLManager.hh"

/**	class BLCMDtrackcuts implements the trackcuts command
 *
 **/
class BLCMDtrackcuts : public BLCommand, public BLManager::SteppingAction,
					public BLManager::StackingAction {
	G4double kineticEnergyCut;
	G4double kineticEnergyMax;
	G4int steppingVerbose;
	G4int killSecondaries;
	G4String kill;
	G4String keep;
	G4double maxTime;
	G4int keepPrimaries;
	std::set<G4int> killPDG;
	std::set<G4int> keepPDG;
	bool stepping;
	BLManager *manager;
public:
	/// Constructor
	BLCMDtrackcuts();

	BLCMDtrackcuts(const BLCMDtrackcuts &r);

	G4String commandName() { return "trackcuts"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	/// UserSteppingAction() from BLManager::SteppingAction.
	void UserSteppingAction(const G4Step *step);

	/// ClassifyNewTrack() from BLManager::StackingAction.
	G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *track);
};
BLCMDtrackcuts defaultTrackCuts;


BLCMDtrackcuts::BLCMDtrackcuts() : BLCommand(), BLManager::SteppingAction(), 
BLManager::StackingAction(), kill(), keep()
{
	registerCommand(BLCMDTYPE_CUTS);
	setSynopsis("Specifies per-track cuts.");
	setDescription("Applied to each track before tracking, and at each "
		"step.\n\n"
		"NOTE: the physics command MUST come before trackcuts.\n\n"
		"This command is not placed into the geometry.");

	kineticEnergyCut = 0.0;
	kineticEnergyMax = DBL_MAX;
	steppingVerbose = 0;
	killSecondaries = 0;
	maxTime = 1.0*millisecond;
	keepPrimaries = 0;
	stepping = false;
	manager = BLManager::getObject();
}

BLCMDtrackcuts::BLCMDtrackcuts(const BLCMDtrackcuts& r) : BLCommand(r), 
BLManager::SteppingAction(), BLManager::StackingAction(), kill(), keep()
{
	kineticEnergyCut = r.kineticEnergyCut;
	kineticEnergyMax = r.kineticEnergyMax;
	steppingVerbose = r.steppingVerbose;
	killSecondaries = r.killSecondaries;
	maxTime = r.maxTime;
	keepPrimaries = r.keepPrimaries;
	stepping = r.stepping;
	manager = BLManager::getObject();
}

int BLCMDtrackcuts::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(BLManager::getObject()->getPhysics() == 0) {
		printError("trackcuts: no physics; the physics command MUST come first.");
		return -1;
	}

	G4ParticleTable *pt = G4ParticleTable::GetParticleTable();

	if(namedArgs.size() == 0) {
		int n = pt->entries();
		for(int i=0; i<n; ++i) {
			G4ParticleDefinition *d = pt->GetParticle(i);
			G4String name = d->GetParticleName();
			G4String type = d->GetParticleType();
			G4String subtype = d->GetParticleSubType();
			G4int pdg = d->GetPDGEncoding();
			printf("%15.15s  %s/%s  PDGEncoding=%d\n",
				name.c_str(),type.c_str(),subtype.c_str(),pdg);
		}
		return 0;
	}

	BLCMDtrackcuts *t = new BLCMDtrackcuts(defaultTrackCuts);

	int retval = t->handleNamedArgs(namedArgs);

	// get PDGid-s for particles to kill into killPDG
	std::vector<G4String> v=splitString(t->kill,",",true);
	for(unsigned i=0; i<v.size(); ++i) {
		if(v[i].size() == 0) continue;
		G4ParticleDefinition *d = pt->FindParticle(v[i]);
		if(!d) {
			printError("Cannot find particle '%s' (physics must come first)",v[i].c_str());
			continue;
		}
		t->killPDG.insert(d->GetPDGEncoding());
	}

	// get PDGid-s for particles to keep into keepPDG
	v=splitString(t->keep,",",true);
	for(unsigned i=0; i<v.size(); ++i) {
		if(v[i].size() == 0) continue;
		G4ParticleDefinition *d = pt->FindParticle(v[i]);
		if(!d) {
			printError("Cannot find particle '%s' (physics must come first)",v[i].c_str());
			continue;
		}
		t->keepPDG.insert(d->GetPDGEncoding());
	}

	BLManager::getObject()->registerSteppingAction(t);
	BLManager::getObject()->registerStackingAction(t);

	BLManager::getObject()->getPhysics()->setMaxTime(t->maxTime);

	t->print("");

	return retval;
}

void BLCMDtrackcuts::defineNamedArgs()
{
	argString(kill,"kill","List of particles to kill (comma separated).");
	argString(keep,"keep","List of particles to keep (kill all others).");
	argInt(killSecondaries,"killSecondaries","Set nonzero to kill all secondaries.");
	argDouble(kineticEnergyCut,"kineticEnergyCut","Minimum K.E. to track (0 MeV).",MeV);
	argDouble(kineticEnergyMax,"kineticEnergyMax","Maximum K.E. to track (infinite MeV).",MeV);
	argDouble(maxTime,"maxTime","Maximum lab time to track (1000000 ns).");
	argInt(keepPrimaries,"keepPrimaries","Set nonzero to keep tracks with ParentID==0 regardless of other tests.");
	argInt(steppingVerbose,"steppingVerbose","Set nonzero to print kills\n"
		"(defaults to parameter value).");
}


void BLCMDtrackcuts::UserSteppingAction(const G4Step *step)
{
	stepping = true;

	G4Track *track = step->GetTrack();
	if(ClassifyNewTrack(track) == fKill)
		((G4Track*)track)->SetTrackStatus(fStopAndKill);

	stepping = false;
}

G4ClassificationOfNewTrack BLCMDtrackcuts::ClassifyNewTrack(
							const G4Track *track)
{
	// handle keepPrimaries independent of other tests
	if(keepPrimaries != 0 && track->GetParentID() == 0)
		return fUrgent;

	if(manager->getSteppingVerbose() > 0) steppingVerbose = 1;

	// Kill particle if Kinetic Energy is below the threshold
	G4double kineticEnergy = track->GetKineticEnergy();
	if(kineticEnergy < kineticEnergyCut)  {
		if(steppingVerbose != 0)
			printf("Track killed: K.E. %.6f < %.6f MeV\n",
				kineticEnergy/MeV, kineticEnergyCut/MeV);
		return fKill;
	}

	// Kill particle if Kinetic Energy is above the Maximum
	if(kineticEnergy > kineticEnergyMax)  {
		if(steppingVerbose != 0)
			printf("Track killed: K.E. %.3f > %.3f MeV\n",
				kineticEnergy/MeV, kineticEnergyMax/MeV);
		return fKill;
	}

	// Kill particle if time is too large
	// Now done in BLCMDphysics via class TrackTimeLimit.

	// shortcut while stepping, avoid unnecessary tests
	if(stepping) return fUrgent;

	// kill secondaries, unless in tune particle mode (the tune particle
	// is frequently re-tracked as a "secondary")
	if(BLManager::getObject()->getState() != TUNE &&
			killSecondaries != 0 && track->GetParentID() != 0) {
		if(steppingVerbose != 0) printf("Secondary Track killed\n");
		return fKill;
	}

	// kill particles in kill
	G4ParticleDefinition *pdg = track->GetDefinition();
	if(killPDG.count(pdg->GetPDGEncoding()) > 0) {
		if(steppingVerbose != 0) 
			printf("Track killed because it's a %s\n",
				pdg->GetParticleName().c_str());
		return fKill;
	}

	// kill particles not in keep
	if(keepPDG.size() > 0 && keepPDG.count(pdg->GetPDGEncoding()) == 0) {
		if(steppingVerbose != 0) 
			printf("Track killed because it's a %s\n",
				pdg->GetParticleName().c_str());
		return fKill;
	}

	return fUrgent;
}
