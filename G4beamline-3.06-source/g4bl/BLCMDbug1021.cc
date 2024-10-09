//	BLCMDbug1021.cc
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

#include "G4VProcess.hh"

#include "BLManager.hh"
#include "BLParam.hh"
#include "BLGlobalField.hh"

/**	class BLCMDbug1021 adds a physics process to all charged particles
 *	to improve the accuracy whenever they turn around in an E field.
 *
 *	This is a workaround for Geant4 bug 1021 -- when that gets fixed
 *	this will no longer be needed.
 **/
class BLCMDbug1021 : public BLCommand, public BLCallback {
	bool alreadyRegistered;
	G4double minStep;
public:
	/// Constructor.
	BLCMDbug1021();

	/// commandName() returns "bug1021".
	virtual G4String commandName() { return "bug1021"; }
	
	/// command() implements the bug1021 command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// callback() from BLCallback.
	void callback(int type);
};
BLCMDbug1021 defaultBug1021;

/** class Bug1021 applies the workaround for Geant4 bug 1021
 **/
class Bug1021 : public G4VProcess {
	G4double minStep;
	bool reflect;
	class ParticleChange : public G4ParticleChange {
	public:
		ParticleChange() : G4ParticleChange() { }
	};
	ParticleChange change;
public:
	Bug1021(G4double _minStep) :
				G4VProcess("Bug1021",fUserDefined), change() {
	    minStep = _minStep;
	    reflect = false;
	    G4ParticleTable::G4PTblDicIterator *myParticleIterator =
			G4ParticleTable::GetParticleTable()->GetIterator();
	    myParticleIterator->reset();
	    while((*myParticleIterator)()) {
		G4ParticleDefinition *pd = myParticleIterator->value();
		G4ProcessManager *pmgr = pd->GetProcessManager();
		if(!pmgr) continue;
		pmgr->AddProcess(this,ordDefault,-1,ordDefault);
	    }
	}

	void StartTracking(G4Track *track) { reflect = false; }

	virtual G4double PostStepGetPhysicalInteractionLength(
			const G4Track& track, G4double   previousStepSize, 
			G4ForceCondition* condition)
		{ *condition = NotForced;
		  G4double q = track.GetDefinition()->GetPDGCharge();
		  if(q == 0.0) return DBL_MAX;
		  const G4ThreeVector &pos = track.GetPosition();
		  G4double point[4], field[6];
		  point[0] = pos[0];
		  point[1] = pos[1];
		  point[2] = pos[2];
		  point[3] = track.GetGlobalTime();
		  BLGlobalField::getObject()->GetFieldValue(point,field);
		  G4ThreeVector E(field[3],field[4],field[5]);
		  G4double Etot = E.mag();
		  if(Etot < 100*volt/meter) return DBL_MAX;
		  G4ThreeVector dir = track.GetMomentumDirection();
		  if((q>0.0?1.0:-1.0)*dir.dot(E.unit()) > -0.99) return DBL_MAX;
		  double distanceToStop = track.GetKineticEnergy()/Etot;
		  double step = distanceToStop / 2.0;
		  if(distanceToStop < minStep) {
		  	step = 0.0;
			reflect = true;
		  }
		  return step;
		}
	virtual G4VParticleChange* PostStepDoIt( const G4Track& track, 
      			const G4Step&  stepData)
		{ change.Initialize(track);
		  if(reflect) {
		  	G4ThreeVector dir=track.GetMomentumDirection();
		  	//@@ need to reflect dir in plane perpendicular to E
			// but dir || E, so this is good enough
			dir[0] = -dir[0];
			dir[1] = -dir[1];
			dir[2] = -dir[2];
		  	//@@ approximation for elapsed time
			G4double time = track.GetGlobalTime();
			time += minStep/track.GetVelocity();
			change.ProposeGlobalTime(time);
		  	change.ProposeMomentumDirection(dir);
			reflect = false;
		  }
		  return &change;
		}
	virtual G4double AlongStepGetPhysicalInteractionLength(
			const G4Track& track, G4double  previousStepSize,
			G4double  currentMinimumStep, G4double& proposedSafety,
			G4GPILSelection* selection)
		{ return -1.0; }
	virtual G4VParticleChange* AlongStepDoIt( const G4Track& track,
			const G4Step& stepData)
		{ change.Initialize(track);
		  return &change;
		}
	virtual G4double AtRestGetPhysicalInteractionLength( 
			const G4Track& track, G4ForceCondition* condition)
		{ return -1.0; }
	virtual G4VParticleChange* AtRestDoIt( const G4Track& track,
			const G4Step& stepData)
		{ change.Initialize(track);
		  return &change;
		}
};
BLCMDbug1021::BLCMDbug1021()
{
	registerCommand(BLCMDTYPE_PHYSICS);
	setSynopsis("Workaround to improve accuracy of bug1021 in E field");
	setDescription("When a charged particle turns around in an E field, "
		"a bug in the Geant4 transportation process can sometimes "
		"give it a wildly-incorrect kinetic energy. "
		"This workaround computes the distance to turn-around, "
		"and limits the step to half that vlue until minStep is "
		"reached; at that point the track is reflected.\n\n"
		"Simulations in which there are no E fields, or no charged "
		"particle ever gets below ~0.001 MeV in an E field, have "
		"no need to apply this workaround.\n\n"
		"This command is not placed into the geometry.");

	alreadyRegistered = false;
	minStep = 0.002 * mm;
}

int BLCMDbug1021::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	int retval = handleNamedArgs(namedArgs);

	if(!alreadyRegistered) {
		BLManager::getObject()->registerCallback(this,0);
		alreadyRegistered = true;
	}

	print("");

	return retval;
}

void BLCMDbug1021::defineNamedArgs()
{
	argDouble(minStep,"minStep","Minimum step in space (mm, default=0.002)",mm);
}

void BLCMDbug1021::callback(int type)
{
	if(type == 0) (void)new Bug1021(minStep);
}
