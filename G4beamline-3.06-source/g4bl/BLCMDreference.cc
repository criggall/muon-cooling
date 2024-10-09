//	BLCMDreference.cc
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

#include "G4RunManager.hh"
#include "G4Track.hh"

#include "BLManager.hh"
#include "BLBeam.hh"
#include "BLParam.hh"
#include "BLGroup.hh"
#include "BLCoordinates.hh"
#include "BLGlobalField.hh"

const G4double UNDEFINED = -3.7e21;
enum State { INIT, TUNING, DONE };
const int MAX_COUNT=10;	// maximum steps while tuning the momentum

/**	class BLCMDreference implements the reference command to track
 *	the tune and reference particles.
 **/
class BLCMDreference : public BLBeam, public BLCommand, 
			public BLManager::TrackingAction,
			public BLManager::ZSteppingAction,
			public BLManager::RunAction {
	G4String particle;
	G4double referenceMomentum;
	G4double beamX;
	G4double beamY;
	G4double beamZ;
	G4double beamT;
	G4double beamXp;
	G4double beamYp;
	G4String rotation;
	G4double tuneZ;
	G4double tuneMomentum;
	G4double tolerance;
	G4int noEfield;
	G4int noBfield;
	G4int noEloss;
	G4int trackID;
	G4RotationMatrix *rotationMatrix;
	G4ThreeVector position;
	G4ParticleGun *particleGun;
	G4ParticleDefinition *particleDefinition;
	G4Track saveTrack;
	State state;
	G4int count;
	static BLCMDreference *lastGen;
	static int referenceID;
public:
	/// Constructor.
	BLCMDreference();

	/// copy constructor
	BLCMDreference(BLCMDreference &r);

	/// commandName() returns "reference".
	virtual G4String commandName() { return "reference"; }
	
	/// command() implements the reference command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// getNEvents() returns the # events to process.
	virtual int getNEvents() const { return nEvents; }

	/// init() will initialize internal variables.
	virtual void init();

	/// generateReferenceParticle() generates the reference particle.
	virtual bool generateReferenceParticle(G4Event *event);

	/// nextBeamEvent() generates the next beam event.
	virtual bool nextBeamEvent(G4Event *event);

	/// summary() will print a summary, if necessary.
	virtual void summary() { }

	// from TrackingAction, ZSteppingAction, and RunAction
	void PreUserTrackingAction(const G4Track *track);
	void PostUserTrackingAction(const G4Track *track);
	void UserZSteppingAction(const G4Track *track);
	void BeginOfRunAction(const G4Run* run) { }
	void EndOfRunAction(const G4Run* run);
};
BLCMDreference *BLCMDreference::lastGen = 0;
int BLCMDreference::referenceID = 0;

BLCMDreference defineReference;

BLCMDreference::BLCMDreference() : BLBeam(), BLCommand(),
			BLManager::TrackingAction(),
			BLManager::ZSteppingAction(),
			BLManager::RunAction(),
			particle(), rotation(), position(), saveTrack()
{
	registerCommand(BLCMDTYPE_BEAM);
	setSynopsis("Define a reference particle.");
	setDescription("The reference particle is nominally headed in the +Z direction.\n"
		"Multiple reference particles can be defined, at different\n"
		"positions, momenta, particle types, etc.\n"
		"All coordinates are centerline coordinates.\n\n"
		"If desired, the referenceMomentum will be tuned to a specific "
		"value at a later z position in the beamline by giving values "
		"for tuneZ, and tuneMomentum; tolerance can be set if desired. "
		"\n\nNormally used in conjunction with a 'beam' command.\n\n"
		"This command is not placed into the geometry.");
	// initial default values:
	particle = "mu+";
	referenceMomentum = UNDEFINED;
	beamX = 0.0;
	beamY = 0.0;
	beamZ = 0.0;
	beamT = 0.0;
	beamXp = 0.0;
	beamYp = 0.0;
	tuneZ = UNDEFINED;
	tuneMomentum = UNDEFINED;
	tolerance = 0.001*MeV;
	noEfield = 0;
	noBfield = 0;
	noEloss = 0;
	trackID = 0;
	rotationMatrix = 0;
	particleDefinition = 0;
	particleGun = 0;
	state = INIT;
	count = 0;
}

BLCMDreference::BLCMDreference(BLCMDreference &r) : BLBeam(), BLCommand(),
			BLManager::TrackingAction(r),
			BLManager::ZSteppingAction(r),
			BLManager::RunAction(r),
			particle(), rotation(), position(), saveTrack()
{
	particle = r.particle;
	referenceMomentum = r.referenceMomentum;
	beamX = r.beamX;
	beamY = r.beamY;
	beamZ = r.beamZ;
	beamT = r.beamT;
	beamXp = r.beamXp;
	beamYp = r.beamYp;
	tuneZ = r.tuneZ;
	tuneMomentum = r.tuneMomentum;
	tolerance = r.tolerance;
	noEfield = r.noEfield;
	noBfield = r.noBfield;
	noEloss = r.noEloss;
	trackID = r.trackID;
	rotationMatrix = 0;
	particleDefinition = 0;
	particleGun = 0;
	state = INIT;
	count = 0;
}

int BLCMDreference::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	BLCMDreference *b = new BLCMDreference(*this);

	b->rotation = "";
	int retval = b->handleNamedArgs(namedArgs);

	b->trackID = ++referenceID;

	if(b->referenceMomentum == UNDEFINED)
		printError("reference: error - need referenceMomentum");

	if(b->rotation != "") {
		b->rotationMatrix = stringToRotationMatrix(b->rotation);
	} else {
		b->rotationMatrix = new G4RotationMatrix();
	}
	// as usual, the order seems backward because (C R C^-1) C = C R.
	*b->rotationMatrix = *BLCoordinates::getCurrentRotation() * 
							*b->rotationMatrix;
	G4ThreeVector local(b->beamX,b->beamY,b->beamZ);
	BLCoordinates::getCurrentGlobal(local,b->position);

	// ensure beam is within the world
	BLGroup::getWorld()->setMinWidth(fabs(b->position[0])*2.0);
	BLGroup::getWorld()->setMinHeight(fabs(b->position[1])*2.0);
	BLGroup::getWorld()->setMinLength(fabs(b->position[2])*2.0);

	BLManager::getObject()->registerReference(b);

	BLManager::getObject()->registerTrackingAction(b);

	// setup for tuning
	if(b->tuneZ != UNDEFINED && b->tuneMomentum != UNDEFINED) {
		BLManager::getObject()->registerZStep(b->tuneZ,b,1);
		BLManager::getObject()->registerRunAction(b,false);
	}

	b->print("");

	return retval;
}

void BLCMDreference::defineNamedArgs()
{
	argString(particle,"particle","Reference particle name");
	argDouble(beamX,"beamX","Reference location in X (mm)");
	argDouble(beamY,"beamY","Reference location in Y (mm)");
	argDouble(beamZ,"beamZ","Reference location in Z (mm)");
	argDouble(beamT,"beamT","Reference time (ns)");
	argString(rotation,"rotation","Rotation of the beam");
	argDouble(referenceMomentum,"referenceMomentum",
				"Reference particle momentum (MeV/c)",MeV);
	argDouble(beamXp,"beamXp","Reference particle Xp (radians)");
	argDouble(beamYp,"beamYp","Reference particle Yp (radians)");
	argDouble(referenceMomentum,"meanMomentum", 
				"Synonymn for referenceMomentum",MeV);
	argDouble(beamXp,"meanXp","Synonym for beamXp.");
	argDouble(beamYp,"meanYp","Synonym for beamYp.");
	argDouble(tuneZ,"tuneZ","Z position for momentum tuning.");
	argDouble(tuneMomentum,"tuneMomentum","Desired momentum for momentum tuning.");
	argDouble(tolerance,"tolerance","tolerance for momentum tuning (0.001 MeV/c).",MeV);
	argInt(noEfield,"noEfield","Set nonzero to make this Tune and Reference"
	    " particle not respond to E fields (ICOOL style)");
	argInt(noBfield,"noBfield","Set nonzero to make this Tune and Reference"
	    " particle not respond to B fields (ICOOL style)");
	argInt(noEloss,"noEloss","Set nonzero to make this Tune and Reference"
	    " particle not respond to ionization energy loss (ICOOL style)");
	argDouble(referenceMomentum,"P","Synonym for referenceMomentum",MeV);
}

void BLCMDreference::init()
{
	if(particleDefinition != 0) return;

	particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle(particle);
	if(!particleDefinition)
		G4Exception("reference command","UnknownParticle",
				FatalException, "Unknown particle type");
	particleGun = new G4ParticleGun(1);
	particleGun->SetParticleDefinition(particleDefinition);
}

bool BLCMDreference::generateReferenceParticle(G4Event *event)
{
	G4double mass = particleDefinition->GetPDGMass();
	G4double ke = sqrt(referenceMomentum*referenceMomentum + mass*mass)
								- mass;
	G4ThreeVector dir;
	dir[0] = beamXp;
	dir[1] = beamYp;
	dir[2] = 1.0/sqrt(1.0 + dir[0]*dir[0] + dir[1]*dir[1]);
	dir[0] *= dir[2];
	dir[1] *= dir[2];
	if(rotationMatrix)
		dir = *rotationMatrix * dir;

	particleGun->SetParticleTime(beamT);
	particleGun->SetParticleEnergy(ke);
	particleGun->SetParticleMomentumDirection(dir);
	particleGun->SetParticlePosition(position);
	particleGun->GeneratePrimaryVertex(event);
	setRandomSeedToTrack(-1);
	lastGen = this;

	BLManager::getObject()->setPrimaryTrackID(trackID,0);
	BLManager::getObject()->clearTrackIDMap();
	BLManager::getObject()->setNextSecondaryTrackID(1000);

	return true;
}

bool BLCMDreference::nextBeamEvent(G4Event *event)
{
	return false;
}

void BLCMDreference::PreUserTrackingAction(const G4Track *track)
{
	if(state == INIT && lastGen == this && 
			tuneZ != UNDEFINED && tuneMomentum != UNDEFINED) {
		saveTrack.CopyTrackInfo(*track);
		saveTrack.SetUserInformation(0);
		state = TUNING;
		count = 0;
		printf("reference: tune begun  referenceMomentum=%.4f\n",
					referenceMomentum);
	}

	BLManagerState mgrState = BLManager::getObject()->getState();
	if((mgrState == TUNE || mgrState == REFERENCE) && lastGen == this) {
	    if(noEfield != 0) {
		BLGlobalField::getObject()->zeroEfield(true);
		printf("***  Tune and Reference Particles ignore E fields.\n");
	    }
	    if(noBfield != 0) {
		BLGlobalField::getObject()->zeroBfield(true);
		printf("***  Tune and Reference Particles ignore B fields.\n");
	    }
	    if(noEloss != 0) {
		BLManager::getObject()->getPhysics()->
						setProcessEnable("Ioni",false);
		printf("***  Tune and Reference Particles ignore Energy loss.\n");
	    }
	}
}

void BLCMDreference::PostUserTrackingAction(const G4Track *track)
{
	BLManagerState mgrState = BLManager::getObject()->getState();
	if((mgrState == TUNE || mgrState == REFERENCE) && lastGen == this) {
		BLGlobalField::getObject()->zeroBfield(false);
		BLGlobalField::getObject()->zeroEfield(false);
		BLManager::getObject()->getPhysics()->
						setProcessEnable("Ioni",true);
	}
}

void BLCMDreference::UserZSteppingAction(const G4Track *track)
{
	if(state != TUNING) return;
	if(++count > MAX_COUNT)
		G4Exception("reference","Tuning Iteration Limit",
							FatalException,"");

	G4double ke = track->GetKineticEnergy();
	G4double mass = particleDefinition->GetPDGMass();
	G4double ptot = sqrt((ke+mass)*(ke+mass) - mass*mass);
	G4double want = sqrt(tuneMomentum*tuneMomentum+mass*mass) - mass;
	printf("reference: tune step %d  got momentum=%.6f (K.E.=%.6f)\n",
							count,ptot,ke);
	if(fabs(ke-want) < tolerance) {
		state = DONE;
		G4double e = saveTrack.GetKineticEnergy() + mass;
		referenceMomentum = sqrt(e*e-mass*mass);
		printf("reference: tune Complete   new referenceMomentum=%.4f\n",
					referenceMomentum);
		goto quit;
	}
	saveTrack.SetKineticEnergy(saveTrack.GetKineticEnergy()+(want-ke));
	// restore saveTrack and kill the current track
	BLManager::getObject()->getSteppingManager()->GetfSecondary()->
			push_back(new G4Track(saveTrack));
	((G4Track *)track)->SetTrackStatus(fStopAndKill);

quit:	;
}

void BLCMDreference::EndOfRunAction(const G4Run* run)
{
	if(state != DONE)
		G4Exception("reference","Tuning failed to converge",
							FatalException,"");

	lastGen = 0;
}
