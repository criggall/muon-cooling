//	BLCMDisotropicsource.cc
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
#define _USE_MATH_DEFINES
#include <math.h>

#include "G4RunManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

#include "BLAssert.hh"
#include "BLBeam.hh"
#include "BLParam.hh"
#include "BLGroup.hh"
#include "BLCoordinates.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

const G4double UNDEFINED = -3.7e21;
const G4int ALL_EVENTS = 0x7FFFFFFF; // used to indicate an unset argument

/**	class BLCMDisotropicsource implements the isotropicsource command.
 **/
class BLCMDisotropicsource : public BLBeam, public BLCommand {
	G4String particle;
	G4int nEvents;
	G4int firstEvent;
	G4double x;
	G4double y;
	G4double z;
	G4double weight;
	G4int secondaryTrackID;
	G4double meanP;
	G4double sigmaP;
	G4double meanT;
	G4double sigmaT;
	// internal variables
	G4int eventsGenerated;
	G4ParticleGun *particleGun;
	G4ParticleDefinition *particleDefinition;
	int eventID;
	int prevEventID;
	// Random number: Gaussian if sigma>0; flat if sigma<0 (-sigma is
	// the halfwidth).
	G4double myrand(G4double mean, G4double sigma) {
		if(sigma >= 0.0) return sigma*CLHEP::RandGauss::shoot() + mean;
		return mean+sigma-2.0*sigma*G4UniformRand();
	}
public:
	/// Constructor.
	BLCMDisotropicsource();

	// (accept the default copy constructor)

	/// commandName() returns "isotropicsource".
	virtual G4String commandName() { return "isotropicsource"; }
	
	/// command() implements the isotropicsource command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// getNEvents() returns the # events to process.
	virtual int getNEvents() const { return nEvents; }

	/// init() will initialize internal variables.
	virtual void init();

	/// generateReferenceParticle() generates the reference particle.
	virtual bool generateReferenceParticle(G4Event *event)
		{ return false; }

	/// nextBeamEvent() generates the next beam event.
	virtual bool nextBeamEvent(G4Event *event);

	/// summary() will print a summary, if necessary.
	virtual void summary() { }
};

BLCMDisotropicsource defineIsotropicsource;

BLCMDisotropicsource::BLCMDisotropicsource() : BLBeam() , BLCommand()
{
	registerCommand(BLCMDTYPE_BEAM);
	setSynopsis("Define an isotropic source.");
	setDescription("Generates tracks isotropically emanating from the "
		"point {x,y,z} (global coordinates).\n\n"
		"If any sigma is < 0, a uniform distribution is generated with "
		"HalfWidth=|sigma|.\n\n"
		"This element must be placed (via the place command)."
	);

	// initial default values:
	particle = "mu+";
	nEvents = 1;
	firstEvent = -1;
	x = y = z = 0.0;
	meanP = 200.0*MeV;
	sigmaP = 0.0;
	meanT = 0.0;
	sigmaT = 0.0;
	weight = 1.0;
	secondaryTrackID = 1001;
	eventsGenerated = 0;
	particleDefinition = 0;
	particleGun = 0;
	eventID = -1;
	prevEventID = -1;
}

int BLCMDisotropicsource::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	BLCMDisotropicsource *b = new BLCMDisotropicsource(*this);

	int retval = b->handleNamedArgs(namedArgs);

	// ensure beam is within the world
	BLGroup::getWorld()->setMinWidth(fabs(b->x)*2.0);
	BLGroup::getWorld()->setMinHeight(fabs(b->y)*2.0);
	BLGroup::getWorld()->setMinLength(fabs(b->z)*2.0);

	BLManager::getObject()->registerBeam(b);

	b->print("");

	return retval;
}

void BLCMDisotropicsource::defineNamedArgs()
{
	argString(particle,"particle","Beam particle name or PDGid");
	argInt(nEvents,"nEvents","Number of events to generate (default=1)");
	argInt(firstEvent,"firstEvent","First event # to generate (default "
		"is the next sequential eventID, 1 if none)");
	argDouble(x,"x","Beam location in X (mm)");
	argDouble(y,"y","Beam location in Y (mm)");
	argDouble(z,"z","Beam location in Z (mm)");
	argDouble(weight,"weight","Weight for events (1.0)");
	argInt(secondaryTrackID,"secondaryTrackID","The next TrackID for secondaries (1001)."); 
	argDouble(meanP,"meanP", "Beam mean momentum (MeV/c)");
	argDouble(meanP,"P","Synonym for meanP.");
	argDouble(sigmaP,"sigmaP","Beam sigma in P (MeV/c)");
	argDouble(meanT,"meanT","Beam mean in T (ns)");
	argDouble(meanT,"t","Synonym for meanT.");
	argDouble(sigmaT,"sigmaT","Beam sigma in T (ns)");
}

void BLCMDisotropicsource::init()
{
	eventsGenerated = 0;
	if(particleDefinition != 0) return;

	if(isdigit(particle(0)) || particle(0) == '-') {
		G4int pdgid = atoi(particle.c_str());
		particleDefinition = G4ParticleTable::GetParticleTable()->
							FindParticle(pdgid);
		if(!particleDefinition && pdgid > 1000000000) { // Ions
		    G4int Z = (pdgid/10000) % 1000;
		    G4int A = (pdgid/10) % 1000;
		    G4double excitationEnergy = 0.0*keV;
		    particleDefinition = G4ParticleTable::GetParticleTable()->
				GetIonTable()->GetIon(Z,A,excitationEnergy);
		} else if(!particleDefinition && pdgid < -1000000000) {
		    G4Exception("isotropicsource command","UnknownParticle",FatalException,
						"Anti-ions not supported");
		}
	} else {
		particleDefinition = G4ParticleTable::GetParticleTable()->
					FindParticle(particle);
	}
	if(!particleDefinition)
		G4Exception("isotropicsource command","UnknownParticle",FatalException,
						"Unknown particle type");
	particleGun = new G4ParticleGun(1);
	particleGun->SetParticleDefinition(particleDefinition);
}

bool BLCMDisotropicsource::nextBeamEvent(G4Event *event)
{
	G4double mass = particleDefinition->GetPDGMass();
	G4ThreeVector pos(x,y,z);
	G4ThreeVector direction;
	G4double time = 0.0;
	G4double momentum = 0.0;
	G4double ke = 0.0;
	G4int PDGid, trackID, parentID;

	BLManager *manager = BLManager::getObject();

	// default eventID -- changed when reading a file unless renumber!=0
	eventID = manager->getEventID();
	if(eventsGenerated == 0 && firstEvent != -1)
		eventID = firstEvent;

	for(;;) {
		if(++eventsGenerated > nEvents) return false;
		if(!manager->skipEvent(eventID)) break;
		++eventID;
	}

	// generate the event
	trackID = 11;
	parentID = 0;
	setRandomSeedToGenerate(eventID);
	double cosTheta = G4UniformRand()*2.0 - 1.0;
	double sinTheta = 1.0-cosTheta*cosTheta;
	if(fabs(sinTheta) < 1E-24) sinTheta = 1E-24;
	sinTheta = sqrt(sinTheta);
	double phi = G4UniformRand()*2.0*M_PI;
	direction[0] = sinTheta * cos(phi);
	direction[1] = sinTheta * sin(phi);
	direction[2] = cosTheta;
	momentum = myrand(meanP,sigmaP);
	time = myrand(meanT,sigmaT);
	trackID = 1;

	ke = sqrt(momentum*momentum + mass*mass) - mass;
	particleGun->SetParticleTime(time);
	particleGun->SetParticlePosition(pos);
	particleGun->SetParticleEnergy(ke);
	particleGun->SetParticleMomentumDirection(direction);
	particleGun->GeneratePrimaryVertex(event);
	event->SetEventID(eventID);
	event->GetPrimaryVertex()->SetWeight(weight);
	if(eventID != prevEventID) {
		setRandomSeedToTrack(eventID);
		manager->clearTrackIDMap();
		manager->setNextSecondaryTrackID(secondaryTrackID);
	}
	prevEventID = eventID;

	manager->setPrimaryTrackID(trackID,parentID);
	if(trackID >= secondaryTrackID) {
		G4Exception("beam command","Large Primary TrackID",JustWarning,
			"Confusion with secondary tracks is likely");
	}
	return true;
}
