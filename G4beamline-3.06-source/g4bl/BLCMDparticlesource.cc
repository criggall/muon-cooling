//	BLCMDparticlesource.cc
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

#include "G4GeneralParticleSource.hh"
#include "G4RunManager.hh"

#include "BLAssert.hh"
#include "BLBeam.hh"
#include "BLParam.hh"
#include "BLGroup.hh"
#include "BLCoordinates.hh"
#include "BLTrackFile.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

#define ALL_EVENTS 0x7FFFFFFF

/**	class BLCMDparticlesource implements the particlesource command.
 **/
class BLCMDparticlesource : public BLBeam, public BLCommand {
	G4int nEvents;
	G4int firstEvent;
	G4int lastEvent;
	G4int secondaryTrackID;
	G4GeneralParticleSource *particleSource;
	G4int eventsGenerated;
	G4int eventID;
	G4int prevEventID;
public:
	/// Constructor.
	BLCMDparticlesource();

	// (accept the default copy constructor)

	/// commandName() returns "particlesource".
	virtual G4String commandName() { return "particlesource"; }
	
	/// command() implements the particlesource command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// init() from BLBeam.
	void init() { }

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// getNEvents() returns the # events to process.
	virtual int getNEvents() const { return nEvents; }

	/// generateReferenceParticle() generates the reference particle.
	virtual bool generateReferenceParticle(G4Event *event);

	/// nextBeamEvent() generates the next beam event.
	virtual bool nextBeamEvent(G4Event *event);

	/// summary() will print a summary, if necessary.
	virtual void summary() { }
};

BLCMDparticlesource defineParticlesource;

BLCMDparticlesource::BLCMDparticlesource() 
{
	registerCommand(BLCMDTYPE_BEAM);
	setSynopsis("Interface to the Geant4 General Particle Source.");
	setDescription(
		"The Geant4 General Particle Source (GPS) is a very flexible "
		"and general way to generate events. It is controlled by "
		"Geant4 commands which should follow this command in "
		"the input.file. If you have a macro, you can include it using "
		"either the G4beamline 'include' command or the Geant4 command "
		"'/control/execute'. Note that G4beamline only recognizes "
		"Geant4 commands when the '/' is in column 1.\n\n"
		"Due to the design of the GPS, only one particlesource "
		"command can be used, but the GPS permits multiple sources "
		"to be combined.\n\n"
		"NOTE: the Geant4 General Particle Source inherently uses "
		"global coordinates, so this is most useful at the "
		"beginning of a beamline when global=centerline. Note also "
		"that it is very easy to generate a beam headed in the -z "
		"direction (this command will rotate to the +z direction: "
		"'/gps/ang/rot1 -1 0 0').\n\n"
		"This command is not placed into the geometry.\n\n"
		"To use this, see the User Manual for the GPS at "
		"http://reat.space.qinetiq.com/gps");
	// initial default values:
	nEvents = ALL_EVENTS;
	firstEvent = -1;
	lastEvent = ALL_EVENTS;
	secondaryTrackID = 1001;
	particleSource = 0;
	eventsGenerated = 0;
	eventID = ALL_EVENTS;
	prevEventID = -9999;
}

int BLCMDparticlesource::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(particleSource != 0) {
		printError("particlesource: only one particlecource command is permitted.");
		return -1;
	}

	int retval = handleNamedArgs(namedArgs);

	if(firstEvent != -1 && lastEvent != ALL_EVENTS)
		nEvents = lastEvent - firstEvent + 1;

	particleSource = new G4GeneralParticleSource();

	BLManager::getObject()->registerBeam(this);

	print("");

	return retval;
}

void BLCMDparticlesource::defineNamedArgs()
{
	argInt(nEvents,"nEvents","Number of events to process (default=1), "
		"set to lastEvent-firstEvent+1 if both are set.");
	argInt(firstEvent,"firstEvent","First event # to process (default "
		"is the next sequential eventID, 1 if none)");
	argInt(lastEvent,"lastEvent","Last  (highest) event # to process");
	argInt(secondaryTrackID,"secondaryTrackID","The next TrackID for secondaries (1001)."); 
}

bool BLCMDparticlesource::generateReferenceParticle(G4Event *event)
{
	return false;
}

bool BLCMDparticlesource::nextBeamEvent(G4Event *event)
{
	BLManager *manager = BLManager::getObject();

	int trackID = 1;
	int parentID = 0;

	BLAssert(particleSource != 0);

	// default eventID 
	eventID = manager->getEventID();
	if(eventsGenerated == 0 && firstEvent != -1)
		eventID = firstEvent;
	--eventID;
	do {
		if(++eventID > lastEvent) return false;
		if(++eventsGenerated > nEvents) return false;
	} while(BLManager::getObject()->skipEvent(eventID));

	setRandomSeedToGenerate(eventID);

	particleSource->GeneratePrimaryVertex(event);
	event->SetEventID(eventID);
	if(eventID != prevEventID) {
		setRandomSeedToTrack(eventID);
		manager->clearTrackIDMap();
		manager->setNextSecondaryTrackID(secondaryTrackID);
		prevEventID = eventID;
	}
	manager->setPrimaryTrackID(trackID,parentID);
	return true;
}
