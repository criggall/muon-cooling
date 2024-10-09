//	BLCMDrandomseed.cc
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
#include <string.h>

#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLBeam.hh"

class BLCMDrandomseed : public BLCommand, public BLManager::EventAction {
	G4int print;
public:
	BLCMDrandomseed();

	G4String commandName() { return "randomseed"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	/// BeginOfEventAction from BLManager::EventAction
	void BeginOfEventAction(const G4Event* event);

	/// EndOfEventAction from BLManager::EventAction
	void EndOfEventAction(const G4Event* event);
};

BLCMDrandomseed defaultRandomseedCommand;

BLCMDrandomseed::BLCMDrandomseed()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("control pseudo random number generator seeds");
	setDescription("This randomseed command controls the pseudo random "
		"number generator seed at the start of each event.\n"
		"The unnamed argument can be any of (case insensitive):\n"
		"   EventNumber\n"
		"   None\n"
		"   Time\n"
		"   Set 12345\n"
		"   Now 12345\n"
		"EventNumber is the default and "
		"permits events to be re-run; None does not re-seed the "
		"PRNG at each event, and Time is like None after seeding "
		"with the time of day in microseconds; Set (Now) seeds the "
		"generator "
		"immediately with the value of the second argument (a long), "
		"and then acts like None.");
}

int BLCMDrandomseed::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() < 1) {
		printError("randomseed: invalid command -- need method");
		return -1;
	}

	handleNamedArgs(namedArgs);

	G4String method = argv[0];
	for(unsigned i=0; i<method.size(); ++i)
		method[i] = tolower(method[i]);

	if(method == "set" || method == "now") {
		if(argv.size() < 2 || argv[1].size() == 0) {
			printError("randomseed: no value to set");
			return 1;
		}
		char *p;
		long seed = strtol(argv[1].c_str(),&p,0);
		if(*p != '\0') {
			printError("randomseed: invalid value to set");
			return 1;
		}
		if(seed < 0L) seed = 0x7FFFFFFEL;
		CLHEP::HepRandom::setTheSeed(seed);
		CLHEP::RandGauss::setFlag(false);
		BLManager::getObject()->setPRNGSeedMethod(NO_SEED);
	} else if(method =="eventnumber")
		BLManager::getObject()->setPRNGSeedMethod(EVENT_NUMBER);
	else if(method =="none")
		BLManager::getObject()->setPRNGSeedMethod(NO_SEED);
	else if(method =="time")
		BLManager::getObject()->setPRNGSeedMethod(TIME_US);
	else
		printError("randomseed: invalid method");

	BLCommand::print("");

	BLManager::getObject()->registerEventAction(this);

	return 0;
}

void BLCMDrandomseed::defineNamedArgs()
{
	argInt(print,"print","Set to 1 to print the next random number at the end of each event.");
}

void BLCMDrandomseed::BeginOfEventAction(const G4Event* event)
{
	if(print)
		printf("randomseed: Event %d Begin next Random: %.8f\n",
				event->GetEventID(),CLHEP::RandGauss::shoot());
}

void BLCMDrandomseed::EndOfEventAction(const G4Event* event)
{
	if(print)
		printf("randomseed: Event %d End next Random: %.8f\n",
				event->GetEventID(),CLHEP::RandGauss::shoot());
}

