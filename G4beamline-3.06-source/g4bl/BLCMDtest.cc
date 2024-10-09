//	BLCMDtest.cc
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
#include "BLParam.hh"

/**	class BLCMDtest implements the test command to test the random
 *	number generator.
 **/
class BLCMDtest : public BLCommand {
	// Random number: Gaussian if sigma>0; flat if sigma<0 (-sigma is
	// the halfwidth).
	G4double myrand(G4double mean, G4double sigma) {
		if(sigma >= 0.0) return sigma*CLHEP::RandGauss::shoot() + mean;
		return mean+sigma-2.0*sigma*G4UniformRand();
	}
public:
	/// Constructor.
	BLCMDtest();

	/// commandName() returns "test".
	virtual G4String commandName() { return "test"; }
	
	/// command() implements the test command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();
};

BLCMDtest defineTest;

BLCMDtest::BLCMDtest()
{
	registerCommand(BLCMDTYPE_OTHER);
	setSynopsis("test random number seeds.");
	setDescription("Test");
}

int BLCMDtest::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("Invalid test command -- need seed");
		return -1;
	}

	G4int seed = atoi(argv[0]);
	if(seed > 0) {
		CLHEP::HepRandom::setTheSeed(seed);
		CLHEP::RandGauss::setFlag(false);
		printf("seed=%6d   %8.4f",seed,myrand(0.0,1.0));
		for(int i=0; i<6; ++i) {
			for(int j=0; j<1000; ++j)
				myrand(0.0,1.0);
			printf(" %8.4f",myrand(0.0,1.0));
		}
		printf("\n");
	} else {
		printf("            %.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
			myrand(0.0,1.0),myrand(0.0,1.0),myrand(0.0,1.0),
			myrand(0.0,1.0),myrand(0.0,1.0),myrand(0.0,1.0),myrand(0.0,1.0));
	}

	return 0;
}

void BLCMDtest::defineNamedArgs()
{
}
