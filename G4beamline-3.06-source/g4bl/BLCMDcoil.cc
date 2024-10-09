//	BLCMDcoil.cc
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

#ifndef G4BL_GSL
int BLCMDcoil_dummy=0;
#else

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

#include "BLCoil.hh"

/**	BLCMDcoil implements the coil command.
 *
 *	Note this class is derived from BLCoil, where the variables are
 *	defined and the real work is done.
 **/
class BLCMDcoil : public BLCommand, public BLCoil {
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDcoil();

	/// Destructor.
	virtual ~BLCMDcoil();

	/// Copy constructor.
	BLCMDcoil(const BLCMDcoil& r);

	/// commandName() returns "coil".
	G4String commandName() { return "coil"; }

	/// command() implements the coil command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the arguments to the command.
	void defineNamedArgs();
};

BLCMDcoil defaultCoil;

BLCMDcoil::BLCMDcoil() : BLCommand(), BLCoil()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("defines a coil (part of a solenoid)");
	setDescription("A coil is a geometrical tube that can carry current\n"
		"when part of a solenoid. The field is computed for a set\n"
		"of nSheets infinitely-thin current sheets evenly spread\n"
		"radially. The solenoid specifies the actual current.\n"
		"For tracking the computation is too slow, so a field map\n"
		"on a grid in r and z is computed and written to filename\n"
		"(defaults to coilname.dat).\n"
		"While there are lots of parameters specifying the field map\n"
		"it is recommended to simply accept the defaults for all\n"
		"but innerRadius, outerRadius, length, material, and\n"
		"possibly tolerance. The other parameters will be determined\n"
		"so that the largest error is less than tolerance times the\n"
		"value of Bz at the center of the coil.\n"
		"If mapFile is given, the file is read in BLFieldMap format.\n"
		"The cache file contains the parameters, and is a field map "
		"in a binary format; it is automatically regenerated if any "
		"parameter changes.\n\n"
		"Note that maxR and maxZ are by default determined to be large "
		"enough so that the field falls below tolerance; this can be "
		"quote large.\n\n"
		"This command is not placed into the geometry.");
}

BLCMDcoil::~BLCMDcoil() 
{
}

BLCMDcoil::BLCMDcoil(const BLCMDcoil& r) : BLCommand(r), BLCoil(r)
{
}

int BLCMDcoil::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("coil: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultCoil.handleNamedArgs(namedArgs);
	}

	BLCMDcoil *t = new BLCMDcoil(defaultCoil);
	t->name = argv[0];
	int retval = t->handleNamedArgs(namedArgs);

	if(t->filename == "")
		t->filename = t->name + ".dat";

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->printCoil();

	// add this coil to the list, and compute its fieldMap
	mapCoil[t->name] = t;
	t->generateFieldMap();

	return retval;
}

void BLCMDcoil::defineNamedArgs()
{
	argDouble(innerRadius,"innerRadius","Inside radius of the coil (mm)",mm,"",false);
	argDouble(outerRadius,"outerRadius","Outside radius of the coil (mm)",mm,"",false);
	argDouble(length,"length","Length of the coil along z (mm)",mm,"",false);
	argString(material,"material","The material of the conductor (default=Cu)",false);
	argDouble(tolerance,"tolerance","The acceptable tolerance of the map.",1,"",false);
	argInt(nSheets,"nSheets","Number of sheets used in the computation",false);
	argDouble(maxR,"maxR","Maxmum r value for the map (automatically determined by default).",mm,"",false);
	argDouble(maxZ,"maxZ","Maxmum z value for the map (automatically determined by default).",mm,"",false);
	argDouble(dR,"dR","R interval between points of the map (automatically determined by default).",mm,"",false);
	argDouble(dZ,"dZ","Z interval between points of the map (automatically determined by default).",mm,"",false);
	argString(filename,"filename","Filename for cache; deaults to name.dat.",false);
	argString(mapFile,"mapFile","Filename for map (e.g. from TOSCA).",false);
	argInt(exactComputation,"exactComputation","Set nonzero to use the exact computation without any field map (0).");
}
#endif // G4BL_GSL
