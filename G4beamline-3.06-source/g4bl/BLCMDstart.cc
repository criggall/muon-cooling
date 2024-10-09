//	BLCMDstart.cc
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

#include "BLCommand.hh"
#include "BLCoordinates.hh"
#include "BLParam.hh"

/**	BLCMDstart sets the starting place and orientation of the centerline 
 *	coordinates.
 *
 **/
class BLCMDstart : public BLCommand {
	G4double x;
	G4double y;
	G4double z;
	G4double initialZ;
	G4String rotation;
	G4double radiusCut;
	G4int ring;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDstart();

	/// Destructor.
	virtual ~BLCMDstart() { }

	/// commandName() returns "start".
	G4String commandName() { return "start"; }

	/// command() implements the start command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();
};

BLCMDstart defaultStart;	// default object

// Default constructor - be sure to use the default constructor BLCommand()
BLCMDstart::BLCMDstart() : BLCommand()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_LAYOUT);
	setSynopsis("Define the initial start of centerline coordinates.");
	setDescription("If used, this command must come before any other "
		"command that puts an element into the world or affects the "
		"centerline coordinates (place, beam, corner, cornerarc, and "
		"reference commands). This command may not always be needed, "
		"but it is needed to eliminate the ambiguities in "
		"the global to centerline coordinate transform, and "
		"when simulating a ring to ensure that "
		"sensible values of the centerline coordinates are used.\n\n"
		"Note that the radiusCut is important to reduce or eliminate "
		"ambiguities in the global to centerline coordinate transform. "
		"It can also be used to 'shield' the beamline to prevent "
		"particles from taking unusual paths around the outside "
		"of beamline elements.");

	// provide initial values for fields
	x = 0.0;
	y = 0.0;
	z = 0.0;
	initialZ = 0.0;
	rotation = "";
	radiusCut = 0.0;
	ring = 0;
}

int BLCMDstart::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	int retval = handleNamedArgs(namedArgs);

	if(radiusCut <= 0.0)
		G4Exception("start command","No radius cut",JustWarning, "");

	G4double point[4];
	point[0] = x;
	point[1] = y;
	point[2] = z;
	point[3] = 0.0;

	G4RotationMatrix *rot = stringToRotationMatrix(rotation);

	BLCoordinates::start(point,initialZ,rot,radiusCut,ring);

	print("");

	return retval;
}

void BLCMDstart::defineNamedArgs()
{
	argDouble(x,"x","The global x position of the start.");
	argDouble(y,"y","The global y position of the start.");
	argDouble(z,"z","The global z position of the start.");
	argDouble(initialZ,"initialZ","The initial centerline z value.");
	argString(rotation,"rotation","The initial rotation of the centerline.");
	argDouble(radiusCut,"radiusCut","The radius cut for the initial segment (mm).");
	argInt(ring,"ring","Set nonzero to indicate a ring is present.");
}

