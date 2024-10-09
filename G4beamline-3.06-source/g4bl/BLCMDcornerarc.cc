//	BLCMDcornerarc.cc
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

#include "G4UserLimits.hh"
#include "G4StepPoint.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"

#include "BLParam.hh"
#include "BLCommand.hh"
#include "BLGroup.hh"
#include "BLCoordinates.hh"
#include "BLCoordinateTransform.hh"

/**	class BLCMDcornerarc - implements a cornerarc in the centerline.
 *	
 *	A cornerarc is an approximation to an arc in the centerline,
 *	implemented by 3 corners such that the overall bend angle and 
 *	centerline path length are the same as the arc would be.
 *
 *	The error (largest distance between the arc and the actual centerline)
 *	increases with angle, but is less than 10% of radius for 90 degrees.
 **/
class BLCMDcornerarc : public BLCommand {
	G4double z;
	G4double centerRadius;
	G4double angle;
	G4String rotation;
	G4double radiusCut;
	G4double finalZ;
public:
	/// Default constructor.
	BLCMDcornerarc();

	/// Destructor.
	~BLCMDcornerarc();

	/// Copy constructor.
	BLCMDcornerarc(const BLCMDcornerarc& r);

	/// commandName() returns "cornerarc".
	G4String commandName() { return "cornerarc"; }

	/// command() implements the cornerarc command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();
};

BLCMDcornerarc defaultCornerArc;

BLCMDcornerarc::BLCMDcornerarc() : BLCommand(), rotation()
{
	registerCommand(BLCMDTYPE_LAYOUT);
	setSynopsis("Implement a cornerarc in the centerline.");
	setDescription("The centerline is bent by a rotation.\n"
		"Three corners are used to approximate an arc; the total\n"
		"angle and path-length are equal to those for the arc.\n"
		"Should be used immediately after an idealsectorbend or\n"
		"a genericbend.\n"
		"The z value is for the front face of the arc.\n"
		"\nNOTE: This command is self-placing, do not use the place\n"
		"command; it also affects all following placements, and\n"
		"it cannot be issued inside a group.\n\n"
		"This command is well matched to a sector bend, but can also "
		"be used with a normal bending magnet -- normally the magnet "
		"is placed before the cornerarc and is rotated by half the "
		"bend angle.\n\n"
		"The only useful rotations are those around the centerline Z.\n"
		"\nz, angle, and centerRadius are required parameters.\n"
		"\nNOTE: all placements before this command must have z\n"
		"values before the corner, and all placements after this\n"
		"command must have z values after the corner. Note also\n"
		"that the angle is limited to 90 degrees.\n\n"
		"Note that the radiusCut is important to reduce or eliminate "
		"ambiguities in the global to centerline coordinate transform. "
		"It can also be used to 'shield' the beamline to prevent "
		"particles from taking unusual paths around the outside "
		"of beamline elements.\n\n"
		"This command places itself into the geometry.\n\n"
		"Note that section 4.6 of the User's Guide has a dimensioned "
		"drawing of a cornerarc.");
	// initial field values
	z = 0.0;
	centerRadius=0.0;
	angle = 0.0;
	rotation = "";
	radiusCut = 0.0;
	finalZ = 0.0;
}

BLCMDcornerarc::~BLCMDcornerarc()
{
}

BLCMDcornerarc::BLCMDcornerarc(const BLCMDcornerarc& r) : BLCommand(r), rotation(r.rotation)
{
	z = r.z;
	centerRadius = r.centerRadius;
	angle = r.angle;
	rotation = r.rotation;
	radiusCut = r.radiusCut;
	finalZ = 0.0;
}
int BLCMDcornerarc::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() >= 1 && argv[0] == "default") {
		return defaultCornerArc.handleNamedArgs(namedArgs);
	}

	if(BLGroup::getWorld() != BLGroup::getCurrent()) {
		printError("cornerarc: cannot be issued inside a group.");
		return -1;
	}

	BLCMDcornerarc *t = new BLCMDcornerarc(defaultCornerArc);
	t->z = -DBL_MAX;
	int retval = t->handleNamedArgs(namedArgs);

	if(t->z == -DBL_MAX) {
		printError("cornerarc: z must be specified");
		return -1;
	}

	G4double angleAbs = fabs(t->angle);

	if(angleAbs < 0.00001 || angleAbs > 3.1416/2.0) {
		printError("cornerarc: invalid angle (0<|angle|<=90 degrees)");
		return -1;
	}
	if(t->centerRadius < 0.001*mm) {
		printError("cornerarc: invalid centerRadius");
		return -1;
	}

	// Compute the angles for the 3 corners that approximate the arc

	G4double beta = acos(sin(angleAbs/2.0)/(angleAbs/2.0));
	G4double phi1 = angleAbs/2.0 - beta;
	G4double phi2 = beta + beta;
	G4double deltaz = t->centerRadius * angleAbs / 2.0;

	// Get the rotationMatrix around z (actually more general)
	char temp[64];
	G4RotationMatrix zrot, invzrot;
	if(t->rotation.size() > 0) {
		sprintf(temp,"%s",t->rotation.c_str());   
		zrot = *stringToRotationMatrix(temp);
		invzrot = zrot.inverse();
        }

	// Get the rotationMatrix-s for the 3 corners
        // When specified, a cornerarc z-rotation (around the beam axis) is
	// applied before and undone after each of the three corners of the
	// cornerarc. This is done to preserve the orientation of the
	// centerline transverse coordinates.
	// This corrected code came from Jean-Francois Ostiguy.
 	sprintf(temp,"Y%.6f",(t->angle>0.0 ? phi1 : -phi1)/deg); 
	G4RotationMatrix m1 = zrot * (*stringToRotationMatrix(temp)) * invzrot;
	sprintf(temp,"Y%.6f",(t->angle>0.0 ? phi2 : -phi2)/deg);
	G4RotationMatrix m2 = zrot * (*stringToRotationMatrix(temp)) * invzrot;
	sprintf(temp,"Y%.6f",(t->angle>0.0 ? phi1 : -phi1)/deg);
	G4RotationMatrix m3 = zrot * (*stringToRotationMatrix(temp)) * invzrot;

	BLCoordinates::corner(t->z,              &m1, t->radiusCut);
	BLCoordinates::corner(t->z+deltaz,       &m2, t->radiusCut);
	BLCoordinates::corner(t->z+deltaz+deltaz,&m3, t->radiusCut);

	t->finalZ = t->z+deltaz+deltaz;

	Param.setParam("Zcl",t->finalZ);

	t->print("");

	return retval;
}

void BLCMDcornerarc::defineNamedArgs()
{
	argDouble(z,"z","The centerline Z of the cornerarc (mm).",mm);
	argDouble(centerRadius,"centerRadius","The radius of the centerline arc (mm).",mm);
	argDouble(angle,"angle","The angle of bend, >0=left, <0=right (degrees).",deg);
	argString(rotation,"rotation","The rotation of the cornerarc (see above).");
	argDouble(radiusCut,"radiusCut","The radius cut for this following segment (mm).",mm);
}
