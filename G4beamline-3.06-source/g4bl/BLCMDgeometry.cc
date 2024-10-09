//	BLCMDgeometry.cc
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
#include <time.h>

#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLGroup.hh"

/**	class BLCMDgeometry defines parameters for the geometry test.
 *
 *	The actual test is in BLGroupElement::testGeometry().
 **/
class BLCMDgeometry : public BLCommand, public BLCallback {
	G4int nPoints;
	G4int printGeometry;
	G4int visual;
	G4double tolerance;
	G4bool registered;
public:
	BLCMDgeometry();

	G4String commandName() { return "geometry"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	void callback(int type);
};

BLCMDgeometry defaultGeometry;

BLCMDgeometry::BLCMDgeometry() : BLCommand(), BLCallback()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("Arranges to perform a geometry test.");
	setDescription("The geometry test checks nPoints points on the surface "
		"of each element, verifying that they are inside the parent\n"
		"element and that they are not inside any sibling element.\n"
		"The default of 100 points is usually sufficient; 0 means\n"
		"omit the geometry test. The first 10-40 points (depending\n"
		"on element) test 'corners', the rest are randomly "
		"distributed on the surface. The default tolerance is "
		"0.002 mm.\n\n"
		"Points which fail are identified with a visible marker. "
		"set visible=1 to put a green marker at all tested points.\n\n"
		"This command is not placed into the geometry.\n\n");

	nPoints = 100;
	printGeometry = 0;
	visual = 0;
	tolerance = 0.002 * mm;
	registered = false;
}

int BLCMDgeometry::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	handleNamedArgs(namedArgs);

	if(!registered) {
		BLManager::getObject()->registerCallback(this,0);
		registered = true;
	}

	print("");

	return 0;
}

void BLCMDgeometry::defineNamedArgs()
{
	argInt(nPoints,"nPoints","The number of surface points to test per element");
	argInt(printGeometry,"printGeometry","Set nonzero to print the entire geometry.");
	argInt(visual,"visual","Set nonzero to display the geometry test points.");
	argDouble(tolerance,"tolerance","Tolerance for inside/outside tests (mm).");
}

void BLCMDgeometry::callback(int type)
{
	if(printGeometry) {
		printf("\nGeometry display:\n");
		BLManager::getObject()->displayGeometry();
	}

	if(nPoints > 0) {
		printf("\nGeometry test nPoints=%d tolerance=%.3f mm:\n",
					nPoints,tolerance);
		time_t start=time(0);
		int err = BLGroup::getWorld()->testGeometry(nPoints,tolerance,
								visual!=0);
		printf("Total geometry errors = %d  %ld seconds\n\n",
					err,time(0)-start);
		// if(err > 0) // exceptions generated during the test
		//     G4Exception("geometry","Geometry Errors",JustWarning,"");
	}
}
