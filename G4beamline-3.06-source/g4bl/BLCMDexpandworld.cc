//	BLCMDexpandworld.cc
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
#include "BLGroup.hh"
#include "BLParam.hh"

/**	class BLCMDexpandworld -- expand the world to include specified points.
 **/
class BLCMDexpandworld : public BLCommand {
public:
	/// Constructor.
	BLCMDexpandworld();

	/// commandName() returns "expandworld".
	G4String commandName() { return "expandworld"; }
	
	/// command() implements the expandworld command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	void defineNamedArgs();
};

BLCMDexpandworld defineExpandworld;

BLCMDexpandworld::BLCMDexpandworld()
{
	registerCommand(BLCMDTYPE_LAYOUT);
	setSynopsis("Expand the world to include specified points.");
	setDescription("Each positional argument is x,y,z of a single point "
		"(global coordinates). The world will be expanded to "
		"include them all.");

}

int BLCMDexpandworld::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() < 1) {
		printError("Invalid expandworld command");
		return -1;
	}

	handleNamedArgs(namedArgs);

	double xMin=0.0, xMax=0.0, yMin=0.0, yMax=0.0, zMin=0.0, zMax=0.0;

	for(unsigned i=0; i<argv.size(); ++i) {
		std::vector<G4double> v = getList(argv[i],",");
		if(v.size() != 3) {
			printError("expandworld: Invalid Point.");
			continue;
		}
		if(xMin < v[0]) xMin = v[0];
		if(xMax > v[0]) xMax = v[0];
		if(yMin < v[1]) yMin = v[1];
		if(yMax > v[1]) yMax = v[1];
		if(zMin < v[2]) zMin = v[2];
		if(zMax > v[2]) zMax = v[2];
	}

	BLGroup::getWorld()->setMinWidth(fabs(xMin)*2.0);
	BLGroup::getWorld()->setMinHeight(fabs(yMin)*2.0);
	BLGroup::getWorld()->setMinLength(fabs(zMin)*2.0);
	BLGroup::getWorld()->setMinWidth(fabs(xMax)*2.0);
	BLGroup::getWorld()->setMinHeight(fabs(yMax)*2.0);
	BLGroup::getWorld()->setMinLength(fabs(zMax)*2.0);

	printf("expandworld  ");
	for(unsigned i=0; i<argv.size(); ++i) {
		printf("%s ",argv[i].c_str());
	}
	printf("\n");

	return 0;
}

void BLCMDexpandworld::defineNamedArgs()
{
}

