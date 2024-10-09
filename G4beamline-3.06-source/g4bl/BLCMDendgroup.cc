//	BLCMDendgroup.cc
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
//	NOTE: the class for the group command is BLGroup, not BLCMDgroup.

#include "BLGroup.hh"

/**	class BLCMDendgroup implements the endgroup command.
 *
 **/
class BLCMDendgroup : public BLCommand {
public:
	BLCMDendgroup();

	G4String commandName() { return "endgroup"; }

	// no arguments.
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);
};

BLCMDendgroup defineEndGroup;

BLCMDendgroup::BLCMDendgroup()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("ends a group definition.");
	setDescription("The group may then be placed as any other element.\n"
		"If the group was not given a length via an argument,\n"
		"endgroup computes the length and adjusts the offsets\n"
		"of all elements in the group so they refer to the center\n"
		"of the group.");
}

int BLCMDendgroup::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 0 || namedArgs.size() != 0) {
		printError("Invalid endgroup command");
		return -1;
	}

	BLGroup *c = BLGroup::getCurrent();
	if(c == BLGroup::getWorld()) {
		printError("Invalid endgroup command");
		return -1;
	}
	c->end();

	printf("endgroup %-6s length=%.1f width=%.1f height=%.1f\n",
		c->getName().c_str(),c->length,c->width,c->height);

	return 0;
}
