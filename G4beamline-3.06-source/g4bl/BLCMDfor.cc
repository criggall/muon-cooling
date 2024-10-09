//	BLCMDfor.cc
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

#include <stdlib.h>
#include "BLCommand.hh"
#include "BLParam.hh"

/**	class BLCMDfor implements the for command for loops
 *
 *	for i v1 v2 v3 ...
 *		... statements of the loop
 *	endfor
 *
 *	Sets i successively to v1, v2, v3, executing the statements of the loop.
 *
 *	NOTE: the endfor command is not a real command in that it is not
 *	a class. It is detected textually by this command class.
 *
 **/
class BLCMDfor : public BLCommand {
public:
	BLCMDfor();
	BLCMDfor(BLCMDfor& r) : BLCommand() { }

	G4String commandName() { return "for"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);
};

BLCMDfor defineFor;

BLCMDfor::BLCMDfor()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("For loop.");
	setDescription("Syntax:\n    for i v1 v2 v3 ...\n        commands ...\n"
		"    endfor\n\nSets the parameter i to "
		"values v1, v2, v3, ..., and executes the commands. Values can "
		"be any strings, including numbers, except they cannot contain "
		"an '=' (named parameter). After completion, i "
		"retains its last value.\n\n"
		"Do-s, for-s, and multi-line if-s can be "
		"nested in any manner.\n\n"
		"Note: the for command will not work from standard-input or "
		"as a command in a single-line if.");
}

int BLCMDfor::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() < 2) {
		printError("Invalid for command");
		return -1;
	}

	// get current pos so we can return to it each time through the loop
	BLCommandPos pos = BLCommand::getPos();
	if(!pos.isValid()) {
		printError("for: cannot implement for loop on this input device.");
		return -1;
	}

	printf("for %s=",argv[0].c_str());
	for(int i=1; i<argv.size(); ++i)
		printf("'%s' ",argv[i].c_str());
	printf("\n");

	// here is the actual for loop
	for(unsigned i=1; i<argv.size(); ++i) {
		BLCommand::setPos(pos);
		Param.setParam(argv[0],argv[i]);
		printf("(for %s='%s')\n",argv[0].c_str(),argv[i].c_str());
		for(;;) {
			G4String *line = BLCommand::getNextCommand(); // trimmed
			if(!line) {
				printError("Unexpected EOF in for loop.");
				return -1;
			}
			if(line->find("endfor") == 0)
				break;
			BLCommand::doCommand(*line);
		}
	}
	printf("endfor\n");
	
	return 0;
}
