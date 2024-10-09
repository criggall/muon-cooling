//	BLCMDdo.cc
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

/**	class BLCMDdo implements the do command for loops
 *
 *	do i 1 10 [1]
 *		... statements of the loop
 *	enddo
 *
 *	Sets i successively from 1 to 10, executing the statements of the loop.
 *	The increment is 1 by default; negative increments work as in FORTRAN.
 *
 *	NOTE: the enddo command is not a real command in that it is not
 *	a class. It is detected textually by this command class.
 *
 **/
class BLCMDdo : public BLCommand {
public:
	BLCMDdo();
	BLCMDdo(BLCMDdo& r) : BLCommand() { }

	G4String commandName() { return "do"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);
};

BLCMDdo defineDo;

BLCMDdo::BLCMDdo()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("Do loop.");
	setDescription("Syntax:\n    do i 1 10 [1]\n        commands ...\n"
		"    enddo\n\nSets the parameter i to "
		"values 1, 2, 3, ... 10, and executes the commands. "
		"The increment is 1 by default, and negative increments are "
		"allowed (with limits reversed). 'do i 1 0' and 'do i 0 1 -1' "
		"will execute no commands. After completion, i retains its "
		"last value. Do-s, for-s, and multi-line if-s can be "
		"nested in any manner.\n\n"
		"Note: the do command will not work from standard-input or "
		"as a command in a single-line if.");
}

int BLCMDdo::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() < 3 || argv.size() > 4) {
		printError("Invalid do command");
		return -1;
	}

	int i = atoi(argv[1].c_str());
	int n = atoi(argv[2].c_str());
	int incr = (argv.size() > 3 ? atoi(argv[3].c_str()) : 1);

	Param.setParam(argv[0],i);

	// special test to execute no commands
	if(incr>0 ? i > n : i < n) {
		printf("do %s=%d,%d,%d -- no commands executed\n",
						argv[0].c_str(),i,n,incr);
		// skip to the matching enddo
		int level = 1;
		for(;;) {
			G4String *line = BLCommand::getNextCommand(); // trimmed
			if(!line) {
				printError("Unexpected EOF in do loop.");
				return -1;
			}
			if(line->find("enddo") == 0) {
				if(--level == 0)
					break;
			} else if(line->find("do") == 0 &&
					line->find_first_of(" \t") == 2) {
				++level;
			}
		}
		return 0;
	}

	// get current pos so we can return to it each time through the loop
	BLCommandPos pos = BLCommand::getPos();
	if(!pos.isValid()) {
		printError("do: cannot implement do loop on this input device.");
		return -1;
	}

	printf("do %s=%d,%d,%d\n",argv[0].c_str(),i,n,incr);

	// here is the actual do loop
	while(incr>0 ? i <= n : i>= n) {
		BLCommand::setPos(pos);
		Param.setParam(argv[0],i);
		printf("(do %s=%d)\n",argv[0].c_str(),i);
		for(;;) {
			G4String *line = BLCommand::getNextCommand(); // trimmed
			if(!line) {
				printError("Unexpected EOF in do loop.");
				return -1;
			}
			if(line->find("enddo") == 0)
				break;
			BLCommand::doCommand(*line);
		}
		i += incr;
	}
	printf("enddo\n");
	
	return 0;
}
