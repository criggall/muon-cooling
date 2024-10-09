//	BLCMDexit.cc
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

/**	class BLCMDexit implements the exit command to exit the input file.
 *
 **/
class BLCMDexit : public BLCommand {
public:
	BLCMDexit();

	G4String commandName() { return "exit"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	// no arguments.
};

BLCMDexit defineExit;

BLCMDexit::BLCMDexit()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("exit a command file.");
	setDescription("The exit command ceases reading the input file,\n"
		"and starts the simulation immediately (ignoring the\n"
		"remainder of the input file).");
}

int BLCMDexit::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	// skip to EOF on current input
	for(;;) {
		if(!BLCommand::getNextCommand()) break;
	}
	return 0;
}
