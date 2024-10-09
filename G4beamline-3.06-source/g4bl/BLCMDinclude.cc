//	BLCMDinclude.cc
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

/**	class BLCMDinclude implements the include command to include a file.
 *
 **/
class BLCMDinclude : public BLCommand {
public:
	BLCMDinclude();

	G4String commandName() { return "include"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	// no arguments.
};

BLCMDinclude defineInclude;

BLCMDinclude::BLCMDinclude()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("includes a command file.");
	setDescription("include requires one argument, the file to include.");
}

int BLCMDinclude::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1 || namedArgs.size() != 0) {
		printError("Invalid include command");
		return -1;
	}

	BLCommand::readFile(argv[0]);

	return 0;
}
