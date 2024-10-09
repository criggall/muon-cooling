//	BLCMDeventcuts.cc
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

#include "BLCommand.hh"
#include "BLManager.hh"

class BLCMDeventcuts : public BLCommand {
	G4String keep;
	G4String skip;
public:
	BLCMDeventcuts();

	G4String commandName() { return "eventcuts"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	void readFile(G4String filename, bool skipflag);
};

BLCMDeventcuts defaultEventcutsCommand;

BLCMDeventcuts::BLCMDeventcuts()
{
	registerCommand(BLCMDTYPE_CUTS);
	setSynopsis("Implements cuts on event number via lists in files.");
	setDescription("The files are ASCII, with one event number per line. "
		"If the keep file is not empty, only event numbers listed in "
		"it will be analyzed (except those listed in the skip file). "
		"The skip file lists event numbers that will be skipped, and "
		"event numbers listed in both files will be skipped. "
		"When reading the files, lines beginning with '#' are ignored; "
		"blank lines are interpreted as event 0.");

	keep = "";
	skip = "";
}

int BLCMDeventcuts::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	int retval = handleNamedArgs(namedArgs);

	print("");

	if(skip.size() > 0) readFile(skip, true);
	if(keep.size() > 0) readFile(keep, false);

	return retval;
}

void BLCMDeventcuts::defineNamedArgs()
{
	argString(keep,"keep","The file containing a list of event numbers to analyze (default is all).");
	argString(skip,"skip","The file containing a list of event numbers to skip.");
	argString(keep,"filename","Synonym for keep.");
	argString(keep,"file","Synonym for keep.");
}

void BLCMDeventcuts::readFile(G4String filename, bool skipflag)
{
	FILE *in = fopen(filename.c_str(),"r");
	if(!in) {
		printError("eventcuts: cannot read '%s'",filename.c_str());
		return;
	}

	char line[128];
	while(fgets(line,sizeof(line),in)) {
		if(line[0] == '#') continue;
		int ev=strtol(line,0,0);
		if(skipflag)
			BLManager::getObject()->setSkipEvent(ev);
		else
			BLManager::getObject()->setKeepEvent(ev);
	}
	fclose(in);
}
