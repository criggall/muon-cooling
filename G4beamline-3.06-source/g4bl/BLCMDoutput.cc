//	BLCMDoutput.cc
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
#include "BLParam.hh"

extern void startupPrint(); // in g4beamline.cc

/**	class BLCMDoutput redirects stdout and stderr to a file.
 *
 **/
class BLCMDoutput : public BLCommand {
public:
	BLCMDoutput();

	G4String commandName() { return "output"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	// no named arguments.
};

BLCMDoutput defineOutput;

BLCMDoutput::BLCMDoutput()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("redirects stdout and stderr to a file");
	setDescription("output requires one argument, the new output file. "
		"Any output generated before this command will not appear in "
		"the file, so this command should come at the start of the "
		"input.file, preceded only by setting parameters. After the "
		"redirection, it will re-print the G4beamline version and the "
		"current parameter values to the file. The most common usage "
		"is to redirect output to a file named like the .root file, "
		"when that is determined by parameter values:\n"
		"    param -unset param1=1.0 param2=3.0\n"
		"    param histofile=$param1,$param2\n"
		"    output $histofile.out");
}

int BLCMDoutput::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1 || namedArgs.size() != 0) {
		printError("Invalid output command");
		return -1;
	}

	FILE *f = fopen(argv[0].c_str(),"w");
	if(!f) {
		printError("output: cannot open '%s'\n",argv[0].c_str());
		return 1;
	}
	fclose(f);

	printf("output redirected to '%s'\n",argv[0].c_str());

	freopen(argv[0].c_str(),"w",stdout);
	freopen(argv[0].c_str(),"w",stderr);

	startupPrint();
	printf("... output redirected to this file '%s'.\n",argv[0].c_str());
	Param.printParam();
	printf("\n");

	return 0;
}
