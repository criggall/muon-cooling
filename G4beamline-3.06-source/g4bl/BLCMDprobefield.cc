//	BLCMDprobefield.cc
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
#include <stdarg.h>
#include <vector>
#ifndef WIN32
#include <unistd.h>
#endif

#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLGlobalField.hh"
#include "BLWriteAsciiFile.hh"
#include "BLMPI.hh"

/**	class BLCMDprobefield is a command to print fields at specified points
 *
 **/
class BLCMDprobefield : public BLCommand, public BLCallback {
	G4String inputFile;
	G4String outputFile;
	FILE *in;
	FILE *out;
	G4int exit;
	std::vector<G4String> args;
public:
	BLCMDprobefield();

	~BLCMDprobefield() { }

	G4String commandName() { return "probefield"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();
	
	/// callback() from BLCallback.
	void callback(int type);

	void handleLine(G4String line);
};

BLCMDprobefield defaultProbeField;	// registers the command, and holds
					// default values of the arguments.
BLCMDprobefield::BLCMDprobefield()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Prints B and E fields at specified points.");
	setDescription("Intended primarily for debugging. Prints Bx,By,Bz "
		"in Tesla, and Ex,Ey,Ez in MegaVolts/meter. "
		"Each input line is x,y,z,t separated by spaces or commas; "
		"omitted values are set to 0.0.\n\n"
		"Field Values are printed after the reference particle is "
		"tracked. Only global coordinates are used.\n\n"
		"Positional arguments are used as input lines before reading "
		"the inputFile.\n\n"
		"This command is not placed into the geometry.");
	inputFile = "-";
	outputFile = "-";
	in = 0;
	out = 0;
	exit = 0;
}

int BLCMDprobefield::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	BLCMDprobefield *p = new BLCMDprobefield(defaultProbeField);
	int retval = p->handleNamedArgs(namedArgs);

	p->args = argv;

	// MPI mode: performed only in rank 1.
	// Non-MPI mode: isRank1() returns true.
	if(BLMPI::isRank1()) {
		// register the object to be called back to do its print.
		// 1 => after reference particle is tracked.
		BLManager::getObject()->registerCallback(p,1);
	} else {
		printf("probefield: performed only in rank 1.\n");
	}

	return retval;
}

void BLCMDprobefield::defineNamedArgs()
{
	argString(inputFile,"inputFile","Filename for reading list of points (- = stdin)");
	argString(outputFile,"outputFile","Filename for output (- = stdout)");
	argInt(exit,"exit","Set nonzero to exit after printing field");
}

void BLCMDprobefield::callback(int callback_type)
{
	if(inputFile == "-") {
		in = stdin;
	} else if(inputFile == "") {
		in = 0;
	} else {
		in = fopen(inputFile.c_str(),"r");
		if(!in)
			printError("probefield cannot open inputFile '%s'",
							inputFile.c_str());
	}

	if(outputFile == "-") {
		out = stdout;
	} else {
		out = BLWriteAsciiFile::fopen(outputFile.c_str());
		if(!out) {
			printError("probefield cannot write '%s'",
							outputFile.c_str());
			return;
		}
	}

	printf("probefield is printing field values:\n");
	fprintf(out,"#x y z t Bx By Bz Ex Ey Ez\n");
	for(unsigned i=0; i<args.size(); ++i)
		handleLine(args[i]);
	while(in != 0) {
		if(in == stdin) {
			fflush(stdout);
			fprintf(stderr,"point: ");
		}
		char line[1024];
		if(!fgets(line,sizeof(line),in)) break;
		handleLine(line);
	}

	// let fclose happen when the program exits
	//if(out != stdout) BLWriteAsciiFile::fclose(out);

	if(exit) {
		G4Exception("probefield","Exiting",FatalException,"");
	}
}

void BLCMDprobefield::handleLine(G4String line)
{
	BLGlobalField *gf = BLGlobalField::getObject();

	for(;;) {
		G4String::size_type pos = line.find(",");
		if(pos == G4String::npos) break;
		line[pos] = ' ';
	}

	G4double pos[4];
	pos[0] = pos[1] = pos[2] = pos[3] = 0.0;
	sscanf(line,"%lf %lf %lf %lf",&pos[0],&pos[1],&pos[2],&pos[3]); 
	pos[0] *= mm; pos[1] *= mm; pos[2] *= mm; pos[3] *= ns;
	G4double field[6];
	gf->GetFieldValue(pos,field);
	fprintf(out,"%.3f %.3f %.3f %.3f  %.6f %.6f %.6f  %.6f %.6f %.6f\n",
		pos[0]/mm,pos[1]/mm,pos[2]/mm,pos[3]/ns,
		field[0]/tesla,field[1]/tesla,field[2]/tesla,
		field[3]/(megavolt/meter),field[4]/(megavolt/meter),
		field[5]/(megavolt/meter));
}
