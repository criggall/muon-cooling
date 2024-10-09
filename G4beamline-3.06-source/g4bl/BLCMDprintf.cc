//	BLCMDprintf.cc
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

#include <vector>
#include <map>
#include <stdio.h>

#include "BLManager.hh"
#include "BLParam.hh"
#include "BLEvaluator.hh"
#include "BLCoordinates.hh"
#include "BLWriteAsciiFile.hh"
#include "BLMPI.hh"

const unsigned MAX_EXPR=16;	// max # expressions in printf command

/**	class BLCMDprintf will print track variables and expressions at a given
 *	z position.
 **/
class BLCMDprintf : public BLCommand, public BLManager::ZSteppingAction,
					public BLManager::EventAction {
	G4String z;
	G4String zloop;
	G4String require;
	G4String file;
	G4String format;
	G4double noNewline;
	G4String coordinates;
	BLCoordinateType coordinateType;
	std::vector<G4String> expr;
	BLEvaluator eval;
public:
	/// Constructor.
	BLCMDprintf();
	BLCMDprintf(BLCMDprintf& r);

	/// commandName() returns "printf".
	virtual G4String commandName() { return "printf"; }
	
	/// command() implements the printf command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// UserZSteppingAction() from BLManager::ZSteppingAction.
	void UserZSteppingAction(const G4Track *track);

	// from BLManager::EventAction
	public:	virtual void BeginOfEventAction(const G4Event* event) { }
	public:	virtual void EndOfEventAction(const G4Event* event)
		{ BLMPI::rank0BufferFlush(); }
};

BLCMDprintf defaultPrintf;

BLCMDprintf::BLCMDprintf() : BLCommand(), BLManager::ZSteppingAction(), 
				expr(), eval()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("prints track variables and expressions");
	setDescription("This is an interface to the C printf() function.\n"
		"The first positional argument is the format, and the "
		"following positional arguments are double expressions "
		"printed with the % fields in the format. Up to 16 "
		"expressions can be printed. The print is performed "
		"only if the 'required' expression is nonzero, as each "
		"track reaches one of the Z positions in the 'z' argument "
		"(centerline coordinates). Multiple printf commands "
		"with the same 'file' will be combined into the file "
		"as tracks reach any of their Z positions. More than 16 "
		"expressions can be broken into multiple printf-s with "
		"noNewline=1 for all but the last.\n\n"
		"The following variables can be used in expressions:\n"
		"   x,y,z,t,Px,Py,Pz\n"
		"   PDGid,EventID,TrackID,ParentId,Weight\n"
		"   Bx,By,Bz, Ex,Ey,Ez\n"
		"   tune (nonzero only for the Tune particle)\n"
		"   reference (nonzero only for the Reference particle)\n"
		"   beam (nonzero for any Beam particle)\n\n"
		"Each value in z and zloop can be an expression using double "
		"constants and the usual C operators and functions.\n\n"
		"Example:\n"
		" printf z=0 'Momentum is %.3f GeV/c' sqrt(Px*Px+Py*Py+Pz*Pz)/1000\n\n"
		"This command is not placed into the geometry.\n\n"
		"NOTE: if format begins 'Ptot=...' the parsing will think it is "
		"a named argument; put a space before the '=' to avoid that error."
	);

	z = "";
	zloop = "";
	require = "1";
	file = "-";
	format = "";
	noNewline = 0.0;
	coordinates = "Centerline";
	coordinateType = BLCOORD_CENTERLINE;
}

BLCMDprintf::BLCMDprintf(BLCMDprintf& r) : BLCommand(r), 
		BLManager::ZSteppingAction(r), expr(), eval()
{
	z = r.z;
	zloop = r.zloop;
	require = r.require;
	file = r.file;
	format = r.format;
	noNewline = r.noNewline;
	coordinates = r.coordinates;
	coordinateType = r.coordinateType;
}

int BLCMDprintf::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	BLCMDprintf *p = new BLCMDprintf(defaultPrintf);

	int retval = p->handleNamedArgs(namedArgs);

	p->coordinateType = BLCoordinates::getCoordinateType(p->coordinates);

	p->format = argv[0];
	for(unsigned i=1; i<argv.size(); ++i)
		p->expr.push_back(argv[i]);

	// register the z position(s) in the string z
	std::vector<G4double> vect = getList(p->z,",");
	if(vect.size() > 0) {
		for(unsigned i=0; i<vect.size(); ++i) {
			BLManager::getObject()->registerZStep(vect[i],p,7);
		}
	} else if(p->z.size() > 0) {
		printError("printf: Syntax error in z");
	}

	// register the z position(s) in the string zloop
	vect = getList(p->zloop,",:");
	if(vect.size() == 3) {
		if(vect[2] < 1.0*mm || vect[0] > vect[1]) {
			printError("printf: invalid zloop '%s'",
							p->zloop.c_str());
		} else {
			while(vect[0] <= vect[1]) {
			    BLManager::getObject()->registerZStep(vect[0],p,7);
			    vect[0] += vect[2];
			}
		}
	} else if(p->zloop.size() > 0) {
		printError("printf: invalid zloop '%s'",p->zloop.c_str());
	}
	BLManager::getObject()->registerEventAction(p,false);

	p->print("");

	return retval;
}

void BLCMDprintf::defineNamedArgs()
{
	argString(z,"z","Comma-separated list of Z positions for printing (mm)");
	argString(zloop,"zloop","Loop in z, first:last:incr (mm)");
	argString(require,"require","logical expression for cutting (default=true)");
	argString(file,"file","Output filename (default=\"-\" => stdout)");
	argString(file,"filename","Synonym for file");
	argDouble(noNewline,"noNewline","set nonzero to omit final newline.");
	argString(coordinates,"coordinates","Coordinates: global, centerline, or reference (default=c).");
}

void BLCMDprintf::UserZSteppingAction(const G4Track *track)
{
	BLManagerState state = BLManager::getObject()->getState();
	if(coordinateType == BLCOORD_REFERENCE && state != BEAM) return;

	// Print Tune and Reference tracks only once.
	if(!BLMPI::isRank1() && state != BEAM) return;

	eval.setTrackVariables(track,coordinateType,"",true);

	eval.setVariable("tune",state==TUNE ? 1.0 : 0.0);
	eval.setVariable("reference",state==REFERENCE ? 1.0 : 0.0);
	eval.setVariable("beam",state==BEAM ? 1.0 : 0.0);

	if(eval.evaluate(require) != 0.0) {
		double a[MAX_EXPR];
		for(unsigned i=0; i<MAX_EXPR; ++i) {
			if(i<expr.size())
				a[i] = eval.evaluate(expr[i]);
			else
				a[i] = 0.0;
		}
		BLMPI::rank0BufferPrintf(file,format.c_str(),
			a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
			a[10],a[11],a[12],a[13],a[14],a[15]);
		if(noNewline == 0.0) BLMPI::rank0BufferPrintf(file,"\n");
	}
}
