//	BLCMDzntuple.cc
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

#include "BLAssert.hh"
#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLParam.hh"
#include "BLCoordinates.hh"
#include "BLTrackNTuple.hh"
#include "BLGlobalField.hh"
#include "BLEvaluator.hh"

#ifdef WIN32
#define snprintf _snprintf
#endif
#ifdef __CYGWIN__
#include "mysnprintf.hh"
#endif


const char TrackFields[] =
    "x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight";
const unsigned NTrackFields = 12;

/**	class BLCMDzntuple will generate an NTuple for each of a list of Z
 *	positions.
 **/
class BLCMDzntuple : public BLCommand, public BLCallback {
	struct Entry : public BLManager::ZSteppingAction {
		BLCMDzntuple *zntuple;
		BLTrackNTuple *ntuple;
		BLCoordinateType coordinateType;
		Entry(BLCMDzntuple *znt, G4double zpos, 
					BLCoordinateType _coordinateType);
		void UserZSteppingAction(const G4Track *track);
	};
	G4String z;
	G4String zloop;
	G4int noSingles;
	G4String format;
	G4String file;
	G4String require;
	G4int referenceParticle;
	G4String coordinates;
	BLCoordinateType coordinateType;
public:
	/// Constructor.
	BLCMDzntuple();
	BLCMDzntuple(BLCMDzntuple& r);

	/// commandName() returns "zntuple".
	virtual G4String commandName() { return "zntuple"; }
	
	/// command() implements the zntuple command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// help() prints help text.
	void help(bool detailed) {
		if(description[description.size()-2] == ':')
			description += BLTrackNTuple::getFormatList(); 
		BLCommand::help(detailed);
	}

	/// callback from BLCallback.
	virtual void callback(int type);
};

BLCMDzntuple defaultZntuple;

BLCMDzntuple::BLCMDzntuple() : BLCommand(), BLCallback() 
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Generate an NTuple for each of a list of Z positions.");
	setDescription("Generates an NTuple like a virtualdetector, but\n"
		"without a physical volume. Tracks are forced to \n"
		"take steps within 2mm surrounding each desired z\n"
		"position, and they are interpolated to the desired\n"
		"z position. Each z position generates a separate\n"
		"NTuple named Z123 (etc.). z accepts a list of z\n"
		"positions, and zloop can generate a set of equally\n"
		"spaced z positions; both can be used.\n\n"
		"Each value in z and zloop can be an expression using "
		"double constants and the usual C operators and "
		"functions.\n\n"
		"This command places itself into the geometry.\n\n"
		"The standard NTuple fields are:\n"
		"    x,y,z (mm)\n"
		"    Px,Py,Pz (MeV/c)\n"
		"    t (ns)\n"
		"    PDGid (11=e-, 13=mu-, 22=gamma, 211=pi+, 2212=proton, ...)\n"
		"    EventID (may be inexact above 16,777,215)\n"
		"    TrackID\n"
		"    ParentID (0 => primary particle)\n"
		"    Weight (defaults to 1.0)\n\n"
		"The following additional fields are appended for "
		"format=Extended, format=asciiExtended, and "
		"format=rootExtended:\n"
		"    Bx, By, Bz (Tesla)\n"
		"    Ex, Ey, Ez (Megavolts/meter)\n"
		"    ProperTime (ns)\n"
		"    PathLength (mm)\n"
		"    PolX, PolY, PolZ (polarization)\n"
		"    InitX, initY, InitZ (initial position, mm)\n"
		"    InitT (initial time, ns)\n"
		"    InitKE (MeV when track was created)\n\n"
		"Valid Formats (ignore case): ");

	z = "";
	zloop = "";
	noSingles = 0;
	format = "";
	file = "";
	require = "";
	referenceParticle = 0;
	coordinates = "Centerline";
	coordinateType = BLCOORD_CENTERLINE;
}

BLCMDzntuple::BLCMDzntuple(BLCMDzntuple& r) : BLCommand(r), BLCallback(r)
{
	z = r.z;
	zloop = r.zloop;
	noSingles = r.noSingles;
	format = r.format;
	file = r.file;
	require = r.require;
	referenceParticle = r.referenceParticle;
	coordinates = r.coordinates;
	coordinateType = r.coordinateType;
}

int BLCMDzntuple::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() >= 1 && argv[0] == "default") {
		return handleNamedArgs(namedArgs);
	}

	BLCMDzntuple *p = new BLCMDzntuple(defaultZntuple);

	int retval = p->handleNamedArgs(namedArgs);

	p->coordinateType = BLCoordinates::getCoordinateType(p->coordinates);

	p->print("");

	// must defer creation of BLTrackNTuple-s
	BLManager::getObject()->registerCallback(p,0);

	// ascii->bltrackfile format, for accuracy and consistency of output
	for(unsigned i=0; i<p->format.size(); ++i)
		p->format[i] = tolower(p->format[i]);
	if(p->format == "ascii")
		p->format = "bltrackfile";

	return retval;
}

void BLCMDzntuple::defineNamedArgs()
{
	argString(z,"z","Comma-separated list of Z positions (mm)");
	argString(zloop,"zloop","Loop in z, first:last:incr (mm)");
	argInt(noSingles,"noSingles","Set to 1 to omit the NTuple for singles.");
	argString(format,"format","NTuple format (see above for list)");
	argString(file,"file","Output filename ('' uses name to determine filename)");
	argString(file,"filename","Synonym for file");
	argString(require,"require","Expression which must be nonzero to include the track (default=1)",false);
	argInt(referenceParticle,"referenceParticle","Set to 1 to include the Reference Particle.");
	argString(coordinates,"coordinates","Coordinates: global, centerline, or reference (default=c).");
}

void BLCMDzntuple::callback(int type)
{
	if(type != 0) return;

	// register the z position(s) in the string z
	std::vector<G4double> vect = getList(z,",");
	if(vect.size() > 0) {
		for(unsigned i=0; i<vect.size(); ++i) {
			Entry *e = new Entry(this,vect[i],coordinateType);
			int when = 4;
			if(referenceParticle != 0 || 
			   (e->ntuple && e->ntuple->needsReference()))
				when |= 2;
			BLManager::getObject()->registerZStep(vect[i],e,when);
		}
	} else if(z.size() > 0) {
		printError("printf: Syntax error in z");
	}

	// register the z position(s) in the string zloop
	vect = getList(zloop,",:");
	if(vect.size() == 3) {
		if(vect[2] < 1.0*mm || vect[0] > vect[1]) {
			printError("printf: invalid zloop '%s'",zloop.c_str());
		} else {
			while(vect[0] <= vect[1]) {
			    Entry *e = new Entry(this,vect[0],coordinateType);
			    int when = 4;
			    if(referenceParticle != 0 || 
			       (e->ntuple && e->ntuple->needsReference()))
			       	when |= 2;
			    BLManager::getObject()->registerZStep(vect[0],
								e,when);
			    vect[0] += vect[2];
			}
		}
	} else if(zloop.size() > 0) {
		printError("printf: invalid zloop '%s'",zloop.c_str());
	}
}

BLCMDzntuple::Entry::Entry(BLCMDzntuple *znt, G4double zpos, 
					BLCoordinateType _coordinateType)
{
	char name[128];
	snprintf(name,sizeof(name),"Z%.0f",floor(zpos+0.5));

	zntuple = znt;

	coordinateType = _coordinateType;

	ntuple = BLTrackNTuple::create(zntuple->format,"NTuple",name,
			zntuple->file,coordinateType,zntuple->require,
			zntuple->noSingles);
}

void BLCMDzntuple::Entry::UserZSteppingAction(const G4Track *track)
{
	// only use reference coordinates when they are valid
	BLManagerState state = BLManager::getObject()->getState();
	if(coordinateType == BLCOORD_REFERENCE && state != BEAM) return;

	ntuple->appendTrack(track);
}
