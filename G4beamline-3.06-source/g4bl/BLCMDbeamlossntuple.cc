//	BLCMDbeamlossntuple.cc
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

#include "BLTrackNTuple.hh"
#include "BLManager.hh"
#include "BLCommand.hh"
#include "BLCoordinates.hh"

const char TrackFields[] =
    "x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight";
const unsigned NTrackFields = 12;

/**	class BLCMDbeamlossntuple - NTuple for tracks where they are lost.
 *	
 **/
class BLCMDbeamlossntuple : public BLCommand, public BLManager::TrackingAction {
	G4String name;
	G4String format;
	G4String filename;
	G4String require;
	G4String coordinates;
	BLCoordinateType coordinateType;
	BLTrackNTuple *ntuple;
public:
	BLCMDbeamlossntuple();

	BLCMDbeamlossntuple(BLCMDbeamlossntuple& r);

	G4String commandName() { return "beamlossntuple"; }

	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	/// help() prints help text.
	void help(bool detailed) {
		if(description[description.size()-2] == ':')
			description += BLTrackNTuple::getFormatList(); 
		BLCommand::help(detailed);
	}

	/// PreUserTrackingAction() from BLManager::TrackingAction.
	void PreUserTrackingAction(const G4Track *track) { }

	/// PostUserTrackingAction() from BLManager::TrackingAction.
	void PostUserTrackingAction(const G4Track *track);
};

BLCMDbeamlossntuple defaultBeamLossNTuple;


BLCMDbeamlossntuple::BLCMDbeamlossntuple() : BLCommand(), BLManager::TrackingAction()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("NTuple containing particle tracks when lost.");
	setDescription("All possible loss mechanisms are included. Tune "
		"particle tracks are omitted, but Reference "
		"track(s) will appear if one or more Reference particles are "
		"tracked.\n\n"
		"This command is not placed into the geometry.\n\n"
		"The NTuple contains the usual track data:\n\n"
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

	name = "";
	format = "";
	filename = "";
	require = "";
	coordinates = "Centerline";
	coordinateType = BLCOORD_CENTERLINE;
	ntuple = 0;
}

BLCMDbeamlossntuple::BLCMDbeamlossntuple(BLCMDbeamlossntuple& r) : BLCommand(r), 
						BLManager::TrackingAction(r)
{
	name = r.name;
	format = r.format;
	filename = r.filename;
	require = r.require;
	coordinates = r.coordinates;
	coordinateType = r.coordinateType;
	ntuple = 0;
}

int BLCMDbeamlossntuple::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("Invalid beamlossntuple command -- need name");
		return -1;
	}

	if(argv[0] == "default") {
		return handleNamedArgs(namedArgs);
	}

	BLCMDbeamlossntuple *t = new BLCMDbeamlossntuple(defaultBeamLossNTuple);
	int retval = t->handleNamedArgs(namedArgs);

	// ascii->bltrackfile format, for accuracy and consistency of output
	for(unsigned i=0; i<t->format.size(); ++i)
		t->format[i] = tolower(t->format[i]);
	if(t->format == "ascii")
		t->format = "bltrackfile";

	t->name = argv[0];
	t->coordinateType = BLCoordinates::getCoordinateType(t->coordinates);

	t->ntuple = BLTrackNTuple::create(t->format,"NTuple",t->name,
				t->filename,t->coordinateType,t->require);

	t->print(argv[0]);

	BLManager::getObject()->registerTrackingAction(t);

	return retval;
}

void BLCMDbeamlossntuple::defineNamedArgs()
{
	argString(format,"format","The NTuple format (see above for list).");
	argString(filename,"filename","The filename for the NTuple.");
	argString(filename,"file","Synonym for filename.");
	argString(require,"require","Expression which must be nonzero to include the track (default=1)",false);
	argString(coordinates,"coordinates","Coordinates: global, centerline, or reference (default=c).");
}

void BLCMDbeamlossntuple::PostUserTrackingAction(const G4Track *track)
{
	// omit tracks that were not tracked or were suspended
	if(track->GetCurrentStepNumber() <= 0 || 
					track->GetTrackStatus() == fSuspend)
		return;

	// only use reference coordinates when they are valid
	BLManagerState state = BLManager::getObject()->getState();
	if(coordinateType == BLCOORD_REFERENCE && state != BEAM) return;

	// omit Tune particle(s)
	if(state == TUNE) return;

	ntuple->appendTrack(track);
}

