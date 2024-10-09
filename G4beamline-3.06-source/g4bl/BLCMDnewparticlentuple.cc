//	BLCMDnewparticlentuple.cc
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

/**	class BLCMDnewparticlentuple - NTuple for tracks where they are created.
 *	
 **/
class BLCMDnewparticlentuple : public BLCommand,
					public BLManager::TrackingAction {
	G4String name;
	G4String format;
	G4String filename;
	G4String require;
	G4String coordinates;
	G4int kill;
	BLCoordinateType coordinateType;
	BLTrackNTuple *ntuple;
public:
	BLCMDnewparticlentuple();

	BLCMDnewparticlentuple(BLCMDnewparticlentuple& r);

	G4String commandName() { return "newparticlentuple"; }

	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	/// help() prints help text.
	void help(bool detailed) {
		if(description[description.size()-2] == ':')
			description += BLTrackNTuple::getFormatList(); 
		BLCommand::help(detailed);
	}

	/// PreUserTrackingAction() from BLManager::TrackingAction.
	void PreUserTrackingAction(const G4Track *track);

	/// PostUserTrackingAction() from BLManager::TrackingAction.
	void PostUserTrackingAction(const G4Track *track) { }
};

BLCMDnewparticlentuple defaultNewParticleNTuple;


BLCMDnewparticlentuple::BLCMDnewparticlentuple() : BLCommand(),
						BLManager::TrackingAction()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("NTuple containing particle tracks when created.");
	setDescription("Note that initial beam particles are included, unless "
		"require is set to 'ParentID>0'.\n\n"
		"This command is not placed into the geometry.\n\n"
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
	kill = 0;
	coordinateType = BLCOORD_CENTERLINE;
	ntuple = 0;
}

BLCMDnewparticlentuple::BLCMDnewparticlentuple(BLCMDnewparticlentuple& r) : BLCommand(r), 
						BLManager::TrackingAction(r)
{
	name = r.name;
	format = r.format;
	filename = r.filename;
	require = r.require;
	coordinates = r.coordinates;
	kill = r.kill;
	coordinateType = r.coordinateType;
	ntuple = 0;
}

int BLCMDnewparticlentuple::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("Invalid newparticlentuple command -- need name");
		return -1;
	}

	if(argv[0] == "default") {
		return handleNamedArgs(namedArgs);
	}

	BLCMDnewparticlentuple *t = new BLCMDnewparticlentuple(defaultNewParticleNTuple);
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

void BLCMDnewparticlentuple::defineNamedArgs()
{
	argString(format,"format","The NTuple format (see above for list).");
	argString(filename,"filename","The filename for the NTuple.");
	argString(filename,"file","Synonym for filename.");
	argString(require,"require","Expression which must be nonzero to include the track (default=1)",false);
	argString(coordinates,"coordinates","Coordinates: global, centerline, or reference (default=c).");
	argInt(kill,"kill","Set nonzero to kill tracks after entering into NTuple; does not kill track if require fails (0).");
}

void BLCMDnewparticlentuple::PreUserTrackingAction(const G4Track *track)
{
	// only use reference coordinates when they are valid
	BLManagerState state = BLManager::getObject()->getState();
	if(coordinateType == BLCOORD_REFERENCE && state != BEAM) return;

	if(ntuple->appendTrack(track)) { 
		if(kill) {
			static BLManager *manager = BLManager::getObject();
			if(manager->getSteppingVerbose() > 0)
			    printf("Track killed by newparticlentuple '%s'\n",
								name.c_str());
			((G4Track*)track)->SetTrackStatus(fStopAndKill);
		}
	}
}

