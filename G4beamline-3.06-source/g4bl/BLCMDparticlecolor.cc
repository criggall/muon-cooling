//	BLCMDparticlecolor.cc
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
//

#include <map>

#include "G4Track.hh"
#include "G4UImanager.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "BLAssert.hh"
#include "BLCommand.hh"
#include "BLCommandAlias.hh"
#include "BLManager.hh"
#include "BLMarkers.hh"

/**	class BLCMDparticlecolor -- set colors to specific particles.
 *	(Visualization only.)
 *
 *	The reference trajectory(s) are displayed as BLMarkers in line mode.
 *	This makes them both persistent and displayed whenever the screen
 *	needs to be refreshed. Their points are collected in TrackingAction
 *	and SteppingAction callbacks.
 **/
class BLCMDparticlecolor : public BLCommand, public BLManager::SteppingAction,
			public BLManager::TrackingAction, public BLCallback {
	BLMarkers *referenceTrajectory;
	std::map<G4String,G4String> name2color;
	G4String plus;
	G4String minus;
	G4String neutral;
	G4String referenceColor;
public:
	/// Constructor.
	BLCMDparticlecolor();

	/// Destructor.
	~BLCMDparticlecolor() { }

	/// commandName() returns "particlecolor"
	G4String commandName()  { return "particlecolor"; }

	/// command() executes the command associated with this element.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command. 
	void defineNamedArgs() { }

	/// UserSteppingAction() from BLManager::SteppingAction
	void UserSteppingAction(const G4Step *step);

	/// from BLManager::TrackingAction
	void PreUserTrackingAction(const G4Track *track);
	void PostUserTrackingAction(const G4Track *track);

	/// callback() from BLCallback.
	void callback(int type);
};
BLCMDparticlecolor defaultParticleColor;

BLCommandAlias aliasParticleColor("trackcolor",defaultParticleColor);

BLCMDparticlecolor::BLCMDparticlecolor() : BLCommand(), 
		BLManager::SteppingAction(), BLManager::TrackingAction(),
		BLCallback(), name2color()
{
	registerCommand(BLCMDTYPE_AUX);
	setSynopsis("Set the colors for tracks by particle name.");
	setDescription("Arguments are of the form 'name=1,1,0', where name is\n"
		"the standard name of a particle, and 1,1,0 is the\n"
		"R,G,B value desired for its color ('' or 'i' for invisible).\n"
		"The special names plus, minus, and neutral will set the "
		"colors for unnamed particles by their charge (default to "
		"blue, red, green).\n"
		"The name reference will apply to the reference\n"
		"track(s) (defaults to invisible).");
	
	referenceTrajectory = 0;
	plus = "0 0 1 1";
	minus = "1 0 0 1";
	neutral = "0 1 0 1";
	referenceColor = "";
}

int BLCMDparticlecolor::command(BLArgumentVector& argv,BLArgumentMap& namedArgs)
{
	int retval = 0;

	if(argv.size() > 0) {
		printError("particlecolor: unnamed arguments are invalid.");
		retval = 1;
	}

	BLArgumentMap::iterator i;
	for(i=namedArgs.begin(); i!=namedArgs.end(); ++i) {
		G4String name = i->first;
		G4String value = i->second; // RGB with commas
		G4String valueRGBA;         // RGBA with spaces
		// Windows VC 2015 has a compiler bug in value[0], use substr().
		if(value == "" || value.substr(1) == "i" ||
						value.substr(1) == "I") {
			valueRGBA = "0 0 0 0"; // invisible
		} else {
			std::vector<G4String> list=splitString(value,',',true);
			if(list.size() != 3) {
				printError("particlecolor: invalid argument "
					"%s='%s'\n",name.c_str(),value.c_str());
				retval = -1;
				continue;
			}
			valueRGBA = list[0]+" "+list[1]+" "+list[2]+" 1";
		}
		if(name == "plus") {
			plus = valueRGBA;
		} else if(name == "minus") {
			minus = valueRGBA;
		} else if(name == "neutral") {
			neutral = valueRGBA;
		} else if(name == "reference") {
			if(valueRGBA != "0 0 0 0") {
				referenceTrajectory = new BLMarkers(value);
				referenceTrajectory->setLine();
				referenceColor = value;
			}
		} else {
			name2color[name] = valueRGBA;
		}
	}

	if(referenceTrajectory != 0) {
		BLManager::getObject()->registerReferenceParticleStep(0,this);
		BLManager::getObject()->registerTrackingAction(this);
	}
	BLManager::getObject()->registerCallback(this,4);

	print("",namedArgs);

	return retval;
}

void BLCMDparticlecolor::UserSteppingAction(const G4Step *step)
{
	if(BLManager::getObject()->getState() != REFERENCE) return;
	if(!referenceTrajectory) return;
	const G4Track *track = step->GetTrack();
	referenceTrajectory->addMarker(track->GetPosition());
}

void BLCMDparticlecolor::PreUserTrackingAction(const G4Track *track)
{
	if(BLManager::getObject()->getState() != REFERENCE) return;
	if(!referenceTrajectory) return;
	if(referenceTrajectory->size() > 0) {
		// keep the old one, so it can be drawn
		referenceTrajectory = new BLMarkers(referenceColor);
		referenceTrajectory->setLine();
	}
	referenceTrajectory->addMarker(track->GetPosition());
}

void BLCMDparticlecolor::PostUserTrackingAction(const G4Track *track)
{
}

void BLCMDparticlecolor::callback(int type)
{
	if(type != 4) return;
	G4UImanager *UI = G4UImanager::GetUIpointer();
	BLAssert(UI != 0);

	UI->ApplyCommand("/vis/modeling/trajectories/create/drawByParticleID "
									"PID");
	G4ParticleTable::G4PTblDicIterator *it =
			G4ParticleTable::GetParticleTable()->GetIterator();
	it->reset();
	while((*it)()) {
		G4ParticleDefinition *pd = it->value();
		G4String name = pd->GetParticleName();
		G4String c;
		if(name2color.count(name) > 0) {
			c = name2color[name];
		} else if(pd->GetPDGCharge() > 0.0) {
			c = plus;
		} else if(pd->GetPDGCharge() < 0.0) {
			c = minus;
		} else {
			c = neutral;
		}
		c = "/vis/modeling/trajectories/PID/setRGBA " + name + " " + c;
		UI->ApplyCommand(c);
	}
	// set default so it will apply to ions
	UI->ApplyCommand("/vis/modeling/trajectories/PID/setDefaultRGBA "+plus);
}

