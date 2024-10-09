//	BLCMDshowmaterial.cc
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

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4VisAttributes.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"
#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "BLManager.hh"
#include "BLCommand.hh"

/**	class BLCMDshowmaterial - command to display a set of materials
 *
 *	Each argument name must be a material name, and name=1,1,0 assigns
 *	a color to that material, overriding any color assigned by the object.
 *	hideOthers=1 causes all other materials to be made invisible.
 *
 *	Applies to ALL logical volumes with each material.
 *
 *	Useful to display Vacuum regions, which are usually invisible and
 *	enclosed in pipes.
 **/
class BLCMDshowmaterial : public BLCommand, public BLCallback {
	static std::map<const G4Material*,const G4VisAttributes*> mat2va;
	int hideOthers;
public:
	/// Constructor.
	BLCMDshowmaterial();

	/// Destructor.
	~BLCMDshowmaterial() { }

	/// commandName() returns "showmaterial"
	G4String commandName()  { return "showmaterial"; }

	/// command() executes the command associated with this element.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// getVisAttributes()  will return the G4VisAttributes for a
	/// material.
	const G4VisAttributes *getVisAttributes(const G4Material *pm) const;

	/// callback() from BLCallback.
	void callback(int type);

	/// handlePV() recursively handles physical volume-s.
	void handlePV(G4VPhysicalVolume *phys);
};
std::map<const G4Material*,const G4VisAttributes*> BLCMDshowmaterial::mat2va;

BLCMDshowmaterial defaultShowMaterial;

BLCMDshowmaterial::BLCMDshowmaterial() : BLCommand()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("Set the colors for selected materials.");
	setDescription("Arguments are of the form 'name=1,1,0', where name is\n"
			"the name of a material, and 1,1,0 is the\n"
			"R,G,B value desired for its color ('' for invisible)\n"
			"Set hideOthers=1 to make all other materials invisible.\n"
			"BEWARE: 'Vacuum' and 'vacuum' are different materials,\n"
			"as are 'Iron' and 'Fe'.");
	hideOthers = 0;
}

int BLCMDshowmaterial::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	int retval = 0;

	if(argv.size() > 0) {
		printError("showmaterial: unnamed arguments are invalid.");
		retval = 1;
	}

	BLArgumentMap::iterator i;
	for(i=namedArgs.begin(); i!=namedArgs.end(); ++i) {
		G4String name = i->first;
		G4String value = i->second;
		if(name == "hideOthers") {
			hideOthers = atoi(value.c_str());
			continue;
		}
		const G4VisAttributes *va = getVisAttrib(value);
		G4Material *pm = getMaterial(name,true);
		if(!pm) {
			printError("showmaterial: material '%s' not found",
								name.c_str());
			retval = 1;
			continue;
		}
		mat2va[pm] = va;
	}

	BLManager::getObject()->registerCallback(this, 0);

	print("",namedArgs);

	return retval;
}

const G4VisAttributes *BLCMDshowmaterial::getVisAttributes(const G4Material *pm) const
{
	if(mat2va.count(pm) > 0)
		return mat2va[pm];
	return 0;
}

void BLCMDshowmaterial::callback(int type)
{
	G4VPhysicalVolume *phys = 
			BLManager::getObject()->getWorldPhysicalVolume();
	handlePV(phys);
}

void BLCMDshowmaterial::handlePV(G4VPhysicalVolume *phys)
{
	G4LogicalVolume *log = phys->GetLogicalVolume();
	const G4VisAttributes *va = getVisAttributes(log->GetMaterial());
	if(!va && hideOthers) va = BLCommand::getVisAttrib("Invisible");
	log->SetVisAttributes(va);

	int n = log->GetNoDaughters();
	for(int i=0; i<n; ++i) {
		G4VPhysicalVolume *p = log->GetDaughter(i);
		if(p)
			handlePV(p);
	}
}
