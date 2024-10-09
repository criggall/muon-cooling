//	BLCMDfieldmap.cc
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

#include <fstream>
#include <string.h>
#include <ctype.h>

#include "G4ElectroMagneticField.hh"
#include "BLFieldMap.hh"
#include "BLCommand.hh"		// command-parsing routines only
#include "BTSpline1D.hh"

#include "BLElement.hh"
#include "BLGlobalField.hh"
#include "BLElementField.hh"
#include "BLTune.hh"

/**	class BLCMDfieldmap implements a FieldMap command
 *
 *	This command has no volumes, logical or physical. it just implements
 *	a field (E and/or B) from an input file based on BLFieldMap.
 **/
class BLCMDfieldmap : public BLElement {
	BLFieldMap *map;
	G4String filename;
	G4double current;
	G4double gradient;
	G4double timeOffset;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDfieldmap();

	/// Destructor.
	virtual ~BLCMDfieldmap() { }

	/// Copy constructor.
	BLCMDfieldmap(const BLCMDfieldmap& r);

	/// clone()
	BLElement *clone() { return new BLCMDfieldmap(*this); }

	/// commandName() returns "fieldmap".
	G4String commandName() { return "fieldmap"; }

	/// command() implements the fieldmap command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() - construct the fieldmap
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns 0 (no physical volume nere)
	G4double getLength() { return 0.0; }

	/// getWidth() returns 0 (no physical volume nere)
	G4double getWidth() { return 0.0; }

	/// getHeight() returns 0 (no physical volume nere)
	G4double getHeight() { return 0.0; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// isOutside() from BLElement. (no volume => every point is "outside")
	bool isOutside(G4ThreeVector &local, G4double tolerance) 
		{ return true; }

	/// generatePoints() from BLElement. (no volume => no generate)
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v) 
		{ v.clear(); }
};

BLCMDfieldmap defaultFieldMapCmd;	// default object

class FieldMapPlacement : public BLElementField {
	BLCoordinateTransform global2local;
	G4RotationMatrix rotation;
	G4double *current;
	G4double *gradient;
	G4double timeOffset;
	BLFieldMap *map;
public:
	FieldMapPlacement(BLFieldMap *_map, G4RotationMatrix *rot,
		G4ThreeVector &pos, G4double& _current, G4double& _gradient,
		G4double _timeOffset);
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};


// Default constructor - be sure to use the default constructor BLElement()
BLCMDfieldmap::BLCMDfieldmap() : BLElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("implements a field map, E and/or B, from a file.");
	setDescription("Reads an input file in BLFieldMap format to define "
		"E and/or B fields, optionally with time dependence. "
		"See the Users Guide for a description of the BLFieldMap "
		"format.\n\n"
		"This field must be placed (via the place command); that "
		"specifies where (x=0,y=0,z=0) is located in the parent. "
		"Note that front=1 can NOT be used to place this field.");

	// provide initial values for fields
	map = 0;
	current = 1.0;
	gradient = 1.0;
	timeOffset = 0.0;
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDfieldmap::BLCMDfieldmap(const BLCMDfieldmap& r) : BLElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	map = r.map;
	BLTune::copyTunableArg(&current,&r.current);
	BLTune::copyTunableArg(&gradient,&r.gradient);
	timeOffset = r.timeOffset;
}

int BLCMDfieldmap::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		BLCommand::printError("fieldmap: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultFieldMapCmd.handleNamedArgs(namedArgs);
	}

	BLCMDfieldmap *t = new BLCMDfieldmap(defaultFieldMapCmd);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	t->map = new BLFieldMap();
	t->map->readFile(t->filename);

	t->print(argv[0]);

	return retval;
}

void BLCMDfieldmap::defineNamedArgs()
{
	argString(filename,"filename","Filename for the Field Map.",false);
	argString(filename,"file","Synonym for filename.",false);
	argTunable(current,"current","Current of the B-field.");
	argTunable(gradient,"gradient","Gradient of the E-field.");
	argDouble(timeOffset,"timeOffset","Time offset (ns).");
}

void BLCMDfieldmap::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)

{
	G4String thisname = parentName+getName();

	// get globalRotation and globalPosition
	G4RotationMatrix *globalRotation = 0;
	if(relativeRotation && parentRotation) {
		globalRotation = 
		    new G4RotationMatrix(*parentRotation * *relativeRotation);
	} else if(relativeRotation) {
		globalRotation = relativeRotation;
	} else if(parentRotation) {
		globalRotation = parentRotation;
	}
	G4ThreeVector globalPosition(relativePosition + parentPosition);
	if(parentRotation)
		globalPosition = *parentRotation * relativePosition +
				parentPosition;

	FieldMapPlacement *p = new FieldMapPlacement(map,globalRotation,
				globalPosition,current,gradient,timeOffset);
	BLGlobalField::getObject()->addElementField(p);

	printf("BLMappedMagnet::Construct %s parent=%s relZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2]);
}

FieldMapPlacement::FieldMapPlacement(BLFieldMap *_map, G4RotationMatrix *rot,
		G4ThreeVector &pos, G4double& _current, G4double& _gradient,
		G4double _timeOffset)
			: global2local(rot,pos), rotation() 
{
	map = _map;
	current = &_current;
	gradient = &_gradient;
	timeOffset = _timeOffset;

	if(global2local.isRotated()) {
		rotation = global2local.getRotation();
		rotation.invert();
	}

	// set global bounding box
	G4double local[4], global[4];
	local[3] = 0.0;
	for(int i=0; i<8; ++i) {
		map->getBoundingPoint(i,local);
		global2local.getGlobal(local,global);
		setGlobalPoint(global);
	}
}

void FieldMapPlacement::addFieldValue(const G4double point[4],
						G4double field[6]) const
{
	G4double local[4], thisField[6];

	global2local.getLocal(local,point);

	local[3] -= timeOffset;  // KBB 25aug10 - shift time back by timeOffset

	map->getFieldValue(local,thisField,*current,*gradient);

	if(map->hasB()) {
		G4ThreeVector B(thisField[0],thisField[1],thisField[2]);
		B = rotation * B;
		field[0] += B[0];
		field[1] += B[1];
		field[2] += B[2];
	}

	if(map->hasE()) {
		G4ThreeVector E(thisField[3],thisField[4],thisField[5]);
		E = rotation * E;
		field[3] += E[0];
		field[4] += E[1];
		field[5] += E[2];
	}
}

