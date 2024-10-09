//	BLCMDgroup.cc
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
//	NOTE: Separating this into a command class and an infrastructure class
//	was not easy, so class BLGroup does both. In keeping with the other
//	classes, BLCMDgroup.cc contains the command functions, and BLGroup.cc
//	contains the infrastructure functions.

#include <algorithm>

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"
#include "BLGroup.hh"
#include "BLParam.hh"
#include "BLCoordinates.hh"

BLGroup defaultGroup;

BLGroup::BLGroup() : BLGroupElement(), offset(0,0,0)
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("begins definition of a group.");
	setDescription("A group is a collection of elements that can be placed\n"
		"together, preserving their relative positions. The group\n"
		"is a LogicalVolume in the geant4 geometry -- this means that "
		"a group cannot overlap any other group or object, even if the "
		"overlapping portion of the group is empty. If you need to "
		"permit overlaps, consider using a macro instead (the define "
		"command).\n\n"
		"If the group is\n"
		"given a length, then children can be placed at specific\n"
		"z offsets relative to the center of the group. If the group\n"
		"is not given a length, then children can only be placed\n"
		"sequentially along z, and the length will be computed by\n"
		"the endgroup command. Width and height are computed\n"
		"from the largest child or the argument.\n"
		"If radius is set to 0, the group will be a cylinder with\n"
		"radius determined by the largest width or height placed\n"
		"into it; if >0 then that is the fixed radius. If radius\n"
		"is not set then a box is used.\n\n"
		"This element must be placed (via the place command).\n\n"
		"Note: when placing objects into a group, if the rename "
		"argument is used, it should begin with a '+' to include the "
		"group's name in the object's name; otherwise there may be "
		"multiple objects with the same name -- this is only a major "
		"problem for virtualdetector-s and other output objects. The "
		"group's name is included by default if no rename is used on "
		"the place command.");

	prevCurrent = 0;
	highZ = 0.0;
	lastZ = 0.0;
	length = 0.0;
	width = 0.0;
	height = 0.0;
	radius = -1.0;
	material = "Vacuum";
	color = "";
	maxStep = -1.0;
	fixedLength = false;
	fixedWidth = false;
	fixedHeight = false;
}

BLGroup::~BLGroup()
{
}

int BLGroup::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("invalid group command -- requires name");
		argv.clear();
		argv.push_back("UnnamedGroup");
	}

	BLGroup *c = new BLGroup(defaultGroup);
	c->setName(argv[0]);
	c->handleNamedArgs(namedArgs);
	if(c->length > 0.0)
		c->lastZ = c->highZ = -c->length/2.0;

	if(c->maxStep < 0.0) c->maxStep = Param.getDouble("maxStep");

	if(c->radius > 0.0)
		c->width = c->height = 2.0*c->radius;
	c->fixedLength = (c->length > 0.0);
	c->fixedWidth = (c->width > 0.0);
	c->fixedHeight = (c->height > 0.0);

	c->prevCurrent = getCurrent();	// (current could be 0)
	current = c;

	// check material exists
	if(c->material.size() > 0) getMaterial(c->material);

	c->print(argv[0]);

	return 0;
}

void BLGroup::defineNamedArgs()
{
	argDouble(length,"length","Overall group length along z (mm)",mm,"",false);
	argDouble(width,"width","Overall group width (mm)",mm,"",false);
	argDouble(height,"height","Overall group height (mm)",mm,"",false);
	argDouble(radius,"radius","Radius for a cylindrical group (mm)",mm,"",false);
	argString(material,"material","Material of the volume outside children");
	argString(color,"color","Color of the volume of the group");
	argDouble(maxStep,"maxStep","The maximum stepsize in the volume (mm)");
}

void BLGroup::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition) 
{
	if(radius == 0.0) radius = std::max<double>(width,height)/2.0;
	bool silent = (getName() == "World");
	if(radius > 0.0 && !silent)
	    printf("BLGroup::Construct %s parent=%s relZ=%.1f globZ=%.1f length=%.1f\n"
			"\tzmin=%.1f zmax=%.1f radius=%.1f\n",
		getName().c_str(),parentName.c_str(),relativePosition[2],
		parentPosition[2],getLength(),
		parentPosition.z()-getLength()/2.0,
		parentPosition.z()+getLength()/2.0,
		radius);
	else if(!silent)
	    printf("BLGroup::Construct %s parent=%s relZ=%.1f globZ=%.1f length=%.1f\n"
			"\tzmin=%.1f zmax=%.1f width=%.1f Height=%.1f\n",
		getName().c_str(),parentName.c_str(),relativePosition[2],
		parentPosition[2],getLength(),
		parentPosition.z()-getLength()/2.0,
		parentPosition.z()+getLength()/2.0,
		getWidth(),getHeight());
	
	if(maxStep <= 0.0) maxStep = Param.getDouble("maxStep");

	// add safety region to the world volume, so steps are guaranteed
	// to remain inside the world volume. Double it to be sure.
	G4double safety = 0.0;
	if(this == world) safety = maxStep * 2.0;

	G4VSolid *solid=0;
	if(radius > 0.0)
		solid = new G4Tubs(getName()+"Tubs", 0.0, radius+safety,
			length/2.0+safety, 0.0, 2.0*pi);
	else
		solid = new G4Box(getName()+"Box",getWidth()/2.0+safety,
				getHeight()/2.0+safety,getLength()/2.0+safety);
	G4Material *mat = getMaterial(material);
	G4LogicalVolume *lv = new G4LogicalVolume(solid,mat,getName()+"LogVol");
	lv->SetVisAttributes(getVisAttrib(color));
	lv->SetUserLimits(new G4UserLimits(maxStep));

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4VPhysicalVolume *phys = new G4PVPlacement(g4rot, relativePosition,lv,
			parentName+getName(),parent, false,0,surfaceCheck);

	G4String xname = parentName + getName();
	if(this == world) {
		worldPhysVol = phys;
		xname = "";	// omit "World" from the start of every name
	}

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

	if(!silent) printf("\tParent pos=%.1f,%.1f,%.1f Relative pos=%.1f,%.1f,%.1f Global pos=%.1f,%.1f,%.1f\n",
		parentPosition[0],parentPosition[1],parentPosition[2],
		relativePosition[0],relativePosition[1],relativePosition[2],
		globalPosition[0],globalPosition[1],globalPosition[2]);

	constructChildren(lv,xname,globalRotation,globalPosition);
}

