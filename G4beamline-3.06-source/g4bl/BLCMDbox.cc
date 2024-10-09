//	BLCMDbox.cc
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

#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLGroupElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLKillTrack.hh"

/**	BLCMDbox implements a box of material.
 *
 **/
class BLCMDbox : public BLGroupElement {
	G4double height;
	G4double width;
	G4double length;
	G4String material;
	G4String color;
	G4int kill;
	G4double maxStep;
	G4Box *box;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDbox();

	/// Destructor.
	virtual ~BLCMDbox() { }

	/// Copy constructor.
	BLCMDbox(const BLCMDbox& r);

	/// clone()
	BLElement *clone() { return new BLCMDbox(*this); }

	/// commandName() returns "box".
	G4String commandName() { return "box"; }

	/// command() implements the box command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// getSolid() returns the solid.
	virtual G4VSolid *getSolid();

	/// getMaterialName() returns the material name.
	virtual G4String getMaterialName() const { return material; }

	/// construct() - construct the box.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the tube/cylinder.
	G4double getLength() { return length; }

	/// getWidth() returns the outer radius of the tube/cylinder.
	G4double getWidth() { return width; }

	/// getHeight() returns the outer radius of the tube/cylinder.
	G4double getHeight() { return height; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// isOutside() from BLElement.
	bool isOutside(G4ThreeVector &local, G4double tolerance);

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isWithin() from BLGroupElement.
	bool isWithin(G4ThreeVector &local, G4double tolerance);
};

BLCMDbox defaultBox;	// default object

// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDbox::BLCMDbox() : BLGroupElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a box.");
	setDescription("This is a direct interface to G4Box.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.\n\n"
	);

	// provide initial values for fields
	height = 0.0;
	width = 0.0;
	length = 0.0;
	material = "";
	color = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	box = 0;
}

// Copy constructor - be sure to use the copy constructor BLGroupElement(r)
BLCMDbox::BLCMDbox(const BLCMDbox& r) : BLGroupElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	height = r.height;
	width = r.width;
	length = r.length;
	material = r.material;
	color = r.color;
	kill = r.kill;
	maxStep = r.maxStep;
	box = r.box;
}

int BLCMDbox::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("box: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultBox.handleNamedArgs(namedArgs);
	}

	BLCMDbox *t = new BLCMDbox(defaultBox);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	return retval;
}

void BLCMDbox::defineNamedArgs()
{
	argDouble(height,"height","The height of the box (mm).");
	argDouble(width,"width","The width of the box (mm).");
	argDouble(length,"length","The length of the box (mm).");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm).");
	argString(material,"material","The material of the box.");
	argString(color,"color","The color of the box (''=invisible).");
	argInt(kill,"kill","Set nonzero to kill every track that enters.");
}

G4VSolid *BLCMDbox::getSolid()
{
	if(!box) box = new G4Box("Box", width/2.0, height/2.0, length/2.0);
	return box;
}

void BLCMDbox::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	getSolid(); // sets box

	G4Material *mat;
	if(material != "")
		mat = getMaterial(material);
	else
		mat = parent->GetMaterial();

	G4LogicalVolume *lv = new G4LogicalVolume(box,mat, thisname+"LogVol");
	lv->SetVisAttributes(getVisAttrib(color));
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");
	lv->SetUserLimits(new G4UserLimits(maxStep));

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4VPhysicalVolume *pv = new G4PVPlacement(g4rot,
					relativePosition,lv,thisname,
					parent,false,0,surfaceCheck);

	// handle kill parameter
	if(kill)
		BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));

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

	constructChildren(lv,thisname,globalRotation,globalPosition);

	printf("BLCMDbox::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2],
		globalPosition[2]-getLength()/2.0,
		globalPosition[2]+getLength()/2.0);
}

G4bool BLCMDbox::isOutside(G4ThreeVector &local, G4double tolerance)
{
	return fabs(local[0]) > width/2.0-tolerance ||
		fabs(local[1]) > height/2.0-tolerance ||
		fabs(local[2]) > length/2.0-tolerance;
}

void BLCMDbox::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	generateBox(npoints,width,height,length,v);
}

G4bool BLCMDbox::isWithin(G4ThreeVector &local, G4double tolerance)
{
	return fabs(local[0]) < width/2.0+tolerance &&
		fabs(local[1]) < height/2.0+tolerance &&
		fabs(local[2]) < length/2.0+tolerance;
}
