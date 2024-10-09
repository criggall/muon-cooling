//	BLCMDtrap.cc
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
#include "G4Trap.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLAssert.hh"
#include "BLGroupElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLKillTrack.hh"

const G4double UNINITIALIZED = -9.99e29;

/**	BLCMDtrap implements a solid trapezoid of material, axis along Z.
 *
 **/
class BLCMDtrap : public BLGroupElement {
	G4double height;
	G4double upperWidth;
	G4double lowerWidth;
	G4double Xul, Xur, Xll, Xlr;
	G4double length;
	G4String material;
	G4String color;
	G4int kill;
	G4double maxStep;
	G4Trap *trap;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDtrap();

	/// Destructor.
	virtual ~BLCMDtrap() { }

	/// Copy constructor.
	BLCMDtrap(const BLCMDtrap& r);

	/// clone()
	BLElement *clone() { return new BLCMDtrap(*this); }

	/// commandName() returns "trap".
	G4String commandName() { return "trap"; }

	/// command() implements the trap command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	virtual G4VSolid *getSolid();

	/// construct() - construct the trapezoid
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the trapezoid
	G4double getLength() { return length; }

	/// getWidth() returns the width of the trapezoid
	G4double getWidth() { 
		return (upperWidth>lowerWidth ? upperWidth : lowerWidth); 
	}

	/// getHeight() returns the height of the tube/cylinder.
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
	bool isOutside(G4ThreeVector &local, G4double tolerance) 
		{ BLAssert(trap!=0);  return trap->Inside(local) != kInside; }

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v)
		{ v.clear();
		  for(int n=0; n<npoints*5; ++n)
		  	v.push_back(trap->GetPointOnSurface());
		}

	/// isWithin() from BLGroupElement.
	bool isWithin(G4ThreeVector &local, G4double tolerance) 
	    { BLAssert(trap!=0);  return trap->Inside(local) != kOutside; }
};

BLCMDtrap defaultTrap;	// default object

// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDtrap::BLCMDtrap() : BLGroupElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a solid trapezoid with axis along z.");
	setDescription("This is a direct interface to G4Trap.\n"
		"The trapezoid is symmetrical left-right, but upper or\n"
		"lower width can be larger or smaller.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.");

	// provide initial values for fields
	height = 0.0;
	upperWidth = UNINITIALIZED;
	lowerWidth = UNINITIALIZED;
	Xul = Xur = Xll = Xlr = 0.0;
	length = 0.0;
	material = "";
	color = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	trap = 0;
}

// Copy constructor - be sure to use the copy constructor BLGroupElement(r)
BLCMDtrap::BLCMDtrap(const BLCMDtrap& r) : BLGroupElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	height = r.height;
	upperWidth = r.upperWidth;
	lowerWidth = r.lowerWidth;
	Xul = r.Xul;
	Xur = r.Xur;
	Xll = r.Xll;
	Xlr = r.Xlr;
	length = r.length;
	material = r.material;
	color = r.color;
	kill = r.kill;
	maxStep = r.maxStep;
	trap = r.trap;
}

int BLCMDtrap::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("trap: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultTrap.handleNamedArgs(namedArgs);
	}

	BLCMDtrap *t = new BLCMDtrap(defaultTrap);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	return retval;
}

void BLCMDtrap::defineNamedArgs()
{
	argDouble(height,"height","The height of the solid trapezoid (mm)");
	argDouble(upperWidth,"upperWidth","The upper width solid trapezoid (mm)");
	argDouble(lowerWidth,"lowerWidth","The lowerWidth of the solid trapezoid (mm)");
	argDouble(Xul,"Xul","X position of upper left corner (mm)");
	argDouble(Xur,"Xur","X position of upper right corner (mm)");
	argDouble(Xll,"Xll","X position of lower left corner (mm)");
	argDouble(Xlr,"Xlr","X position of lower right corner (mm)");
	argDouble(length,"length","The length of the solid trapezoid (mm)");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)");
	argString(material,"material","The material of the trapezoid");
	argString(color,"color","The color of the trapezoid (''=invisible)");
	argInt(kill,"kill","Set nonzero to kill every track that enters.");
}

G4VSolid *BLCMDtrap::getSolid()
{
	if(lowerWidth != UNINITIALIZED) {
		Xll = -lowerWidth/2.0;
		Xlr = lowerWidth/2.0;
	}
	if(upperWidth != UNINITIALIZED) {
		Xul = -upperWidth/2.0;
		Xur = upperWidth/2.0;
	}

	if(!trap) {
		G4ThreeVector p[8];
		p[0][2] = p[1][2] = p[2][2] = p[3][2] = -length/2.0;
		p[4][2] = p[5][2] = p[6][2] = p[7][2] = length/2.0;
		p[0][1] = p[1][1] = p[4][1] = p[5][1] = -height/2.0;
		p[2][1] = p[3][1] = p[6][1] = p[7][1] = height/2.0;
		p[0][0] = p[4][0] = Xll;
		p[1][0] = p[5][0] = Xlr;
		p[2][0] = p[6][0] = Xul;
		p[3][0] = p[7][0] = Xur;
		trap = new G4Trap("Trap", p);
	}
	return trap;
}

void BLCMDtrap::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	getSolid(); // sets trap

	G4Material *mat;
	if(material != "")
		mat = getMaterial(material);
	else
		mat = parent->GetMaterial();

	G4LogicalVolume *lv = new G4LogicalVolume(trap,mat, thisname+"LogVol");
	lv->SetVisAttributes(getVisAttrib(color));
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");
	lv->SetUserLimits(new G4UserLimits(maxStep));

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4VPhysicalVolume *pv = new G4PVPlacement(g4rot, relativePosition,lv,
					thisname, parent,false,0,surfaceCheck);

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

	printf("BLCMDtrap::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f kill=%d\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2],
		globalPosition[2]-getLength()/2.0,
		globalPosition[2]+getLength()/2.0,
		kill);

}

