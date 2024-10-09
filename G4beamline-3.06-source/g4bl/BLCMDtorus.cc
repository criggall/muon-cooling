//	BLCMDtorus.cc
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

#include "G4VisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Torus.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4Circle.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLAssert.hh"
#include "BLGroupElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLKillTrack.hh"

/**	BLCMDtorus implements a torus of material.
 *
 **/
class BLCMDtorus : public BLGroupElement {
	G4double innerRadius;
	G4double outerRadius;
	G4double majorRadius;
	G4double initialPhi;
	G4double finalPhi;
	G4String material;
	G4String color;
	G4int kill;
	G4double maxStep;
	G4Torus *torus;
	friend class Surface;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDtorus();

	/// Destructor.
	virtual ~BLCMDtorus() { }

	/// Copy constructor.
	BLCMDtorus(const BLCMDtorus& r);

	/// clone()
	BLElement *clone() { return new BLCMDtorus(*this); }

	/// commandName() returns "torus".
	G4String commandName() { return "torus"; }

	/// command() implements the torus command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	G4VSolid *getSolid();

	/// construct() - construct the torus.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the torus, assuming
	/// full 2pi around the major radius.
	G4double getLength() { return outerRadius*2; }

	/// getWidth() returns the width of the torus, assuming
	/// full 2pi around the major radius.
	G4double getWidth() { return (majorRadius+outerRadius)*2; }

	/// getHeight() returns the height of the tube/cylinder.
	G4double getHeight() { return (majorRadius+outerRadius)*2; }

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
		{ BLAssert(torus!=0);  return torus->Inside(local) != kInside; }

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v)
		{ v.clear();
		  for(int n=0; n<npoints*5; ++n)
		  	v.push_back(torus->GetPointOnSurface());
		}

	/// isWithin() from BLGroupElement.
	bool isWithin(G4ThreeVector &local, G4double tolerance) 
		{ BLAssert(torus!=0);  return torus->Inside(local) != kOutside; }
};
BLCMDtorus defaultTorus;	// default object

/** class Surface will display points on the surface of a torus
 **/
class Surface : public BLCallback {
	BLCMDtorus *torus;
	G4RotationMatrix *rot;
	G4ThreeVector position;
public:
	Surface(BLCMDtorus *_torus, G4RotationMatrix *globalRotation,
					G4ThreeVector &globalPosition) 
		{ torus = _torus; rot = globalRotation; position = globalPosition; 
		  BLManager::getObject()->registerCallback(this,4);
		}

	void callback(int type);
};


// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDtorus::BLCMDtorus() : BLGroupElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a torus.");
	setDescription("This is a direct interface to G4Torus. "
		"The major radius is in the X-Y plane, with phi=0 "
		"along X.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.");

	// provide initial values for fields
	innerRadius = 0.0;
	outerRadius = 0.0;
	majorRadius = 0.0;
	initialPhi = 0.0;
	finalPhi = 360*deg;
	material = "";
	color = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	torus = 0;
}

// Copy constructor - be sure to use the copy constructor BLGroupElement(r)
BLCMDtorus::BLCMDtorus(const BLCMDtorus& r) : BLGroupElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	innerRadius = r.innerRadius;
	outerRadius = r.outerRadius;
	majorRadius = r.majorRadius;
	initialPhi = r.initialPhi;
	finalPhi = r.finalPhi;
	material = r.material;
	color = r.color;
	kill = r.kill;
	maxStep = r.maxStep;
	torus = r.torus;
}

int BLCMDtorus::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("torus: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultTorus.handleNamedArgs(namedArgs);
	}

	BLCMDtorus *t = new BLCMDtorus(defaultTorus);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	return retval;
}

void BLCMDtorus::defineNamedArgs()
{
	argDouble(innerRadius,"innerRadius","The inner radius of the torus (mm)");
	argDouble(outerRadius,"outerRadius","The outer radius of the torus (mm)");
	argDouble(majorRadius,"majorRadius","The major radius of the torus (mm)");
	argDouble(initialPhi,"initialPhi","The initial phi around major radius (0 degrees).",deg);
	argDouble(finalPhi,"finalPhi","The final phi around major radius (360 degrees).",deg);
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)");
	argString(material,"material","The material of the torus");
	argString(color,"color","The color of the torus (''=invisible)");
	argInt(kill,"kill","Set nonzero to kill every track that enters.");
}

G4VSolid *BLCMDtorus::getSolid()
{
	if(!torus)
		torus = new G4Torus("Torus",innerRadius,outerRadius,
			majorRadius,initialPhi,finalPhi-initialPhi);
	return torus;
}

void BLCMDtorus::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	getSolid(); // sets torus

	G4Material *mat;
	if(material != "")
		mat = getMaterial(material);
	else
		mat = parent->GetMaterial();

	G4LogicalVolume *lv = new G4LogicalVolume(torus,mat, thisname+"LogVol");
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

///@@	new Surface(this,globalRotation,globalPosition);

	printf("BLCMDtorus::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2],
		globalPosition[2]-getLength()/2.0,
		globalPosition[2]+getLength()/2.0);
}

void Surface::callback(int type)
{
	G4VVisManager* visManager = G4VVisManager::GetConcreteInstance(); 
	if(visManager) { 
		printf("Surface::callback(%d)\n",type);
		for(int i=0; i<5000; ++i) {
			G4ThreeVector pos = torus->torus->GetPointOnSurface();
			if(rot) pos = *rot * pos;
			pos += position;
			G4Circle circle(pos);
			circle.SetScreenDiameter(1.0);
			circle.SetFillStyle(G4Circle::filled);
			G4Colour colour(1.,0.,0.);
			G4VisAttributes attribs(colour);
			circle.SetVisAttributes(attribs);
			visManager->Draw(circle);
		}
	}
}
