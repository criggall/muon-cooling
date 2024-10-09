//	BLCMDsphere.cc
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
#include "G4Sphere.hh"
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

/**	BLCMDsphere implements a sphere of material.
 *
 **/
class BLCMDsphere : public BLGroupElement {
	G4double innerRadius;
	G4double outerRadius;
	G4double initialPhi;
	G4double finalPhi;
	G4double initialTheta;
	G4double finalTheta;
	G4String material;
	G4String color;
	G4int kill;
	G4double maxStep;
	G4Sphere *sphere;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDsphere();

	/// Destructor.
	virtual ~BLCMDsphere() { }

	/// Copy constructor.
	BLCMDsphere(const BLCMDsphere& r);

	/// clone()
	BLElement *clone() { return new BLCMDsphere(*this); }

	/// commandName() returns "sphere".
	G4String commandName() { return "sphere"; }

	/// command() implements the sphere command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	virtual G4VSolid *getSolid();

	/// construct() - construct the sphere
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the sphere
	G4double getLength() { return 2.0*outerRadius; }

	/// getWidth() returns the width of the sphere
	G4double getWidth() { return 2.0*outerRadius; }

	/// getHeight() returns the height of the sphere
	G4double getHeight() { return 2.0*outerRadius; }

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
		{ BLAssert(sphere!=0);  return sphere->Inside(local) != kInside; }

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v)
		{ v.clear();
		  for(int n=0; n<npoints*5; ++n)
		  	v.push_back(sphere->GetPointOnSurface());
		}

	/// isWithin() from BLGroupElement.
	bool isWithin(G4ThreeVector &local, G4double tolerance) 
	    { BLAssert(sphere!=0);  return sphere->Inside(local) != kOutside; }
};

BLCMDsphere defaultSphere;	// default object

// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDsphere::BLCMDsphere() : BLGroupElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a sphere (or section of one)");
	setDescription("This is a direct interface to G4Sphere.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.");

	// provide initial values for fields
	innerRadius = 0.0;
	outerRadius = 0.0;
	initialPhi = 0.0;
	finalPhi = 360.0*deg;
	initialTheta = 0.0;
	finalTheta = 360.0*deg;
	material = "";
	color = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	sphere = 0;
}

// Copy constructor - be sure to use the copy constructor BLGroupElement(r)
BLCMDsphere::BLCMDsphere(const BLCMDsphere& r) : BLGroupElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	innerRadius = r.innerRadius;
	outerRadius = r.outerRadius;
	initialPhi = r.initialPhi;
	finalPhi = r.finalPhi;
	initialTheta = r.initialTheta;
	finalTheta = r.finalTheta;
	material = r.material;
	color = r.color;
	kill = r.kill;
	maxStep = r.maxStep;
	sphere = r.sphere;
}

int BLCMDsphere::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("sphere: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultSphere.handleNamedArgs(namedArgs);
	}

	BLCMDsphere *t = new BLCMDsphere(defaultSphere);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	return retval;
}

void BLCMDsphere::defineNamedArgs()
{
	argDouble(innerRadius,"innerRadius","The inside radius of the sphere (mm)");
	argDouble(outerRadius,"outerRadius","The outer radius of the sphere (mm)");
	argDouble(initialPhi,"initialPhi","The initial Phi value (deg; 0 for all)",deg);
	argDouble(finalPhi,"finalPhi","The final Phi value (deg; 360 for all)",deg);
	argDouble(initialTheta,"initialTheta","The initialTheta of the sphere (deg, 0 for all)",deg);
	argDouble(finalTheta,"finalTheta","The finalTheta of the sphere (deg, 180 for all)",deg);
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)");
	argString(material,"material","The material of the sphere");
	argString(color,"color","The color of the sphere (''=invisible)");
	argInt(kill,"kill","Set nonzero to kill every track that enters.");
}

G4VSolid *BLCMDsphere::getSolid()
{
	if(!sphere)
		sphere = new G4Sphere("Sphere", innerRadius, outerRadius,
				initialPhi,(finalPhi-initialPhi),
				initialTheta,(finalTheta-initialTheta));
	return sphere;
}

void BLCMDsphere::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	getSolid(); // sets sphere

	G4Material *mat;
	if(material != "")
		mat = getMaterial(material);
	else
		mat = parent->GetMaterial();

	G4LogicalVolume *lv = new G4LogicalVolume(sphere,mat, thisname+"LogVol");
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

	printf("BLCMDsphere::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f kill=%d\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2],
		globalPosition[2]-getLength()/2.0,
		globalPosition[2]+getLength()/2.0,
		kill);
}

