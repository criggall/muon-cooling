//	BLCMDtubs.cc
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
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLCommandAlias.hh"
#include "BLGroupElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLKillTrack.hh"

/**	BLCMDtubs implements a tube or cylinder of material, axis along Z.
 *
 **/
class BLCMDtubs : public BLGroupElement {
	G4double innerRadius;
	G4double outerRadius;
	G4double initialPhi;
	G4double finalPhi;
	G4double length;
	G4String material;
	G4String color;
	G4int kill;
	G4double maxStep;
	G4Tubs *tubs;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDtubs();

	/// Destructor.
	virtual ~BLCMDtubs() { }

	/// Copy constructor.
	BLCMDtubs(const BLCMDtubs& r);

	/// clone()
	BLElement *clone() { return new BLCMDtubs(*this); }

	/// commandName() returns "tubs".
	G4String commandName() { return "tubs"; }

	/// command() implements the tubs command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// getSolid() return the solid.
	G4VSolid *getSolid();

	G4String getMaterialName() const { return material; }

	/// construct() - construct the tube/cylinder.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the tube/cylinder.
	G4double getLength() { return length; }

	/// getWidth() returns the width of the tube/cylinder.
	G4double getWidth() { return 2.0*outerRadius; }

	/// getHeight() returns the height of the tube/cylinder.
	G4double getHeight() { return 2.0*outerRadius; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// generatePoints() from BLElement
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isOutside() from BLElement
	G4bool isOutside(G4ThreeVector &local, G4double tolerance);

	/// isWithin() from BLGroupElement
	bool isWithin(G4ThreeVector &local, G4double tolerance);
};

BLCMDtubs defaultTubs;	// default object

BLCommandAlias aliasTubs("cylinder",defaultTubs);	// Alias "cylinder"
BLCommandAlias aliasTubs2("tube",defaultTubs);		// Alias "tube"

// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDtubs::BLCMDtubs() : BLGroupElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a tube or cylinder with axis along z.");
	setDescription("This is a direct interface to G4Tubs, which can\n"
		"implement a tube or cylinder; either can subtend less "
		"than 360 degrees in phi.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.");

	// provide initial values for fields
	innerRadius = 0.0;
	outerRadius = 0.0;
	initialPhi = 0.0;
	finalPhi = 360.0*deg;
	length = 0.0;
	material = "";
	color = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	tubs = 0;
}

// Copy constructor - be sure to use the copy constructor BLGroupElement(r)
BLCMDtubs::BLCMDtubs(const BLCMDtubs& r) : BLGroupElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	innerRadius = r.innerRadius;
	outerRadius = r.outerRadius;
	initialPhi = r.initialPhi;
	finalPhi = r.finalPhi;
	length = r.length;
	material = r.material;
	color = r.color;
	kill = r.kill;
	maxStep = r.maxStep;
	tubs = r.tubs;
}

int BLCMDtubs::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("tubs: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultTubs.handleNamedArgs(namedArgs);
	}

	BLCMDtubs *t = new BLCMDtubs(defaultTubs);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	return retval;
}

void BLCMDtubs::defineNamedArgs()
{
	argDouble(innerRadius,"innerRadius","The inside of the tube, 0.0 for cylinder (mm)",
			mm);
	argDouble(outerRadius,"outerRadius","The outer radius of the tube or cylinder (mm)",
			mm);
	argDouble(initialPhi,"initialPhi","The initial Phi value (deg; 0 for all)",
			deg);
	argDouble(finalPhi,"finalPhi","The final Phi value (deg; 360 for all)",
			deg);
	argDouble(length,"length","The length of the tube or cylinder (mm)",
			mm);
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)",
			mm);
	argString(material,"material","The material of the tube or cylinder");
	argString(color,"color","The color of the tube or cylinder (''=invisible)");
	argInt(kill,"kill","Set nonzero to kill every track that enters.");
	argDouble(outerRadius,"radius","Synonym for outerRadius (mm)",mm);
}

G4VSolid *BLCMDtubs::getSolid()
{
	if(!tubs) tubs = new G4Tubs("Tubs", innerRadius, outerRadius,
			length/2.0, initialPhi,finalPhi-initialPhi);
	return tubs;
}

void BLCMDtubs::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	getSolid(); // sets tubs

	G4Material *mat;
	if(material != "")
		mat = getMaterial(material);
	else
		mat = parent->GetMaterial();

	G4LogicalVolume *lv = new G4LogicalVolume(tubs,mat, thisname+"LogVol");
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

	printf("BLCMDtubs::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f kill=%d\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2],
		globalPosition[2]-getLength()/2.0,
		globalPosition[2]+getLength()/2.0,
		kill);
}

void BLCMDtubs::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	generateTubs(npoints, innerRadius, outerRadius, initialPhi, finalPhi,
			length, v);
}

G4bool BLCMDtubs::isOutside(G4ThreeVector &local, G4double tolerance)
{
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	G4double phi = atan2(local[1],local[0]);
	bool phiInside = phi > initialPhi+0.001 && phi < finalPhi-0.001;
	if(!phiInside && phi < 0.0) {
		phi += pi+pi;
		phiInside = phi > initialPhi+0.001 && phi < finalPhi-0.001;
	}
	return r < innerRadius+tolerance || r > outerRadius-tolerance ||
		fabs(local[2]) > length/2.0-tolerance || !phiInside;
}

bool BLCMDtubs::isWithin(G4ThreeVector &local, G4double tolerance)
{
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	G4double phi = atan2(local[1],local[0]);
	if(phi < 0.0 && phi < initialPhi-0.001) phi += pi+pi;
	return r > innerRadius-tolerance && r < outerRadius+tolerance &&
		fabs(local[2]) < length/2.0+tolerance &&
		phi >= initialPhi-0.001 && phi < finalPhi+0.001;
}
