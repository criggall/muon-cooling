//	BLCMDboolean.cc
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
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

#include "BLGroupElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLKillTrack.hh"

/**	BLCMDboolean implements a boolean of material.
 *
 **/
class BLCMDboolean : public BLGroupElement {
	G4double height;
	G4double width;
	G4double length;
	G4String material;
	G4String color;
	G4int kill;
	G4double maxStep;
	G4String op;
	G4double x;
	G4double y;
	G4double z;
	G4String rotation;
	G4String element1;
	G4String element2;
	G4VSolid *boolean;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDboolean();

	/// Destructor.
	virtual ~BLCMDboolean() { }

	/// Copy constructor.
	BLCMDboolean(const BLCMDboolean& r);

	/// clone()
	BLElement *clone() { return new BLCMDboolean(*this); }

	/// commandName() returns "boolean".
	G4String commandName() { return "boolean"; }

	/// command() implements the boolean command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// setup() will setup the boolean solid.
	int setup();

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// getSolid() returns the solid.
	virtual G4VSolid *getSolid() { return boolean; }

	virtual G4String getMaterialName() const { return material; }

	/// construct() - construct the boolean.
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

BLCMDboolean defaultBoolean;	// default object

// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDboolean::BLCMDboolean() : BLGroupElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("Construct an element from a boolean operation between elements.");
	setDescription("This command takes a pair of simple elements and "
	"creates a new element that is a boolean operation between the solids "
	"of the initial elements. Boolean operations are: union, subtraction, "
	"and intersection. Simple elements are those consisting of a "
	"single solid, with no EM field: boolean, box, extrusion, polycone "
	"sphere, tessellatedsolid, torus, trap, tubs.\n\n"
	"The command is:\n"
	"    boolean op=subtraction e3 e1 e2 [args...]\n"
	"This creates element e3 from the subtraction of e2 from e1. e1, e2, "
	"and e3 can be used in further boolean operations, or can be placed "
	"into the world.\n\n"
	"The local coordinates of the new element are those of the first "
	"element. If material is not specified, the material of the "
	"first element is used.\n\n"
	"The two input solids should intersect, and the resulting solid should "
	"be a single piece; if not the results are undefined.\n\n"
	"Note that for op=subtraction, the two elements should not share "
	"any common faces; to create a hole, the subtracted (second) element "
	"should extend at least a micron outside the first element.\n\n"
	"The values of width, height, and length default to those of e1. If "
	"taht is not correct, provide the correct values as arguments.\n\n"
	"Note also that visualization of boolean solids may not work properly, "
	"except for the RayTracer viewer. Tracking is fine.\n\n"
	"This element must be placed (via the place command).");

	// provide initial values for fields
	height = 0.0;
	width = 0.0;
	length = 0.0;
	material = "";
	color = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	op = "";
	x = y = z = 0.0;
	rotation = "";
	boolean = 0;
}

// Copy constructor - be sure to use the copy constructor BLGroupElement(r)
BLCMDboolean::BLCMDboolean(const BLCMDboolean& r) : BLGroupElement(r)
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
	op = r.op;
	x = r.x;
	y = r.y;
	z = r.z;
	rotation = r.rotation;
	boolean = r.boolean;
}

int BLCMDboolean::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 3) {
		printError("boolean: Invalid command, must have name, element1, and element2");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultBoolean.handleNamedArgs(namedArgs);
	}

	BLCMDboolean *t = new BLCMDboolean(defaultBoolean);
	t->setName(argv[0]);
	t->element1 = argv[1];
	t->element2 = argv[2];
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	retval = t->setup();

	t->print(argv[0]);

	return retval;
}

int BLCMDboolean::setup()
{
	BLElement *e1 = BLElement::find(element1);
	BLElement *e2 = BLElement::find(element2);
	G4VSolid *s1 = (e1 ? e1->getSolid() : 0);
	G4VSolid *s2 = (e2 ? e2->getSolid() : 0);
	if(!e1 || !s1) {
	    printError("boolean: cannot find element1 '%s', or it is unusable",
							e1->getName().c_str());
	    return -1;
	}
	if(!e2 || !s2) {
	    printError("boolean: cannot find element2 '%s', or it is unusable",
							e2->getName().c_str());
	    return -1;
	}

	// get material, check it exists
	if(material.size() == 0) material = e1->getMaterialName();
	if(material.size() == 0) material = "Vacuum";
	getMaterial(material);

	// construct the solid
	G4ThreeVector offset(x,y,z);
	G4RotationMatrix *rot = stringToRotationMatrix(rotation);
	if(op == "union") {
		boolean = new G4UnionSolid(getName(),s1,s2,rot,offset);
	} else if(op == "subtraction") {
		boolean = new G4SubtractionSolid(getName(),s1,s2,rot,offset);
	} else if(op == "intersection") {
		boolean = new G4IntersectionSolid(getName(),s1,s2,rot,offset);
	} else {
		printError("boolean: invalid op '%s'\n",op.c_str());
		return -1;
	}

	if(height <= 0.0) height = e1->getHeight();
	if(width <= 0.0) width = e1->getWidth();
	if(length <= 0.0) length = e1->getLength();

	return 0;
}

void BLCMDboolean::defineNamedArgs()
{
	argDouble(maxStep,"maxStep","The maximum stepsize (mm).");
	argString(material,"material","The material of the boolean.");
	argString(color,"color","The color of the boolean (white).");
	argInt(kill,"kill","Set nonzero to kill every track that enters.");
	argString(op,"op","The boolean operation (union, subtraction, or intersection).");
	argDouble(x,"x","The x offset of element2 relative to element1 (mm).");
	argDouble(y,"y","The y offset of element2 relative to element1 (mm).");
	argDouble(z,"z","The z offset of element2 relative to element1 (mm).");
	argString(rotation,"rotation","The rotation of element2 ('').");
}

void BLCMDboolean::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	G4Material *mat = getMaterial(material);
	G4LogicalVolume *lv = new G4LogicalVolume(boolean,mat, thisname+"LogVol");
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

	printf("BLCMDboolean::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2],
		globalPosition[2]-getLength()/2.0,
		globalPosition[2]+getLength()/2.0);
}

G4bool BLCMDboolean::isOutside(G4ThreeVector &local, G4double tolerance)
{
	return (boolean->Inside(local) == kOutside);
}

void BLCMDboolean::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	v.clear();
	for(int n=0; n<npoints*5; ++n)
		v.push_back(boolean->GetPointOnSurface());
}

G4bool BLCMDboolean::isWithin(G4ThreeVector &local, G4double tolerance)
{
	return (boolean->Inside(local) != kOutside);
}
