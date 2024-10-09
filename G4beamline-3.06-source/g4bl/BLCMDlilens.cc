//	BLCMDlilens.cc
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
#include "BLElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLKillTrack.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

/**	BLCMDlilens implements a Lithium lens.
 *
 *	This is only the current-carrying Li cylinder, and its magnetic field.
 *	The field exists only between the end planes of the cylinder, out to
 *	radial infinity.
 *
 **/
class BLCMDlilens : public BLElement {
	G4double radius;
	G4double length;
	G4double current;
	G4String material;
	G4String color;
	G4double maxStep;
	G4Tubs *lilens;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDlilens();

	/// Destructor.
	virtual ~BLCMDlilens() { }

	/// Copy constructor.
	BLCMDlilens(const BLCMDlilens& r);

	/// clone()
	BLElement *clone() { return new BLCMDlilens(*this); }

	/// commandName() returns "lilens".
	G4String commandName() { return "lilens"; }

	/// command() implements the lilens command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() - construct the Li lens
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the Li lens
	G4double getLength() { return length; }

	/// getWidth() returns the width of the Li lens
	G4double getWidth() { return 2.0*radius; }

	/// getHeight() returns the height of the Li lens
	G4double getHeight() { return 2.0*radius; }

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
};

BLCMDlilens defaultLiLens;	// default object

class LiLensField : public BLElementField {
	BLCoordinateTransform global2local;
	G4RotationMatrix rotation;
	G4double radius;
	G4double length;
	G4double current;
public:
	LiLensField(G4double _radius, G4double _length, G4double _current,
				G4RotationMatrix *rot, G4ThreeVector &pos);
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};

BLCMDlilens::BLCMDlilens() : BLElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a simple Lithium lens.");
	setDescription(" This element consists of a current-carrying cylinder "
		"and its field. The field exists only between the end planes "
		"of the cylinder, out to radial infinity.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.");
	radius = 5.0;
	length = 100.0;
	current = 100000.0;
	material = "Li";
	color = "1,1,1";
	maxStep = -1.0;
	lilens = 0;
}

BLCMDlilens::BLCMDlilens(const BLCMDlilens& r) : BLElement(r)
{
	radius = r.radius;
	length = r.length;
	current = r.current;
	material = r.material;
	color = r.color;
	maxStep = r.maxStep;
	lilens = r.lilens;
}

int BLCMDlilens::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("lilens: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultLiLens.handleNamedArgs(namedArgs);
	}

	BLCMDlilens *t = new BLCMDlilens(defaultLiLens);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	return retval;
}

void BLCMDlilens::defineNamedArgs()
{
	argDouble(radius,"radius","The radius of the current-carrying cylinder (5 mm)", mm);
	argDouble(length,"length","The length of the cylinder (100 mm).", mm);
	argDouble(current,"current","The current in the cylinder (100000 Amp).");
	argString(material,"material","The material of the cylinder (Li).");
	argString(color,"color","The color of the tube or cylinder (''=invisible)");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm).",
			mm);
}

void BLCMDlilens::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	if(!lilens)
		lilens = new G4Tubs(thisname+"Tubs", 0.0, radius,
						length/2.0, 0.0,2.0*M_PI);
	G4Material *mat = getMaterial(material);
	G4LogicalVolume *lv = new G4LogicalVolume(lilens,mat,thisname+"LogVol");
	lv->SetVisAttributes(getVisAttrib(color));
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");
	lv->SetUserLimits(new G4UserLimits(maxStep));

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4VPhysicalVolume *pv = new G4PVPlacement(g4rot, relativePosition,lv,
					thisname, parent,false,0,surfaceCheck);

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

	new LiLensField(radius,length,current,globalRotation,globalPosition);

	printf("BLCMDlilens::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2],
		globalPosition[2]-getLength()/2.0,
		globalPosition[2]+getLength()/2.0);
}

void BLCMDlilens::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	generateTubs(npoints, 0.0, radius, 0.0, 2.0*M_PI, length, v);
}

G4bool BLCMDlilens::isOutside(G4ThreeVector &local, G4double tolerance)
{
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	return r > radius-tolerance || fabs(local[2]) > length/2.0-tolerance;
}

LiLensField::LiLensField(G4double _radius, G4double _length, G4double _current,
				G4RotationMatrix *rot, G4ThreeVector &pos)
					: global2local(rot,pos), rotation(*rot)
{
	radius = _radius;
	length = _length;
	current = _current;

	if(global2local.isRotated()) {
		rotation = global2local.getRotation();
		rotation.invert();
	}

	// keep default infinite bounding box

	BLGlobalField::getObject()->addElementField(this);
}

void LiLensField::addFieldValue(const G4double point[4],G4double field[6]) const
{
	G4double local[4], thisField[6];

	global2local.getLocal(local,point);
	if(fabs(local[2]) > length/2.0) return;

	double r=sqrt(local[0]*local[0]+local[1]*local[1]);
	double phi=atan2(local[1],local[0]);

	// A current of 1 A at a radius 1 mm has a field 0.00019999 T.
	// For current along +z and x>0, B is in the +y direction.
	double Bmax = 0.00019999*tesla*current/(radius/mm);
	G4ThreeVector B(0.0,0.0,0.0);
	if(r < radius)
		B[1] = Bmax * r / radius;
	else
		B[1] = Bmax * radius / r;

	B[0] = -B[1] * sin(phi);
	B[1] = B[1] * cos(phi);

	if(!rotation.isIdentity()) B = rotation * B;

	field[0] += B[0];
	field[1] += B[1];
	field[2] += B[2];
}
