//	BLCMDidealsectorbend.cc
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

#include <algorithm>

#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLElement.hh"
#include "BLElementField.hh"
#include "BLGlobalField.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLTune.hh"
#include "BLKillTrack.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

/**	BLCMDidealsectorbend implements an ideal bending magnet with a sector
 *	field. Ideal means a block field with no fringes.
 *
 *	The sector axis is along Y, offset by fieldCenterRadius in the +X
 *	direction from the centerline, for angle>0. angle>0 means a left bend,
 *	while angle<0 means a right bend (toward -X).
 *
 *	The magnetic field is uniform throught the sector, along the local Y 
 *	axis.
 *
 *	NOTE: the width of this object is not correct (too hard to compute).
 *	That's OK as long as some element is placed both before and after it.
 **/
class BLCMDidealsectorbend : public BLElement {
	G4double angle;
	G4double fieldCenterRadius;
	G4double fieldInnerRadius;
	G4double fieldOuterRadius;
	G4double fieldHeight;
	G4double ironInnerRadius;
	G4double ironOuterRadius;
	G4double ironHeight;
	G4double By;
	G4String fieldMaterial;
	G4String fieldColor;
	G4String ironMaterial;
	G4String ironColor;
	G4int kill;
	G4int openAperture;
	G4double maxStep;
	G4Tubs *fieldTubs;
	G4Tubs *ironTubs;
	G4VSolid *ironSolid;
	friend class IdealSectorBendField;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDidealsectorbend();

	/// Destructor.
	virtual ~BLCMDidealsectorbend() { }

	/// Copy constructor.
	BLCMDidealsectorbend(const BLCMDidealsectorbend& r);

	/// clone()
	BLElement *clone() { return new BLCMDidealsectorbend(*this); }

	/// commandName() returns "idealsectorbend".
	G4String commandName() { return "idealsectorbend"; }

	/// command() implements the idealsectorbend command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() - construct the ideal sector bending magnet
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length along the Z axis
	G4double getLength() { return ironOuterRadius*sin(fabs(angle))*2.0; }

	/// getWidth() returns the width of the iron region
	/// (this is an overestimate on the inside)
	G4double getWidth() { return std::max<double>(fieldCenterRadius,
				ironOuterRadius-fieldCenterRadius)*2.0; }

	/// getHeight() returns the height of the iron region
	G4double getHeight() { return ironHeight; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,0.0);
		if(index == 1) {
			double z = fieldCenterRadius * sin(angle);
			double x = fieldCenterRadius * (1.0 - cos(angle));
			return G4ThreeVector(x,0.0,z);
		}
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// isOutside() from BLElement.
	bool isOutside(G4ThreeVector &local, G4double tolerance)
		{ return true; }

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v)
		{ v.clear(); printf("BLCMDidealsectorbend does not participate in the geometry test\n"); }
};

BLCMDidealsectorbend defaultIdealSectorBend;	// default object

/**	IdealSectorBendField represents one placement of a ideal sector
 *	bending magnet.
 *
 **/
class IdealSectorBendField : public BLElementField {
	BLCoordinateTransform global2local;
	G4double *By;
	G4RotationMatrix rotation;
	G4double r2min;
	G4double r2max;
	G4double halfheight;
	G4double tanangle;
public:
	/// constructor. 
	IdealSectorBendField(BLCoordinateTransform& _global2local,
						BLCMDidealsectorbend *bend);

	/// addFieldValue() adds the field for this solenoid into field[].
	/// point[] is in global coordinates.
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};


// Default constructor - be sure to use the default constructor BLElement()
BLCMDidealsectorbend::BLCMDidealsectorbend() : BLElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct an ideal sector bending magnet.");
	setDescription("The field region is a sector with By specified.\n"
		"Unlike most Elements, the position of the idealsectorbend\n"
		"is the center of the front face of its field (aperture).\n"
		"angle>0 bends to the left around Y; angle<0 bends right.\n"
		"The only useful rotations are around the centerline Z.\n"
		"This element should normally be followed immediately by a "
		"cornerarc. Note that -90<=angle<=90 degrees.\n\n"
		"No fringe fields are implemented.\n\n"
		"if ironHeight<=0, no iron is generated.\n\n"
		"Note that there is no field inside the 'iron'; this can "
		"result in gross tracking errors for particles in the iron, "
		"and implies that kill=1 is desirable.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.\n\n"
		);

	// provide initial values for fields
	angle = 0.0;
	fieldCenterRadius = 0.0;
	fieldInnerRadius = 0.0;
	fieldOuterRadius = 0.0;
	fieldHeight = 0.0;
	ironInnerRadius = 0.0;
	ironOuterRadius = 0.0;
	ironHeight = 0.0;
	By = 0.0;
	fieldMaterial = "Vacuum";
	fieldColor = "";
	ironMaterial = "Fe";
	ironColor = "";
	kill = 0;
	openAperture = 0;
	maxStep = -1.0;
	fieldTubs = 0;
	ironTubs = 0;
	ironSolid = 0;
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDidealsectorbend::BLCMDidealsectorbend(const BLCMDidealsectorbend& r) : BLElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	angle = r.angle;
	fieldCenterRadius = r.fieldCenterRadius;
	fieldInnerRadius = r.fieldInnerRadius;
	fieldOuterRadius = r.fieldOuterRadius;
	fieldHeight = r.fieldHeight;
	ironInnerRadius = r.ironInnerRadius;
	ironOuterRadius = r.ironOuterRadius;
	ironHeight = r.ironHeight;
	BLTune::copyTunableArg(&By,&r.By);
	fieldMaterial = r.fieldMaterial;
	fieldColor = r.fieldColor;
	ironMaterial = r.ironMaterial;
	ironColor = r.ironColor;
	kill = r.kill;
	openAperture = r.openAperture;
	maxStep = r.maxStep;
	fieldTubs = r.fieldTubs;
	ironTubs = r.ironTubs;
	ironSolid = r.ironSolid;
}

int BLCMDidealsectorbend::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("idealsectorbend: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultIdealSectorBend.handleNamedArgs(namedArgs);
	}

	BLCMDidealsectorbend *t = new BLCMDidealsectorbend(defaultIdealSectorBend);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	if(fabs(t->angle) > 90.0)
	    printError("idealsectorbend: must have -90<=angle<=90 degrees.");

	// check material exists
	getMaterial(t->fieldMaterial);
	getMaterial(t->ironMaterial);

	t->print(argv[0]);

	return retval;
}

void BLCMDidealsectorbend::defineNamedArgs()
{
	argDouble(angle,"angle","Angle of bend (degrees).",deg);
	argDouble(fieldCenterRadius,"fieldCenterRadius","Center radius of field (mm).");
	argDouble(fieldInnerRadius,"fieldInnerRadius","Inner radius of field (mm).");
	argDouble(fieldOuterRadius,"fieldOuterRadius","Outer radius of field (mm).");
	argDouble(fieldHeight,"fieldHeight","Height of field (mm).");
	argDouble(ironInnerRadius,"ironInnerRadius","Inner radius of iron (mm).");
	argDouble(ironOuterRadius,"ironOuterRadius","Outer radius of iron (mm).");
	argDouble(ironHeight,"ironHeight","Height of iron (mm).");
	argTunable(By,"By","Magnetic field (Tesla).",tesla);
	argString(fieldMaterial,"fieldMaterial","Material of field.");
	argString(fieldColor,"fieldColor","Color of field.");
	argString(ironMaterial,"ironMaterial","Material of iron.");
	argString(ironColor,"ironColor","Color of iron.");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)");
	argInt(kill,"kill","Set nonzero to kill particles hitting the iron.");
	argInt(openAperture,"openAperture","Set nonzero to keep the aperture open (no physical volume).");
}

void BLCMDidealsectorbend::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)

{
	G4String thisname = parentName+getName();

	// values for right bend (angle <0)
	G4double phi = angle;	// start angle for iron Tubs
	G4double dphi = -angle;	// end-start angle for iron Tubs
	G4double xOffset = -fieldCenterRadius;
	G4RotationMatrix *g4rot = stringToRotationMatrix("X90");

	// angle>0 bends to the left
	if(angle >= 0) {
		phi = pi;
		dphi = angle;
		xOffset = fieldCenterRadius;
	}

	if(!fieldTubs) {
		G4double extra = 0.1*mm/fieldCenterRadius;
		fieldTubs = new G4Tubs(thisname+"IdealSectorBendField",
				fieldInnerRadius,fieldOuterRadius,
				fieldHeight/2.0,phi-extra,dphi+extra*2.0);
	}
	G4Material *mat = getMaterial(fieldMaterial);
	G4LogicalVolume *lvfield = new G4LogicalVolume(fieldTubs,mat,
					thisname+"LogVol");
	lvfield->SetVisAttributes(getVisAttrib(fieldColor));
	if(maxStep <= 0.0) maxStep = Param.getDouble("maxStep");
	lvfield->SetUserLimits(new G4UserLimits(maxStep));

	G4LogicalVolume *lviron = 0;

	if(ironHeight > 0.0) {
	    if(!ironTubs)
		ironTubs = new G4Tubs(thisname+"IdealSectorBendIron",
					ironInnerRadius,ironOuterRadius,
					ironHeight/2.0,phi,dphi);
	    if(!ironSolid)
		ironSolid = new G4SubtractionSolid(thisname+"IdealSectorBendSolid",
					ironTubs, fieldTubs);
	    mat = getMaterial(ironMaterial);
	    lviron = new G4LogicalVolume(ironSolid,mat, thisname+"LogVol");
	    lviron->SetVisAttributes(getVisAttrib(ironColor));
	    lviron->SetUserLimits(new G4UserLimits(maxStep));
	}

	// (geant4 rotation convention is backwards from g4beamline)
	if(relativeRotation)
		*g4rot = *g4rot * relativeRotation->inverse();

	G4ThreeVector pos(xOffset,0.0,0.0);
	if(relativeRotation)
		pos = *relativeRotation * pos;

	pos += relativePosition;

	if(openAperture == 0)
		new G4PVPlacement(g4rot,pos,lvfield,thisname,
					parent,false,0,surfaceCheck);
	if(ironHeight > 0.0) {
	    G4PVPlacement *pv = new G4PVPlacement(g4rot,pos,lviron,
					thisname, parent,false,0,surfaceCheck);
	    if(kill)
		BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
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
	G4ThreeVector globalPosition(pos + parentPosition);
	if(parentRotation)
		globalPosition = *parentRotation * relativePosition +
				parentPosition;

	BLCoordinateTransform global2local(globalRotation,globalPosition);

	IdealSectorBendField *p = new IdealSectorBendField(global2local,this);
	BLGlobalField::getObject()->addElementField(p);

	printf("BLCMDidealsectorbend::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2]);
}

IdealSectorBendField::IdealSectorBendField(BLCoordinateTransform& _global2local,
			BLCMDidealsectorbend *bend)
		       : BLElementField(), global2local(), rotation()
{
	global2local = _global2local;
	By = &bend->By;
	rotation = global2local.getRotation().inverse();
	r2min = bend->fieldInnerRadius;
	r2max = bend->fieldOuterRadius;
	r2min *= r2min;
	r2max *= r2max;
	halfheight = bend->fieldHeight/2.0;
	tanangle = tan(bend->angle);

	// set global bounding box
	G4double local[4], global[4];
	local[3] = 0.0;
	double r1 = bend->fieldInnerRadius;
	double r2 = bend->fieldOuterRadius;
	double phi0 = 0.0;	 // values for angle < 0 (to the right)
	double phi1 = -bend->angle;
	if(bend->angle >= 0.0)  {
		phi0 = M_PI-bend->angle; // values for angle >= 0 (to the left)
		phi1 = M_PI;
	}
	phi0 -= 0.05; // expand bounding box slightly (hard edges)
	phi1 += 0.05;
	r1 *= 0.95;
	r2 *= 1.05;
	double dphi = (phi1-phi0)/100.0;
	BLAssert(dphi > 0.0);
	for(double phi=phi0; phi<=phi1; phi+=dphi) {
		local[0] = r1*cos(phi);
		local[1] = halfheight;
		local[2] = r1*sin(phi);
		global2local.getGlobal(local,global);
		setGlobalPoint(global);
		local[1] = -halfheight;
		global2local.getGlobal(local,global);
		setGlobalPoint(global);
		local[0] = r2*cos(phi);
		local[1] = halfheight;
		local[2] = r2*sin(phi);
		global2local.getGlobal(local,global);
		setGlobalPoint(global);
		local[1] = -halfheight;
		global2local.getGlobal(local,global);
		setGlobalPoint(global);
	}
}

void IdealSectorBendField::addFieldValue(const G4double point[4], G4double field[6]) const
{
	G4ThreeVector global(point[0],point[1],point[2]);
	G4ThreeVector local;

	global2local.getLocal(local,global);
	G4double r2 = local[0]*local[0] + local[2]*local[2];
	if(fabs(local[1]) > halfheight || r2 < r2min || r2 > r2max)
		return;
	if(tanangle > 0.0) {
		if(local[0] > 0.0 || local[2] < 0.0 || 
		   local[2] > -local[0]*tanangle)
			return;
	} else {
		if(local[0] < 0.0 || local[2] < 0.0 || 
		   local[2] > -local[0]*tanangle)
			return;
	}

	G4ThreeVector B(0.0,*By,0.0);
	if(global2local.isRotated())
		B = rotation * B;
	field[0] += B[0];
	field[1] += B[1];
	field[2] += B[2];
}
