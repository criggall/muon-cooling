//	BLCMDgenericsectorbend.cc
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
 
 
 This class is written by Yu Bao, 03/21/2016, baoyubaoyu@gmail.com
 
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
#include "BLEngeFunction.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

const G4double FRINGE_ACCURACY=1.0e-4;

/**	BLCMDgenericsectorbend implements an generic bending magnet with a sector
 *	field. The magnet can be defined by the user as a multi-function magnet, up to octuple.
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
class BLCMDgenericsectorbend : public BLElement {
	G4double angle;
	G4double fieldCenterRadius;
	G4double fieldInnerRadius;
	G4double fieldOuterRadius;
	G4double fieldHeight;
	G4double ironInnerRadius;
	G4double ironOuterRadius;
	G4double ironHeight;
	G4String fieldMaterial;
	G4String fieldColor;
	G4String ironMaterial;
	G4String ironColor;
	G4int kill;
	G4double maxStep;
	G4Tubs *fieldTubs;
	G4Tubs *ironTubs;
	G4VSolid *ironSolid;
	////// Add arg for solenoidal, dipole, quadrupole, sextupole...
	////// Bao 3/17/16
	G4double SolenoidalField;
	G4double DipoleField;
	G4double SkewQuadrupoleField;
	G4double QuadrupoleField;
	G4double SextupoleField;
	G4double SkewSextupoleField;
	G4double OctupoleField;
	G4double SkewOctupoleField;
	
	//////// Add fringe field. Bao 4/18/2016 /////////////////////
	G4double fringe;
	G4double fringeFactor;
	
	friend class GenericSectorBendField;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDgenericsectorbend();
	
	/// Destructor.
	virtual ~BLCMDgenericsectorbend() { }
	
	/// Copy constructor.
	BLCMDgenericsectorbend(const BLCMDgenericsectorbend& r);
	
	/// clone()
	BLElement *clone() { return new BLCMDgenericsectorbend(*this); }
	
	/// commandName() returns "genericsectorbend".
	G4String commandName() { return "genericsectorbend"; }
	
	/// command() implements the genericsectorbend command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);
	
	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();
	
	/// construct() - construct the generic sector bending magnet
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
	{ v.clear(); printf("BLCMDgenericsectorbend does not participate in the geometry test\n"); }
};

BLCMDgenericsectorbend defaultGenericSectorBend;	// default object

/**	GenericSectorBendField represents one placement of a generic sector
 *	bending magnet.
 *
 **/
class GenericSectorBendField : public BLElementField {
	BLCoordinateTransform global2local;
	G4RotationMatrix rotation;
	G4double fieldInnerRadius;
	G4double fieldOuterRadius;
	G4double halfheight;
	G4double tanangle;
	
	G4double BendAngle;
	G4double *SolenoidalField;
	G4double *DipoleField;
	G4double *SkewQuadrupoleField;
	G4double *QuadrupoleField;
	G4double *SextupoleField;
	G4double *SkewSextupoleField;
	G4double *OctupoleField;
	G4double *SkewOctupoleField;
	G4double fieldRadius;
	G4double fieldHeight;
	
	G4double fringeMaxZ;
	G4double fringeDepth;
	
	BLEngeFunction enge_0;//No Fringe
	BLEngeFunction enge_Dipole;//Dipole Fringe
	BLEngeFunction enge_Quad;//Quad/SkewQ Fringe
	BLEngeFunction enge_Sext;//Sextupole or higher field fringe
	
public:
	/// constructor.
	GenericSectorBendField(BLCoordinateTransform& _global2local,
						   BLCMDgenericsectorbend *bend);
	
	/// addFieldValue() adds the field for this solenoid into field[].
	/// point[] is in global coordinates.
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};


// Default constructor - be sure to use the default constructor BLElement()
BLCMDgenericsectorbend::BLCMDgenericsectorbend() : BLElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct an generic sector bending magnet.");
	setDescription("The field region is a sector with multi-function field specified. For straight magnets set fieldRadius=-1. \n"
				   "Unlike most Elements, the position of the genericsectorbend\n"
				   "is the center of the front face of its field (aperture).\n"
				   "angle>0 bends to the left around Y; angle<0 bends right.\n"
				   "The only useful rotations are around the centerline Z.\n"
				   "This element should normally be followed immediately by a "
				   "cornerarc. Note that -90<=angle<=90 degrees.\n\n"
				   "Note that there is no field inside the 'iron'; this can "
				   "result in gross tracking errors for particles in the iron, "
				   "and implies that kill=1 is desirable.\n\n"
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
	fieldMaterial = "Vacuum";
	fieldColor = "";
	ironMaterial = "Fe";
	ironColor = "";
	kill = 0;
	maxStep = -1.0;
	fieldTubs = 0;
	ironTubs = 0;
	ironSolid = 0;
	
	//// Added /////////////////////////////
	SolenoidalField = 0;
	DipoleField = 0;
	SkewQuadrupoleField = 0;
	QuadrupoleField = 0;
	SextupoleField = 0;
	SkewSextupoleField = 0;
	OctupoleField = 0;
	SkewOctupoleField = 0;
	
	fringe = 0;
	fringeFactor = 1.0;
	
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDgenericsectorbend::BLCMDgenericsectorbend(const BLCMDgenericsectorbend& r) : BLElement(r)
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
	fieldMaterial = r.fieldMaterial;
	fieldColor = r.fieldColor;
	ironMaterial = r.ironMaterial;
	ironColor = r.ironColor;
	kill = r.kill;
	maxStep = r.maxStep;
	fieldTubs = r.fieldTubs;
	ironTubs = r.ironTubs;
	ironSolid = r.ironSolid;
	
	SolenoidalField = r.SolenoidalField;
	DipoleField = r.DipoleField;
	SkewQuadrupoleField = r.SkewQuadrupoleField;
	QuadrupoleField = r.QuadrupoleField;
	SextupoleField = r.SextupoleField;
	SkewSextupoleField = r.SkewSextupoleField;
	OctupoleField = r.OctupoleField;
	SkewOctupoleField = r.SkewOctupoleField;
	
	fringe = r.fringe;
	fringeFactor = r.fringeFactor;
	
}

int BLCMDgenericsectorbend::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("genericsectorbend: Invalid command, must have name");
		return -1;
	}
	
	if(argv[0] == "default") {
		return defaultGenericSectorBend.handleNamedArgs(namedArgs);
	}
	
	BLCMDgenericsectorbend *t = new BLCMDgenericsectorbend(defaultGenericSectorBend);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);
	
	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");
	
	if(fabs(t->angle) > 90.0)
	    printError("genericsectorbend: must have -90<=angle<=90 degrees.");
	
	// check material exists
	getMaterial(t->fieldMaterial);
	getMaterial(t->ironMaterial);
	
	t->print(argv[0]);
	
	return retval;
}

void BLCMDgenericsectorbend::defineNamedArgs()
{
	argDouble(angle,"angle","Angle of bend (degrees).",deg);
	argDouble(fieldCenterRadius,"fieldCenterRadius","Center radius of field (mm).");
	argDouble(fieldInnerRadius,"fieldInnerRadius","Inner radius of field (mm).");
	argDouble(fieldOuterRadius,"fieldOuterRadius","Outer radius of field (mm).");
	argDouble(fieldHeight,"fieldHeight","Height of field (mm).");
	argDouble(ironInnerRadius,"ironInnerRadius","Inner radius of iron (mm).");
	argDouble(ironOuterRadius,"ironOuterRadius","Outer radius of iron (mm).");
	argDouble(ironHeight,"ironHeight","Height of iron (mm).");
	argString(fieldMaterial,"fieldMaterial","Material of field.");
	argString(fieldColor,"fieldColor","Color of field.");
	argString(ironMaterial,"ironMaterial","Material of iron.");
	argString(ironColor,"ironColor","Color of iron.");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)");
	argInt(kill,"kill","Set nonzero to kill particles hitting the iron.");
	
	argDouble(SolenoidalField,"SolenoidalField","Solenoidal field (Tesla).",tesla);
	argDouble(DipoleField,"DipoleField","Dipole field (Tesla).",tesla);
	argDouble(QuadrupoleField,"QuadrupoleField","Quadrupole field (Tesla/meter).",tesla/meter);
	argDouble(SkewQuadrupoleField,"SkewQuadrupoleField","Skew Quadrupole field (Tesla/meter).",tesla/meter);
	argDouble(SextupoleField,"SextupoleField","Sextupole Field (T/m^2).",tesla/meter/meter);
	argDouble(SkewSextupoleField,"SkewSextupoleField","Skew Sextupole Field (T/m^2).",tesla/meter/meter);
	argDouble(OctupoleField,"OctupoleField","Octupole Field (T/m^3).",tesla/meter/meter/meter);
	argDouble(SkewOctupoleField,"SkewOctupoleField","Skew Octupole field (T/m^3).",tesla/meter/meter/meter);
	
	argDouble(fringe,"fringe","Fringe field computation, set to non-zero to enable.");
	argDouble(fringeFactor,"fringeFactor","Fringe depth factor (1.0).");
	
}

void BLCMDgenericsectorbend::construct(G4RotationMatrix *relativeRotation,
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
		fieldTubs = new G4Tubs(thisname+"GenericSectorBendField",
							   fieldInnerRadius,fieldOuterRadius,
							   fieldHeight/2.0,phi-extra,dphi+extra*2.0);
	}
	G4Material *mat = getMaterial(fieldMaterial);
	G4LogicalVolume *lvfield = new G4LogicalVolume(fieldTubs,mat,
												   thisname+"LogVol");
	lvfield->SetVisAttributes(getVisAttrib(fieldColor));
	if(maxStep <= 0.0) maxStep = Param.getDouble("maxStep");
	lvfield->SetUserLimits(new G4UserLimits(maxStep));
	
	if(!ironTubs)
		ironTubs = new G4Tubs(thisname+"GenericSectorBendIron",
							  ironInnerRadius,ironOuterRadius,
							  ironHeight/2.0,phi,dphi);
	if(!ironSolid)
		ironSolid = new G4SubtractionSolid(thisname+"GenericSectorBendSolid",
										   ironTubs, fieldTubs);
	mat = getMaterial(ironMaterial);
	G4LogicalVolume *lviron = new G4LogicalVolume(ironSolid,mat,
												  thisname+"LogVol");
	lviron->SetVisAttributes(getVisAttrib(ironColor));
	lviron->SetUserLimits(new G4UserLimits(maxStep));
	
	// (geant4 rotation convention is backwards from g4beamline)
	if(relativeRotation)
		*g4rot = *g4rot * relativeRotation->inverse();
	
	G4ThreeVector pos(xOffset,0.0,0.0);
	if(relativeRotation)
		pos = *relativeRotation * pos;
	
	pos += relativePosition;
	
	new G4PVPlacement(g4rot,pos,lvfield,thisname,
					  parent,false,0,surfaceCheck);
	G4PVPlacement *pv = new G4PVPlacement(g4rot,pos,lviron,
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
	G4ThreeVector globalPosition(pos + parentPosition);
	if(parentRotation)
		globalPosition = *parentRotation * relativePosition +
		parentPosition;
	
	BLCoordinateTransform global2local(globalRotation,globalPosition);
	
	GenericSectorBendField *p = new GenericSectorBendField(global2local,this);
	BLGlobalField::getObject()->addElementField(p);
	
	printf("BLCMDgenericsectorbend::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		   thisname.c_str(),parentName.c_str(),relativePosition[2],
		   globalPosition[2]);
}

GenericSectorBendField::GenericSectorBendField(BLCoordinateTransform& _global2local,
											   BLCMDgenericsectorbend *bend)
: BLElementField(), global2local(), rotation(), enge_0(ENGE_BLOCK), enge_Dipole(ENGE_BEND),enge_Quad(ENGE_QUAD), enge_Sext(ENGE_QUAD)
// Use ENGE_QUAD for sextupole, for lack of anything better.
{
	global2local = _global2local;
	
	//// Added by Bao/////////////////////////////
	SolenoidalField = &bend->SolenoidalField;
	DipoleField = &bend->DipoleField;
	QuadrupoleField = &bend->QuadrupoleField;
	SkewQuadrupoleField = &bend->SkewQuadrupoleField;
	SextupoleField = &bend->SextupoleField;
	SkewSextupoleField = &bend->SkewSextupoleField;
	OctupoleField = &bend->OctupoleField;
	SkewOctupoleField = &bend->SkewOctupoleField;
	
	fieldRadius = bend->fieldCenterRadius;
	fieldInnerRadius = bend->fieldInnerRadius;
	fieldOuterRadius = bend->fieldOuterRadius;
	fieldHeight = bend->fieldHeight;
	fringeDepth = bend->fringeFactor * fieldHeight;
	
	// Get fringe field region, decide fringeMaxZ . Only consider Dipole fringe as it propagates the furtherest//////////////
	if(bend->fringe == 0) fringeMaxZ=0;
	else{
		for(int i=0; i<1000; ++i) {
			fringeMaxZ = i*fieldHeight/10.0;
			if(enge_Dipole(fringeMaxZ/fringeDepth) < FRINGE_ACCURACY)
				break;
		}
	}
	
	////////////////////////////////////////
	
	rotation = global2local.getRotation().inverse();
	tanangle = tan(bend->angle);
	BendAngle = bend->angle;
	
	// set global bounding box
/*	G4double local[4], global[4];
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
		local[1] = fieldHeight/2.;
		local[2] = r1*sin(phi);
		global2local.getGlobal(local,global);
		setGlobalPoint(global);
		local[1] = -fieldHeight/2.;
		global2local.getGlobal(local,global);
		setGlobalPoint(global);
		local[0] = r2*cos(phi);
		local[1] = fieldHeight/2.;
		local[2] = r2*sin(phi);
		global2local.getGlobal(local,global);
		setGlobalPoint(global);
		local[1] = -fieldHeight/2.;
		global2local.getGlobal(local,global);
		setGlobalPoint(global);
	}
*/
}

void GenericSectorBendField::addFieldValue(const G4double point[4], G4double field[6]) const
{
	G4ThreeVector global(point[0],point[1],point[2]);
	G4ThreeVector local;
	global2local.getLocal(local,global);
	G4double r2 = local[0]*local[0] + local[2]*local[2];
	
	G4ThreeVector B(0.0,0.0,0.0);
	
	G4double bs=*SolenoidalField;	//Solenoid
	G4double b0=*DipoleField;	//Dipole
	G4double b1=*QuadrupoleField;	//Quadrupole
	G4double a1=*SkewQuadrupoleField;	//Skew Quadrupole
	G4double b2=*SextupoleField;	//Sextupole
	G4double a2=*SkewSextupoleField;	//Sextupole rotated
	G4double b3=*OctupoleField;	//Octupole
	G4double a3=*SkewOctupoleField;	//Octupole rotated
	
	G4double bs_prime=0, bs_second=0, bs_third=0;
	G4double b0_prime=0, b0_second=0, b0_third=0;
	G4double a1_prime=0, a1_second=0, a1_third=0;
	G4double b1_prime=0, b1_second=0, b1_third=0;
	G4double a2_prime=0, a2_second=0, a2_third=0;
	G4double b2_prime=0, b2_second=0, b2_third=0;
	G4double a3_prime=0, a3_second=0, a3_third=0;
	G4double b3_prime=0, b3_second=0, b3_third=0;
	
	/// Get positions in centerline coordinate /////////
	
	G4double K=1./fieldRadius;
	//	if(K==-1) K=0;
	
	G4double alpha=0,x=0,y=0,s=0;
	
	// alpha is the angular position in the field region ////
	if (tanangle>0) {
		alpha=-atan(local[2]/local[0]);
		
	}else{
		alpha=atan(local[2]/local[0]);
	}
	
	// define x and s in three different region; x+ is towards the bending center ///
	if(alpha>=0 && tan(alpha)<=fabs(tanangle)){
		x=fieldRadius-sqrt(local[0]*local[0]+local[2]*local[2]);
		s=alpha*fieldRadius;
	}
	if(alpha<0) {
		x=fieldRadius-fabs(local[0]);
		s=local[2];
	}
	if(alpha>fabs(BendAngle)) {
		x=fieldRadius-fabs(sin(BendAngle))*(local[2]-fabs(local[0]*tanangle))-fabs(local[0]/cos(BendAngle));
		s=(local[2]-fabs(local[0]*tanangle))*cos(BendAngle)+fieldRadius*BendAngle;
	}
	y=local[1];
	
	/// Boundary of field ////////////
	if(alpha>=0 && alpha<=fabs(BendAngle) && (sqrt(local[0]*local[0]+local[2]*local[2])<fieldInnerRadius||sqrt(local[0]*local[0]+local[2]*local[2])>fieldOuterRadius)) return;
	
	if(alpha<0 && local[2]>-fringeMaxZ && (fabs(local[0])<fieldInnerRadius ||fabs(local[0])>fieldOuterRadius)) return;
	
	if(alpha<0 && local[2]<-fringeMaxZ) return;
	
	if(alpha>fabs(BendAngle) && (local[2]-fabs(local[0]*tanangle))*cos(BendAngle)<fringeMaxZ && (x>fieldRadius-fieldInnerRadius || x<fieldRadius-fieldOuterRadius || fabs(y)>fieldHeight)) return;
	
	if(alpha>fabs(BendAngle) && (local[2]-fabs(local[0]*tanangle))*cos(BendAngle)>fringeMaxZ ) return;

	/// Enable Enge factor ///////////
	G4double FringZD = 0;
	if(fringeMaxZ > 0.){

		if(alpha>fabs(BendAngle/2.)) FringZD = (s-fieldRadius*BendAngle)/fringeDepth;
		else FringZD = -s/fringeDepth;

		b0 = b0*enge_Dipole(FringZD);
		b0_prime = b0*enge_Dipole.prime(FringZD)/fringeDepth;
		b0_second = b0*enge_Dipole.second(FringZD)/fringeDepth/fringeDepth;
		b0_third = b0*enge_Dipole.third(FringZD)/fringeDepth/fringeDepth/fringeDepth;
		b1 = b1*enge_Quad(FringZD);
		b1_prime = b1*enge_Quad.prime(FringZD)/fringeDepth;
		b1_second = b1*enge_Quad.second(FringZD)/fringeDepth/fringeDepth;
		b1_third = b1*enge_Quad.third(FringZD)/fringeDepth/fringeDepth/fringeDepth;
		b2 = b2*enge_Sext(FringZD);
		b2_prime = b2*enge_Sext.prime(FringZD)/fringeDepth;
		b2_second = b2*enge_Sext.second(FringZD)/fringeDepth/fringeDepth;
		b2_third = b2*enge_Sext.third(FringZD)/fringeDepth/fringeDepth/fringeDepth;
		b3 = b3*enge_Sext(FringZD);
		b3_prime = b3*enge_Sext.prime(FringZD)/fringeDepth;
		b3_second = b3*enge_Sext.second(FringZD)/fringeDepth/fringeDepth;
		b3_third = b3*enge_Sext.third(FringZD)/fringeDepth/fringeDepth/fringeDepth;
		a1 = a1*enge_Quad(FringZD);
		a1_prime = a1*enge_Quad.prime(FringZD)/fringeDepth;
		a1_second = a1*enge_Quad.second(FringZD)/fringeDepth/fringeDepth;
		a1_third = a1*enge_Quad.third(FringZD)/fringeDepth/fringeDepth/fringeDepth;
		a2 = a2*enge_Sext(FringZD);
		a2_prime = a2*enge_Sext.prime(FringZD)/fringeDepth;
		a2_second = a2*enge_Sext.second(FringZD)/fringeDepth/fringeDepth;
		a2_third = a2*enge_Sext.third(FringZD)/fringeDepth/fringeDepth/fringeDepth;
		a3 = a3*enge_Sext(FringZD);
		a3_prime = a3*enge_Sext.prime(FringZD)/fringeDepth;
		a3_second = a3*enge_Sext.second(FringZD)/fringeDepth/fringeDepth;
		a3_third = a3*enge_Sext.third(FringZD)/fringeDepth/fringeDepth/fringeDepth;
	}
	G4double Bx=a1*x+b1*y+a2*x*x+2*b2*x*y-0.5*(2*a2+K*(a1-2*bs_prime))*y*y+a3*x*x*x+3*b3*x*x*y-0.5*(6*a3+a1_second+2*K*a2-2*K*K*(a1-3*bs_prime))*x*y*y-1./6.*(6*b3+b1_second+2*K*(b2-b0_second)-K*K*b1)*y*y*y;
	G4double By=b0+b1*x-(a1+bs_prime)*y+b2*x*x-(2*a2+K*(a1-2*bs_prime))*x*y-0.5*(2*b2+b0_second+K*b1)*y*y+b3*x*x*x-0.5*(6*a3+a1_second+2*K*a2-2*K*K*(a1-3*bs_prime))*x*x*y-0.5*(6*b3+b1_second+2*K*(b2-b0_second)-K*K*b1)*x*y*y+1./6.*(6*a3+2*a1_second+bs_third+K*4*a2-K*K*(a1-4*bs_prime))*y*y*y;
	G4double Bs=bs-K*bs*x+b0_prime*y+0.5*(a1_prime+2*K*K*bs)*x*x+(b1_prime-K*b0_prime)*x*y-0.5*(a1_prime+bs_second)*y*y+1./6.*(2*a2_prime-3*K*a1_prime-6*K*K*K*bs)*x*x*x+(b2_prime-K*b1_prime+K*K*b0_prime)*x*x*y-0.5*(2*a2_prime-3*K*bs_second)*x*y*y-1./6.*(2*b2_prime+b0_third+K*b1_prime)*y*y*y;

	if(tanangle>0)	B[0]=Bs*sin(alpha)+Bx*cos(alpha);
	else B[0]=-Bs*sin(alpha)-Bx*cos(alpha);
	B[1]=By;
	B[2]=Bs*cos(alpha)-Bx*sin(alpha);
	if(alpha<fabs(BendAngle/2.)) B[2] = -B[2];

//	if(fmod(s,1)<0.01)
	//	G4cout<<alpha<<";  l0:"<<local[0]<<";  l1:"<<local[1]<<";  l2:"<<local[2]<<";  x:"<<x<<";  y:"<<y<<";  s:"<<s<<";   Bx:"<<Bx<<";  By:"<<By<<";  Bz:"<<Bs<<G4endl;
	//	G4cout<<alpha<<";  l0:"<<local[0]<<";  l1:"<<local[1]<<";  l2:"<<local[2]<<";   x:"<<x<<";  y:"<<y<<";  s:"<<s<<"; By:"<<B[1]<<";  Fzd:"<<FringZD<<";  g0:"<<global[0]<<";  g1:"<<global[1]<<";  g2:"<<global[2]<<G4endl;
	
	//	G4cout<<x<<" "<<y<<" "<<s<<" "<<alpha<<" "<<Bx<<" "<<By<<" "<<Bs<<G4endl;
	
	if(global2local.isRotated())
		B = rotation * B;
	field[0] += B[0];
	field[1] += B[1];
	field[2] += B[2];
}
