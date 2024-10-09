//	BLCMDgenericquad.cc
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
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLElement.hh"
#include "BLElementField.hh"
#include "BLGlobalField.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLEngeFunction.hh"
#include "BLKillTrack.hh"

const G4double FRINGE_ACCURACY=1.0e-4;


/**	BLCMDgenericquad implements an generic quadrupole magnet with a cylindrical
 *	field.
 *
 *	The magnetic field is quadrupole in x and y, uniform along the z axis,
 *	modified by the fringe field computation described below.
 *	There's no field in the iron, so setting kill=1 is a good idea.
 *
 *	Set apertureRadius to get a circular aperture.
 *	Set poleTipRadius, coilRadius, and coilHalfwidth to get a "rounded
 *	+ sign" aperture, using circles for the pole tips.
 **/
class BLCMDgenericquad : public BLElement {
	G4double fieldLength;
	G4double ironLength;
	G4double ironRadius;
	G4double apertureRadius;
	G4double poleTipRadius;
	G4double coilRadius;
	G4double coilHalfwidth;
	G4double gradient;
	G4String ironMaterial;
	G4String fieldMaterial;
	G4String ironColor;
	G4int kill;
	G4double maxStep;
	G4double fieldRadius;
	G4String fringe;
	G4double fringeFactor;
	G4int openAperture;
	friend class GenericQuadField;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDgenericquad();

	/// Destructor.
	virtual ~BLCMDgenericquad() { }

	/// Copy constructor.
	BLCMDgenericquad(const BLCMDgenericquad& r);

	/// clone()
	BLElement *clone() { return new BLCMDgenericquad(*this); }

	/// commandName() returns "genericquad".
	G4String commandName() { return "genericquad"; }

	/// command() implements the genericquad command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() - construct the generic quadrupole magnet
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the ironLength of the quad
	G4double getLength() 
		{ return ironLength>0.0 ? ironLength : fieldLength; }

	/// getWidth() returns the outer radius of the quad
	G4double getWidth() { return ironRadius*2.0; }

	/// getHeight() returns the outer radius of the quad
	G4double getHeight() { return ironRadius*2.0; }

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

BLCMDgenericquad defaultGenericQuad;	// default object


/**	GenericQuadField represents one placement of a generic quadrupole magnet.
 *
 **/
class GenericQuadField : public BLElementField {
	G4double fieldRadius;
	G4double halflength;
	G4double gradient;
	BLCoordinateTransform global2local;
	G4RotationMatrix rotation;
	G4double fringeMaxZ;
	G4double fringeDepth;
	BLEngeFunction enge;
public:
	/// constructor. 
	GenericQuadField(BLCoordinateTransform& _global2local, BLCMDgenericquad *quad);

	/// addFieldValue() adds the field for this solenoid into field[].
	/// point[] is in global coordinates.
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};


// Default constructor - be sure to use the default constructor BLElement()
BLCMDgenericquad::BLCMDgenericquad() : BLElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a generic quadrupole magnet.");
	setDescription("The field region is a tubs with gradient specified.\n"
		"A positive gradient yields a horizontally-focusing\n"
		"quad for positive particles.\n"
		"If apertureRadius>0 the quad has a circular aperture.\n"
		"For a 'rounded +' aperture using circles for the poles,\n"
		"set poleTipRadius, coilRadius, coilHalfwidth.\n"
		"Due to visualization bugs, in the latter case you cannot\n"
		"see through the aperture; it is solid black.\n"
		"A fringe field computation based on the method of COSY\n"
		"INFINITY is included by default, extending the field region.\n"
		"This is first order only, "
		"and the fringe field extends outside of "
		"the magnet aperture only in a cylinder extending the "
		"aperture straight along local z. As the fringe field is first "
		"order only, it is slightly non-Maxwellian. It is computed using "
		"Enge functions.\n\n"
		"Note that there is no field inside the 'iron'; this can "
		"result in gross tracking errors for particles in the iron, "
		"and implies that kill=1 is desirable.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.\n\n"
		"If ironLength <= 0, no iron is constructed.\n\n"
		"Note that section 4.7 of the user's Guide has a dimensioned "
		"drawing of the genericquad aperture.");

	// provide initial values for fields
	fieldLength = 0.0;
	ironLength = 0.0;
	ironRadius = 0.0;
	apertureRadius = 0.0;
	poleTipRadius = 0.0;
	coilRadius = 0.0;
	coilHalfwidth = 0.0;
	gradient = 0.0;
	ironMaterial = "Fe";
	fieldMaterial = "Vacuum";
	ironColor = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	fieldRadius = 0.0;
	fringe = "";
	fringeFactor = 1.0;
	openAperture = 0;
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDgenericquad::BLCMDgenericquad(const BLCMDgenericquad& r) : BLElement(r) 
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	fieldLength = r.fieldLength;
	ironLength = r.ironLength;
	ironRadius = r.ironRadius;
	apertureRadius = r.apertureRadius;
	poleTipRadius = r.poleTipRadius;
	coilRadius = r.coilRadius;
	coilHalfwidth = r.coilHalfwidth;
	gradient = r.gradient;
	ironMaterial = r.ironMaterial;
	fieldMaterial = r.fieldMaterial;
	ironColor = r.ironColor;
	kill = r.kill;
	maxStep = r.maxStep;
	fieldRadius = r.fieldRadius;
	fringe = r.fringe;
	fringeFactor = r.fringeFactor;
	openAperture = r.openAperture;
}

int BLCMDgenericquad::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("genericquad: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultGenericQuad.handleNamedArgs(namedArgs);
	}

	BLCMDgenericquad *t = new BLCMDgenericquad(defaultGenericQuad);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	getMaterial(t->fieldMaterial);
	getMaterial(t->ironMaterial);

	t->print(argv[0]);

	return retval;
}

void BLCMDgenericquad::defineNamedArgs()
{
	argDouble(fieldLength,"fieldLength","The length of the field region (mm)");
	argDouble(ironLength,"ironLength","The length of the iron (mm)");
	argDouble(ironRadius,"ironRadius","The outer radius of the iron (mm)");
	argDouble(apertureRadius,"apertureRadius","The radius of the aperture (mm)");
	argDouble(poleTipRadius,"poleTipRadius","The inner radius of the pole tips (mm).");
	argDouble(coilRadius,"coilRadius","The radius of the inside of the coil (mm).");
	argDouble(coilHalfwidth,"coilHalfwidth","The halfwidth of the coil (mm).");
	argDouble(coilHalfwidth,"coilHalfWidth","Synonym for coilHalfwidth.");
	argDouble(gradient,"gradient","The magnetic field gradient, dBy/dx (Tesla/meter)",tesla/meter);
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)");
	argString(ironMaterial,"ironMaterial","The material of the iron region.");
	argString(fieldMaterial,"fieldMaterial","The material of the field region.");
	argString(ironColor,"ironColor","The color of the iron region.");
	argInt(kill,"kill","Set nonzero to kill tracks hitting the iron.");
	argString(fringe,"fringe","Fringe field computation, set to 0 to disable,"
		" or a comma-separated list of 6 Enge function parameters.");
	argDouble(fringeFactor,"fringeFactor","Fringe depth factor (1.0).");
	argInt(openAperture,"openAperture","Set nonzero to omit the aperture volume.",false);
}

void BLCMDgenericquad::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)

{
	G4String thisname = parentName+getName();

	G4LogicalVolume *lviron=0;
	G4LogicalVolume *lvfield=0;

	if(ironLength > 0.0) {
	    G4Material *ironMat = getMaterial(ironMaterial);
	    G4Material *fieldMat = getMaterial(fieldMaterial);
	    if(apertureRadius > 0.0) {
		G4VSolid *ironSolid = new G4Tubs(thisname, apertureRadius,
				ironRadius,ironLength/2.0,0.0,360.0*deg);
		lviron = new G4LogicalVolume(ironSolid,ironMat,thisname);
		lviron->SetVisAttributes(getVisAttrib(ironColor));
		lviron->SetUserLimits(new G4UserLimits(kill ? 0.1 :maxStep));
		if(openAperture == 0) {
			G4VSolid *fieldSolid = new G4Tubs(thisname, 0,
				apertureRadius,ironLength/2.0,0.0,360.0*deg);
			lvfield = new G4LogicalVolume(fieldSolid,fieldMat,
								thisname);
			lvfield->SetVisAttributes(
					BLCommand::getVisAttrib("Invisible"));
			lvfield->SetUserLimits(new G4UserLimits(maxStep));
		}
		fieldRadius = apertureRadius;
	    } else {
		G4VSolid *box=0;
		G4VSolid *ironSolid = new G4Tubs(thisname, 0.0,
				ironRadius,ironLength/2.0,0.0,360.0*deg);
		box = new G4Box(thisname,coilRadius,
					coilRadius, ironLength/2.0+0.02*mm);
		// compute circle approximating the hyperbolic pole tip --
		// center is at (x,x) on the diagonal; (x0,y0) is the pole tip,
		// (x1,y1) is the intersection of pole face with coil.
		G4double x0=poleTipRadius/sqrt(2.0);
		G4double y0=x0;
		G4double x1=coilRadius;
		G4double y1=coilHalfwidth;
		G4double x=0.5*(x1*x1+y1*y1-x0*x0-y0*y0)/(x1+y1-x0-y0);
		G4double r=sqrt(2.0)*x-poleTipRadius;
		G4Tubs *t = new G4Tubs("xxx",0.0,r,ironLength/2.0+0.1,
						0.0,360.0*deg);
		// Subtract 4 copies of t from box, to get the "rounded +"
		G4ThreeVector offset(x,x,0.0);
		G4RotationMatrix norot;
		box = new G4SubtractionSolid(thisname,box,
				t,&norot,offset);
		offset[0] = -x;
		box = new G4SubtractionSolid(thisname,box,
				t,&norot,offset);
		offset[1] = -x;
		box = new G4SubtractionSolid(thisname,box,
				t,&norot,offset);
		offset[0] = x;
		box = new G4SubtractionSolid(thisname,box,
				t,&norot,offset);
		lviron = new G4LogicalVolume(ironSolid,ironMat,thisname);
		lviron->SetVisAttributes(getVisAttrib(ironColor));
		lviron->SetUserLimits(new G4UserLimits(kill ? 0.1 : maxStep));
		// place the "rounded +" as daughter in the iron Tubs; black.
		G4LogicalVolume *lvbox = new G4LogicalVolume(box,fieldMat,
					thisname);
		lvbox->SetVisAttributes(getVisAttrib("0,0,0"));
		lvbox->SetUserLimits(new G4UserLimits(maxStep));
		G4ThreeVector no_offset(0.0,0.0,0.0);
		new G4PVPlacement(0,no_offset,lvbox,thisname,lviron,false,0,
								surfaceCheck);
		fieldRadius = 
			sqrt(coilRadius*coilRadius+coilHalfwidth*coilHalfwidth);
		lvfield = 0;
		// set apertureRadius for fringe field computation
		apertureRadius = poleTipRadius;
	    }
	} else {
		fieldRadius = apertureRadius;
	}

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	if(lviron) {
	    G4PVPlacement *pv = new G4PVPlacement(g4rot,relativePosition,lviron,
					thisname,parent,false,0,surfaceCheck);
	    if(lvfield) new G4PVPlacement(g4rot,relativePosition,lvfield,
					thisname,parent,false,0,surfaceCheck);

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
	G4ThreeVector globalPosition(relativePosition + parentPosition);
	if(parentRotation)
		globalPosition = *parentRotation * relativePosition +
				parentPosition;

	BLCoordinateTransform global2local(globalRotation,globalPosition);

	GenericQuadField *p = new GenericQuadField(global2local,this);
	BLGlobalField::getObject()->addElementField(p);

	printf("BLCMDgenericquad::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2]);
}

GenericQuadField::GenericQuadField(BLCoordinateTransform& _global2local,
	BLCMDgenericquad *quad) : BLElementField(), rotation(), enge(ENGE_QUAD)
{
	fieldRadius = quad->fieldRadius;
	halflength = quad->fieldLength/2.0;
	gradient = quad->gradient;
	global2local = _global2local;
	rotation = global2local.getRotation().inverse();
	fringeDepth = quad->fringeFactor * quad->apertureRadius * 2.0;

	if(quad->fringe != "") {
		std::vector<G4double> v = BLCommand::getList(quad->fringe,',');
		if(v.size() == 1 && v[0] == 0.0)
			enge.set(0,0,0,0,0,0);
		else if(v.size() == 6)
			enge.set(v[0],v[1],v[2],v[3],v[4],v[5]);
		else
			BLCommand::printError("Invalid fringe value '%s'\n",
						quad->fringe.c_str());
	}

	if(enge.getType() == ENGE_BLOCK) {
		fringeMaxZ = halflength;
	} else {
		for(int i=0; i<1000; ++i) {
			fringeMaxZ = i*fieldRadius/10.0 + halflength;
			if(enge((fringeMaxZ-halflength)/fringeDepth) <
							FRINGE_ACCURACY)
				break;
		}
	}

	// set global bounding box
	G4double local[4], global[4];
	local[3] = 0.0;
	for(int i=0; i<2; ++i) {
		local[0] = (i==0 ? -1.0 : 1.0) * fieldRadius;
		for(int j=0; j<2; ++j) {
			local[1] = (j==0 ? -1.0 : 1.0) * fieldRadius;
			for(int k=0; k<2; ++k) {
				local[2] = (k==0 ? -1.0 : 1.0) * fringeMaxZ;
				global2local.getGlobal(local,global);
				setGlobalPoint(global);
			}
		}
	}
}

void GenericQuadField::addFieldValue(const G4double point[4], G4double field[6])
								const
{
	G4ThreeVector global(point[0],point[1],point[2]);
	G4ThreeVector local;

	global2local.getLocal(local,global);
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	if(r > fieldRadius || fabs(local[2]) > fringeMaxZ)
		return;

	/* apply enge() to the scalar potential phi=-G0*x*y*enge(z);
	   B is minus its gradient. Handle both edges properly. */
	G4double G0 = gradient;
	double fringeZ = (fabs(local[2])-halflength)/fringeDepth;
	G4double f = enge(fringeZ);
	G4double fp = enge.prime(fringeZ)/fringeDepth;
	G4ThreeVector B(G0*f*local[1],G0*f*local[0],G0*fp*local[0]*local[1]);
	if(local[2] < 0.0) B[2] = -B[2];

	if(global2local.isRotated())
		B = rotation * B;

	field[0] += B[0];
	field[1] += B[1];
	field[2] += B[2];
}

void BLCMDgenericquad::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	v.clear();
	if(ironLength <= 0.0) return;
	generateTubs(npoints, apertureRadius, ironRadius, 0.0, 360.0*deg,
			ironLength, v);
}

G4bool BLCMDgenericquad::isOutside(G4ThreeVector &local, G4double tolerance)
{
	if(ironLength <= 0.0) return true;
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	return r < apertureRadius+tolerance || r > ironRadius-tolerance ||
		fabs(local[2]) > ironLength/2.0-tolerance;
}
