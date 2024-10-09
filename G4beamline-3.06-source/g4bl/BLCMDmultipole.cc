//	BLCMDmultipole.cc
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


/**	BLCMDmultipole implements a generic multipole magnet with a cylindrical
 *	field.
 *
 *	The magnetic field is multipole in x and y, uniform along the z axis,
 *	modified by the fringe field computation described below.
 *	There's no field in the iron, so setting kill=1 is a good idea.
 **/
class BLCMDmultipole : public BLElement {
	G4double fieldLength;
	G4double ironLength;
	G4double ironRadius;
	G4double apertureRadius;
	G4String ironMaterial;
	G4String fieldMaterial;
	G4double dipole;
	G4double quadrupole;
	G4double sextupole;
	G4double octopole;
	G4double decapole;
	G4double dodecapole;
	G4String ironColor;
	G4int kill;
	G4double maxStep;
	G4String fringe;
	G4double fringeFactor;
	G4int openAperture;
	friend class MultipoleField;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDmultipole();

	/// Destructor.
	virtual ~BLCMDmultipole() { }

	/// Copy constructor.
	BLCMDmultipole(const BLCMDmultipole& r);

	/// clone()
	BLElement *clone() { return new BLCMDmultipole(*this); }

	/// commandName() returns "multipole".
	G4String commandName() { return "multipole"; }

	/// command() implements the multipole command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() - construct the multipole magnet
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the ironLength of the multiple magnet
	G4double getLength() 
		{ return (ironLength>fieldLength ? ironLength : fieldLength); }

	/// getWidth() returns the outer radius of the multiple magnet
	G4double getWidth() { return ironRadius*2.0; }

	/// getHeight() returns the outer radius of the multiple magnet
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

BLCMDmultipole defaultMultipole;	// default object

/**	MultipoleField represents one placement of a multipole magnet.
 *
 **/
class MultipoleField : public BLElementField {
	G4double fieldRadius;
	G4double halflength;
	G4double gradient;
	BLCoordinateTransform global2local;
	G4RotationMatrix rotation;
	G4double fringeMaxZ;
	G4double fringeDepth;
	BLEngeFunction enge;
	BLCMDmultipole *mp;
public:
	/// constructor. 
	MultipoleField(BLCoordinateTransform& _global2local, BLCMDmultipole *mp);

	/// addFieldValue() adds the field for this solenoid into field[].
	/// point[] is in global coordinates.
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};


// Default constructor - be sure to use the default constructor BLElement()
BLCMDmultipole::BLCMDmultipole() : BLElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a generic multipole magnet.");
	setDescription("Multipole magnetic fields from dipole through dodecapole "
		"are implemented with a cylindrical field region and "
		"optional surrounding iron (ironLength=0 or ironRadius=0 "
		"omits it). "
		"All fields with positive strengths are oriented so in the "
		"X-Z plane for X>0 (beam left) the field is purely By. "
		"Negative strengths are allowed and reverse the field. "
		"The fringe field computation is not implemented.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.");

	// provide initial values for fields
	fieldLength = 0.0;
	ironLength = 0.0;
	ironRadius = 0.0;
	apertureRadius = 0.0;
	ironMaterial = "Fe";
	fieldMaterial = "Vacuum";
	dipole = 0.0;
	quadrupole = 0.0;
	sextupole = 0.0;
	octopole = 0.0;
	decapole = 0.0;
	dodecapole = 0.0;
	ironColor = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	fringe = "";
	fringeFactor = 1.0;
	openAperture = 0;
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDmultipole::BLCMDmultipole(const BLCMDmultipole& r) : BLElement(r) 
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	fieldLength = r.fieldLength;
	ironLength = r.ironLength;
	ironRadius = r.ironRadius;
	apertureRadius = r.apertureRadius;
	ironMaterial = r.ironMaterial;
	fieldMaterial = r.fieldMaterial;
	dipole = r.dipole;
	quadrupole = r.quadrupole;
	sextupole = r.sextupole;
	octopole = r.octopole;
	decapole = r.decapole;
	dodecapole = r.dodecapole;
	ironColor = r.ironColor;
	kill = r.kill;
	maxStep = r.maxStep;
	fringe = r.fringe;
	fringeFactor = r.fringeFactor;
	openAperture = r.openAperture;
}

int BLCMDmultipole::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("multipole: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultMultipole.handleNamedArgs(namedArgs);
	}

	BLCMDmultipole *t = new BLCMDmultipole(defaultMultipole);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	getMaterial(t->fieldMaterial);
	getMaterial(t->ironMaterial);

	t->print(argv[0]);

	return retval;
}

void BLCMDmultipole::defineNamedArgs()
{
	argDouble(fieldLength,"fieldLength","The length of the field region (mm)",mm);
	argDouble(ironLength,"ironLength","The length of the iron (mm)",mm);
	argDouble(ironRadius,"ironRadius","The outer radius of the iron (mm)",mm);
	argDouble(apertureRadius,"apertureRadius","The radius of the aperture (mm)",mm);
	argString(ironMaterial,"ironMaterial","The material of the iron region.");
	argString(fieldMaterial,"fieldMaterial","The material of the field region.");
	argDouble(dipole,"dipole","Strength of dipole (Tesla)",tesla);
	argDouble(quadrupole,"quadrupole","Strength of quadrupole (T/m)",tesla/meter);
	argDouble(sextupole,"sextupole","Strength of sextupole (T/m^2)",tesla/meter/meter);
	argDouble(octopole,"octopole","Strength of octopole (T/m^3)",tesla/meter/meter/meter);
	argDouble(decapole,"decapole","Strength of decapole (T/m^4)",tesla/meter/meter/meter/meter);
	argDouble(dodecapole,"dodecapole","Strength of dodecapole (T/m^5)",tesla/meter/meter/meter/meter/meter);
	argString(ironColor,"ironColor","The color of the iron region.");
	argInt(kill,"kill","Set nonzero to kill tracks hitting the iron.");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)");
	argString(fringe,"fringe","Fringe field computation, set to 0 to disable");
	argDouble(fringeFactor,"fringeFactor","Fringe depth factor (1.0).");
	argInt(openAperture,"openAperture","Set nonzero to omit the aperture volume.",false);
}

void BLCMDmultipole::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)

{
	G4String thisname = parentName+getName();

	G4LogicalVolume *lviron=0;
	G4LogicalVolume *lvfield=0;

	if(ironLength > 0.0 && ironRadius > 0.0) {
		G4Material *ironMat = getMaterial(ironMaterial);
		G4Material *fieldMat = getMaterial(fieldMaterial);
		G4VSolid *ironSolid = new G4Tubs(thisname, apertureRadius,
				ironRadius,ironLength/2.0,0.0,360.0*deg);
		G4VSolid *fieldSolid = new G4Tubs(thisname, 0,
				apertureRadius,ironLength/2.0,0.0,360.0*deg);
		lviron = new G4LogicalVolume(ironSolid,ironMat,thisname);
		lviron->SetVisAttributes(getVisAttrib(ironColor));
		lviron->SetUserLimits(new G4UserLimits(kill ? 0.1 :maxStep));
		if(openAperture == 0) {
		    lvfield = new G4LogicalVolume(fieldSolid,fieldMat,thisname);
		    lvfield->SetVisAttributes(
		    			BLCommand::getVisAttrib("Invisible"));
		    lvfield->SetUserLimits(new G4UserLimits(maxStep));
		}
	}

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4PVPlacement *pv=0;
	if(lviron) pv = new G4PVPlacement(g4rot,relativePosition,lviron,
					thisname,parent,false,0,surfaceCheck);
	if(lvfield) new G4PVPlacement(g4rot,relativePosition,lvfield,
					thisname,parent,false,0,surfaceCheck);

	if(kill && pv)
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

	BLCoordinateTransform global2local(globalRotation,globalPosition);

	MultipoleField *p = new MultipoleField(global2local,this);
	BLGlobalField::getObject()->addElementField(p);

	printf("BLCMDmultipole::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2]);
}


MultipoleField::MultipoleField(BLCoordinateTransform& _global2local,
	BLCMDmultipole *_mp) : BLElementField(), rotation(), enge(ENGE_QUAD)
{
	mp = _mp;
	fieldRadius = mp->apertureRadius;
	halflength = mp->fieldLength/2.0;
	global2local = _global2local;
	rotation = global2local.getRotation().inverse();
	fringeDepth = mp->fringeFactor * fieldRadius * 2.0;

	// @@@ fringe field not supported
	if(mp->fringe != "" && mp->fringe != "0") {
		BLCommand::printError("Invalid fringe value\n");
	}
	enge.set(0,0,0,0,0,0);
	fringeMaxZ = halflength;

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

void MultipoleField::addFieldValue(const G4double point[4], G4double field[6])
								const
{
	G4ThreeVector global(point[0],point[1],point[2]);
	G4ThreeVector local;

	global2local.getLocal(local,global);
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	if(r > fieldRadius || fabs(local[2]) > fringeMaxZ)
		return;

	G4ThreeVector B(0.0,0.0,0.0);
	G4double phi = atan2(local[1],local[0]);

	if(mp->dipole != 0.0) {
		B[1] += mp->dipole;
	}
	if(mp->quadrupole != 0.0) {
		B[0] += mp->quadrupole * r * sin(phi);
		B[1] += mp->quadrupole * r * cos(phi);
	}
	if(mp->sextupole != 0.0) {
		B[0] += mp->sextupole * r*r * sin(2.0*phi);
		B[1] += mp->sextupole * r*r * cos(2.0*phi);
	}
	if(mp->octopole != 0.0) {
		B[0] += mp->octopole * r*r*r * sin(3.0*phi);
		B[1] += mp->octopole * r*r*r * cos(3.0*phi);
	}
	if(mp->decapole != 0.0) {
		B[0] += mp->decapole * r*r*r*r * sin(4.0*phi);
		B[1] += mp->decapole * r*r*r*r * cos(4.0*phi);
	}
	if(mp->dodecapole != 0.0) {
		B[0] += mp->dodecapole * r*r*r*r*r * sin(5.0*phi);
		B[1] += mp->dodecapole * r*r*r*r*r * cos(5.0*phi);
	}

	if(global2local.isRotated())
		B = rotation * B;

	field[0] += B[0];
	field[1] += B[1];
	field[2] += B[2];
}

void BLCMDmultipole::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	generateTubs(npoints, apertureRadius, ironRadius, 0.0, 360.0*deg,
			ironLength, v);
}

G4bool BLCMDmultipole::isOutside(G4ThreeVector &local, G4double tolerance)
{
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	return r < apertureRadius+tolerance || r > ironRadius-tolerance ||
		fabs(local[2]) > ironLength/2.0-tolerance;
}
