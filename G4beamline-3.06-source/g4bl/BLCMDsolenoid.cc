//	BLCMDsolenoid.cc
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

#ifndef G4BL_GSL
int BLCMDsolenoid_dummy=0;
#else

#include <vector>

#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLElement.hh"
#include "BLManager.hh"
#include "BLElementField.hh"
#include "BLGlobalField.hh"
#include "BLCoil.hh"
#include "BLKillTrack.hh"

/**	SolenoidField represents one placement of a solenoid.
 *
 **/
class SolenoidField : public BLElementField {
	BLCoil *coil;
	G4double current;
	BLCoordinateTransform global2local;
	G4RotationMatrix *rotation;
public:
	/// constructor. _zmin/max are global coordinates.
	SolenoidField(BLCoordinateTransform& _global2local, BLCoil *_coil,
		G4double _current);

	/// addFieldValue() adds the field for this solenoid into field[].
	/// Calls coil->addField() after converting to relative coords.
	/// point[] is in global coordinates.
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};

/**	BLCMDsolenoid implements a solenoid.
 *
 *	A solenoid consists of a BLCoil and a current value. It is a
 *	BLElement and can be placed geometrically (multiple times). The
 *	geometrical parameters are in the BLCoil, not the BLCMDsolenoid;
 *	this class uses the variables of its BLCoil to determine the
 *	geometry of the solenoid.
 **/
class BLCMDsolenoid : public BLElement {
	G4String coilName;
	G4double current;
	G4String color;
	G4int alternate;
	G4int kill;
	BLCoil *coil;
	std::vector<SolenoidField *> solenoidField;
	friend class SolenoidField;
public:
	/// Default constructor. Defines the command, etc.
	BLCMDsolenoid();

	/// Destructor.
	virtual ~BLCMDsolenoid() { }

	/// Copy constructor.
	BLCMDsolenoid(const BLCMDsolenoid& r);

	/// clone()
	BLElement *clone() { return new BLCMDsolenoid(*this); }

	/// commandName() returns "solenoid".
	G4String commandName() { return "solenoid"; }

	/// command() implements the solenoid command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments of the command.
	void defineNamedArgs();

	/// construct() - construct the solenoid.
	/// Creates a new SolenoidField and adds it to BLGlobalField.
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the solenoid
	G4double getLength() { return coil->length; }

	/// getWidth() returns the outer radius of the solenoid
	G4double getWidth() { return coil->outerRadius*2.0; }

	/// getHeight() returns the outer radius of the solenoid
	G4double getHeight() { return coil->outerRadius*2.0; }

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

BLCMDsolenoid defaultSolenoid;	// default object


// Default constructor - be sure to use the default constructor BLElement()
BLCMDsolenoid::BLCMDsolenoid() : BLElement(), solenoidField()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("defines a solenoid (a coil and current)");
	setDescription("A solenoid is a coil and a current. If alternate is\n"
		"nonzero, then each placement of the solenoid (or an\n"
		"enclosing group) will flip the sign of current.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.");

	// provide initial values for fields
	coilName = "";
	current = 0.0;
	color = "1,1,1";
	alternate = 0;
	kill = 0;
	coil = 0;
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDsolenoid::BLCMDsolenoid(const BLCMDsolenoid& r) : BLElement(r), solenoidField()
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	coilName = r.coilName;
	current = r.current;
	color = r.color;
	alternate = r.alternate;
	kill = r.kill;
	coil = r.coil;
}

int BLCMDsolenoid::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("solenoid: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultSolenoid.handleNamedArgs(namedArgs);
	}

	BLCMDsolenoid *t = new BLCMDsolenoid(defaultSolenoid);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	t->coil = BLCoil::find(t->coilName);
	if(t->coil == 0) {
		printError("solenoid: coil '%s' not found.",
							t->coilName.c_str());
		retval = -1;
	}

	t->print(argv[0]);

	return retval;
}

void BLCMDsolenoid::defineNamedArgs()
{
	argString(coilName,"coilName","The name of the coil (must exist)",false);
	// the map was generated for 1 Amp/mm^2, so this factor is actually unitless
	argDouble(current,"current","The current density in the conductor (Amp/mm^2)",1.0);
	argString(color,"color","The color of the solenoid (''=invisible).");
	argInt(alternate,"alternate","Set nonzero to alternate sign each placement.");
	argInt(kill,"kill","Set nonzero to kill all tracks that hit the coil.");
	argString(coilName,"coil","Synonym for coilName.",false);
}

void BLCMDsolenoid::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

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

	// add this solenoid to the GlobalField
	SolenoidField *sf = new SolenoidField(global2local,coil, current);
	BLGlobalField::getObject()->addElementField(sf);
	solenoidField.push_back(sf);

	if(alternate)
		current = -current;

	if(coil->length == 0.0)
		return;

	// construct the coil
	G4String thisname = parentName+getName();
	G4Tubs *tubs = new G4Tubs(thisname+"Tubs", coil->innerRadius, 
			coil->outerRadius, coil->length/2.0, 0.0,2.0*pi);
	G4Material *mat = getMaterial(coil->material);
	G4LogicalVolume *lv = new G4LogicalVolume(tubs,mat, thisname+"LogVol");
	lv->SetVisAttributes(getVisAttrib(color));

	G4PVPlacement *pv = new G4PVPlacement(g4rot,relativePosition,lv,
					thisname, parent,false,0,surfaceCheck);

	if(kill)
		BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));

	printf("BLCMDsolenoid::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		global2local.getPosition()[2]);
	printf("\tglobal pos=%.1f,%.1f,%.1f  ",globalPosition[0],
				globalPosition[1],globalPosition[2]);
}

SolenoidField::SolenoidField( BLCoordinateTransform& _global2local, 
		BLCoil *_coil, G4double _current) : BLElementField() 
{
	coil = _coil;
	current = _current;
	global2local = _global2local;
	rotation = 0;

	if(global2local.isRotated()) {
		rotation = new G4RotationMatrix(global2local.getRotation());
		rotation->invert();
	}

	// set global bounding box
	G4double local[4], global[4];
	local[3] = 0.0;
	for(int i=0; i<2; ++i) {
		local[0] = (i==0 ? -1.0 : 1.0) * coil->maxR;
		for(int j=0; j<2; ++j) {
			local[1] = (j==0 ? -1.0 : 1.0) * coil->maxR;
			for(int k=0; k<2; ++k) {
				local[2] = (k==0 ? coil->minZ : coil->maxZ);
				global2local.getGlobal(local,global);
				setGlobalPoint(global);
			}
		}
	}
}

void SolenoidField::addFieldValue(const G4double point[4], G4double field[6]) const
{
	G4double relPoint[4];
	global2local.getLocal(relPoint,point);
	if(rotation) {
		G4double f[6];
		f[0] = f[1] = f[2] = f[3] = f[4] = f[5] = 0.0;
		coil->addField(relPoint,f,current);
		G4ThreeVector B(f[0],f[1],f[2]);
		B = *rotation * B;
		field[0] += B[0];
		field[1] += B[1];
		field[2] += B[2];
	} else {
		coil->addField(relPoint,field,current);
	}
}

void BLCMDsolenoid::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	v.clear();
	if(coil->length == 0) return;
	generateTubs(npoints, coil->innerRadius, coil->outerRadius, 0.0, 
			360.0*deg, coil->length, v);
}

G4bool BLCMDsolenoid::isOutside(G4ThreeVector &local, G4double tolerance)
{
	if(coil->length == 0) return true;
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	return r < coil->innerRadius+tolerance ||
	        r > coil->outerRadius-tolerance ||
		fabs(local[2]) > coil->length/2.0-tolerance;
}
#endif // G4BL_GSL
