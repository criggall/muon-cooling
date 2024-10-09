//	BLCMDabsorber.cc
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
#include "G4VSolid.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLGroupElement.hh"
#include "BLParam.hh"
#include "BLWindowShape.hh"

/**	BLCMDabsorber implements an absorber using windows defined by 
 *	BLWindowShape. This is optimized for a Liquid Hydrogen absorber,
 *	including optional safety windows.
 *
 *	Geometry test:
 *	For intersections with siblings and parent, behaves as if it were
 *	the enclosing cylinder of the outer windows. For placing children
 *	inside, behaves as the largest cylinder inside the LH2 (i.e. it
 *	does not extend into the windows).
 **/
class BLCMDabsorber : public BLGroupElement {
	G4String absWindow;
	G4String safetyWindow;
	G4double insideLength;
	G4String absMaterial;
	G4String windowMaterial;
	G4String safetyMaterial;
	G4double safetyDistance;
	G4String color;
	G4double maxStep;
	G4double totalLength;
	G4double outerRadius;
	G4double absTubsLength;
	G4double absTubsRadius;
	BLWindowShape *absWindowShape;
	BLWindowShape *safetyWindowShape;
	G4VSolid *safetyWinSolid;
	G4VSolid *safetyVolSolid;
	G4VSolid *safetyPipeSolid;
	G4VSolid *absWinSolid;
	G4VSolid *absVolSolid;
	G4VSolid *absPipeSolid;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDabsorber();

	/// Destructor.
	virtual ~BLCMDabsorber() { }

	/// Copy constructor.
	BLCMDabsorber(const BLCMDabsorber& r);

	/// clone()
	BLElement *clone() { return new BLCMDabsorber(*this); }

	/// commandName() returns "absorber".
	G4String commandName() { return "absorber"; }

	/// command() implements the absorber command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() - construct the tube/cylinder.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the totalLength of the absorber.
	G4double getLength() { return totalLength; }

	/// getWidth() returns the outer radius of the tube/cylinder.
	G4double getWidth() { return outerRadius*2.0; }

	/// getHeight() returns the outer radius of the tube/cylinder.
	G4double getHeight() { return outerRadius*2.0; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// constructSolids() constructs all necessary G4VSolid-s
	void constructSolids();

	/// generatePoints() from BLElement
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isOutside() from BLElement
	G4bool isOutside(G4ThreeVector &local, G4double tolerance);

	/// isWithin() from BLGroupElement
	bool isWithin(G4ThreeVector &local, G4double tolerance);
};

BLCMDabsorber defaultAbsorber;	// default object

// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDabsorber::BLCMDabsorber() : BLGroupElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct an absorber");
	setDescription("The absorber has two windows with beampipe and an\n"
		"absorber material. Optionally it has an additional two safety windows with\n"
		"beampipe. The WindowShape(s) are read from a file, and they\n"
		"determine the thickness and length of the beampipe(s).\n"
		"For geometry testing, acts like a cylinder enclosing the\n"
		"windows. For placing children, acts like a cylinder inside\n"
		"the central absorber.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.\n\n"
		"Note that section 4.5 of the User's Guide has a dimensioned "
		"drawing of an absorber.");

	// provide initial values for fields
	absWindow = "";
	safetyWindow = "";
	insideLength = 0.0;
	absMaterial = "LH2";
	windowMaterial = "Al";
	safetyMaterial = "Vacuum";
	safetyDistance = 12.0*cm;
	color = "1,1,1";
	maxStep = -1.0;
	totalLength = 0.0;
	outerRadius = 0.0;
	absTubsLength = 0.0;
	absTubsRadius = 0.0;
	absWindowShape = 0;
	safetyWindowShape = 0;
	safetyWinSolid = 0;
	safetyVolSolid = 0;
	safetyPipeSolid = 0;
	absWinSolid = 0;
	absVolSolid = 0;
	absPipeSolid = 0;
}

// Copy constructor - be sure to use the copy constructor BLGroupElement(r)
BLCMDabsorber::BLCMDabsorber(const BLCMDabsorber& r) : BLGroupElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	absWindow = r.absWindow;
	safetyWindow = r.safetyWindow;
	insideLength = r.insideLength;
	absMaterial = r.absMaterial;
	windowMaterial = r.windowMaterial;
	safetyMaterial = r.safetyMaterial;
	safetyDistance = r.safetyDistance;
	color = r.color;
	maxStep = r.maxStep;
	totalLength = r.totalLength;
	outerRadius = r.outerRadius;
	absTubsLength = r.absTubsLength;
	absTubsRadius = r.absTubsRadius;
	absWindowShape = r.absWindowShape;
	safetyWindowShape = r.safetyWindowShape;
	safetyWinSolid = r.safetyWinSolid;
	safetyVolSolid = r.safetyVolSolid;
	safetyPipeSolid = r.safetyPipeSolid;
	absWinSolid = r.absWinSolid;
	absVolSolid = r.absVolSolid;
	absPipeSolid = r.absPipeSolid;
}

int BLCMDabsorber::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("absorber: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultAbsorber.handleNamedArgs(namedArgs);
	}

	BLCMDabsorber *t = new BLCMDabsorber(defaultAbsorber);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// get the WindowShape-s
	if(t->absWindowShape == 0)
		t->absWindowShape = new BLWindowShape(t->absWindow);
	if(t->safetyWindowShape == 0 && t->safetyWindow != "")
		t->safetyWindowShape = new BLWindowShape(t->safetyWindow);

	t->constructSolids();

	// check that materials exist
	getMaterial(t->windowMaterial);
	getMaterial(t->safetyMaterial);
	getMaterial(t->absMaterial);

	t->print(argv[0]);

	return retval;
}

void BLCMDabsorber::defineNamedArgs()
{
	argString(absWindow,"absWindow","The name of the absorber window.",false);
	argString(safetyWindow,"safetyWindow","The name of the safety window.",false);
	argDouble(insideLength,"insideLength","Absorber length inside windows (mm)",mm,"",false);
	argString(absMaterial,"absMaterial","The material of the absorber",false);
	argString(windowMaterial,"windowMaterial","The material of the window(s)",false);
	argString(safetyMaterial,"safetyMaterial","The material inside the safety windows.",false);
	argDouble(safetyDistance,"safetyDistance","Distance between absorber and safety windows(mm)",mm,"",false);
	argString(color,"color","The color of the absorber (''=invisible)");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)");
}

void BLCMDabsorber::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4Material *winMat = getMaterial(windowMaterial);
	G4Material *sftyMat = getMaterial(safetyMaterial);
	G4Material *absMat = getMaterial(absMaterial);
	if(!winMat || !winMat || !absMat) {
		G4Exception("absorber","Missing material",FatalException, "");
	}

	G4UserLimits *limits = new G4UserLimits(maxStep);

	G4LogicalVolume *swL = 0;
	G4LogicalVolume *svL = 0;
	G4LogicalVolume *spL = 0;
	if(safetyWinSolid && safetyVolSolid) {
		swL = new G4LogicalVolume(safetyWinSolid,winMat,
						getName()+"SftyWinLV");
		svL = new G4LogicalVolume(safetyVolSolid,sftyMat,
						getName()+"SftyVolLV");
		spL = new G4LogicalVolume(safetyPipeSolid,winMat,
						getName()+"SftyPipeLV");
		swL->SetUserLimits(limits);
		svL->SetUserLimits(limits);
		spL->SetUserLimits(limits);
	}
	G4LogicalVolume *awL = new G4LogicalVolume(absWinSolid,winMat,
						getName()+"AbsWinLV");
	G4LogicalVolume *avL = new G4LogicalVolume(absVolSolid,absMat,
						getName()+"AbsVolLV");
	G4LogicalVolume *apL = new G4LogicalVolume(absPipeSolid,winMat,
						getName()+"AbsPipeLV");
	awL->SetUserLimits(limits);
	avL->SetUserLimits(limits);
	apL->SetUserLimits(limits);
#ifdef G4BL_VISUAL
	const G4VisAttributes *visible = getVisAttrib(color);
	const G4VisAttributes *vis2 = visible;
	if(color == "0,1,0") vis2 = getVisAttrib("0,0.6,0");
	const G4VisAttributes *invisible = BLCommand::getVisAttrib("Invisible");
	if(swL) swL->SetVisAttributes(color!="" ? vis2 : invisible);
	if(svL) svL->SetVisAttributes(color!="" ? vis2 : invisible);
	if(spL) spL->SetVisAttributes(color!="" ? vis2 : invisible);
	awL->SetVisAttributes(invisible);
	avL->SetVisAttributes(color!="" ? visible : invisible);
	apL->SetVisAttributes(color!="" ? visible : invisible);
#endif

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4ThreeVector loc(relativePosition);
	G4LogicalVolume *pL = parent;
	if(swL) {
		new G4PVPlacement(g4rot,loc,swL,parentName+getName()+"SftyWin",
						pL,false,0,surfaceCheck);
		new G4PVPlacement(g4rot,loc,spL,parentName+getName()+"SftyPipe",
						pL,false,0,surfaceCheck);
		loc = G4ThreeVector(0.0,0.0,0.0);
		g4rot = 0;
		new G4PVPlacement(g4rot,loc,svL,parentName+getName()+"SftyVol",
						swL,false,0,surfaceCheck);
		pL = svL;
	}
	new G4PVPlacement(g4rot,loc,awL,parentName+getName()+"AbsWin",
						pL,false,0,surfaceCheck);
	new G4PVPlacement(g4rot,loc,apL,parentName+getName()+"AbsPipe",
						pL,false,0,surfaceCheck);
	loc = G4ThreeVector(0.0,0.0,0.0);
	new G4PVPlacement(0,loc,avL,parentName+getName()+"AbsVol",
						awL,false,0,surfaceCheck);

	// Childrren need global position and rotation
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

	constructChildren(avL,parentName+getName(),globalRotation,globalPosition);
}

void BLCMDabsorber::constructSolids()
{
	if(absWindowShape && absWinSolid == 0) {
		// allocate arrays for 2 G4Polycone-s (window and volume)
		int n = absWindowShape->r.size();
		G4double z0 = absWindowShape->z[0];
		G4double *r = new G4double[n+n];
		G4double *zw = new G4double[n+n];
		G4double *zv = new G4double[n+n];
		G4double *zero = new G4double[n+n];
		for(int j=0; j<n; ++j) {
			int k = n + n - j - 1;
			zero[j] = zero[k] = 0.0;
			// small flat at end (i.e. no corner at center)
			r[j] = r[k] = (j==0 ? absWindowShape->r[1]/4.0 :
						absWindowShape->r[j]);
			zv[k] = absWindowShape->z[j] - z0 + 0.5*insideLength;
			zv[j] = -zv[k];
			zw[k] = zv[k] + absWindowShape->t[j];
			zw[j] = -zw[k];
		}
		absWinSolid = new G4Polycone(getName()+"AbsWin",
						0.0,360.0*deg,n+n,zw,zero,r);
		absVolSolid = new G4Polycone(getName()+"AbsVol",
						0.0,360.0*deg,n+n,zv,zero,r);
		// flange dimensions become dimensions of the absorber pipe
		G4double pipeLength = insideLength + 
					2.0*(absWindowShape->flangeOutsideZ-z0);
		absPipeSolid = new G4Tubs(getName()+"AbsPipe", 
					absWindowShape->flangeInnerRadius, 
					absWindowShape->flangeOuterRadius,
					pipeLength/2.0, 0.0,360.0*deg);

		totalLength = insideLength + 2.0*absWindowShape->t[0];
		if(totalLength < pipeLength) totalLength = pipeLength;
		outerRadius = absWindowShape->r[n-1];
		if(outerRadius < absWindowShape->flangeOuterRadius)
			outerRadius = absWindowShape->flangeOuterRadius;
		// size of absTubs (for children)
		absTubsRadius = r[n-1];
		absTubsLength = insideLength - 2.0*z0;
	}
	if(safetyWindowShape && safetyWinSolid == 0) {
		// allocate arrays for 2 G4Polycone-s (window and volume)
		int n = safetyWindowShape->r.size();
		G4double z0 = safetyWindowShape->z[0];
		G4double *r = new G4double[n+n];
		G4double *zw = new G4double[n+n];
		G4double *zv = new G4double[n+n];
		G4double *zero = new G4double[n+n];
		for(int j=0; j<n; ++j) {
			int k = n + n - j - 1;
			zero[j] = zero[k] = 0.0;
			// small flat at end (i.e. no corner at center)
			r[j] = r[k] = (j==0 ? safetyWindowShape->r[1]/4.0 :
						safetyWindowShape->r[j]);
			zv[k] = safetyWindowShape->z[j] - z0 + 
					0.5*insideLength + safetyDistance;
			zv[j] = -zv[k];
			zw[k] = zv[k] + safetyWindowShape->t[j];
			zw[j] = -zw[k];
		}
		safetyWinSolid = new G4Polycone(getName()+"SftyWin",
						0.0,360.0*deg,n+n,zw,zero,r);
		safetyVolSolid = new G4Polycone(getName()+"SftyVol",
						0.0,360.0*deg,n+n,zv,zero,r);
		// flange dimensions become dimensions of the absorber pipe
		G4double pipeLength = insideLength + 2.0*safetyDistance +
				2.0*(safetyWindowShape->flangeOutsideZ-z0);
		safetyPipeSolid = new G4Tubs(getName()+"SftyPipe", 
					safetyWindowShape->flangeInnerRadius, 
					safetyWindowShape->flangeOuterRadius,
					pipeLength/2.0, 0.0,360.0*deg);

		totalLength = insideLength + 2.0*safetyWindowShape->t[0];
		if(totalLength < pipeLength) totalLength = pipeLength;
		outerRadius = safetyWindowShape->r[n-1];
		if(outerRadius < safetyWindowShape->flangeOuterRadius)
			outerRadius = safetyWindowShape->flangeOuterRadius;
	}
}

void BLCMDabsorber::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	generateTubs(npoints, 0.0, outerRadius, 0.0, 360*deg, totalLength, v);
}

G4bool BLCMDabsorber::isOutside(G4ThreeVector &local, G4double tolerance)
{
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	return r > outerRadius-tolerance ||
		fabs(local[2]) > totalLength/2.0-tolerance;
}

bool BLCMDabsorber::isWithin(G4ThreeVector &local, G4double tolerance)
{
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	return r < absTubsRadius+tolerance &&
		fabs(local[2]) < absTubsLength/2.0+tolerance;
}
