//	BLCMDpillbox.cc
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

#ifdef G4BL_GSL

#include <vector>
#include <stdio.h>

#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SteppingManager.hh"
#include "G4Material.hh"

#include "gsl/gsl_sf_bessel.h"

#include "BLElement.hh"
#include "BLElementField.hh"
#include "BLGlobalField.hh"
#include "BLManager.hh"
#include "BLParam.hh"
#include "BLFieldMap.hh"
#include "BLTune.hh"
#include "BLKillTrack.hh"

extern void g4bl_exit(int);

const double UNDETERMINED = -4.7e21;	// large, negative, unlikely value
inline bool isUndetermined(G4double v) { return v == UNDETERMINED; }

const G4int ITERATION_LIMIT = 10;	// limit on iterations, for 
					// setting timeOffset

class BLCMDpillbox; // forward reference

/**	PillboxField represents one placement of a pillbox.
 *
 *	frequency=0.0 is accepted and yields a pillbox with a constant
 *	Ez = maxGradient (useful to verify units are correct).
 **/
class PillboxField : public BLElementField, public BLManager::SteppingAction {
	enum State {TIMING_UNKNOWN, SETTING_TIMING, TIMING_COMPLETE};
	G4String name;
	BLCMDpillbox *pillbox;
	G4VPhysicalVolume *timingPV;
	BLCoordinateTransform global2local;
	G4RotationMatrix *rotation;
	G4double timeOffset;
	G4Track saveTrack;
	G4int timeCount;
	G4bool validSaveTrack;
	State state;
	BLFieldMap *fieldMap;
	friend class BLCMDpillbox;
public:
	/// constructor. _zmin/max are global coordinates.
	PillboxField(G4String& _name,BLCoordinateTransform& _global2local,
			BLCMDpillbox *_pillbox, G4VPhysicalVolume *_timingPV);

	/// getName() returns the name of this placement.
	G4String getName() { return name; }

	/// addFieldValue() adds the field for this pillbox into field[].
	/// point[] is in global coordinates.
	void addFieldValue(const G4double point[4], G4double field[6]) const;

	/// doStep() handles a step of the tune particle. 
	void UserSteppingAction(const G4Step *step);
};

/**	BLCMDpillbox implements a pillbox.
 *
 *	Each placement of a BLCMDpillbox creates a PillboxField that is linked
 *	into BLGlobalField. Each PillboxField sets its timeOffset via the
 *	tune particle (re-stepping through its TimingVol if necessary).
 *
 *	If fieldMapFile is non-null, it is read as a BLFieldMap and determines
 *	the peak field; both B and E are multiplied by maxGradient of the 
 *	cavity.
 **/
class BLCMDpillbox : public BLElement {
	G4double maxGradient;
	G4String color;
	G4double frequency;
	G4double innerLength;
	G4double innerRadius;
	G4double pipeThick;
	G4double wallThick;
	G4double irisRadius;
	G4double collarRadialThick;
	G4double collarThick;
	G4double win1Thick;
	G4double win1OuterRadius;
	G4double win2Thick;
	G4String winMat;
	G4double phaseAcc;
	G4double skinDepth;
	G4double timingTolerance;
	G4double maxStep;
	G4String cavityMaterial;
	G4double timeOffset;
	G4double timeIncrement;
	G4String fieldMapFile;
	G4int kill;
	// non-arguments:
	G4double overallLength;
	G4double outerRadius;
	G4double omega;
	G4double rFactor;
	G4double E2Bfactor;
	BLFieldMap *fieldMap;
	std::vector<PillboxField *> pillboxField;
	friend class PillboxField;
public:
	/// Default constructor. Defines the command, etc.
	BLCMDpillbox();

	/// Destructor.
	virtual ~BLCMDpillbox() { }

	/// Copy constructor.
	BLCMDpillbox(const BLCMDpillbox& r);

	/// clone()
	BLElement *clone() { return new BLCMDpillbox(*this); }

	/// commandName() returns "pillbox".
	G4String commandName() { return "pillbox"; }

	/// command() implements the pillbox command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments of the command.
	void defineNamedArgs();

	/// argChanged() does internal computations after some arg changed
	void argChanged();

	/// construct() - construct the pillbox.
	/// Creates a new PillboxField and adds it to BLGlobalField.
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume* parent,
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the overallLength of the pillbox
	G4double getLength() { return overallLength; }

	/// getWidth() returns the outer radius of the pillbox
	G4double getWidth() { return outerRadius*2.0; }

	/// getHeight() returns the outer radius of the pillbox
	G4double getHeight() { return outerRadius*2.0; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true only if all placed pillboxes of this
	/// type have had their timing set, and if the tuning (if
	/// used) converged.
	G4bool isOK();

	/// generatePoints() from BLElement
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isOutside() from BLElement
	G4bool isOutside(G4ThreeVector &local, G4double tolerance);
};

BLCMDpillbox defaultPillbox;	// default object

// Default constructor - be sure to use the default constructor BLElement()
BLCMDpillbox::BLCMDpillbox() : BLElement(), pillboxField()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("Defines a pillbox RF cavity");
	setDescription("A Pillbox RF cavity is the basic RF element used to\n"
		"construct a linac. The phaseAcc parameter sets the\n"
		"phase of the tune particle at the center of the\n"
		"cavity, and the timing offset of the cavity is determined\n"
		"from that the first time that the Tune particle is tracked "
		"through the cavity. Zero degrees is the rising zero-crossing "
		"of the Ez field. If timeOffset is specified, it is used "
		"rather than setting it from the Tune particle.\n\n"
		"The Pipe, walls, and collars are always made of copper. "
		"Pipe, wall, collar, win1, and win2 can be omitted by setting\n"
		"their thickness to 0. Common usage is to set the collar "
		"values so by placing multiple pillboxes sequentially the "
		"collars form a beam pipe between them.\n\n"
		"Note that section 4.4 of the User's "
		"This element must be placed (via the place command).\n\n"
		"Guide has a dimensioned drawing of a pillbox.");

	// provide initial values for fields
	maxGradient = UNDETERMINED; // will be converted to MV/m
	color = "1.0,0.0,0.0";
	frequency = UNDETERMINED;
	innerLength = UNDETERMINED;
	innerRadius = UNDETERMINED;
	pipeThick = 3.0*mm;
	wallThick = 3.0*mm;
	irisRadius = 11.0*cm;
	collarRadialThick = 5.0*mm;
	collarThick = 2.5*mm;
	win1Thick = 0.200*mm;
	win1OuterRadius = 5.0*cm;
	win2Thick = 0.500*mm;
	winMat = "Be";
	phaseAcc = 40.0*deg;
	skinDepth = 0.002*mm;
	timingTolerance = 0.001*ns;
	maxStep = -1.0;
	cavityMaterial = "Vacuum";
	timeOffset = UNDETERMINED;
	timeIncrement = 0.0;
	fieldMapFile = "";
	kill = 0;
	// non-arguments:
	overallLength = UNDETERMINED;
	outerRadius = UNDETERMINED;
	omega = UNDETERMINED;
	rFactor = UNDETERMINED;
	E2Bfactor = UNDETERMINED;
	fieldMap = 0;
}

// Copy constructor - be sure to use the copy constructor BLElement(r)
BLCMDpillbox::BLCMDpillbox(const BLCMDpillbox& r) : BLElement(r), pillboxField()
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	BLTune::copyTunableArg(&maxGradient,&r.maxGradient);
	color = r.color;
	frequency = r.frequency;
	innerLength = r.innerLength;
	innerRadius = r.innerRadius;
	pipeThick = r.pipeThick;
	wallThick = r.wallThick;
	irisRadius = r.irisRadius;
	collarRadialThick = r.collarRadialThick;
	collarThick = r.collarThick;
	win1Thick = r.win1Thick;
	win1OuterRadius = r.win1OuterRadius;
	win2Thick = r.win2Thick;
	winMat = r.winMat;
	phaseAcc = r.phaseAcc;
	skinDepth = r.skinDepth;
	timingTolerance = r.timingTolerance;
	maxStep = r.maxStep;
	cavityMaterial = r.cavityMaterial;
	timeOffset = r.timeOffset;
	timeIncrement = r.timeIncrement;
	fieldMapFile = r.fieldMapFile;
	kill = r.kill;

	// non-arguments:
	overallLength = r.overallLength;
	outerRadius = r.outerRadius;
	omega = r.omega;
	rFactor = r.rFactor;
	E2Bfactor = r.E2Bfactor;
	fieldMap = r.fieldMap;
}

int BLCMDpillbox::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	// If the name of the cavity is "testUnits", perform the test and exit.
	if(argv[0] == "testUnits") {
		// set the test case values, from Al Moretti's SuperFish example
		maxGradient = 1.0 * megavolt/meter;
		frequency = 0.956188/ns;
		double Hcorrect = 1545.0; // A/m
		double BcorrectTesla = 4.0*pi*1E-7 * Hcorrect; // MKS
		argChanged();
		// the first maximum of J1() has a value of 0.5819
		double Bmax = maxGradient*E2Bfactor*0.5819;
		printf("pillobx testUnits:\n           "         
			"Freq(GHz)  Emax(MV/m)  Bmax(Tesla)  Radius(mm)\n");
		printf("Values:  %10.4f %10.4f %12.4e %11.4f\n",
			frequency/(1e9*hertz), maxGradient/(megavolt/meter),
			Bmax/tesla, innerRadius/mm);
		// these values already are in the correct units:
		printf("Correct: %10.4f %10.4f %12.4e %8.1f\n",
			0.956188,1.0,BcorrectTesla, 120.0);
		printf("pillbox testUnits: exiting\n");
		g4bl_exit(0);
	}

	if(argv.size() != 1) {
		printError("pillbox: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultPillbox.handleNamedArgs(namedArgs);
	}

	BLCMDpillbox *p = new BLCMDpillbox(defaultPillbox);
	p->setName(argv[0]);
	int retval = p->handleNamedArgs(namedArgs);

	if(p->fieldMapFile != "") {
		p->fieldMap = new BLFieldMap();
		p->fieldMap->readFile(p->fieldMapFile);
	}

	// check material exists
	if(p->cavityMaterial.size() > 0) getMaterial(p->cavityMaterial);

	p->print(argv[0]);

	return retval;
}

void BLCMDpillbox::defineNamedArgs()
{
	argTunable(maxGradient,"maxGradient","The peak gradient of the cavity (MV/m)",megavolt/meter);
	argString(color,"color","The color of the cavity");
	argDouble(frequency,"frequency","The frequency of the cavity (GHz)",1/ns,"",false);
	argDouble(innerLength,"innerLength","The inside length of the cavity (mm)",mm,"",false);
	argDouble(innerRadius,"innerRadius","The inside radius of the cavity (mm)",mm,"",false);
	argDouble(pipeThick,"pipeThick","The thickness of the pipe wall (mm)",mm,"",false);
	argDouble(wallThick,"wallThick","The thickness of the cavity walls (mm)",mm,"",false);
	argDouble(irisRadius,"irisRadius","The radius of the iris (mm)",mm,"",false);
	argDouble(collarRadialThick,"collarRadialThick","The radial thickness of the collar (mm)",mm,"",false);
	argDouble(collarThick,"collarThick","The thickness of the collar along z(mm)",mm,"",false);
	argDouble(win1Thick,"win1Thick","The thickness of the central portion of the\n"
		"windows; zero for no window (mm)",mm,"",false);
	argDouble(win1OuterRadius,"win1OuterRadius","The radius of the central portion of\n"
		"the windows (mm)",mm,"",false);
	argDouble(win2Thick,"win2Thick","The thickness of the outer portion of the\n"
		"windows; zero for no window (mm)",mm,"",false);
	argString(winMat,"winMat","The material of the windows");
	argDouble(phaseAcc,"phaseAcc","The reference phase of the cavity (degrees)",deg);
	argDouble(skinDepth,"skinDepth","The skin depth (mm)",mm,"",false);
	argDouble(timingTolerance,"timingTolerance","Tolerance for timing tuning (ns)",ns);
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm).",mm);
	argString(cavityMaterial,"cavityMaterial","Material of cavity volume (Vacuum).");
	argDouble(timeOffset,"timeOffset","Time offset for cavity (default: tuned by tune particle) (ns).",ns);
	argDouble(timeIncrement,"timeIncrement","Increment to timeOffset, applied AFTER tuning. (ns).",ns);
	argString(fieldMapFile,"fieldMapFile","Filename for BLFieldMap (pillbox if null).",false);
	argInt(kill,"kill","Set nonzero to kill tracks that hit the pipe, walls, or collars (0).");
}

void BLCMDpillbox::argChanged()
{
	// compute non-argument variables
	omega = 2.0 * pi * frequency;
	if(innerRadius == UNDETERMINED) {
		if(frequency > 0.0) {
			// the first zero of J0(r) is at r=2.405
			innerRadius =  2.405 / (omega/c_light);
		} else {
			printError("pillbox: innerRadius is undefined");
		}
	}
	double wall=collarThick;
	if(wallThick > wall) wall = wallThick;
	if(win1Thick > wall) wall = win1Thick;
	if(win2Thick > wall) wall = win2Thick;
	overallLength = innerLength + 2.0 * wall;
	outerRadius = innerRadius + pipeThick;
	rFactor = omega/c_light;
	// Previous versions had the incorrect value from the BeamTools:
	//     E2Bfactor = 1.0e-9 / c_light;   // MeV/m, MKS units
	// C.f. Eq. 8.77 in J.D.Jackson, _Classical_Electrodynamics_,
	//	 first edition (1962). His value is 1/c in vacuum.
	E2Bfactor = 1.0/c_light;   // Geant4 units
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");
}

void BLCMDpillbox::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume* parent,
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	G4Material *cu = getMaterial("Cu");
	G4Material *wm = getMaterial(winMat);
	G4Material *cavityMat = getMaterial(cavityMaterial);
	G4UserLimits *userLimits = new G4UserLimits(maxStep);

	// find enclosing Tubs
	G4double halflen = innerLength/2+wallThick;
	halflen = std::max(halflen,innerLength/2+collarThick);
	halflen = std::max(halflen,innerLength/2+win1Thick);
	halflen = std::max(halflen,innerLength/2+win2Thick);
	G4Tubs *tubs = new G4Tubs(getName(),0.0,innerRadius+pipeThick,
					halflen,0.0,2.0*pi);
	G4LogicalVolume *pillbox = new G4LogicalVolume(tubs,cavityMat,
					parentName+getName(),0,0,userLimits);

	G4ThreeVector loc;
	G4VPhysicalVolume *pv;

	// The outer pipe
	G4LogicalVolume *pipe=0;
	if(pipeThick > 0.001*mm) {
		tubs = new G4Tubs(getName()+"pipe",innerRadius,
					innerRadius+pipeThick,
					innerLength/2+wallThick,0.0,2.0*pi);
		pipe = new G4LogicalVolume(tubs,cu,
					getName()+"pipeLog",0,0,userLimits);
		pv = new G4PVPlacement(0,loc,pipe,parentName+getName()+"Pipe",
					pillbox,false,0,surfaceCheck);
		if(kill) BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
	}

	// the wall (us and ds) - radially: collarOuterRadius to innerRadius
	G4LogicalVolume *wall=0;
	if(wallThick > 0.001*mm) {
		tubs = new G4Tubs(getName()+"wall",irisRadius+collarRadialThick,
					innerRadius,wallThick/2,0.0,2.0*pi);
		wall = new G4LogicalVolume(tubs,cu,
					getName()+"wallLog",0,0,userLimits);
		G4double wallCenterZ = innerLength/2.0 + wallThick/2.0;
		loc.setZ(-wallCenterZ);
		pv = new G4PVPlacement(0,loc,wall,parentName+getName()+"UsWall",
					pillbox,false,0,surfaceCheck);
		if(kill) BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
		loc.setZ(+wallCenterZ);
		pv = new G4PVPlacement(0,loc,wall,parentName+getName()+"DsWall",
					pillbox,false,0,surfaceCheck);
		if(kill) BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
	}
	
	// the collar (us and ds) - includes wall thickness (i.e. the inner
	// edge of the collar is the same as the inner edge of the wall)
	G4LogicalVolume *collar=0;
	if(collarThick > 0.001*mm) {
		tubs = new G4Tubs(getName()+"collar",irisRadius,
					irisRadius+collarRadialThick,
					collarThick/2,0.0,2.0*pi);
		collar = new G4LogicalVolume(tubs,cu,
					getName()+"collarLog",0,0,userLimits);
		G4double collarCenterZ = innerLength/2.0 + collarThick/2.0;
		loc.setZ(-collarCenterZ);
		pv = new G4PVPlacement(0,loc,collar,parentName+getName()+"UsCollar",
					pillbox,false,0,surfaceCheck);
		if(kill) BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
		loc.setZ(+collarCenterZ);
		pv = new G4PVPlacement(0,loc,collar,parentName+getName()+"DsCollar",
					pillbox,false,0,surfaceCheck);
		if(kill) BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
	}
	
	// win1 (the inner window)
	G4LogicalVolume *win1=0;
	if(win1Thick > 0.001*mm) {
		tubs = new G4Tubs(getName()+"win1",0.0,win1OuterRadius,
					win1Thick/2,0.0,2.0*pi);
		win1 = new G4LogicalVolume(tubs,wm,
					getName()+"win1Log",0,0,userLimits);
		G4double win1CenterZ = innerLength/2.0 + win1Thick/2.0;
		loc.setZ(-win1CenterZ);
		new G4PVPlacement(0,loc,win1,parentName+getName()+"UsWin1",
					pillbox,false,0,surfaceCheck);
		loc.setZ(+win1CenterZ);
		new G4PVPlacement(0,loc,win1,parentName+getName()+"DsWin1",
					pillbox,false,0,surfaceCheck);
	}

	// win2 (the outer window)
	G4LogicalVolume *win2=0;
	if(win2Thick > 0.001*mm) {
		tubs = new G4Tubs(getName()+"win2",win1OuterRadius,irisRadius,
					win2Thick/2,0.0,2.0*pi);
		win2 = new G4LogicalVolume(tubs,wm,
					getName()+"win2Log",0,0,userLimits);
		G4double win2CenterZ = innerLength/2.0 + win2Thick/2.0;
		loc.setZ(-win2CenterZ);
		new G4PVPlacement(0,loc,win2,parentName+getName()+"UsWin2",
					pillbox,false,0,surfaceCheck);
		loc.setZ(+win2CenterZ);
		new G4PVPlacement(0,loc,win2,parentName+getName()+"DsWin2",
					pillbox,false,0,surfaceCheck);
	}

	// half of the interior
	tubs = new G4Tubs(getName()+"half",0.0,innerRadius,
					innerLength/4,0.0,2.0*pi);
	G4LogicalVolume *half = new G4LogicalVolume(tubs,cavityMat,getName()+"halfLog",
	 				0,0,userLimits);
	G4double halfCenterZ = innerLength/4.0;
	loc.setZ(+halfCenterZ);
	new G4PVPlacement(0,loc,half,parentName+getName()+"DsHalf",
						pillbox,false,0,surfaceCheck);
	loc.setZ(-halfCenterZ);
	G4VPhysicalVolume *timingPV = new G4PVPlacement(0,loc,half,
		parentName+getName()+"TimingVol",pillbox,false,0,surfaceCheck);

#ifdef G4BL_VISUAL
	G4VisAttributes *c1 = new G4VisAttributes(true,G4Color(0.5,0.0,0.0));
	G4VisAttributes *c2 = new G4VisAttributes(true,G4Color(0.7,0.0,0.0));
	G4VisAttributes *c3 = new G4VisAttributes(true,G4Color(1.0,0.0,0.0));
	const G4VisAttributes *invisible = BLCommand::getVisAttrib("Invisible");
	if(pipe) pipe->SetVisAttributes(c1);
	if(wall) wall->SetVisAttributes(c1);
	if(collar) collar->SetVisAttributes(c1);
	if(win1) win1->SetVisAttributes(c3);
	if(win2) win2->SetVisAttributes(c2);
	half->SetVisAttributes(invisible);
	pillbox->SetVisAttributes(invisible);
#endif

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	// place the pillbox into the parent
	new G4PVPlacement(g4rot,relativePosition,pillbox,
			parentName+getName(),parent,false,0,surfaceCheck);

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

	// add this pillbox field to the GlobalField, and to BLManager
	G4double dz = innerLength/2.0 + wallThick;
	// if rotated, make an overestimate of the field occupancy along z
	if(global2local.isRotated())
		dz += innerRadius;
	PillboxField *pf = 
		new PillboxField(thisname,global2local,this,timingPV);
	BLGlobalField::getObject()->addElementField(pf);
	pillboxField.push_back(pf);
	BLManager::getObject()->registerTuneParticleStep(timingPV,pf);

	printf("BLCMDpillbox::construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		global2local.getPosition()[2], global2local.getPosition()[2]-dz,
		global2local.getPosition()[2]+dz);
}

G4bool BLCMDpillbox::isOK()
{
	G4bool retval = true;

	// verify all placed pillboxField-s have set their timeOffset 
	std::vector<PillboxField *>::iterator i;
	for(i=pillboxField.begin(); i!=pillboxField.end(); ++i) {
		if((*i)->state != PillboxField::TIMING_COMPLETE) {
			printf("*** BLCMDpillbox %s has not determined timing\n",
				(*i)->getName().c_str());
			retval = false;
		}
	}

	return retval;
}

PillboxField::PillboxField(G4String& _name, BLCoordinateTransform& _global2local,
			BLCMDpillbox *_pillbox,
			G4VPhysicalVolume *_timingPV) :
			BLElementField(), BLManager::SteppingAction()
{
	name = _name;
	pillbox = _pillbox;
	timingPV = _timingPV;
	global2local = _global2local;
	rotation = 0;
	timeOffset = pillbox->timeOffset;
	timeCount = 0;
	validSaveTrack = false;
	state = TIMING_UNKNOWN;
	fieldMap = pillbox->fieldMap;

	if(global2local.isRotated()) {
		rotation = new G4RotationMatrix(global2local.getRotation());
		rotation->invert();
	}

	if(pillbox->frequency == 0.0 || pillbox->maxGradient == 0.0)
		timeOffset = 0.0;

	// don't tune timeOffset if it is known
	if(!isUndetermined(timeOffset))
		state = TIMING_COMPLETE;
	else
		timeOffset = 0.0; // don't start tuning at UNDETERMINED

	// set global bounding box
	G4double local[4], global[4];
	local[3] = 0.0;
	G4double dz = pillbox->innerLength/2.0 + pillbox->wallThick;
	for(int i=0; i<2; ++i) {
		local[0] = (i==0 ? -1.0 : 1.0) * pillbox->innerRadius;
		for(int j=0; j<2; ++j) {
			local[1] = (j==0 ? -1.0 : 1.0) * pillbox->innerRadius;
			for(int k=0; k<2; ++k) {
				local[2] = (k==0 ? -1.0 : 1.0) * dz;
				global2local.getGlobal(local,global);
				setGlobalPoint(global);
			}
		}
	}
}

void PillboxField::addFieldValue(const G4double point[4], G4double field[6]) const
{
	G4double local[4];
	global2local.getLocal(local,point);

	if(fieldMap) {
		G4double myField[6];
		fieldMap->getFieldValue(local,myField,pillbox->maxGradient,
							pillbox->maxGradient);
		// frequency == 0.0 means a d.c. electric field
		G4double timeE=1.0, timeB=0.0;
		if(pillbox->frequency != 0) {
			double tArg = pillbox->omega * (point[3] - timeOffset);
			timeE = cos(tArg);
			timeB = sin(tArg);
		}
		if(rotation) {
			G4ThreeVector B(myField[0],myField[1],myField[2]);
			G4ThreeVector E(myField[3],myField[4],myField[5]);
			B = *rotation * B;
			E = *rotation * E;
			field[0] += B[0] * timeB;
			field[1] += B[1] * timeB;
			field[2] += B[2] * timeB;
			field[3] += E[0] * timeE;
			field[4] += E[1] * timeE;
			field[5] += E[2] * timeE;
		} else {
			field[0] += myField[0] * timeB;
			field[1] += myField[1] * timeB;
			field[2] += myField[2] * timeB;
			field[3] += myField[3] * timeE;
			field[4] += myField[4] * timeE;
			field[5] += myField[5] * timeE;
		}
		return;
	}

	G4double x = local[0];
	G4double y = local[1];
	G4double z = local[2];
	G4double r = sqrt(x*x + y*y);
	G4double dz = pillbox->innerLength/2.0 + pillbox->wallThick;
	if(z < -dz || z > dz || r > pillbox->innerRadius) return;

	// frequency == 0.0 means a d.c. electric field
	if(pillbox->frequency == 0.0) {
		// zmin and zmax include the walls, so re-test z position
		if(z >= -pillbox->innerLength/2.0 && 
		   z <=  pillbox->innerLength/2.0) {
			if(rotation) {
				G4ThreeVector E(0.0, 0.0, pillbox->maxGradient);
				E = *rotation * E;
				field[3] += E[0];
				field[4] += E[1];
				field[5] += E[2];
			} else {
				field[5] += pillbox->maxGradient;
			}
		}
		return;
	}

	// compute normal RF pillbox field
	double arg = pillbox->rFactor*r;
	if(fabs(arg) < 1.0e-20) arg = 0.0;
	double j0 = gsl_sf_bessel_J0(arg);
	double j1 = gsl_sf_bessel_J1(arg);
	double tArg = pillbox->omega * (point[3] - timeOffset);
	double ez = pillbox->maxGradient * j0 * cos(tArg);
	double bphi = -pillbox->maxGradient*pillbox->E2Bfactor*j1*sin(tArg);
	// handle skin depth if in wall or window
	double f = 1.0;
	if(z > pillbox->innerLength/2.0)
		f=exp(-(z-pillbox->innerLength/2.0)/pillbox->skinDepth);
	else if(z < -pillbox->innerLength/2.0)
		f=exp(-(-z-pillbox->innerLength/2.0)/pillbox->skinDepth);
	//@ double phi = atan2(point[1],point[0]);
	//@ TJR 3/29/2011 - this surely should be local, not point
	double phi = atan2(local[1],local[0]);
	if(rotation) {
		G4ThreeVector B(-bphi*sin(phi)*f, bphi*cos(phi)*f, 0.0);
		G4ThreeVector E(0.0, 0.0, ez*f);
		B = *rotation * B;
		E = *rotation * E;
		field[0] += B[0];
		field[1] += B[1];
		field[2] += B[2];
		field[3] += E[0];
		field[4] += E[1];
		field[5] += E[2];
	} else {
		field[0] += -bphi * sin(phi) * f;
		field[1] += bphi * cos(phi) * f;
		field[5] += ez * f;
	}
}

// called ONLY in tune particle mode, for this pillbox's timingVol
void  PillboxField::UserSteppingAction(const G4Step *step)
{
    // do nothing if d.c. or no field or timing is complete
    if(pillbox->frequency <= 0.0 || pillbox->maxGradient == 0.0) 
    	state = TIMING_COMPLETE;
    if(state == TIMING_COMPLETE) return;

    G4Track *track = step->GetTrack();
    G4StepPoint *prePoint = step->GetPreStepPoint();
    if(!prePoint) return;
    G4VPhysicalVolume *prePV = prePoint->GetPhysicalVolume();
    G4StepPoint *postPoint = step->GetPostStepPoint();
    if(!postPoint) return;
    G4VPhysicalVolume *postPV = postPoint->GetPhysicalVolume();
    G4SteppingManager *steppingMgr =
				BLManager::getObject()->getSteppingManager();

    if(prePV == postPV) return;     // neither entering nor leaving

    if(postPV == timingPV) {	    // entering timingPV
        // save track for setting timeOffset
        saveTrack.CopyTrackInfo(*step->GetTrack());
	saveTrack.SetUserInformation(0);
        validSaveTrack = true;
	state = SETTING_TIMING;
        // estimate timingOffset by neglecting acceleration
        double v = saveTrack.GetVelocity();
        double dist = pillbox->innerLength/2.0;
        double arrival = saveTrack.GetGlobalTime() + dist/v;
        // pi/2 because E is a cos() and we want phaseAcc=0
        // to be the rising zero-crossing
        timeOffset = arrival - (pillbox->phaseAcc-pi/2.0)/pillbox->omega;
    } else if(prePV == timingPV) {	    // leaving timingPV
        if(++timeCount > ITERATION_LIMIT)
		G4Exception("pillbox","Iteration Limit",FatalException,"");
        // check timing offset
        double t = track->GetGlobalTime();
        double phase = pillbox->omega * (t-timeOffset) + pi/2.0;
        double dt = (phase-pillbox->phaseAcc)/pillbox->omega;
        timeOffset = timeOffset + dt;
        if(fabs(dt) <= pillbox->timingTolerance) {
            // Success!
            printf("pillbox %s: Time OK  timeOffset=%.4f ns, incremented to "
	    	"%.4f ns\n",getName().c_str(),timeOffset,
					timeOffset+pillbox->timeIncrement);
	    timeOffset += pillbox->timeIncrement;
            state = TIMING_COMPLETE;
            return;
        }
        // restore saved data (i.e. jump back to when the track 
        // enteered TimingVol)
        if(!validSaveTrack)
		G4Exception("pillbox","Invalid Step",FatalException,"");
        steppingMgr->GetfSecondary()->push_back(new G4Track(saveTrack));
        track->SetTrackStatus(fStopAndKill);
    }
}

void BLCMDpillbox::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	generateTubs(npoints, 0.0, outerRadius, 0.0, 360.0*deg,
			overallLength, v);
}

G4bool BLCMDpillbox::isOutside(G4ThreeVector &local, G4double tolerance)
{
	G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
	return r > outerRadius-tolerance ||
		fabs(local[2]) > overallLength/2.0-tolerance;
}

#else // G4BL_GSL
int dummypillbox=0;
#endif // G4BL_GSL
