//	BLCMDdetector.cc
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

#include <map>
#include <vector>

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4UserLimits.hh"
#include "G4StepPoint.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4EmSaturation.hh"

#include "BLElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLTrackNTuple.hh"
#include "BLCoordinates.hh"
#include "BLGlobalField.hh"
#include "BLEvaluator.hh"
#include "BLKillTrack.hh"

const char DetectorFields[] =
    "x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight:Edep:VisibleEdep:Ntracks";
const unsigned NDetectorFields = 15;
static const char *fieldNames[] = {
	"x", "y", "z", "Px", "Py", "Pz", "t", "PDGid", "EventID", "TrackID",
	"ParentID", "Weight", "Edep", "VisibleEdep", "Ntracks"
};

/**	class BLCMDdetector - implements a Detector
 *	Each placement of this class generates an NTuple of the beam as it
 *	enters the physical volume of the Detector. This is therefore a
 *	"perfect" detector in that it does not perturb the beam at all, and
 *	intrinsically has the resolution of a float (it measures position
 *	and momentum and time).
 *
 *	The NTuple for a detector can be added to a BLCMDntuple by
 *	including a pattern that matches its name in the 'detectors'
 *	argument to the ntuple command.
 *
 *	Note that if a BLCoordinates instance is linked into the track,
 *	its centerline coordinates are used; otherwise global coordinates
 *	are used.
 **/
class BLCMDdetector : public BLElement {
	G4double radius;
	G4double innerRadius;
	G4double height;
	G4double width;
	G4double length;
	G4String material;
	G4double birks;
	G4String color;
	G4String solid;
	G4VSolid *detSolid;
	G4double maxStep;
	G4int noSingles;
	G4String format;
	G4String filename;
	G4String require;
	G4int referenceParticle;
	G4String coordinates;
	G4int kill;
	G4int perTrack;
	BLCoordinateType coordinateType;
public:
	/// Default constructor.
	BLCMDdetector();

	/// Destructor.
	 ~BLCMDdetector();

	/// clone()
	BLElement *clone() { return new BLCMDdetector(*this); }

	/// Copy constructor.
	BLCMDdetector(const BLCMDdetector& r);

	/// commandName() returns "detector".
	G4String commandName() { return "detector"; }

	/// command() implements the detector command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// argChanged() does internal computations after some arg changed
	void argChanged();

	/// help() prints help text.
	void help(bool detailed) {
		if(description[description.size()-2] == ':')
			description += BLTrackNTuple::getFormatList(); 
		BLCommand::help(detailed);
	}

	/// construct() will construct the detector.
	/// Used for normal placements of a Detector object.
	 void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns this element's Length along the Z axis.
	 G4double getLength() { return length; }

	/// getWidth() returns this element's Width along the X axis.
	 G4double getWidth() { return width; }

	/// getHeight() returns this element's height along the Y axis.
	 G4double getHeight() { return height; }

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

BLCMDdetector defaultDetector;

/**	class BLDetectorNTuple implements an NTuple for a BLCMDdetector.
 **/
class BLDetectorNTuple : public BLManager::SteppingAction, 
			public BLManager::EventAction {
	G4String name;
	G4VPhysicalVolume *thisVol;
	BLCoordinateType coordinateType;
	int noSingles;
	G4String require;
	int kill;
	int perTrack;
	BLNTuple *ntuple;
	bool verbose;
	double data[NDetectorFields];
	G4EmSaturation *emSaturation;
	BLEvaluator *eval;
	int prevTrackID;
	friend class BLCMDdetector;
public:
	/// constructor.
	BLDetectorNTuple(G4String type, G4String category, 
		G4String _name, G4VPhysicalVolume *pv, int _noSingles,
		G4String filename, G4String _require, 
		BLCoordinateType _coordinateType, int _kill, int _perTrack);

	bool needsReference() 
		{ return ntuple ? ntuple->needsReference() : false; }

	/// from BLManager::EventAction:
	void BeginOfEventAction(const G4Event* event);
	void EndOfEventAction(const G4Event* event);

	/// from BLManager::SteppingAction:
	void UserSteppingAction(const G4Step *step);

	void clear();
	void setFirstTrackFields(const G4Step *step, const G4Track *track);
	void accumulate(const G4Step *step, const G4Track *track);
	void finish();
};

BLCMDdetector::BLCMDdetector() : BLElement()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Construct a Detector that generates an NTuple.");
	setDescription("A Detector generates an NTuple of tracks when they\n"
		"enter the physical volume of the Detector, summing up the "
		"total energy deposited until it leaves the volume. "
		"By default an entry in the NTuple is generated for each "
		"event, with the track variables set by the first track as it "
		"enters the detector. Set perTrack=1 to generate "
		"an entry for each track.\n\n"
		"A detector may\n"
		"be placed via multiple place commands (usually with a "
		"'rename=det#' argument to distinguish the different "
		"placements). Each placement creates an individual NTuple.\n\n"
		"If material is not specified, it uses Scintillator.\n\n"
		"There are three ways to specify the geometry:\n"
		" * set solid to the name of an existing solid.\n"
		" * set height, width, and length to get a rectangular box.\n"
		" * set radius and length (and innerRadius) to get a cylinder\n"
		"   (tube).\n"
		"By using the solid argument you can make a detector out of "
		"any solid, including sphere, boolean, and extrusion.\n\n"
		"The NTuple by default uses local coordinates.\n\n"
		"The NTuple of the detector can be included in an "
		"ntuple command by including a pattern that matches its name "
		"in the 'detectors' argument to the ntuple command. Note that "
		"must match the name as placed (i.e. includes rename=), not "
		"the name given to this command. The noSingles argument may be "
		"useful in this case to avoid a huge NTuple of singles (an "
		"empty NTuple may be created).\n\n"
		"Note that secondary particles created within the detector "
		"will not get an entry until they have taken one step. "
		"They are guaranteed to do so.\n\n"
		"For VisibleEdep, this command knows the Birks constants for "
		"Scintillator, POLYSTYRENE, BGO, and lAr. For other materials "
		"you must supply a value, or 0.0 will be used.\n\n"
		"This element must be placed (via the place command).\n\n"
		"The NTuple fields are:\n"
		"    x,y,z (mm)\n"
		"    Px,Py,Pz (MeV/c)\n"
		"    t (ns)\n"
		"    PDGid (11=e-, 13=mu-, 22=gamma, 211=pi+, 2212=proton...)\n"
		"    EventID (may be inexact above 16,777,215)\n"
		"    TrackID\n"
		"    ParentID (0 => primary particle)\n"
		"    Weight (defaults to 1.0)\n"
		"    Edep (Energy deposited, MeV)\n"
		"    VisibleEdep (Energy deposited, MeV, reduced by Birks effect)\n"
		"    Ntracks (# tracks of this event that hit the detector)\n"
		"Valid formats: root, ascii (no extended formats).");
	
	// initial field values
	radius = 0.0;
	innerRadius = 0.0;
	height = 0.0;
	width = 0.0;
	length = 1.0*mm;
	material = "Scintillator";
	birks = 0.0*mm/MeV;
	color = "1,1,1";
	solid = "";
	detSolid = 0;
	maxStep = -99.0;
	noSingles = 0;
	format = "";
	filename = "";
	require = "";
	referenceParticle = 0;
	coordinates = "Local";
	kill = 0;
	perTrack = 0;
	coordinateType = BLCOORD_LOCAL;
}

BLCMDdetector::~BLCMDdetector()
{
	if(detSolid) delete detSolid;
}

BLCMDdetector::BLCMDdetector(const BLCMDdetector& r) : BLElement(r)
{
	radius = r.radius;
	innerRadius = r.innerRadius;
	height = r.height;
	width = r.width;
	length = r.length;
	material = r.material;
	birks = r.birks;
	color = r.color;
	solid = r.solid;
	detSolid = r.detSolid;
	maxStep = r.maxStep;
	noSingles = r.noSingles;
	format = r.format;
	filename = r.filename;
	referenceParticle = r.referenceParticle;
	coordinates = r.coordinates;
	kill = r.kill;
	perTrack = r.perTrack;
	require = r.require;
	coordinateType = r.coordinateType;
}

int BLCMDdetector::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("detector: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return handleNamedArgs(namedArgs);
	}

	BLCMDdetector *t = new BLCMDdetector(defaultDetector);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	t->coordinateType = BLCoordinates::getCoordinateType(t->coordinates);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists, handle Birks constant
	G4Material *mat = getMaterial(t->material);
	BLAssert(mat != 0);
	if(t->birks == 0.0)
		t->birks = mat->GetIonisation()->GetBirksConstant();
	if(t->birks == 0.0) {
		G4String s = G4String("G4_")+t->material;
		G4Material *m = getMaterial(s,true);
		if(m) t->birks = m->GetIonisation()->GetBirksConstant();
	}
	if(t->birks == 0.0 && t->material == "Scintillator") {
		G4Material *m = getMaterial("G4_POLYSTYRENE",true);
		if(m) t->birks = m->GetIonisation()->GetBirksConstant();
	}
	if(t->birks == 0.0) {
		// values copied from Geant4.9.2.p02 G4EmSaturation.cc
		if(t->material == "Scintillator" ||
		   t->material == "G4_POLYSTYRENE" ||
		   t->material == "POLYSTYRENE") {
			t->birks = 0.07943*mm/MeV;
		}
		if(t->material == "G4_BGO" || t->material == "BGO") {
			t->birks = 0.008415*mm/MeV;
		}
		if(t->material == "G4_lAr" || t->material == "lAr") {
			t->birks = 0.1576*mm/MeV;
		}
	}
	mat->GetIonisation()->SetBirksConstant(t->birks);

	t->print(argv[0]);

	return retval;
}

void BLCMDdetector::defineNamedArgs()
{
	argDouble(radius,"radius","The radius of the circular Detector (mm).");
	argDouble(innerRadius,"innerRadius","The inner radius of the circular Detector (0 mm).");
	argDouble(height,"height","The height of the rectangular Detector (mm).");
	argDouble(width,"width","The width of the rectangular Detector (mm).");
	argDouble(length,"length","The length of the Detector (mm).");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm).");
	argString(material,"material","The material of the Detector (Scintillator).");
	argDouble(birks,"birks","The Birks constant (0 mm/MeV, unless known)",mm/MeV);
	argString(color,"color","The color of the Detector (white, ''=invisible).");
	argString(solid,"solid","The name of the solid to use (overrides height,width,radius,length).");
	argInt(noSingles,"noSingles","Set to 1 to omit the NTuple for singles.");
	argString(format,"format","NTuple format: (see above for list).");
	argString(filename,"filename","filename ('' uses name to determine filename)");
	argString(filename,"file","alias for filename");
	argString(require,"require","Expression which must be nonzero to include the track (default=1)",false);
	argInt(referenceParticle,"referenceParticle","Set to 1 to include the Reference Particle.");
	argString(coordinates,"coordinates","Coordinates: global, centerline, reference, or local (default=local).");
	argString(coordinates,"coord","Alias for coordinates.");
	argInt(kill,"kill","Set to 1 kill all tracks after entering them into NTuple(s).");
	argInt(perTrack,"perTrack","Set to 1 generate an entry for each track (default is per event).");
}

void BLCMDdetector::argChanged()
{
	if(radius > 0.0) width = height = radius*2.0;
}

void BLCMDdetector::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4Material *mat = getMaterial(material);

	G4String thisname = parentName+getName();

	if(!detSolid) {
		if(solid != "") {
			BLElement *e = BLElement::find(solid);
			if(!e) {
				printError("detector solid '%s' not found",
								solid.c_str());
				return;
			}
			detSolid = e->getSolid();
			if(!detSolid) {
				printError("detector solid '%s' is invalid",
								solid.c_str());
				return;
			}
		} else if(radius > 0.0) {
			detSolid = new G4Tubs(thisname+"Tubs", innerRadius,
					radius, length/2.0, 0.0, 2.0*pi);
		} else if(height > 0.0 && width > 0.0) {
			detSolid = new G4Box(thisname+"Box",width/2.0,
					height/2.0,length/2.0);
		} else {
			printError("detector::construct %s INVALID - no "
				"radius or height&width",thisname.c_str());
			return;
		}
	}
	G4LogicalVolume *lv = new G4LogicalVolume(detSolid,mat,
							thisname+"LogVol");
	lv->SetVisAttributes(getVisAttrib(color));
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");
	lv->SetUserLimits(new G4UserLimits(maxStep));

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4VPhysicalVolume *pv = new G4PVPlacement(g4rot, relativePosition,
				lv,thisname,parent,false,0,surfaceCheck);

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
	if(globalRotation)
		globalPosition = *globalRotation * relativePosition +
				parentPosition;

	for(unsigned i=0; i<format.size(); ++i)
		format[i] = tolower(format[i]);
	if(format != "ascii" && format != "root")
		format = "";
	BLDetectorNTuple *nt = new BLDetectorNTuple(format,
			"Detector",thisname,pv,noSingles,filename,
			require,coordinateType,kill,perTrack);

	BLManager::getObject()->registerBeamStep(pv,nt);
	if(referenceParticle != 0 || nt->needsReference())
		BLManager::getObject()->registerReferenceParticleStep(pv,nt);
	if(!perTrack)
		BLManager::getObject()->registerEventAction(nt,false);

	// special case: if kill!=0 and viewer != "none", then register a
	// BLKillTrack. (in viewer mode we don't accumulate NTuples.)
	if(kill != 0 && Param.getString("viewer") != "none") {
		BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
	}

	printf("BLCMDdetector::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2]);
}

BLDetectorNTuple::BLDetectorNTuple(G4String type,
	G4String category, G4String _name, G4VPhysicalVolume *pv, int _noSingles,
	G4String filename, G4String _require,BLCoordinateType _coordinateType,
	int _kill, int _perTrack)

{
	name = _name;
	thisVol = pv;
	noSingles = _noSingles;
	coordinateType = _coordinateType;
	require = _require;
	eval = 0;
	if(require != "") eval = new BLEvaluator();
	ntuple = BLNTuple::create(type,category,name,DetectorFields,filename);
	kill = _kill;
	perTrack = _perTrack;
	verbose = (BLManager::getObject()->getSteppingVerbose() > 0);
	emSaturation = new G4EmSaturation(0);
	if(verbose) emSaturation->SetVerbose(verbose);
	prevTrackID = -1;
}

void BLDetectorNTuple::BeginOfEventAction(const G4Event* event)
{
	clear();
}

void BLDetectorNTuple::EndOfEventAction(const G4Event* event)
{
	if(!perTrack) finish();
}

void BLDetectorNTuple::UserSteppingAction(const G4Step *step)
{
	if(BLManager::getObject()->getState() == SPECIAL) return;

	// only use reference coordinates when they are valid
	BLManagerState state = BLManager::getObject()->getState();
	if(coordinateType == BLCOORD_REFERENCE && state != BEAM) return;

	// get basic physical-volume info
	G4Track *track = step->GetTrack();
	G4StepPoint *prePoint = step->GetPreStepPoint();
	if(!prePoint) return;
	G4VPhysicalVolume *preVol = prePoint->GetPhysicalVolume();
	if(!preVol) return;
	G4StepPoint *postPoint = step->GetPostStepPoint();
	if(!postPoint) return;
	G4VPhysicalVolume *postVol = postPoint->GetPhysicalVolume();
	if(!postVol) return;
	
	// determine if entering, inside, or leaving
	if(preVol != thisVol && postVol == thisVol) {		// enter
		if(perTrack) clear();
		setFirstTrackFields(step,track);
	} else if(preVol == thisVol && (postVol != thisVol ||	// leave
			track->GetTrackStatus() != fAlive)) {
		accumulate(step,track);
		if(perTrack) finish();
		if(kill) {
		    const_cast<G4Track*>(track)->SetTrackStatus(fStopAndKill);
		    if(verbose) printf("Track killed by '%s' with kill=1\n",
				name.c_str());
		}
	} else if(preVol == thisVol) {				// inside
		accumulate(step,track);
	}
}

void BLDetectorNTuple::clear()
{
	prevTrackID = -1;
	for(int i=0; i<NDetectorFields; ++i)
		data[i] = 0.0;
	data[9] = -1.0;
}

void BLDetectorNTuple::setFirstTrackFields(const G4Step *step,
							const G4Track *track)
{
	if(data[9] < 0.0) {
		G4ThreeVector position = track->GetPosition();
		G4ThreeVector momentum = track->GetMomentum();
		G4double time = track->GetGlobalTime();
		G4RunManager* runmgr = G4RunManager::GetRunManager();
		const G4Event* event = runmgr->GetCurrentEvent();
		int evId = event->GetEventID();

		// transform to desired coordinates, if available
		BLCoordinates *coord = (BLCoordinates *)track->GetUserInformation();
		if(coord && coord->isValid()) {
			if(coordinateType == BLCOORD_LOCAL) {
				G4StepPoint *prePoint = step->GetPreStepPoint();
				BLAssert(prePoint != 0);
				G4StepPoint *postPoint = 
						step->GetPostStepPoint();
				BLAssert(postPoint != 0);
				G4StepPoint *point = postPoint;
				if(prePoint->GetPhysicalVolume() == thisVol)
					point = prePoint;
				BLAssert(point->GetPhysicalVolume() == thisVol);
				const G4AffineTransform &trans = point->
					GetTouchableHandle()->GetHistory()->
					GetTopTransform();
				position = trans.TransformPoint(position);
				momentum = trans.TransformAxis(momentum);
			} else {
				coord->getCoords(coordinateType,position);
				momentum = coord->getRotation() * momentum;
			}
		} else {
			printf("BLCMDdetector::SteppingAction: track has no "
				"BLCoordinates object\n");
		}

		data[0] = position[0]/mm;		// x (mm)
		data[1] = position[1]/mm;		// y (mm)
		data[2] = position[2]/mm;		// z (mm)
		data[3] = momentum[0]/MeV;		// Px (MeV/c)
		data[4] = momentum[1]/MeV;		// Py (MeV/c)
		data[5] = momentum[2]/MeV;		// Pz (MeV/c)
		data[6] = time/ns;			// t (ns)
		data[7] = track->GetDefinition()->GetPDGEncoding();
		data[8] = evId;				// Event ID
		data[9] = BLManager::getObject()->getExternalTrackID(track);
		data[10] = BLManager::getObject()->getExternalParentID(track);
		data[11] = track->GetWeight();		// Weight
		data[12] = 0.0;				// Edep
		data[13] = 0.0;				// VisibleEdep
		data[14] = 0.0;				// Ntracks
	}

	if(prevTrackID != track->GetTrackID())  {
		data[14] += 1.0;
		prevTrackID = track->GetTrackID();
	}
}

void BLDetectorNTuple::accumulate(const G4Step *step, const G4Track *track)
{
	setFirstTrackFields(step,track);
	data[12] += step->GetTotalEnergyDeposit();
	data[13] += emSaturation->VisibleEnergyDepositionAtAStep(step);
}

void BLDetectorNTuple::finish()
{
	if(data[9] < 0.0) return;

	bool ignore=false;

	// implement require
	if(eval) {
		for(int i=0; i<NDetectorFields; ++i)
			eval->setVariable(fieldNames[i],data[i]);
		double v=eval->evaluate(require.c_str());
		if(!eval->isOK())
			G4Exception("BLTrackNTuple",
				"Invalid require expression",
				FatalException,require);
		if(fabs(v) < 1.0e-12) {
			ignore = true;
		}
	}

	// add entry into the NTuple
	if(!ignore) {
	    if(noSingles) 
		ntuple->doCallbacks(data,NDetectorFields);
	    else
		ntuple->appendRow(data,NDetectorFields); // calls doCallbacks()
	}

	if(verbose) 
		printf("detector %s: EDep=%.3f VisibleEdep=%.3f Ntracks=%.0f\n",
			name.c_str(),data[12]/MeV,data[13]/MeV,data[14]);

	clear();
}


void BLCMDdetector::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	if(radius > 0.0)
		generateTubs(npoints, innerRadius, radius, 0.0, 360.0*deg, length, v);
	else
		generateBox(npoints,width,height,length,v);
}

G4bool BLCMDdetector::isOutside(G4ThreeVector &local, G4double tolerance)
{
	if(radius > 0.0) {
		G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
		return r > radius-tolerance || r < innerRadius+tolerance ||
			fabs(local[2]) > length/2.0-tolerance;
	} else {
		return fabs(local[0]) > width/2.0-tolerance ||
			fabs(local[1]) > height/2.0-tolerance ||
			fabs(local[2]) > length/2.0-tolerance;
	}
}

