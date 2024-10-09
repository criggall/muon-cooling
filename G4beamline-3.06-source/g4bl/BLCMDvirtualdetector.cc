//	BLCMDvirtualdetector.cc
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

#include "BLElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLTrackNTuple.hh"
#include "BLCoordinates.hh"
#include "BLCoordinateTransform.hh"
#include "BLGlobalField.hh"
#include "BLKillTrack.hh"

const char TrackFields[] =
    "x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight";
const unsigned NTrackFields = 12;

/**	class BLCMDvirtualdetector - implements a VirtualDetector
 *	Each placement of this class generates an NTuple of the beam as it
 *	enters the physical volume of the VirtualDetector. This is therefore a
 *	"perfect" detector in that it does not perturb the beam at all, and
 *	intrinsically has the resolution of a float (it measures position
 *	and momentum and time).
 *
 *	The NTuple for a virtualdetector can be added to a BLCMDntuple by
 *	including a pattern that matches its name in the 'detectors'
 *	argument to the ntuple command.
 *
 *	Note that if a BLCoordinates instance is linked into the track,
 *	its centerline coordinates are used; otherwise global coordinates
 *	are used.
 **/
class BLCMDvirtualdetector : public BLElement {
	G4double radius;
	G4double innerRadius;
	G4double height;
	G4double width;
	G4double length;
	G4String material;
	G4String color;
	G4VSolid *solid;
	G4double maxStep;
	G4int noSingles;
	G4String format;
	G4String filename;
	G4String require;
	G4int referenceParticle;
	G4String coordinates;
	G4int kill;
	BLCoordinateType coordinateType;
public:
	/// Default constructor.
	BLCMDvirtualdetector();

	/// Destructor.
	virtual ~BLCMDvirtualdetector();

	/// clone()
	BLElement *clone() { return new BLCMDvirtualdetector(*this); }

	/// Copy constructor.
	BLCMDvirtualdetector(const BLCMDvirtualdetector& r);

	/// commandName() returns "virtualdetector".
	G4String commandName() { return "virtualdetector"; }

	/// command() implements the virtualdetector command.
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

	/// construct() will construct the virtualdetector.
	/// Used for normal placements of a VirtualDetector object.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns this element's Length along the Z axis.
	virtual G4double getLength() { return length; }

	/// getWidth() returns this element's Width along the X axis.
	virtual G4double getWidth() { return width; }

	/// getHeight() returns this element's height along the Y axis.
	virtual G4double getHeight() { return height; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	virtual G4bool isOK() { return true; }

	/// generatePoints() from BLElement
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isOutside() from BLElement
	G4bool isOutside(G4ThreeVector &local, G4double tolerance);
};

BLCMDvirtualdetector defaultVirtualDetector;

/**	class BLVirtualDetectorNTuple implements an NTuple for a BLCMDvirtualdetector.
 **/
class BLVirtualDetectorNTuple : public BLManager::SteppingAction {
	G4String name;
	G4VPhysicalVolume *thisVol;
	BLCoordinateType coordinateType;
	int noSingles;
	G4String require;
	int kill;
	BLTrackNTuple *ntuple;
	bool verbose;
	int prevEventID;
	int prevTrackID;
	BLCoordinateTransform *global2local;
	long nHitsForThisTrack;
	friend class BLCMDvirtualdetector;
public:
	/// constructor.
	BLVirtualDetectorNTuple(G4String type, G4String category, 
		G4String _name, G4VPhysicalVolume *pv, int _noSingles,
		G4String filename, G4String _require, 
		BLCoordinateType _coordinateType, int _kill);

	bool needsReference() 
		{ return ntuple ? ntuple->needsReference() : false; }

	/// UserSteppingAction() from BLManager::SteppingAction.
	void UserSteppingAction(const G4Step *step);
};

BLCMDvirtualdetector::BLCMDvirtualdetector() : BLElement()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Construct a VirtualDetector that generates an NTuple.");
	setDescription("A VirtualDetector generates an NTuple of any track when it\n"
		"enters the physical volume of the VirtualDetector. It may\n"
		"be placed via multiple place commands (usually with a "
		"'rename=det#' argument to distinguish the different placements).\n"
		"If material is not specified, it uses the material of the\n"
		"enclosing element. Every placement creates an individual NTuple.\n"
		"For a circular VirtualDetector give radius; for a\n"
		"rectangular one give height and width; length is usually\n"
		"left at 1 mm, but can be set to correspond to the length "
		"of a physical detector.\n"
		"The NTuple by default uses centerline coordinates.\n"
		"The NTuple of the virtualdetector can be included in an "
		"ntuple command by including a pattern that matches its name "
		"in the 'detectors' argument to the ntuple command. Note that "
		"must match the name as placed (i.e. includes rename=), not "
		"the name given to this command. The noSingles argument may be "
		"useful in this case to avoid a huge NTuple of singles (an "
		"empty NTuple may be created).\n\n"
		"This element must be placed (via the place command).\n\n"
		"Note that secondary particles created within the virtualdetector "
		"will not get an entry until they have taken one step. "
		"They are guaranteed to do so.\n\n"
		"The standard NTuple fields are:\n"
		"    x,y,z (mm)\n"
		"    Px,Py,Pz (MeV/c)\n"
		"    t (ns)\n"
		"    PDGid (11=e-, 13=mu-, 22=gamma, 211=pi+, 2212=proton, ...)\n"
		"    EventID (may be inexact above 16,777,215)\n"
		"    TrackID\n"
		"    ParentID (0 => primary particle)\n"
		"    Weight (defaults to 1.0)\n\n"
		"The following additional fields are appended for "
		"format=Extended, format=asciiExtended, and "
		"format=rootExtended:\n"
		"    Bx, By, Bz (Tesla)\n"
		"    Ex, Ey, Ez (Megavolts/meter)\n"
		"    ProperTime (ns)\n"
		"    PathLength (mm)\n"
		"    PolX, PolY, PolZ (polarization)\n"
		"    InitX, initY, InitZ (initial position, mm)\n"
		"    InitT (initial time, ns)\n"
		"    InitKE (MeV when track was created)\n\n"
		"Valid Formats (ignore case): ");
	
	// initial field values
	radius = 0.0;
	innerRadius = 0.0;
	height = 0.0;
	width = 0.0;
	length = 1.0*mm;
	material = "";
	color = "1,1,1";
	solid = 0;
	maxStep = -99.0;
	noSingles = 0;
	format = "";
	filename = "";
	require = "";
	referenceParticle = 0;
	coordinates = "Centerline";
	kill = 0;
	coordinateType = BLCOORD_CENTERLINE;
}

BLCMDvirtualdetector::~BLCMDvirtualdetector()
{
	if(solid) delete solid;
}

BLCMDvirtualdetector::BLCMDvirtualdetector(const BLCMDvirtualdetector& r) : BLElement(r)
{
	radius = r.radius;
	innerRadius = r.innerRadius;
	height = r.height;
	width = r.width;
	length = r.length;
	material = r.material;
	color = r.color;
	solid = r.solid;
	maxStep = r.maxStep;
	noSingles = r.noSingles;
	format = r.format;
	filename = r.filename;
	referenceParticle = r.referenceParticle;
	coordinates = r.coordinates;
	kill = r.kill;
	require = r.require;
	coordinateType = r.coordinateType;
}

int BLCMDvirtualdetector::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("virtualdetector: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return handleNamedArgs(namedArgs);
	}

	BLCMDvirtualdetector *t = new BLCMDvirtualdetector(defaultVirtualDetector);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	t->coordinateType = BLCoordinates::getCoordinateType(t->coordinates);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	// ascii->bltrackfile format, for accuracy and consistency of output
	for(unsigned i=0; i<t->format.size(); ++i)
		t->format[i] = tolower(t->format[i]);
	if(t->format == "ascii")
		t->format = "bltrackfile";

	t->print(argv[0]);

	return retval;
}

void BLCMDvirtualdetector::defineNamedArgs()
{
	argDouble(radius,"radius","The radius of the circular VirtualDetector (mm).");
	argDouble(innerRadius,"innerRadius","The inner radius of the circular VirtualDetector (0 mm, solid).");
	argDouble(height,"height","The height of the rectangular VirtualDetector (mm).");
	argDouble(width,"width","The width of the rectangular VirtualDetector (mm).");
	argDouble(length,"length","The length of the VirtualDetector (mm).");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm).");
	argString(material,"material","The material of the VirtualDetector.");
	argString(color,"color","The color of the VirtualDetector (''=invisible).");
	argInt(noSingles,"noSingles","Set to 1 to omit the NTuple for singles.");
	argString(format,"format","NTuple format: (see above for list).");
	argString(filename,"filename","filename ('' uses name to determine filename)");
	argString(filename,"file","alias for filename");
	argString(require,"require","Expression which must be nonzero to include the track (default=1)",false);
	argInt(referenceParticle,"referenceParticle","Set to 1 to include the Reference Particle.");
	argString(coordinates,"coordinates","Coordinates: global, centerline, or reference (default=c).");
	argInt(kill,"kill","Set to 1 kill all tracks after entering them into NTuple(s).");
}

void BLCMDvirtualdetector::argChanged()
{
	if(radius > 0.0) width = height = radius*2.0;
}

void BLCMDvirtualdetector::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4Material *mat;
	if(material != "")
		mat = getMaterial(material);
	else
		mat = parent->GetMaterial();

	G4String thisname = parentName+getName();

	if(!solid) {
		if(radius > 0.0) {
			solid = new G4Tubs(thisname+"Tubs", innerRadius, radius,
					length/2.0, 0.0, 2.0*pi);
		} else if(height > 0.0 && width > 0.0) {
			solid = new G4Box(thisname+"Box",width/2.0,
					height/2.0,length/2.0);
		} else {
			printError("virtualdetector::construct %s INVALID - no "
				"radius or height&width",thisname.c_str());
			return;
		}
	}
	G4LogicalVolume *lv = new G4LogicalVolume(solid,mat, thisname+"LogVol");
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

	BLVirtualDetectorNTuple *nt = new BLVirtualDetectorNTuple(format,
			"VirtualDetector",thisname,pv,noSingles,filename,
			require,coordinateType,kill);

	BLManager::getObject()->registerBeamStep(pv,nt);
	if(referenceParticle != 0 || nt->needsReference())
		BLManager::getObject()->registerReferenceParticleStep(pv,nt);

	// special case: if kill!=0 and viewer != "none", then register a
	// BLKillTrack. (in viewer mode we don't accumulate NTuples.)
	if(kill != 0 && Param.getString("viewer") != "none") {
		BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));
	}

	printf("BLCMDvirtualdetector::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2]);
}

BLVirtualDetectorNTuple::BLVirtualDetectorNTuple(G4String type,
	G4String category, G4String _name, G4VPhysicalVolume *pv, int _noSingles,
	G4String filename, G4String _require,BLCoordinateType _coordinateType,
	int _kill)

{
	name = _name;
	thisVol = pv;
	noSingles = _noSingles;
	coordinateType = _coordinateType;
	require = _require;
	ntuple = BLTrackNTuple::create(type,category,name,filename,
			coordinateType,require,noSingles);
	kill = _kill;
	verbose = BLManager::getObject()->getSteppingVerbose() > 0;
	prevEventID = prevTrackID = -999;
	nHitsForThisTrack = 0;
	global2local = 0;
	if(coordinateType == BLCOORD_LOCAL)
		global2local = new BLCoordinateTransform(pv->GetRotation(),
							pv->GetTranslation());
}

void BLVirtualDetectorNTuple::UserSteppingAction(const G4Step *step)
{
	if(BLManager::getObject()->getState() == SPECIAL) return;

	// only use reference coordinates when they are valid
	BLManagerState state = BLManager::getObject()->getState();
	if(coordinateType == BLCOORD_REFERENCE && 
				(state != BEAM && state != VISUAL)) return;

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
	
	// return if not entering thisVol
	if(preVol == postVol || postVol != thisVol) return;

	BLCoordinates *coord = dynamic_cast<BLCoordinates*>(
						track->GetUserInformation());
	if(coord) {
		coord->setLocalTransform(global2local);
	}

	// add empty line between tracks, if >= 4 hits for this track
	if(BLManager::getObject()->getExternalTrackID(track) == prevTrackID && 
			BLManager::getObject()->getEventID() == prevEventID) {
		++nHitsForThisTrack;
	} else {
		if(nHitsForThisTrack >= 4) ntuple->annotate("");
		nHitsForThisTrack = 0;
		prevTrackID = BLManager::getObject()->getExternalTrackID(track);
		prevEventID = BLManager::getObject()->getEventID();
	}

	ntuple->appendTrack(track);
	if(coord) coord->setLocalTransform(0);

	if(kill) {
		track->SetTrackStatus(fStopAndKill);
		if(verbose) printf("Track killed by '%s' with kill=1\n",
			name.c_str());
	}
}

void BLCMDvirtualdetector::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	if(radius > 0.0)
		generateTubs(npoints, innerRadius, radius, 0.0, 360.0*deg, length, v);
	else
		generateBox(npoints,width,height,length,v);
}

G4bool BLCMDvirtualdetector::isOutside(G4ThreeVector &local, G4double tolerance)
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

