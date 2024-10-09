//	BLCMDsample.cc
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

#include "Polygon2D.h"

const char TrackFields[] =
    "x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight";
const unsigned NTrackFields = 12;

/**	class BLCMDsample - implements a sample (samples the beam to an NTuple)
 *	Each placement of this class generates an NTuple of the beam as it
 *	crosses the plane of the Sample object. This is therefore a
 *	"perfect" detector in that it does not perturb the beam at all, and
 *	intrinsically has the resolution of a float. Tracking is not affected
 *	in any way, and the two steps that bracket the sample plane are
 *	linearly interpolated.
 *
 *	The NTuple for a sample can be added to a BLCMDntuple by
 *	including a pattern that matches its name in the 'detectors'
 *	argument to the ntuple command.
 **/
class BLCMDsample : public BLElement {
	G4double radius;
	G4String polygon;
	G4int dir;
	G4double tolerance;
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
	BLCMDsample();

	/// Destructor.
	virtual ~BLCMDsample() { }

	/// clone()
	BLElement *clone() { return new BLCMDsample(*this); }

	/// Copy constructor.
	BLCMDsample(const BLCMDsample& r);

	/// commandName() returns "sample".
	G4String commandName() { return "sample"; }

	/// command() implements the sample command.
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

	/// construct() will construct the sample.
	/// Used for normal placements of a Sample object.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// This element has no volume
	virtual G4double getLength() { return 0.0; }
	virtual G4double getWidth() { return 0.0; }
	virtual G4double getHeight() { return 0.0; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		return G4ThreeVector(0.0,0.0,0.0);
	}

	/// isOK() returns true.
	virtual G4bool isOK() { return true; }

	/// generatePoints() from BLElement
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v)
		{ v.clear(); }

	/// isOutside() from BLElement
	G4bool isOutside(G4ThreeVector &local, G4double tolerance);
};

BLCMDsample defaultSample;

/**	class BLSampleNTuple implements an NTuple for a BLCMDsample.
 **/
class BLSampleNTuple : public BLManager::SteppingAction {
	G4String name;
	BLCoordinateType coordinateType;
	int noSingles;
	G4String require;
	int kill;
	BLTrackNTuple *ntuple;
	bool verbose;
	int previousEventID;
	int previousTrackID;
	long nHitsForThisTrack;
	G4ThreeVector point;
	G4ThreeVector normal;
	G4double radius;
	G4double tolerance;
	G4int dir;
	Polygon2D poly;
	BLCoordinateTransform *global2local;
	G4ThreeVector prevPos;
	friend class BLCMDsample;
public:
	/// constructor.
	BLSampleNTuple(G4String type, G4String category, 
		G4String _name, G4ThreeVector pos, G4RotationMatrix *rot,
		int _noSingles, G4String filename, G4String _require, 
		BLCoordinateType _coordinateType, int _kill, 
		G4double _tolerance, G4int _dir, G4double _radius,
		G4String polygon);

	bool needsReference() 
		{ return ntuple ? ntuple->needsReference() : false; }

	/// UserSteppingAction() from BLManager::SteppingAction.
	void UserSteppingAction(const G4Step *step);

#ifdef STUB
	/// copy a track, including info Geant4 does not copy.
	void copyTrack(G4Track &dest, const G4Track &src);
#endif //STUB
};

BLCMDsample::BLCMDsample() : BLElement()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Sample tracks into an NTuple.");
	setDescription("A sample generates an NTuple of any track when it\n"
		"crosses the local X-Y plane. By default either direction is "
		"sampled, but if dir=1 only Pz>0 is sampled, and if "
		"dir=-1 only Pz<0 is sampled. These are all in local "
		"coordinates. While this element does not define a volume, "
		"it must be placed via the 'place' command; the local origin "
		"is placed (and the local Z axis is rotated with the sample "
		"plane).\n\n"
		"The sample is taken for a disk defined by the specifiied "
		"radius from (x=0,y=0). If polygon is not empty, then it must "
		"be simple (edges do not intersect each other, no holes); only "
		"tracks inside it are sampled, and radius is set for the "
		"vertex furthest from the origin (the polygon need not include "
		"the origin, but sampling is most efficient if the local "
		"origin is near the center of the polygon). Remember that the "
		"+x direction is beam LEFT, +y is up, +z is along the beam.\n\n"
		"Note this element linearly interpolates the track to its "
		"plane, using pairs of steps which bracket the plane. So "
		"if tracks are not straight (EM fields, multiple scattering), "
		"be sure that maxStep in surrounding volumes is small enough "
		"for good results.\n\n"
		"This element does not affect tracking, but it can be placed "
		"along the side of a solid that does; tolerance is used to "
		"ensure that 0-length steps when entering a daughter volume "
		"will not generate multiple samples.\n\n"
		"This element must be placed (via the place command).\n\n"
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
	polygon = "";
	dir = 0;
	tolerance = 0.002*mm;
	noSingles = 0;
	format = "";
	filename = "";
	require = "";
	referenceParticle = 0;
	coordinates = "Centerline";
	kill = 0;
	require = "";
	coordinateType = BLCOORD_CENTERLINE;
}

BLCMDsample::BLCMDsample(const BLCMDsample& r) : BLElement(r)
{
	radius = r.radius;
	polygon = r.polygon;
	dir = r.dir;
	tolerance = r.tolerance;
	noSingles = r.noSingles;
	format = r.format;
	filename = r.filename;
	referenceParticle = r.referenceParticle;
	coordinates = r.coordinates;
	kill = r.kill;
	require = r.require;
	coordinateType = r.coordinateType;
}

int BLCMDsample::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("sample: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return handleNamedArgs(namedArgs);
	}

	BLCMDsample *t = new BLCMDsample(defaultSample);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	t->coordinateType = BLCoordinates::getCoordinateType(t->coordinates);

	// ascii->bltrackfile format, for accuracy and consistency of output
	for(unsigned i=0; i<t->format.size(); ++i)
		t->format[i] = tolower(t->format[i]);
	if(t->format == "ascii")
		t->format = "bltrackfile";

	t->print(argv[0]);

	return retval;
}

void BLCMDsample::defineNamedArgs()
{
	argDouble(radius,"radius","The radius of the circular Sample (mm).");
	argString(polygon,"polygon","The vertices of an X-Y polygon outlining the sample region (x1,y1;x2,y2;x3,y3;... mm).");
	argInt(dir,"dir","Required direction: 1=+Z, -1=-Z, 0=either (0).");
	argDouble(tolerance,"tolerance","Tracking tolerance (2 microns).");
	argInt(noSingles,"noSingles","Set to 1 to omit this NTuple "
			"(still can be used in the ntuple command) (0).");
	argString(format,"format","NTuple format: (see above for list).");
	argString(filename,"filename","filename ('' uses name to determine filename)");
	argString(filename,"file","alias for filename");
	argString(require,"require","Expression which must be nonzero to include the track (default=1)",false);
	argInt(referenceParticle,"referenceParticle","Set to 1 to include the Reference Particle.");
	argString(coordinates,"coordinates","Coordinates: global, local, centerline, or reference (default=c).");
	argInt(kill,"kill","Set to 1 kill all tracks after entering them into NTuple(s).");
}

void BLCMDsample::argChanged()
{
}

void BLCMDsample::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

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

	BLSampleNTuple *nt = new BLSampleNTuple(format,
			"Sample",thisname,globalPosition,globalRotation,
			noSingles,filename, require,coordinateType,kill,
			tolerance,dir,radius,polygon);

	BLManager::getObject()->registerBeamStep(0,nt);
	if(referenceParticle != 0 || nt->needsReference())
		BLManager::getObject()->registerReferenceParticleStep(0,nt);

	printf("BLCMDsample::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2]);
}

BLSampleNTuple::BLSampleNTuple(G4String type,
	G4String category, G4String _name, G4ThreeVector pos, 
	G4RotationMatrix *rot, int _noSingles,
	G4String filename, G4String _require,BLCoordinateType _coordinateType,
	int _kill, G4double _tolerance, G4int _dir, G4double _radius,
	G4String polygon) :
	poly(), prevPos(-DBL_MAX,-DBL_MAX,-DBL_MAX), global2local(0)

{
	name = _name;
	coordinateType = _coordinateType;
	global2local = new BLCoordinateTransform(rot,pos);
	noSingles = _noSingles;
	require = _require;
	kill = _kill;
	ntuple = BLTrackNTuple::create(type,category,name,filename,
			coordinateType,require,noSingles);
	verbose = BLManager::getObject()->getSteppingVerbose() > 0;
	previousEventID = previousTrackID = -999;
	nHitsForThisTrack = 0;
	point = pos;
	normal = G4ThreeVector(0.0,0.0,1.0); // z axis
	if(rot) normal = *rot * normal;
	tolerance = _tolerance;
	dir = _dir;
	radius = _radius;
	if(polygon.size() > 0) {
		radius = 0.0;
		std::vector<G4String> vs = BLCommand::splitString(polygon,
								";",true);
		std::vector<Point2D> vtx; // Point2D from Polygon2D.h
		for(unsigned i=0; i<vs.size(); ++i) {
			float x,y;
			if(sscanf(vs[i].c_str(),"%f,%f",&x,&y) != 2) continue;
			vtx.push_back(Point2D(x,y));
			float r = x*x+y*y;
			if(radius < r) radius = r;
		}
		if(vtx.size() < 3 || vtx.size() != vs.size())
			G4Exception("sample","invalid polygon",FatalException,
								polygon);
		poly.setVertices(vtx);
		radius += tolerance;
	}
}

void BLSampleNTuple::UserSteppingAction(const G4Step *step)
{
	if(BLManager::getObject()->getState() == SPECIAL) return;

	// only use reference coordinates when they are valid
	BLManagerState state = BLManager::getObject()->getState();
	if(coordinateType == BLCOORD_REFERENCE && state != BEAM) return;

	G4Track *track = step->GetTrack();
	G4StepPoint *prePoint = step->GetPreStepPoint();
	G4StepPoint *postPoint = step->GetPostStepPoint();

	// check if new track
	if(previousTrackID != track->GetTrackID() || 
	   previousEventID != BLManager::getObject()->getEventID()) {
		prevPos = G4ThreeVector(-DBL_MAX,-DBL_MAX,-DBL_MAX);
		previousTrackID = track->GetTrackID();
		previousEventID = BLManager::getObject()->getEventID();
	}

	// determine if this step brackets the plane, get interpolation factors
	double preDot = normal.dot(prePoint->GetPosition()-point);
	double postDot = normal.dot(postPoint->GetPosition()-point);
	double f=1.0; // interpolation factor (0 => pre, 1 => post)
	if(fabs(postDot) <= tolerance) {
		// this step ended right on the plane (well, within tolerance)
		f = 1.0; // this interpolation factor
	} else if(preDot*postDot >= 0.0) {
		// this step does not bracket the plane, not on the plane
		prevPos = G4ThreeVector(-DBL_MAX,-DBL_MAX,-DBL_MAX);
		return;
	} else {
		// this step brackets the plane
		double deltaZ = fabs(postDot - preDot);
		f = (deltaZ>tolerance ? fabs(preDot)/deltaZ : 1.0);
	}
	double g = 1.0 - f; // the other interpolation factor

	// check track direction
	if(dir != 0) {
		double d = normal.dot(track->GetMomentumDirection());
		if(dir < 0 && d > -1E-6 || dir > 0 && d < 1E-6) {
			prevPos = G4ThreeVector(-DBL_MAX,-DBL_MAX,-DBL_MAX);
			return;
		}
	}

	// get position on our plane
	G4ThreeVector pos(track->GetPosition() - g*step->GetDeltaPosition());
	BLAssert(fabs(normal.dot(pos-point)) < tolerance); // verify on plane

	// check if inside disk
	double r = (pos - point).mag();
	if(r > radius) {
		prevPos = G4ThreeVector(-DBL_MAX,-DBL_MAX,-DBL_MAX);
		return;
	}

	// check if inside polygon
	if(poly.getNvertices() > 0) {
		G4ThreeVector local, global = 
			track->GetPosition() - g*step->GetDeltaPosition();
		global2local->getLocal(local,global);
		BLAssert(fabs(local.z()) < tolerance);
		if(!poly.isInside(local.x(),local.y())) {
			prevPos = G4ThreeVector(-DBL_MAX,-DBL_MAX,-DBL_MAX);
			return;
		}
	}

	// check if duplicate of previous step
	if((pos-prevPos).mag() < tolerance) return;
	prevPos = pos;

	// handle local coordinates
	BLCoordinates *coord = dynamic_cast<BLCoordinates*>(
						track->GetUserInformation());
	if(coord) coord->setLocalTransform(global2local);

	// linear interpolation between step points
	// (Creating a new G4Track has issues, so use *track, then restore it.)
	if(f != 1.0) {
		// get current track values
		G4ThreeVector position = track->GetPosition();
		G4ThreeVector mom = track->GetMomentum();
		G4double time = track->GetGlobalTime();
		G4double ke = track->GetKineticEnergy();
		G4double properTime = track->GetProperTime();
		G4double trackLength = track->GetTrackLength();
		G4double localTime = track->GetLocalTime();
		G4ThreeVector pol = track->GetPolarization();
		// interpolate linearly to z, global coords
		G4ThreeVector deltaPos=step->GetDeltaPosition();
		G4ThreeVector deltaMom=(postPoint->GetMomentum() -
						prePoint->GetMomentum());
		G4double deltaTime=step->GetDeltaTime();
		G4double deltaE=(postPoint->GetKineticEnergy() -
						prePoint->GetKineticEnergy());
		G4double deltaProperTime=(postPoint->GetProperTime() -
						prePoint->GetProperTime());
		G4double deltatrackLength=step->GetStepLength();
		G4ThreeVector deltaPol=(postPoint->GetPolarization() -
						prePoint->GetPolarization());
		track->SetPosition(position-g*deltaPos);
		track->SetMomentumDirection((mom-g*deltaMom).unit());
		track->SetGlobalTime(time-g*deltaTime);
		track->SetKineticEnergy(ke-g*deltaE);
		track->SetProperTime(properTime-g*deltaProperTime);
		track->AddTrackLength(-g*deltatrackLength);
		track->SetLocalTime(localTime-g*deltaTime);
		track->SetPolarization(pol-g*deltaPol);
		BLCoordinates::update(track);
		// add track to the NTuple
		ntuple->appendTrack(track);
		// restore current values to track
		track->SetPosition(position);
		track->SetMomentumDirection(mom.unit());
		track->SetGlobalTime(time);
		track->SetKineticEnergy(ke);
		track->SetProperTime(properTime);
		track->AddTrackLength(g*deltatrackLength);
		track->SetLocalTime(localTime);
		track->SetPolarization(pol);
		BLCoordinates::update(track);
	} else {
		// add track to the NTuple
		ntuple->appendTrack(track);
	}

	// remove local coordinates
	if(coord) coord->setLocalTransform(0);

	if(verbose) {
		printf("sample '%s' sampled the track at (%.3f,%.3f,%.3f)\n",
			name.c_str(),pos.x()/mm,pos.y()/mm,pos.z()/mm);
	}

	if(kill) {
		track->SetTrackStatus(fStopAndKill);
		if(verbose) printf("Track killed by '%s' with kill=1\n",
								name.c_str());
	}
}

G4bool BLCMDsample::isOutside(G4ThreeVector &local, G4double tolerance)
{
	return true;
}

#ifdef STUB
void BLSampleNTuple::copyTrack(G4Track &dest, const G4Track &src)
{
	dest.CopyTrackInfo(src);
	dest.SetTrackID(src.GetTrackID());
	dest.SetParentID(src.GetParentID());
	dest.SetCreatorProcess(src.GetCreatorProcess());
	dest.SetUserInformation(src.GetUserInformation());
}
#endif //STUB

