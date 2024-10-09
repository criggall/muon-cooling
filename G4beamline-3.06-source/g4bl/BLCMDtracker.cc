//	BLCMDtracker.cc
//	NOTE: This file includes three command classes:
//	BLCMDtracker, BLCMDtrackerplane, and BLCMDtrackermode.
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
int BLCMDtracker_dummy=0;
#else

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

#include "BLAssert.hh"
#include "BLElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLRunManager.hh"
#include "BLNTuple.hh"
#include "BLCoordinates.hh"
#include "BLEvaluator.hh"
#include "BLGlobalField.hh"
#include "BLMinimize.hh"

static const G4double UNDEFINED = -3.7e21;

class TrackerPlaneInstance;
class BLCMDtracker;
static const int NO_HIT=0x80000000;

static const char *HIT_FIELDS="true_x:true_y:true_z:true_Px:true_Py:true_Pz:true_t:true_PDGid:true_EventID:true_TrackID:true_ParentID:true_Weight:truereport_x:truereport_y:truereport_z:truereport_Px:truereport_Py:truereport_Pz:truereport_t";
#define  N_HIT_FIELDS 19

static const char *FIT_FIELDS="x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight:ChisqPerDF:nDF:nHit:nIter:true_x:true_y:true_z:true_Px:true_Py:true_Pz:true_t";
#define  N_FIT_FIELDS 23

static const char *FOR009_FIELDS="x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight:Bx:By:Bz:Ex:Ey:Ez";
#define  N_FOR009_FIELDS 18


/** enum BLTrackerMode specifies the mode of the tacker and all its planes
 **/
enum BLTrackerMode {BLTRACKER_TRUE, BLTRACKER_FIT, BLTRACKER_IGNORE};

/**	class BLCMDtracker implements a tracker that can fit tracks to
 *	wire hits and timing in its trackerplanes.
 *
 *	Modes:	true	tracks the "true" track -- the trackerplane-s
 *			with wires report which wire was hit, and the
 *			trackerplane-s with timing report the time.
 *		fit	fits a track to the "true" track -- the
 *			trackerplane-s report their conrtibution to chisq.
 *		ignore	the track is ignored
 *
 **/
class BLCMDtracker : public BLCommand, public BLManager::TrackingAction,
				public BLManager::ZSteppingAction,
				public BLCallback, public BLMinimizeFunction {
	struct FitParam {
		std::vector<G4String> name;
		std::vector<double> scale;
		std::vector<double> sumDX;
		std::vector<double> sumDX2;
		int nDelta;
		void define(G4String _name, double _scale) {
			name.push_back(_name);
			scale.push_back(_scale);
			sumDX.push_back(0.0);
			sumDX2.push_back(0.0);
			nDelta = 0; // (doing this multiple times is OK)
		}
		int size() { return name.size(); }
		void setDelta(G4String _name, double dx) {
			for(unsigned i=0; i<name.size(); ++i) {
				if(_name != name[i]) continue;
				sumDX[i] += dx;
				sumDX2[i] += dx*dx;
				if(i == 0) ++nDelta;
				break;
			}
		}
		void printSummary() {
			printf("\nSummary of track fitting parameters (fit-true) for %d tracks:\n"
			  "  Name     Mean      Sigma      Scale    Flag\n"
			  "  ----   --------   --------   -------- -------\n",
			  nDelta);
			for(unsigned i=0; i<name.size(); ++i) {
				double mean = sumDX[i]/nDelta;
				double sigma = sqrt(fabs(sumDX2[i]/nDelta
					- mean*mean));
				const char *p = "";
				if(scale[i] != 0 && (sigma/scale[i]<0.2 ||
				   sigma/scale[i]>5.0)) p = "RESCALE";
				printf("  %4s %10.6f %10.6f %10.6f %s\n",
					name[i].c_str(),mean,sigma,scale[i],p);
			}
			printf("\n");
		}
		G4String getNameList() {
			G4String list=name[0];
			for(unsigned i=1; i<name.size(); ++i) {
				list += ":";
				list += name[i];
			}
			return list;
		}
	};
	FitParam fitParam;
	G4String name;
	G4double trackerZ;
	G4double reportZ;
	bool reportOnly;
	G4double scaleX;
	G4double scaleXp;
	G4double scalePtot;
	G4double scaleT;
	G4double minPz;
	G4double tolerance;
	G4int maxIter;
	G4int verbose;
	G4String format;
	G4String filename;
	G4int for009;
	BLTrackerMode mode;
	std::vector<TrackerPlaneInstance*> plane;
	std::vector<int> hitWire;
	std::vector<double> hitTime;
	double chisq;
	int ndf;
	G4ThreeVector truePosition;
	G4ThreeVector trueMomentum;
	G4double trueTime;
	G4int trueTrackID;
	G4int trueParentID;
	G4ParticleDefinition *trueDefinition;
	G4ThreeVector reportPosition;
	G4ThreeVector reportMomentum;
	G4double reportTime;
	G4int trueHits;
	G4int fitHits;
	int minHits;
	G4RotationMatrix rotationMatrix;
	G4ThreeVector beamPosition;
	BLRunManager *runmgr;
	BLManager *mgr;
	BLNTuple *ntuple_hit;
	BLNTuple *ntuple_fit;
	BLNTuple *for009_fit;
	G4bool chisqFailure;
	G4String trackFields_hit;
	G4String trackFields_fit;
	friend class TrackerPlaneInstance;
public:
	static std::vector<BLCMDtracker*> list;
public:
	/// Constructors.
	BLCMDtracker();
	BLCMDtracker(BLCMDtracker &r);

	/// commandName() returns "tracker".
	virtual G4String commandName() { return "tracker"; }
	
	/// command() implements the tracker command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	virtual G4String getName() const { return name; }
	virtual void setName(G4String _name) { name = _name; }

	/// findTracker() returns a pointer to the named tracker.
	/// returns NULL if not found.
	static BLCMDtracker *findTracker(G4String name);

	/// registerTrackerPlane() registers a BLCMDtrackerplane with this
	/// tracker. Returns the plane id #.
	int registerTrackerPlane(TrackerPlaneInstance *plane,
							const G4String name);

	/// planeHit() tells the tracker that a plane was hit by the current
	/// G4Track being tracked. Called only when mode=BLTRACKER_TRUE.
	/// wire and/or time should be zero if this plane does not measure them.
	void planeHit(int planeID, int wire, double time);

	/// delta() tells the tracker of this track's offset from the hit wire,
	/// divided by sigma. Called only when mode=BLTRACKER_FIT.
	/// Returns true to continue tracking, false to have caller kill the track.
	bool delta(int planeID, double deltaPos, double deltaTime, 
					bool validWire, bool validTime);

	/// fitTrack() fits a track to the true track.
	void fitTrack();

	/// setMode() sets the tracker mode.
	void setMode(BLTrackerMode _mode) { mode = _mode; }

	/// virtual void UserZSteppingAction() from BLManager::ZSteppingAction.
	/// Handles the true track in mode=true, and handles the reference
	/// particle in all modes.
	virtual void UserZSteppingAction(const G4Track *track);

	/// PreUserTrackingAction() from BLManager::TrackingAction.
	/// Initializes stuff.
	virtual void PreUserTrackingAction(const G4Track *track);

	/// PostUserTrackingAction() from BLManager::TrackingAction.
	/// Writes the true track in mode=true, only if it hit all planes.
	virtual void PostUserTrackingAction(const G4Track *track);

	/// callback() from BLCallback.
	/// Constructs NTuples --  must wait until all planes are defined
	/// and mode is known (trackermode could follow this command)
	void callback(int type);

	/// constructTrack() will construct a track from the parameters.
	/// The caller must delete the track;
	G4Track *constructTrack(const std::vector<double> &x);

	/// operator() from BLMinimizeFunction.
	/// Computes the chisq by tracking the trake constructed from the
	/// current parameters.
	virtual double operator()(const std::vector<double> &x);

	/// handlePreviousTracks() reads a previously-written NTuple contining
	/// 'true' tracks, and fits tracks to them.
	void handlePreviousTracks(const G4String filename);
};
BLCMDtracker defaultTracker;
std::vector<BLCMDtracker*> BLCMDtracker::list;


/**	class BLCMDtrackerplane - implements one plane of a tracker.
 *
 *	A tracker plane is a detector consisting of either:
 *	  A) a large number of parallel wires arranged in a plane
 *	  B) a counter that gives timing information
 *	  C) both (A) and (B)
 *	For a given track hitting the plane, one and only one wire is hit.
 *	Arguments to the command specify the wire
 *	spacing, orientation, and offset; the number of wires is implicit
 *	from their spacing and orientation and the size of the plane.
 *	if wireSpacing=0, no wires are present.
 *
 *	If sigmaTime>0 this plane measures the time of the track. A
 *	Gaussian random number with that sigma is added to the true value.
 *
 *	This class works with the BLCMDtracker class to permit the fitting of
 *	tracks to the wire hits. modes:
 *	true	reports the hit wire # to the BLCMDtracker instance
 *		(this is G4beamline tracking "true" tracks)
 *	fit	computes the plane's contribution to Chisq for the
 *		track and reports it to the BLCMDtracker instance
 *	ignore	track is ignored
 **/
class BLCMDtrackerplane : public BLElement {
	G4String tracker;
	G4double radius;
	G4double innerRadius;
	G4double height;
	G4double width;
	G4double length;
	G4String material;
	G4String color;
	G4VSolid *solid;
	G4double maxStep;
	G4double theta;
	G4double wireSpacing;
	G4double wireOffset;
	G4String errType;
	G4double errTheta;
	G4double errSpacing;
	G4double errOffset;
	G4double sigmaTime;
	friend class TrackerPlaneInstance;
public:
	/// Default constructor.
	BLCMDtrackerplane();

	/// Destructor.
	virtual ~BLCMDtrackerplane();

	/// clone()
	BLElement *clone() { return new BLCMDtrackerplane(*this); }

	/// Copy constructor.
	BLCMDtrackerplane(const BLCMDtrackerplane& r);

	/// commandName() returns "trackerplane".
	G4String commandName() { return "trackerplane"; }

	/// command() implements the trackerplane command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// argChanged() does internal computations after some arg changed
	void argChanged();

	/// construct() will construct the trackerplane.
	/// Used for normal placements of a tracker plane object.
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
BLCMDtrackerplane defaultTrackerPlane;

/**	class TrackerPlaneInstance implements one placement of a tracker 
 *	plane.
 **/
class TrackerPlaneInstance : public BLManager::SteppingAction {
	G4String name;
	G4VPhysicalVolume *thisVol;
	BLTrackerMode mode;
	BLCMDtracker *tracker;
	int planeID;
	int hitWire;
	double hitTime;
	bool haveBeenHit;
	double sinTheta, cosTheta;
	double wireSpacing;
	double wireOffset;
	double sigmaWire;
	double sigmaTime;
	double sinThetaErr, cosThetaErr;
	double wireSpacingErr;
	double wireOffsetErr;
public:
	/// constructor.
	TrackerPlaneInstance(G4String _name, G4VPhysicalVolume *pv,
				BLCMDtrackerplane *plane);

	G4String getName() { return name; }

	/// UserSteppingAction() from BLManager::SteppingAction.
	void UserSteppingAction(const G4Step *step);

	void setTrueMode() { mode = BLTRACKER_TRUE; setNotHit(); }
	void setIgnoreMode() { mode = BLTRACKER_IGNORE; setNotHit(); }
	void setFitMode(int w, double t) 
	  { mode = BLTRACKER_FIT; hitWire = w; hitTime = t; haveBeenHit=false; }
	void setNotHit() { hitWire=NO_HIT; hitTime=0.0; haveBeenHit=false; }
	bool isHit() { return haveBeenHit; }
};

/**	class BLCMDtrackermode sets the mode for all trackers, and
 *	manages track fitting.
 *
 *	modes:
 *		true	normal g4beamline operation; the "true" track is
 *			tracked. Same as if no trackermode command was given.
 *		fit	special BLRunManager mode in which the Root
 *			file from an earlier run is read for the TrackerHits
 *			NTuple, and a track is fit to each such track.
 *		both	special BLRunManager mode in which a single event
 *			is simulated, and then a track is fitted to the
 *			first track in each tracker. Primarily for testing.
 **/
class BLCMDtrackermode : public BLCommand, public BLCallback {
	G4String filename;
	G4String mode;
	G4bool registered;
public:
	BLCMDtrackermode();

	G4String commandName() { return "trackermode"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	/// callback() from BLCallback.
	/// Implements the event loop for mode=fit and mode=both.
	void callback(int type);

	const G4String& getMode() const { return mode; }
};
BLCMDtrackermode trackermodeInstance;


BLCMDtracker::BLCMDtracker() : BLCommand(), BLManager::TrackingAction(),
				BLManager::ZSteppingAction(),
				BLCallback(), BLMinimizeFunction(),
				plane(), hitWire(), truePosition(), 
				trueMomentum(), reportPosition(),
				reportMomentum(), rotationMatrix(),
				beamPosition()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Defines a tracker.");
	setDescription("A tracker consists of several trackerplane-s and can "
		"fit a track to wire hits and times in the trackerplanes. "
		"This is a simple algorithm that does not handle backgrounds "
		"or multiple hits. It assumes "
		"that every track hits each trackerplane at most once. "
		"It is intended "
		"to be used to explore resolutions and the effects of survey "
		"errors. A tracker is a logical combination of its "
		"trackerplanes -- the tracker cannot be placed, but its "
		"trackerplanes must be placed into the system.\n\n"
		"The fitting algorithm used requires that all of its "
		"parameters have comparable scales, so the 'scaleX', "
		"'scaleXp', 'scalePtot', and 'scaleT' arguments should be set "
		"to the approximate sigmas "
		"of the tracker. They should be within a factor of 10 of "
		"the actual values, but closer is better. At the end of "
		"fitting tracks a summary is printed that flags each "
		"parameter with 'RESCALE' if its scale is too different from "
		"its sigma (factor of 5 or more).\n\n"
		"NOTE: if the tracker cannot measure Ptot, then 'scalePtot' "
		"MUST be set to zero. If the tracker cannot measure T, then "
		"scaleT MUST be set to zero. Parameters with zero scales are "
		"held fixed at their true-track values.\n\n"
		"NOTE: the trackerplane-s of a tracker MUST be placed in the "
		"order that particles will hit them; the code does not sort "
		"them. Usually this means that each place command of a "
		"trackerplane must have a larger z value than the previous "
		"place command. They must also come after the trackerZ value "
		"of the tracker command.\n\n"
		"A tracker has 3 modes of operation:\n"
		"   true      Each track of each event is written to a\n"
		"             TrackerHits NTuple if it hits all trackerplanes;\n"
		"             the NTuple includes both the true track values\n"
		"             and the individual wire hits for each trackerplane.\n"
		"   fit       A track is fit to the wire hits from a previous\n"
		"             run, and the fit track is written to a TrackerFit\n"
		"             NTuple (includes true values).\n"
		"   ignore    Any track is ignored.\n\n"
		"See the 'trackermode' command to control the mode of trackers.\n\n"
		"The TrackerHits NTuple written in 'true' mode contains:\n"
		"  true_x, true_y, true_z, true_Px, true_Py, true_Pz, true_t, \n"
		"  true_PDGid, true_EventID, true_TrackID, true_ParentID, \n"
		"  true_Weight, ... plus 1 hit and 1 time per trackerplane.\n"
		"  (the first 12 are the same as a BLTrackFile.)\n\n"
		"The TrackerFit Ntuple written in 'fit' mode contains:\n"
		"  x, y, z, Px, Py, Pz, t, PDGid, EventID, TrackID, ParentID,\n"
		"  Weight, ChisqPerDF, nDF, nHit, nIter, true_x, true_y, true_z, true_Px,\n"
		"  true_Py, true_Pz, true_t.\n"
		"(the first 12 are from the fit and are the same as a "
		"BLTrackFile.)\n\n"
		"The parameters of the fit are: x, y, dxdz, dydz, Ptot, time. "
		"You must ensure that there are at least as many data points "
		"as free parameters (scaleT=0 fixes time; scalePtot=0 fixes "
		"Ptot). Each trackerplane with nonzero wireSpacing provides "
		"a data point; each trackerplane with nonzero sigmaT provides "
		"a data point; trackerplanes that measure both provide two "
		"data points. You must have enough trackerplanes to meet this "
		"requirement, and must set minHits large enough to meet it. "
		"The TrackerFit NTuple has a field nDF that gives the number "
		"of degrees of freedom for the fit, which is defined as "
		"(#DataPoints)-(#FreeParameters); it also has nHit which gives "
		"the number of trackerplane-s hit.\n\n"
		"Note the tracker can simulate survey errors -- see the "
		"'trackerplane' command for details (each plane can have "
		"different errors).\n\n"
		"Both the true and the fit tracks are reported at reportZ, "
		"which defaults to the Z position of the tracker.\n\n"
		"Note that beamlossntuple and newparticlentuple will get many "
		"entries per track when used in mode=fit -- the fit runs the "
		"track many times through the tracker (30-100, up to "
		"maxIter).\n\n"
		"NOTE: the trackermode command must preceed all tracker "
		"commands in the input file, and each tracker command must "
		"preceed all of its trackerplane commands. The trackerZ value "
		"must also preceed all trackerplane-s, but reportZ can be "
		"equal to or anywhere after trackerZ.\n\n"
		"Note that each trackerplane must have a unique name. This "
		"means you should either have a separate trackerplane command "
		"for each one (with unique name), or use the rename= argument "
		"to the place command (again with unique name). If you use "
		"groups for trackerplane-s, use rename=+ in the group.\n\n"
		"NOTE: This command does not work properly in collective "
		"tracking mode."
	);

	name = "default";
	trackerZ = UNDEFINED;
	reportZ = UNDEFINED;
	reportOnly = false;
	scaleX = 1.0*mm;
	scaleXp = 0.001;
	scalePtot = 0.0;
	scaleT = 0.0;
	minPz = 10.0*MeV;
	tolerance = 0.01;
	maxIter = 200;
	verbose = 0;
	format = "";
	filename = "";
	for009 = 0;
	chisq = 0.0;
	ndf = 0;
	trueTime = 0.0;
	trueTrackID = -9999;
	trueParentID = -9999;
	trueDefinition = 0;
	reportTime = UNDEFINED;
	trueHits = -9999;
	fitHits = -9999;
	minHits = -1;
	runmgr = 0;
	mgr = 0;
	ntuple_hit = 0;
	ntuple_fit = 0;
	for009_fit = 0;
	chisqFailure = false;
	trackFields_hit = HIT_FIELDS;
	trackFields_fit = FIT_FIELDS;
}

BLCMDtracker::BLCMDtracker(BLCMDtracker &r) : BLCommand(r), 
	BLManager::TrackingAction(r), BLManager::ZSteppingAction(),
	BLCallback(r), BLMinimizeFunction(r),
	plane(), hitWire(), truePosition(), trueMomentum(), reportPosition(),
	reportMomentum(), rotationMatrix(), beamPosition()
{
	name = "";
	trackerZ = r.trackerZ;
	reportZ = r.reportZ;
	reportOnly = r.reportOnly;
	scaleX = r.scaleX;
	scaleXp = r.scaleXp;
	scalePtot = r.scalePtot;
	scaleT = r.scaleT;
	minPz = r.minPz;
	tolerance = r.tolerance;
	maxIter = r.maxIter;
	verbose = r.verbose;
	format = r.format;
	filename = r.filename;
	for009 = r.for009;
	mode = BLTRACKER_IGNORE;
	chisq = 0.0;
	ndf = 0;
	trueTime = 0.0;
	trueTrackID = -9999;
	trueParentID = -9999;
	trueDefinition = 0;
	reportTime = UNDEFINED;
	trueHits = -9999;
	fitHits = -9999;
	minHits = -1;
	runmgr = 0;
	mgr = 0;
	ntuple_hit = 0;
	ntuple_fit = 0;
	for009_fit = 0;
	chisqFailure = false;
	trackFields_hit = r.trackFields_hit;
	trackFields_fit = r.trackFields_fit;
}

int BLCMDtracker::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("tracker: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return handleNamedArgs(namedArgs);
	}

	BLCMDtracker *t = new BLCMDtracker(defaultTracker);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->trackerZ == UNDEFINED)
		printError("tracker: trackerZ must be defined");
	if(t->reportZ == UNDEFINED)
		t->reportZ = t->trackerZ;
	if(t->reportZ < t->trackerZ)
		printError("tracker: reportZ must be >= trackerZ");

	// define the fit parameters
	t->fitParam.define("x",t->scaleX);
	t->fitParam.define("y",t->scaleX);
	t->fitParam.define("Xp",t->scaleXp);
	t->fitParam.define("Yp",t->scaleXp);
	t->fitParam.define("Ptot",t->scalePtot);
	t->fitParam.define("t",t->scaleT);

	t->rotationMatrix = *BLCoordinates::getCurrentRotation();
	G4ThreeVector tmp(0.0,0.0,t->trackerZ);
	BLCoordinates::getCurrentGlobal(tmp,t->beamPosition);

	list.push_back(t);

	t->print(t->name);

	BLManager::getObject()->registerTrackingAction(t);
	BLManager::getObject()->registerCallback(t,0);
	BLManager::getObject()->registerZStep(t->trackerZ,t,6);
	if(t->trackerZ != t->reportZ)
		BLManager::getObject()->registerZStep(t->reportZ,t,6);

	return retval;
}

void BLCMDtracker::defineNamedArgs()
{
	argDouble(trackerZ,"trackerZ","The Z position of the tracker (Centerline, mm).");
	argDouble(reportZ,"reportZ","The Z position at which fit tracks are reported (Centerline, mm, default=trackerZ).");
	argDouble(scaleX,"scaleX","Scale for X and Y (mm); default=1mm.");
	argDouble(scaleXp,"scaleXp","Scale for dxdz and dydz (radians), default=0.001.");
	argDouble(scalePtot,"scalePtot","Scale for Ptot (MeV), default=0. Set to 0.0 if tracker cannot determine momentum.",MeV);
	argDouble(scaleT,"scaleT","Scale for t (ns), default=0. Set to 0.0 if tracker cannot determine time.",ns);
	argDouble(minPz,"minPz","Minimum Pz for valid tracks (default=10 MeV/c)",MeV);
	argInt(minHits,"minHits","Minimum number of hits for fitting (# planes).");
	argDouble(tolerance,"tolerance","Track fitting tolerance (mm) (default=0.01 mm)");
	argInt(maxIter,"maxIter","Maximum iterations during fitting (default=200).");
	argInt(verbose,"verbose","0=none, 1=result, 2=iterations, 3=detail, default=0.");
	argString(format,"format","Format of output NTuple.");
	argString(filename,"filename","Filename of output NTuple.");
	argInt(for009,"for009","Set nonzero to also output TrackerFit as FOR009.Dat.");
	argString(filename,"file","Synonym for filename.");
}

BLCMDtracker *BLCMDtracker::findTracker(G4String name)
{
	for(unsigned i=0; i<list.size(); ++i) {
		if(name == list[i]->getName()) return list[i];
	}
	return 0;
}

int BLCMDtracker::registerTrackerPlane(TrackerPlaneInstance *_plane, 
							const G4String name)
{
	for(unsigned i=0; i<plane.size(); ++i) {
		if(name == plane[i]->getName()) {
			char tmp[256];
			sprintf(tmp,"trackerplane %s has duplicate name",
				name.c_str());
			G4Exception("tracker","Duplicate trackerplane-s",
				FatalException, tmp);
		}
	}

	plane.push_back(_plane);
	hitWire.push_back(NO_HIT);
	hitTime.push_back(0.0);
	trackFields_hit += ":";
	trackFields_hit += name + "_w:" + name + "_t";
	return plane.size()-1;
}

void BLCMDtracker::planeHit(int planeID, int wire, double time)
{
	BLAssert(mode==BLTRACKER_TRUE);
	if(verbose >= 3)
		printf("Tracker %s Plane %d: hit wire %d time %.4f\n",
						name.c_str(),planeID,wire,time);
	BLAssert(planeID>=0 && planeID<(int)hitWire.size());
	hitWire[planeID] = wire;
	hitTime[planeID] = time;
	++trueHits;
}


bool BLCMDtracker::delta(int planeID, double deltaWire, double deltaTime,
					bool validWire, bool validTime)
{
	BLAssert(mode==BLTRACKER_FIT);
	if(verbose >= 3)
		printf("Tracker %s Plane %d: deltaWire %.3f deltaTime %.4f\n",
				name.c_str(),planeID,deltaWire,deltaTime);
	BLAssert(planeID>=0 && planeID<(int)hitWire.size());
	if(validWire) {
		chisq += deltaWire*deltaWire;
		++ndf;
	}
	if(validTime) {
		chisq += deltaTime*deltaTime;
		++ndf;
	}

	if(validWire || validTime) ++fitHits;

	if(planeID == (int)hitWire.size()-1) {
		// kill after last plane, unless reportOnly
		if(verbose >= 2 && !reportOnly)
			printf("TrackerPlane %d: last plane, kill track\n",
								planeID);
		return reportOnly;
	}
	return true;
}

G4Track *BLCMDtracker::constructTrack(const std::vector<double> &x)
{
	G4ThreeVector pos(x[0],x[1],0.0);
	G4double dxdz = x[2];
	G4double dydz = x[3];
	G4double Ptot = x[4];
	G4double time = x[5];
	G4double mass = trueDefinition->GetPDGMass();
	G4double ke = sqrt(Ptot*Ptot+mass*mass) - mass;
	G4double norm = sqrt(1.0+dxdz*dxdz+dydz*dydz);
	G4ThreeVector dir(dxdz/norm,dydz/norm,1.0/norm);

	// convert to centerline coordinates
	dir = rotationMatrix * dir;
	pos = rotationMatrix * pos + beamPosition;

	G4DynamicParticle *particle =
				new G4DynamicParticle(trueDefinition,dir,ke);
	G4Track *track = new G4Track(particle,time,pos);
	track->SetTrackID(trueTrackID+10000);
	track->SetParentID(trueParentID);
	track->SetGlobalTime(time);

	return track;
}

double BLCMDtracker::operator()(const std::vector<double> &x)
{
	if(chisqFailure) return 0.0; // short-circuit iterations in BLMinimize

	mode = BLTRACKER_FIT;

	G4Track *track = constructTrack(x);

	BLManager::getObject()->setPrimaryTrackID(trueTrackID+10000,
								trueParentID);
	BLRunManager::getObject()->processOneTrack(track);

	delete track;

	if(verbose >= 3)
		printf("%s tracking chisq=%.4f ndf=%d  fitHits=%d trueHits=%d\n",
					name.c_str(),chisq,ndf,fitHits,trueHits);

	if(fitHits < minHits || fitHits < trueHits) {
		chisqFailure = true;
		if(verbose >= 3)
		    printf("%s too few hits -- minimization short-circuited\n",
		    			name.c_str());
		return 0.0; // will short-circuit iteration in BLMinimize
	}

	return (ndf>0 ? chisq/ndf : chisq);
}

void BLCMDtracker::fitTrack()
{
	mode = BLTRACKER_FIT;
	if(verbose >= 3)
		printf("Tracker %s fitTrack() entered\n",name.c_str());

	if(!runmgr) runmgr = BLRunManager::getObject();
	if(!mgr) mgr = BLManager::getObject();

	if(mgr->getState() != BEAM) return;
	if(trueHits < minHits) {
		if(verbose >= 3)
			printf("Tracker %s fitTrack FAILED trueHits=%d < %d\n",
				name.c_str(),trueHits,minHits);
		return;
	}

	mgr->getPhysics()->setDoStochastics(FORCE_OFF,0);
	
	G4ThreeVector trueReportPosition = reportPosition;
	G4ThreeVector trueReportMomentum = reportMomentum;
	G4double trueReportTime = reportTime;

	BLMinimize min((verbose>=2 ? 2 : 0));
	min.setFunc(*this,fitParam.scale,fitParam.getNameList());
	min.setChisqMin(0.05);

	// this code must match the parameter definitions in command()
	// initial values come from the true track.
	std::vector<double> x;
	x.push_back(truePosition[0]);
	x.push_back(truePosition[1]);
	x.push_back(trueMomentum[0]/trueMomentum[2]);
	x.push_back(trueMomentum[1]/trueMomentum[2]);
	x.push_back(trueMomentum.mag());
	x.push_back(trueTime);

	chisqFailure = false;

	int status = min.minimize(x,(scaleX>0.0? tolerance/scaleX : tolerance),
								maxIter);

	if(chisqFailure)
		status = GSL_EFAILED;

	// try again if large chisq
	if(status==0 && min.value()>3.0 && min.getIterations()<maxIter/2) {
		if(verbose>0) printf("Tracker %s fitTrack Try Again\n",
								name.c_str());
		status = min.minimize(x,
		 	(scaleX>0.0 ? tolerance/scaleX : tolerance),maxIter);
		if(chisqFailure)
			status = GSL_EFAILED;
	}

	if(verbose > 0 && status == 0) {
		printf("%s fitTrack chisq=%.3f size=%.3f iter=%d\n",
			name.c_str(),min.value(),min.getSize(),
			min.getIterations());
		printf("   true x,y,dxdz,dydz,Ptot,t: %7.1f %7.1f %7.4f %.4f %7.1f %9.3f\n",
			truePosition[0],truePosition[1], trueMomentum[0]/trueMomentum[2],
			trueMomentum[1]/trueMomentum[2],trueMomentum.mag(),
			trueTime);
		printf("   fit  x,y,dxdz,dydz,Ptot,t: %7.1f %7.1f %7.4f %.4f %7.1f %9.3f\n",
			x[0],x[1],x[2],x[3],x[4],x[5]);
	} else if(verbose > 0) {
		printf("%s fitTrack FAILED\n",name.c_str());
	}

	if(status == 0 && ntuple_fit != 0) {
		// update the summary -- code must match param defs in command()
		fitParam.setDelta("x",x[0]-truePosition[0]);
		fitParam.setDelta("y",x[1]-truePosition[1]);
		fitParam.setDelta("Xp",x[2]-trueMomentum[0]/trueMomentum[2]);
		fitParam.setDelta("Yp",x[3]-trueMomentum[1]/trueMomentum[2]);
		fitParam.setDelta("Ptot",x[4]-trueMomentum.mag());
		fitParam.setDelta("t",x[5]-trueTime);

//printf("fitTrack() is tracking the fit track\n");
//verbose=3;
//mgr->setSteppingVerbose(3);
		// need to re-track the fitted track to z=reportZ
		reportOnly = true;
		G4Track *track = constructTrack(x);
		trueTrackID += 10000;
		BLManager::getObject()->setPrimaryTrackID(trueTrackID,
								trueParentID);
		BLRunManager::getObject()->processOneTrack(track);
		trueTrackID -= 10000;
		delete track;
		reportOnly = false;
//mgr->setSteppingVerbose(0);
//verbose=0;

		if(reportTime != UNDEFINED) {
//printf("fitTrack() reported OK\n");
			// this code must agree with the field names in 
			// trackFields_fit and FOR009_FIELDS
			double data[N_FIT_FIELDS];
#if N_FIT_FIELDS!=23
#error N_FIT_FIELDS!=23
#endif

			data[0] = reportPosition[0]/mm;
			data[1] = reportPosition[1]/mm;
			data[2] = reportPosition[2]/mm;
			data[3] = reportMomentum[0]/MeV;
			data[4] = reportMomentum[1]/MeV;
			data[5] = reportMomentum[2]/MeV;
			data[6] = reportTime/ns;
			data[7] = trueDefinition->GetPDGEncoding();
			data[8] = runmgr->GetCurrentEvent()->GetEventID();
			data[9] = trueTrackID+10000;
			data[10] = trueParentID;
			data[11] = 1.0;
			data[12] = min.value();
			data[13] = ndf;
			data[14] = fitHits;
			data[15] = min.getIterations();
			data[16] = trueReportPosition[0]/mm;
			data[17] = trueReportPosition[1]/mm;
			data[18] = trueReportPosition[2]/mm;
			data[19] = trueReportMomentum[0]/MeV;
			data[20] = trueReportMomentum[1]/MeV;
			data[21] = trueReportMomentum[2]/MeV;
			data[22] = trueReportTime/ns;
			ntuple_fit->appendRow(data,N_FIT_FIELDS);
			if(for009_fit != 0) {
#if N_FOR009_FIELDS!=18
#error N_FOR009_FIELDS!=18
#endif
				G4double point[4], field[6];
				point[0] = reportPosition[0];
				point[1] = reportPosition[1];
				point[2] = reportPosition[2];
				point[3] = reportTime;
				BLGlobalField::getObject()->GetFieldValue(point,
									field);
				data[12] = field[0]/tesla;		// Bx
				data[13] = field[1]/tesla;		// By
				data[14] = field[2]/tesla;		// Bz
				data[15] = field[3]/(megavolt/meter);	// Ex
				data[16] = field[4]/(megavolt/meter);	// Ey
				data[17] = field[5]/(megavolt/meter);	// Ez
				for009_fit->appendRow(data,N_FOR009_FIELDS);
			}
		} else {
//printf("fitTrack() -- track has reportTime==UNDEFINED\n");
		}
	}

	mgr->getPhysics()->setDoStochastics(NORMAL,0);

	mode = BLTRACKER_IGNORE;
	trueHits = -9999;

	if(verbose >= 3)
		printf("%s fitTrack() returns; status=%d\n",name.c_str(),status);
}

void BLCMDtracker::UserZSteppingAction(const G4Track *track)
{
	// get Centerline position and momentum
	BLCoordinates *coord = (BLCoordinates *)track->GetUserInformation();
	if(coord && !coord->isValid()) coord = 0;
	G4ThreeVector pos;
	G4ThreeVector momentum = track->GetMomentum();
	if(coord) {
		coord->getCoords(BLCOORD_CENTERLINE,pos);
		momentum = coord->getRotation() * momentum;
	}

	if(BLManager::getObject()->getState() == REFERENCE) {
		if(fabs(pos[2]-trackerZ) > 0.010*mm &&
		   fabs(pos[2]-reportZ) > 0.010*mm) {
			if(verbose >= 3) printf("tracker %s ZStep ignored\n",
							name.c_str());
		   	return;
		}
	} else if(reportOnly) {
		if(fabs(pos[2]-reportZ) > 0.010*mm) {
			if(verbose >= 3) printf("tracker %s ZStep ignored\n",
							name.c_str());
		   	return;
		}
	} else if(mode == BLTRACKER_TRUE) {
		if(fabs(pos[2]-trackerZ) > 0.010*mm &&
		   fabs(pos[2]-reportZ) > 0.010*mm) {
			if(verbose >= 3) printf("tracker %s ZStep ignored\n",
							name.c_str());
		   	return;
		}
	} else {
		if(verbose >= 3) printf("tracker %s ZStep ignored\n",
						name.c_str());
		return;
	}

	// here if Reference particle at either trackerZ or reportZ,
	// reportOnly at reportZ, or true particle at either trackerZ or reportZ

	if(momentum[2] < minPz) {
		if(verbose >= 3)
			printf("tracker %s setTrueTrack() -- Pz < %.3f\n",
					name.c_str(),minPz);
		return;
	}

	if(fabs(pos[2]-trackerZ) < 0.010*mm) {
		if(verbose >= 3)
			printf("tracker %s setTrueTrack() z=%.3f minHits=%d\n",
						name.c_str(),pos[2],minHits);
		truePosition = pos;
		trueMomentum = momentum;
		trueTime = track->GetGlobalTime();
		trueTrackID = BLManager::getObject()->getExternalTrackID(track);
		trueParentID = BLManager::getObject()->getExternalParentID(track);
		trueDefinition = track->GetDefinition();
		trueHits = 0;
	}
	if(fabs(pos[2]-reportZ) < 0.010*mm) {
		if(trueTrackID != BLManager::getObject()->getExternalTrackID(track) ||
		   trueParentID != BLManager::getObject()->getExternalParentID(track) ||
		   trueDefinition != track->GetDefinition()) {
			if(verbose >= 3) 
				printf("tracker %s ZStep wrong track: "
				"%d,%d %d,%d %08lX,%08lX\n",
				name.c_str(),trueTrackID,track->GetTrackID(),
				trueParentID,track->GetParentID(),
				(long)trueDefinition,(long)track->GetDefinition());
			return;
		}
		if(verbose >= 3)
			printf("tracker %s setReportTrack() z=%.3f\n",
							name.c_str(),pos[2]);
		reportPosition = pos;
		reportMomentum = momentum;
		reportTime = track->GetGlobalTime();
		if(reportOnly) {
			if(verbose >= 3)
				printf("tracker %s reportOnly killing track\n",
							name.c_str());
			((G4Track *)track)->SetTrackStatus(fStopAndKill);
		}
	}

	if(BLManager::getObject()->getState() != REFERENCE) return;
	if(for009_fit == 0 || reportTime == UNDEFINED) return;

	// here to append the reference to for009_fit

	// this code must agree with the field names in for009Track_fit
	G4double point[4], field[6];
	point[0] = reportPosition[0];
	point[1] = reportPosition[1];
	point[2] = reportPosition[2];
	point[3] = reportTime;
	BLGlobalField::getObject()->GetFieldValue(point,field);
	double data[18];
	data[0] = reportPosition[0]/mm;
	data[1] = reportPosition[1]/mm;
	data[2] = reportPosition[2]/mm;
	data[3] = reportMomentum[0]/MeV;
	data[4] = reportMomentum[1]/MeV;
	data[5] = reportMomentum[2]/MeV;
	data[6] = reportTime/ns;
	data[7] = track->GetDefinition()->GetPDGEncoding();
	data[8] = -1;		// EventID
	data[9] = 1;		// TrackID
	data[10] = 0;		// ParentID
	data[11] = 1.0;		// Weight
	data[12] = field[0]/tesla;		// Bx
	data[13] = field[1]/tesla;		// By
	data[14] = field[2]/tesla;		// Bz
	data[15] = field[3]/(megavolt/meter);	// Ex
	data[16] = field[4]/(megavolt/meter);	// Ey
	data[17] = field[5]/(megavolt/meter);	// Ez
	for009_fit->appendRow(data,18);
}

void BLCMDtracker::PreUserTrackingAction(const G4Track *track)
{
	BLManagerState state = BLManager::getObject()->getState();
	if(state != BEAM) return;

	if(BLRunManager::getObject()->getCollectiveMode())
		G4Exception("tracker command","Collective mode error",
				FatalException,
				"Does not work in collective tracking mode");

	for(unsigned i=0; i<plane.size(); ++i) {
		switch(mode) {
		case BLTRACKER_TRUE:
			plane[i]->setTrueMode();
			hitWire[i] = NO_HIT;
			hitTime[i] = 0.0;
			break;
		case BLTRACKER_FIT:
			plane[i]->setFitMode(hitWire[i],hitTime[i]);
			break;
		case BLTRACKER_IGNORE:
			plane[i]->setIgnoreMode();
			break;
		}
	}

	// reportTime is used as a flag to indicate that the report variables
	// are set for this track
	reportTime = UNDEFINED;

	switch(mode) {
	case BLTRACKER_TRUE:
		trueHits = 0;
		break;
	case BLTRACKER_FIT:
		fitHits = 0;
		chisq = 0.0;
		ndf = 0;
		break;
	case BLTRACKER_IGNORE:
		break;
	}
}

void BLCMDtracker::PostUserTrackingAction(const G4Track *track)
{
	BLManagerState state = BLManager::getObject()->getState();
	if(state != BEAM) return;

	if(mode == BLTRACKER_TRUE && reportTime != UNDEFINED &&
	   ntuple_hit != 0 && trueHits >= minHits) {
		if(!runmgr) runmgr = BLRunManager::getObject();
		static double *data=0;
		static unsigned ndata=0;
		if(ndata < N_HIT_FIELDS+2*hitWire.size()) {
			if(data) free(data);
			ndata = N_HIT_FIELDS+2*hitWire.size();
			data = new double[ndata];
		}
		// this code must agree with the field names in trackFields_hit
		data[0] = truePosition[0]/mm;
		data[1] = truePosition[1]/mm;
		data[2] = truePosition[2]/mm;
		data[3] = trueMomentum[0]/MeV;
		data[4] = trueMomentum[1]/MeV;
		data[5] = trueMomentum[2]/MeV;
		data[6] = trueTime/ns;
		data[7] = trueDefinition->GetPDGEncoding();
		data[8] = runmgr->GetCurrentEvent()->GetEventID();
		data[9] = trueTrackID;
		data[10] = trueParentID;
		data[11] = 1.0;
		data[12] = reportPosition[0]/mm;
		data[13] = reportPosition[1]/mm;
		data[14] = reportPosition[2]/mm;
		data[15] = reportMomentum[0]/MeV;
		data[16] = reportMomentum[1]/MeV;
		data[17] = reportMomentum[2]/MeV;
		data[18] = reportTime/ns;
		int j=N_HIT_FIELDS;
		for(unsigned i=0; i<hitWire.size(); ++i) {
			data[j++] = hitWire[i];
			data[j++] = hitTime[i];
		}
		ntuple_hit->appendRow(data,N_HIT_FIELDS+2*hitWire.size());
	}

	if(mode == BLTRACKER_FIT) {
		// update ndf for fitted parameters
		for(unsigned i=0; i<fitParam.scale.size(); ++i) {
			if(fitParam.scale[i] > 0.0)
				--ndf;
		}
		if(ndf < 0) ndf = 0;
	}

	for(unsigned i=0; i<plane.size(); ++i) {
		plane[i]->setIgnoreMode();
	}
}

void BLCMDtracker::callback(int type)
{
	// had to wait until all trackerplane-s are registered
	// and mode is known
	G4String m = trackermodeInstance.getMode();
	if(m == "true" || m == "both") {
	    ntuple_hit = BLNTuple::create(format,"TrackerHits",name,
	    					trackFields_hit,filename);
	} else if(m == "fit" || m == "both") {
	    ntuple_fit = BLNTuple::create(format,"TrackerFit",name,
	    					trackFields_fit,filename);
	   if(for009) for009_fit = BLNTuple::create("for009","TrackerFit",name,
	    					FOR009_FIELDS,filename);
	} else {
		G4Exception("tracker","Invalid mode",FatalException,m.c_str());
	}

	if(minHits < 0) minHits = hitWire.size();
}

void BLCMDtracker::handlePreviousTracks(const G4String file)
{
	printf("==================== Begin Fitting Tracker %s ======================\n",name.c_str());

	BLNTuple *in = BLNTuple::read(format, "TrackerHits", name,
							trackFields_hit, file);
	if(!in) {
		printError("tracker::handlePreviousTracks cannot open NTuple");
		return;
	}
	if(in->getNData() != N_HIT_FIELDS+2*(int)hitWire.size()) {
		in->close();
		printError("tracker::handlePreviousTracks NTuple wrong # fields");
		return;
	}

	if(!runmgr) runmgr = BLRunManager::getObject();
	if(!mgr) mgr = BLManager::getObject();
	runmgr->beginRun();
	mgr->setState(BEAM);

	// this code must agree with the field names in trackFields_hit
	static double *data=0;
	static unsigned ndata=0;
	if(ndata < N_HIT_FIELDS+2*hitWire.size()) {
		if(data) free(data);
		ndata = N_HIT_FIELDS+2*hitWire.size();
		data = new double[ndata];
	}
	while(in->readRow(data,ndata)) {
		int eventID = (int)data[8];
		if(eventID < 0) continue; // skip fitting reference particle
		truePosition[0] = data[0]*mm;
		truePosition[1] = data[1]*mm;
		truePosition[2] = data[2]*mm;
		trueMomentum[0] = data[3]*MeV;
		trueMomentum[1] = data[4]*MeV;
		trueMomentum[2] = data[5]*MeV;
		trueTime = data[6]*ns;
		trueDefinition = G4ParticleTable::GetParticleTable()->FindParticle((int)data[7]);
		trueTrackID = (int)data[9];
		trueParentID = (int)data[10];
		reportPosition[0] = data[12]*mm;
		reportPosition[1] = data[13]*mm;
		reportPosition[2] = data[14]*mm;
		reportMomentum[0] = data[15]*MeV;
		reportMomentum[1] = data[16]*MeV;
		reportMomentum[2] = data[17]*MeV;
		reportTime = data[18]*ns;
		trueHits = 0;
		int j=N_HIT_FIELDS;
		for(unsigned i=0; i<hitWire.size(); ++i) {
			hitWire[i] = (int)data[j++];
			hitTime[i] = data[j++];
			if(hitWire[i] != NO_HIT) ++trueHits;
		}
		if(trueHits < minHits) continue;
		mgr->setEventID(eventID);
		runmgr->beginEvent(eventID);
		fitTrack();
		runmgr->endEvent();
	}

	runmgr->endRun();
	in->close();

	fitParam.printSummary();
}


BLCMDtrackerplane::BLCMDtrackerplane() : BLElement()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Construct a tracker plane.");
	setDescription("A trackerplane belongs to a specific tracker, and "
		"represents one measuring element of the tracker. "
		"While the term 'wire' is used, a trackerplane can model "
		"any planar measuring device that measures one dimension "
		"using equally spaced detectors that are either on or off "
		"(hit or not hit). A trackerplane can be "
		"circular (specify radius and possibly innerRadius) or "
		"rectangular (specify height and width).\n\n"
		"The wires are at angle theta from the vertical, so theta=0 "
		"means vertical wires that measure x; theta=90 means "
		"horizontal wires that measure y; theta=180 also "
		"measures x, but increasing x means decreasing wire #.\n\n"
		"Setting wireSpacing=0 means there are no wires, which is "
		"useful for a trackerplane that models a scintillator used "
		"for track timing (set sigmaTime >= 0 to indicate that).\n\n"
		"sigmaTime>0 means this plane can measure the time of the "
		"track with that resolution. A Gaussian random number is added "
		"to the true track's time at this plane when reading the "
		"TrackHit NTuple.\n\n"
		"Survey errors can be modeled using the err-arguments -- "
		"the values they specify are applied during trackfitting "
		"(but not for true tracks). "
		"errType: 'fixed' means error values given are the actual values, "
		"'gaussian' means error values are the sigma of a Gaussian random "
		"number, 'rect' means error values are the half-width of a uniform "
		"random number. The random number is picked before the run "
		"begins; the random-number seed is set from the clock so "
		"every run will have different random errors.\n\n"
		"trackerplane has 3 modes:\n"
		"  true    The trackerplane reports the hit wire and time to\n"
		"          the tracker;\n"
		"  fit     The trackerplane reports the Chisq contribution of\n"
		"          the fit track distance to the true track's hit wire\n"
		"          center, plus survey errors (if any). The track time\n"
		"          also contributes to the chisq.\n"
		"  ignore  Any track is ignored.\n\n"
		"NOTE: the trackerplane-s of a tracker MUST be placed in the "
		"order that particles will hit them; the code does not sort "
		"them. Usually this means that each place command of a "
		"trackerplane must have a larger z value than the previous "
		"place command. They must also come after the trackerZ value "
		"of the tracker command."
	);
	
	// initial field values
	tracker = "";
	radius = 0.0;
	innerRadius = 0.0;
	height = 0.0;
	width = 0.0;
	length = 1.0*mm;
	material = "";
	color = "1,1,1";
	solid = 0;
	maxStep = -99.0;
	theta = 0.0;
	wireSpacing = -1.0;
	wireOffset = 0.0;
	errType = "fixed";
	errTheta = 0.0;
	errSpacing = 0.0;
	errOffset = 0.0;
	sigmaTime = -1.0;
}

BLCMDtrackerplane::~BLCMDtrackerplane()
{
	if(solid) delete solid;
}

BLCMDtrackerplane::BLCMDtrackerplane(const BLCMDtrackerplane& r) : BLElement(r)
{
	tracker = r.tracker;
	radius = r.radius;
	innerRadius = r.innerRadius;
	height = r.height;
	width = r.width;
	length = r.length;
	material = r.material;
	color = r.color;
	solid = r.solid;
	maxStep = r.maxStep;
	theta = r.theta;
	wireSpacing = r.wireSpacing;
	wireOffset = r.wireOffset;
	errType = r.errType;
	errTheta = r.errTheta;
	errSpacing = r.errSpacing;
	errOffset = r.errOffset;
	sigmaTime = r.sigmaTime;
}

int BLCMDtrackerplane::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("trackerplane: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultTrackerPlane.handleNamedArgs(namedArgs);
	}

	BLCMDtrackerplane *t = new BLCMDtrackerplane(defaultTrackerPlane);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	if(BLCMDtracker::findTracker(t->tracker) == 0)
		printError("trackerplane: invalid tracker");

	// handle errType
	if(t->errType.at(0) == 'g' || t->errType.at(0) == 'G') {
		t->errTheta = t->errTheta*CLHEP::RandGauss::shoot();
		t->errSpacing = t->errSpacing*CLHEP::RandGauss::shoot();
		t->errOffset = t->errOffset*CLHEP::RandGauss::shoot();
	} else if(t->errType.at(0) == 'r' || t->errType.at(0) == 'R') {
		t->errTheta = t->errTheta - t->errTheta*2.0*G4UniformRand();
		t->errSpacing = t->errSpacing - t->errSpacing*2.0*G4UniformRand();
		t->errOffset = t->errOffset - t->errOffset*2.0*G4UniformRand();
	}

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	return retval;
}

void BLCMDtrackerplane::defineNamedArgs()
{
	argString(tracker,"tracker","The tracker to which this plane belongs; REQUIRED.");
	argDouble(radius,"radius","The radius of the circular tracker plane (mm).");
	argDouble(innerRadius,"innerRadius","The inner radius of the circular tracker plane (0 mm).");
	argDouble(height,"height","The height of the rectangular tracker plane (mm).");
	argDouble(width,"width","The width of the rectangular tracker plane (mm).");
	argDouble(length,"length","The length of the tracker plane (mm).");
	argDouble(theta,"theta","Wire angle in X-Y plane (deg). 0=>x, 90=>y...",deg);
	argDouble(wireSpacing,"wireSpacing","Wire spacing (mm).");
	argDouble(wireOffset,"wireOffset","Wire # 0 offset (mm).");
	argString(errType,"errType","Error type: fixed, gaussian, rect.");
	argDouble(errTheta,"errTheta","Error in wire angle (deg).",deg);
	argDouble(errSpacing,"errSpacing","Error in wire spacing (mm).");
	argDouble(errOffset,"errOffset","Error in wire offset (mm).");
	argDouble(sigmaTime,"sigmaTime","Sigma for timing plane (ns), default=-1.");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm).");
	argString(material,"material","The material of the tracker plane.");
	argString(color,"color","The color of the tracker plane (''=invisible).");
}

void BLCMDtrackerplane::argChanged()
{
	if(radius > 0.0) width = height = radius*2.0;
}

void BLCMDtrackerplane::construct(G4RotationMatrix *relativeRotation,
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
			printError("trackerplane::construct %s INVALID - no "
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

	TrackerPlaneInstance *tpi = new TrackerPlaneInstance(thisname,pv,this);

	BLManager::getObject()->registerBeamStep(pv,tpi);
	// I'm not sure if this needs reference step any more...
	BLManager::getObject()->registerReferenceParticleStep(pv,tpi);

	printf("BLCMDtrackerplane::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2]);
}

void BLCMDtrackerplane::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	if(radius > 0.0)
		generateTubs(npoints, innerRadius, radius, 0.0, 360.0*deg,
								length, v);
	else
		generateBox(npoints,width,height,length,v);
}

G4bool BLCMDtrackerplane::isOutside(G4ThreeVector &local, G4double tolerance)
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

TrackerPlaneInstance::TrackerPlaneInstance(G4String _name,G4VPhysicalVolume *pv,
				BLCMDtrackerplane *plane)

{
	name = _name;
	thisVol = pv;
	mode = BLTRACKER_TRUE;
	tracker = BLCMDtracker::findTracker(plane->tracker);
	if(!tracker) {
		BLCommand::printError("TrackerPlaneInstance: invalid tracker");
		planeID = -1;
	} else {
		planeID = tracker->registerTrackerPlane(this,name);
	}
	hitWire = NO_HIT;
	hitTime = 0.0;
	haveBeenHit = false;
	sinTheta = sin(plane->theta);
	cosTheta = cos(plane->theta);
	wireSpacing = plane->wireSpacing;
	wireOffset = plane->wireOffset;
	sigmaWire = plane->wireSpacing/sqrt(12.0);
	sigmaTime = plane->sigmaTime;
	sinThetaErr = sin(plane->theta+plane->errTheta);
	cosThetaErr = cos(plane->theta+plane->errTheta);
	wireSpacingErr = plane->wireSpacing+plane->errSpacing;
	wireOffsetErr = plane->wireOffset+plane->errOffset;

	setNotHit();
}

void TrackerPlaneInstance::UserSteppingAction(const G4Step *step)
{
	BLManagerState state = BLManager::getObject()->getState();
	if(state != BEAM || mode == BLTRACKER_IGNORE) return;

	if(tracker->verbose >= 3)
		printf("trackerplane %s UserSteppingAction\n", name.c_str());
	if(isHit()) {
		if(tracker->verbose >= 3)
			printf("trackerplane %s already hit\n", name.c_str());
		return;
	}

	// get basic physical-volume info
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

	G4Track *track = step->GetTrack();

	G4ThreeVector position = track->GetPosition();
	G4double time = track->GetGlobalTime();
	G4ThreeVector momentum = track->GetMomentum();

	// transform to centerline coordinates, if available
	BLCoordinates *c = (BLCoordinates *)track->GetUserInformation();
	if(c && c->isValid()) {
		c->getCoords(BLCOORD_CENTERLINE,position);
		momentum = c->getRotation() * momentum;
	}

	// Remember that isHit() is known to be false.
	haveBeenHit = true;

	if(mode == BLTRACKER_TRUE) {
		hitWire = 0; // fake wire hit if no wires
		if(wireSpacing > 0.0) {
			double x = cosTheta*position[0] + sinTheta*position[1];
			hitWire = (int)floor((x-wireOffset)/wireSpacing+0.5);
		}
		hitTime = 0.0; // fake time if no timing
		if(sigmaTime > 0.0)
			hitTime = time + CLHEP::RandGauss::shoot(0.0,sigmaTime);
		tracker->planeHit(planeID,hitWire,hitTime);
		if(tracker->verbose >= 3)
			printf("trackerplane %s hitWire=%d hitTIme=%.3f\n",
						name.c_str(),hitWire,hitTime);
	} else if(mode == BLTRACKER_FIT) {
		if(hitWire == NO_HIT) {
			if(tracker->verbose >= 3)
				printf("trackerplane %s no true hit\n",
								name.c_str());
			return;
		}
		double dw=0.0, dt=0.0;
		bool validWire=false, validTime=false;
		if(wireSpacing > 0.0) {
			double x = cosThetaErr*position[0] +
							sinThetaErr*position[1];
			dw = (x-wireSpacingErr*hitWire-wireOffsetErr)/sigmaWire;
			validWire = true;
		}
		if(sigmaTime > 0.0) {
			dt = (time-hitTime)/sigmaTime;
			validTime = true;
		}
		if(tracker->verbose >= 3)
			printf("trackerplane %s hit dw=%.3f dt=%.3f\n",
							name.c_str(),dw,dt);
		if(!tracker->delta(planeID,dw,dt,validWire,validTime)) {
			track->SetTrackStatus(fStopAndKill);
			if(tracker->verbose >= 3) 
				printf("trackerplane %s kills fitting track\n",
								name.c_str());
		}
	}
}

BLCMDtrackermode::BLCMDtrackermode() : BLCommand(), BLCallback()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("Sets mode for all trackers, manages track fitting.");
	setDescription("USAGE: trackermode mode [file=...]\n"
		" mode can be any of:\n"
		"    true     tracks true tracks (normal operation)\n"
		"    fit      fits tracks to previous 'true' output\n"
		"    both     does both true and fit at once\n"
		"'fit' requires the filename argument to be the output of "
		"a previous 'true' run (filename is ignored in other modes); "
		"each tracker processes all of its tracks in the file. "
		"'true' mode simply denotes the standard G4beamline "
		"operation, and the simulated tracks are taken to be the "
		"'true' tracks of the system; the response of the tracker(s) "
		"to these tracks is then simulated in 'fit' mode. 'both' "
		"tracks a 'true' event, and then a fit is peerformed in each "
		"tracker for which its first track hit all trackerplane-s.\n\n"
		"Note that in 'true' "
		"mode every track that hits all trackerplane-s is considered, "
		"but in 'both' mode only the first track of the event can be "
		"considered.\n\n"
		"In fit mode, the filename argument MUST be different from "
		"the parameter 'histoFile', because this run must not "
		"overwrite the Root file from the previous (true) run.\n\n"
		"One trackermode command controls the mode of all trackers. "
		"'true' mode is normal G4beamline operation, and is the same "
		"as if no trackermode command was present.\n\n"
		"Note that the geometry of the system must not change "
		"between a 'true' run and a 'fit' run. You can, however, "
		"make small variations in fields to explore how errors in "
		"setting them will affect the fitted tracks. The trackerplane "
		"command can simulate survey errors.\n\n"
		"NOTE: the trackermode command must preceed all tracker "
		"commands in the input file."
	);

	filename = "";
	mode = "true";
	registered = false;
}

int BLCMDtrackermode::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(registered) {
		printError("trackermode: Multiple commands not allowed");
		return -1;
	}

	if(BLCMDtracker::list.size() > 0) {
		// Ordering is required for the callbacks to work properly.
		printError("trackermode: must preceed all tracker commands.");
		return -1;
	}

	handleNamedArgs(namedArgs);

	if(argv.size() >= 1) mode = argv[0];

	if(mode == "true") {
		;
	} else if(mode == "fit") {
		if(filename == "") 
			printError("trackermode: mode 'fit' requires a file\n");
		if(filename == Param.getString("histoFile") ||
		   filename == (G4String)(Param.getString("histoFile")+".root"))
			// Immediate fatal exception, trying not to clobber
			// the Root file from the previous (true) run
			G4Exception("trackermode","Overwriting input file",
				FatalException, 
				"filename must be different from histoFile");
	} else if(mode == "both") {
		;
	} else {
		printError("trackermode: invalid mode '%s'\n",mode.c_str());
		mode = "true";
	}

	if(!registered) {
		if(mode == "true")
			BLManager::getObject()->registerCallback(this,0);
		else
			BLManager::getObject()->registerCallback(this,3);
		registered = true;
	}

	print(mode);

	return 0;
}

void BLCMDtrackermode::defineNamedArgs()
{
	argString(filename,"filename","Filename to read for fitting tracks.");
	argString(filename,"file","Synonym for filename.");
}

void BLCMDtrackermode::callback(int type)
{
	if(BLCMDtracker::list.size() == 0) {
		G4Exception("trackermode","No trackers",FatalException, "");
	}

	if(mode == "true") {
		BLAssert(type == 0);
		for(unsigned i=0; i<BLCMDtracker::list.size(); ++i) {
			BLCMDtracker::list[i]->setMode(BLTRACKER_TRUE);
		}
		return;
	} else if(mode == "fit") {
		BLAssert(type==3);
		if(filename == Param.getString("histoFile") ||
		   filename == (G4String)(Param.getString("histoFile")+".root"))
			// probably already clobbered the Root file (user
			// could change histoFile between command() and here)
			G4Exception("trackermode","Overwriting input file",
				FatalException, 
				"filename must be different from histoFile");
		for(unsigned i=0; i<BLCMDtracker::list.size(); ++i) {
			BLCMDtracker::list[i]->handlePreviousTracks(filename);
		}
		return;
	} else if(mode != "both") {
		printError("trackermode::callback: invalid mode '%s'",
			mode.c_str());
		return;
	}
	
	BLAssert(type==3);

	BLRunManager *runmgr = BLRunManager::getObject();
	G4EventManager *evmgr = runmgr->getEventManager();
	BLManager *mgr = BLManager::getObject();

	printf("================== Prepare Tracking Beam for Tracker Fit ==================\n");
	mgr->getPhysics()->setDoStochastics(NORMAL);
	runmgr->beginRun();

	printf("================== Begin Tracking Beam for Tracker Fit ===============\n");
	mgr->setState(BEAM);

	for(;;) {
		runmgr->beginEvent(0);
		G4Event *currentEvent = (G4Event *)runmgr->GetCurrentEvent();
		mgr->GeneratePrimaries(currentEvent);
		if(runmgr->getRunAborted()) break;
		for(unsigned i=0; i<BLCMDtracker::list.size(); ++i)
			BLCMDtracker::list[i]->setMode(BLTRACKER_TRUE);
		evmgr->ProcessOneEvent(currentEvent);
		if(runmgr->getRunAborted()) break;
		for(unsigned i=0; i<BLCMDtracker::list.size(); ++i) {
			BLCMDtracker::list[i]->setMode(BLTRACKER_FIT);
			BLCMDtracker::list[i]->fitTrack();
		}
		runmgr->endEvent();
		delete currentEvent;
		currentEvent = 0;
	} 

	runmgr->endRun();
}

#endif // G4BL_GSL
