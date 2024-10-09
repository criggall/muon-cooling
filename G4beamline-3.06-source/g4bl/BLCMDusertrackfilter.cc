//	BLCMDusertrackfilter.cc
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
#include "G4EventManager.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Material.hh"

#include "BLElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLCoordinates.hh"

#include "BLUserCode.hh"

/**	class BLCMDusertrackfilter - implements a track filter using user code
 *	Each placement of this class implements a filter on all tracks
 *	entering its physical volume. The filter is implemented in user
 *	code dynamically loaded at run time. The user code can kill or
 *	modify the input track, and can create any number of secondary tracks.
 **/
class BLCMDusertrackfilter : public BLElement, BLCallback,
						BLManager::SteppingAction {
	G4double radius;
	G4double innerRadius;
	G4double height;
	G4double width;
	G4double length;
	G4String material;
	G4String color;
	G4double maxStep;
	G4String filterName;
	G4String init;
	G4VPhysicalVolume *thisVol;
	int steppingVerbose;
	BLUserTrackFilter *filter;
	G4VSolid *solid;
public:
	/// Default constructor.
	BLCMDusertrackfilter();

	/// Destructor.
	virtual ~BLCMDusertrackfilter();

	/// clone()
	BLElement *clone() { return new BLCMDusertrackfilter(*this); }

	/// Copy constructor.
	BLCMDusertrackfilter(const BLCMDusertrackfilter& r);

	/// commandName() returns "usertrackfilter".
	G4String commandName() { return "usertrackfilter"; }

	/// command() implements the usertrackfilter command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() will construct the usertrackfilter.
	/// Used for normal placements of a track filter object.
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

	/// callback() from BLCallback.
	void callback(int type);

	/// UserSteppingAction() from BLManager::SteppingAction.
	void UserSteppingAction(const G4Step *step); 
};

BLCMDusertrackfilter defaultUserTrackFilter;

BLCMDusertrackfilter::BLCMDusertrackfilter() : 
			BLElement(), BLCallback(), BLManager::SteppingAction()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Construct a usertrackfilter that filters tracks via user code.");
	setDescription("...");
	
	// initial field values
	radius = 0.0;
	innerRadius = 0.0;
	height = 0.0;
	width = 0.0;
	length = 1.0*mm;
	material = "";
	color = "1,1,1";
	maxStep = -99.0;
	filterName = "";
	init = "";
	solid = 0;
	filter = 0;
	thisVol = 0;
	steppingVerbose = 0;
}

BLCMDusertrackfilter::~BLCMDusertrackfilter()
{
	if(solid) delete solid;
}

BLCMDusertrackfilter::BLCMDusertrackfilter(const BLCMDusertrackfilter& r) : BLElement(r), BLCallback(), BLManager::SteppingAction()
{
	radius = r.radius;
	innerRadius = r.innerRadius;
	height = r.height;
	width = r.width;
	length = r.length;
	material = r.material;
	color = r.color;
	maxStep = r.maxStep;
	filterName = r.filterName;
	init = r.init;
	solid = r.solid;
	filter = r.filter;
	thisVol = 0;
	steppingVerbose = 0;
}

int BLCMDusertrackfilter::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("usertrackfilter: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultUserTrackFilter.handleNamedArgs(namedArgs);
	}

	BLCMDusertrackfilter *t = 
			new BLCMDusertrackfilter(defaultUserTrackFilter);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	// find the BLUserTrackFilter
	G4String knownFilters;
	std::vector<BLUserCode*> list = BLManager::getObject()->
					getUserCodeInstances("usertrackfilter");
	for(unsigned i=0; i<list.size(); ++i) {
		if(t->filterName == list[i]->getName()) {
			t->filter = dynamic_cast<BLUserTrackFilter*>(list[i]);
			break;
		}
		knownFilters += list[i]->getName();
		knownFilters += " ";
	}
	if(!t->filter) {
		printError("usertrackfilter: cannot find UserTrackFilter '%s'",
			t->filterName.c_str());
		printError("       known Filters: %s\n",knownFilters.c_str());
	}

	// call the user's setup
	if(t->filter) t->filter->setup(t->init.c_str());

	// register callback for after tracking
	BLManager::getObject()->registerCallback(t,2);

	return retval;
}

void BLCMDusertrackfilter::defineNamedArgs()
{
	argDouble(radius,"radius","The radius of the circular element (mm).");
	argDouble(innerRadius,"innerRadius","The inner radius of the circular element (0 mm, solid).");
	argDouble(height,"height","The height of the rectangular element (mm).");
	argDouble(width,"width","The width of the rectangular element (mm).");
	argDouble(length,"length","The length of the element (mm).");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm).");
	argString(material,"material","The material of the element.");
	argString(color,"color","The color of the element (''=invisible).");
	argString(filterName,"filterName","Name of the UserTrackFilter.");
	argString(filterName,"filter","Synonym for filterName.");
	argString(init,"init","Initialization string passed to user setup().");
}

void BLCMDusertrackfilter::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	steppingVerbose = Param.getInt("steppingVerbose");

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
			printError("usertrackfilter::construct %s INVALID - no "
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
	thisVol = pv;

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

	BLManager::getObject()->registerSteppingAction(pv,this);

	printf("usertrackfilter::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2]);
}

void BLCMDusertrackfilter::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	if(radius > 0.0)
		generateTubs(npoints, innerRadius, radius, 0.0, 360.0*deg, length, v);
	else
		generateBox(npoints,width,height,length,v);
}

G4bool BLCMDusertrackfilter::isOutside(G4ThreeVector &local, G4double tolerance)
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

void BLCMDusertrackfilter::callback(int type)
{
	if(type == 2) {
		if(filter) filter->complete(init.c_str());
	}
}

void BLCMDusertrackfilter::UserSteppingAction(const G4Step *step)
{
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
	int trackID = track->GetTrackID();
	int eventID = BLManager::getObject()->getEventID();
	std::vector<G4Track*> secondaries;

	if(filter) filter->filter(track,eventID,secondaries,steppingVerbose);

	// (Geant4 handles weight differently here)
	postPoint->SetWeight(track->GetWeight());

	// prevent user from screwing up trackID
	track->SetTrackID(trackID);

	for(unsigned i=0; i<secondaries.size(); ++i) {
		G4Track *t = secondaries[i];
		t->SetTrackID(0);
		t->SetParentID(track->GetTrackID());
		G4EventManager::GetEventManager()->GetTrackingManager()->
			GimmeSecondaries()->push_back(t);
	}
}
