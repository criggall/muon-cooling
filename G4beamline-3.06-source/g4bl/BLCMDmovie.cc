//	BLCMDmovie.cc
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

#include <stdio.h>
#include <vector>

#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4Polymarker.hh"
#include "G4VVisManager.hh"

#include "BLAssert.hh"
#include "BLCommand.hh"
#include "BLNTuple.hh"
#include "BLParam.hh"
#include "BLCoordinates.hh"

const int N_POINTS_ON_SURFACE = 200;

const char OutlineFields[] =
 "type:Xcl:Ycl:Zcl:Xg:Yg:Zg:Height:Width:Length:XZAngle:YZangle:Red:Green:Blue";
const int NOutlineFields = 15;

const char ReferenceFields[] = "T:Zref:Xcl:Ycl:Zcl:Xg:Yg:Zg:Ptot";
const int NReferenceFields = 9;

const char TrackFields[] =
    "x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight";
const unsigned NTrackFields = 12;

const char ElementFields[] = "ID:Xmin:Xmax:Ymin:Ymax:Zmin:Zmax:R:G:B";
const unsigned NElementFields = 10;

/**	class BLCMDmovie - output NTuple-s suitable for a movie
 *
 **/
class BLCMDmovie : public BLCommand, public BLCallback, 
					public BLManager::SteppingAction,
					public BLManager::TrackingAction {
	G4String coordinates;
	BLCoordinateType coordinateType;
	BLNTuple *outline;	// NTuple giving outlines of beamline elements
	BLNTuple *reference;	// NTuple giving reference particle coords
	BLNTuple *trace;	// NTuple tracing all beam tracks
	BLNTuple *elements;	// NTuple containing all BLElements
	BLManager *manager;
	std::vector<double> timeV;	// for interpolating in time the values
	std::vector<double> ZrefV;	// of Zref and Ptot for the
	std::vector<double> PtotV;	// reference particle.
public:
	BLCMDmovie();

	G4String commandName() { return "movie"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	// interpolate t in the interpolation vectors, returning Zref and Ptot.
	void interpolate(double t, double &Zref, double &Ptot);

	// dumpGeometry() dumps the geometry to the Movie/Elements NTuple.
	void dumpGeometry(G4VPhysicalVolume *phys, 
		G4ThreeVector &parent_offset, G4RotationMatrix &parent_rot);
	void dumpGeometry(G4VPhysicalVolume *phys);

	/// callback() from BLCallback.
	void callback(int type);

	/// UserSteppingAction() from BLManager::SteppingAction.
	void UserSteppingAction(const G4Step *step);

	/// from BLManager::TrackingAction
	virtual void PreUserTrackingAction(const G4Track *track);
	virtual void PostUserTrackingAction(const G4Track *track);

	void referenceStep(const G4Track *track);
	void beamStep(const G4Track *track);
	void generateFakeReference();
};

BLCMDmovie defaultMovieCommand;

BLCMDmovie::BLCMDmovie() : BLCallback(), BLManager::SteppingAction(),
						timeV(), ZrefV(), PtotV()
{
	registerCommand(BLCMDTYPE_OTHER);
	setSynopsis("Generate movie NTuple.");
	setDescription("This command outputs a set of NTuples suitable for "
		"generating a movie.");

	// initialize class variables here
	coordinates = "Reference";
	outline = 0;
	reference = 0;
	trace = 0;
	elements = 0;
	manager = 0;
}

int BLCMDmovie::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	manager = BLManager::getObject();

	handleNamedArgs(namedArgs);

	print("");

	coordinateType = BLCoordinates::getCoordinateType(coordinates);

	// register reference coords before registering ourselves
	if(coordinateType == BLCOORD_REFERENCE)
		BLCoordinates::useReferenceCoordinates();

	manager->registerTrackingAction(this);
	manager->registerReferenceParticleStep(0,this);
	manager->registerBeamStep(0,this);
	manager->registerCallback(this,0);
	manager->registerCallback(this,1);
	manager->registerCallback(this,4);

	//outline = BLNTuple::create("root", "Movie", "Outline",
	//						OutlineFields, "");
	reference = BLNTuple::create("root", "Movie", "Reference", 
							ReferenceFields, "");
	trace = BLNTuple::create("root", "Movie", "Tracks", 
							TrackFields, "");
	elements = BLNTuple::create("root","Movie","Elements",ElementFields,"");

	if(coordinateType != BLCOORD_REFERENCE)
		generateFakeReference();

	print("");

	return 0;
}

void BLCMDmovie::defineNamedArgs()
{
	argString(coordinates,"coordinates","Coordinates: global, centerline, or reference (default=r).");
}	

void BLCMDmovie::interpolate(double t, double &Zref, double &Ptot)
{
	if(timeV.size() < 2) {
		Zref = 0.0;
		Ptot = 0.0;
		return;
	}

	// simple linear search
	unsigned i;
	for(i=1; i<timeV.size(); ++i) {
		if(t < timeV[i]) break;
	}
	if(i >= timeV.size()) i = timeV.size() - 1;
	// now i points to the second entry to use for linear interpolation
	// (extrapolates the first or last pair if t is too small or too big).
	BLAssert(i > 0 && i < timeV.size());
	double f = (t - timeV[i-1])/(timeV[i] - timeV[i-1]);
	Zref = (1.0-f)*ZrefV[i-1] + f*ZrefV[i];
	Ptot = (1.0-f)*PtotV[i-1] + f*PtotV[i];
}

#ifdef DISPLAY_VISIBLE_MARKERS
static void display_dot(G4ThreeVector *p)
{
	static G4VVisManager* pVM=0;
	static G4Polymarker *marker=0;
	if(!p) {
		if(!!pVM && !!marker) {
			marker->SetMarkerType(G4Polymarker::circles);
			marker->SetScreenSize(5);
			marker->SetFillStyle(G4VMarker::filled);
			G4VisAttributes va(G4Colour(1.,1.,1.));  // white
			marker->SetVisAttributes(&va);
			pVM->Draw(*marker);
		}
		return;
	} else if(!pVM) {
		pVM = G4VVisManager::GetConcreteInstance();
		if(!pVM) return;
		marker = new G4Polymarker();
	}
	marker->push_back(*p);
}
#endif //DISPLAY_VISIBLE_MARKERS

void BLCMDmovie::dumpGeometry(G4VPhysicalVolume *phys, 
		G4ThreeVector &parent_offset, G4RotationMatrix &parent_rot)
{
	static BLCoordinates coord;
	static int number=0;

	G4ThreeVector my_offset = parent_offset + 
				parent_rot * phys->GetObjectTranslation();
	G4RotationMatrix my_rot = 
			parent_rot * phys->GetObjectRotationValue();
	G4LogicalVolume *log = phys->GetLogicalVolume();
	const G4VisAttributes *va = log->GetVisAttributes();
	if(va->IsVisible()) {
		G4VSolid *solid = log->GetSolid();
		double Xmin=DBL_MAX, Xmax=-DBL_MAX, Ymin=DBL_MAX, Ymax=-DBL_MAX;
		double Zmin=DBL_MAX, Zmax=-DBL_MAX;
		// find the element's extent in the desired coordinates by
		// generating points on its surface -- simple but slow.
		for(int i=0; i<N_POINTS_ON_SURFACE; ++i) {
			G4ThreeVector p = solid->GetPointOnSurface();
			p = my_rot * p + my_offset;
#ifdef DISPLAY_VISIBLE_MARKERS
display_dot(&p);
#endif
			coord.setGlobal(p,0.0);
			coord.getCoords(coordinateType,p);
			if(Xmin > p[0]) Xmin = p[0];
			if(Xmax < p[0]) Xmax = p[0];
			if(Ymin > p[1]) Ymin = p[1];
			if(Ymax < p[1]) Ymax = p[1];
			if(Zmin > p[2]) Zmin = p[2];
			if(Zmax < p[2]) Zmax = p[2];
		}
		BLAssert(NElementFields == 10);
		double data[NElementFields];
		data[0] = ++number;
		data[1] = Xmin;
		data[2] = Xmax;
		data[3] = Ymin;
		data[4] = Ymax;
		data[5] = Zmin;
		data[6] = Zmax;
		data[7] = va->GetColor().GetRed();
		data[8] = va->GetColor().GetGreen();
		data[9] = va->GetColor().GetBlue();
		elements->appendRow(data,NElementFields);
	}
	int n = log->GetNoDaughters();
	for(int i=0; i<n; ++i) {
		G4VPhysicalVolume *p = log->GetDaughter(i);
		if(!p) continue;
		dumpGeometry(p,my_offset,my_rot);
	}
}

void BLCMDmovie::callback(int type)
{
	if(type == 0)  {	// pre-tune callback
	} else if(type == 1) {	// post-reference callback
		G4ThreeVector offset;
		G4RotationMatrix rot;
		printf("Generating Movie/Elements NTuple... "); fflush(stdout);
		dumpGeometry(manager->getWorldPhysicalVolume(),offset,rot);
		printf("Done.\n"); fflush(stdout);
	} else if(type == 4) {	// visualization
#ifdef DISPLAY_VISIBLE_MARKERS
		G4ThreeVector offset;
		G4RotationMatrix rot;
		printf("Generating Movie/Elements NTuple... "); fflush(stdout);
		dumpGeometry(manager->getWorldPhysicalVolume(),offset,rot);
		printf("Done.\n"); fflush(stdout);
		display_dot(0);
#endif
	}
}

void BLCMDmovie::UserSteppingAction(const G4Step *step)
{
	switch(manager->getState()) {
	case IDLE:
		break;
	case VISUAL:
		break;
	case TUNE:
		break;
	case REFERENCE:
		if(coordinateType == BLCOORD_REFERENCE)
			referenceStep(step->GetTrack());
		break;
	case BEAM:
		beamStep(step->GetTrack());
		break;
	case SPECIAL:
		break;
	case SOURCE:
		break;
	}
}

void BLCMDmovie::PreUserTrackingAction(const G4Track *track)
{
	if(manager->getState() == REFERENCE && 
	   coordinateType == BLCOORD_REFERENCE)
		referenceStep(track);
	//@ if(manager->getState() == BEAM) beamStep(track);
}

void BLCMDmovie::PostUserTrackingAction(const G4Track *track)
{
	if(manager->getState() == REFERENCE && 
	   coordinateType == BLCOORD_REFERENCE)
		referenceStep(track);
	if(manager->getState() == BEAM) beamStep(track);
}

void BLCMDmovie::referenceStep(const G4Track *track)
{
	G4ThreeVector posGlobal = track->GetPosition();
	G4ThreeVector posCL = posGlobal;
	G4double Zref = BLCoordinates::getMostRecentReferenceZ();
	G4double time = track->GetGlobalTime();
	G4ThreeVector momentum = track->GetMomentum();

	// transform to centerline coordinates
	BLCoordinates *c = (BLCoordinates *)track->GetUserInformation();
	if(c && c->isValid()) {
		c->getCoords(BLCOORD_CENTERLINE,posCL);
	}

	// enter this step into interpolation vectors.
	if(timeV.size() == 0 || time-timeV.back() > 1.0*picosecond) {
		timeV.push_back(time);
		ZrefV.push_back(Zref);
		PtotV.push_back(momentum.mag());
	}

	BLAssert(NReferenceFields == 9);
	double data[NReferenceFields];
	data[0] = time/ns;
	data[1] = Zref/mm;
	data[2] = posCL[0]/mm;
	data[3] = posCL[1]/mm;
	data[4] = posCL[2]/mm;
	data[5] = posGlobal[0]/mm;
	data[6] = posGlobal[1]/mm;
	data[7] = posGlobal[2]/mm;
	data[8] = momentum.mag()/MeV;

	reference->appendRow(data,NReferenceFields);
}

void BLCMDmovie::beamStep(const G4Track *track)
{
	G4RunManager* runmgr = G4RunManager::GetRunManager();
	const G4Event* event = runmgr->GetCurrentEvent();
	int evNum = event->GetEventID();
	G4ThreeVector position = track->GetPosition();
	G4double time = track->GetGlobalTime();
	G4ThreeVector momentum = track->GetMomentum();

	// transform to desired coordinates
	BLCoordinates *c = (BLCoordinates *)track->GetUserInformation();
	if(c && c->isValid()) {
		c->getCoords(coordinateType,position);
		momentum = c->getRotation() * momentum;
	}

	// subtract interpolated Zref and Pref for reference coords
	if(coordinateType == BLCOORD_REFERENCE) {
		double Zref, Pref;
		interpolate(time,Zref,Pref);
		position[2] -= Zref;
		momentum[2] -= Pref;
	}

	BLAssert(NTrackFields == 12);
	double data[NTrackFields];
	data[0] = position[0]/mm;		// x (mm)
	data[1] = position[1]/mm;		// y (mm)
	data[2] = position[2]/mm;		// z (mm)
	data[3] = momentum[0]/MeV;		// Px (MeV/c)
	data[4] = momentum[1]/MeV;		// Py (MeV/c)
	data[5] = momentum[2]/MeV;		// Pz (MeV/c)
	data[6] = time/ns;			// t (ns)
	data[7] = track->GetDefinition()->GetPDGEncoding();
	data[8] = evNum;				// Event ID
	data[9] = 
	    BLManager::getObject()->getExternalTrackID(track);
	data[10] = 
	    BLManager::getObject()->getExternalParentID(track);
	data[11] = track->GetWeight();		// Weight

	trace->appendRow(data,NTrackFields);
}

void BLCMDmovie::generateFakeReference()
{
	BLAssert(NReferenceFields == 9);
	double data[NReferenceFields];
	data[0] = -(DBL_MAX/4.0);
	data[1] = -(DBL_MAX/4.0);
	data[2] = 0.0;
	data[3] = 0.0;
	data[4] = 0.0;
	data[5] = 0.0;
	data[6] = 0.0;
	data[7] = 0.0;
	data[8] = 0.0;
	reference->appendRow(data,NReferenceFields);

	data[0] = (DBL_MAX/4.0);
	data[1] = (DBL_MAX/4.0);
	reference->appendRow(data,NReferenceFields);
}

