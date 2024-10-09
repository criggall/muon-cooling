//	BLCoordinates.cc
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

#include "BLAssert.hh"
#include "BLCoordinates.hh"
#include "BLCommand.hh"
#include "BLParam.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

class BLCoordinates::ReferenceCoordinates : public BLManager::TrackingAction, 
					public BLManager::SteppingAction {
	BLManager *manager;
	G4ThreeVector prevDir;
	G4double prevZ;
	G4ThreeVector prevPos;
	G4double stepZ;
	G4ThreeVector stepPos;
	G4ThreeVector prevStepPos;
	static G4ThreeVector unitZ;
public:
	ReferenceCoordinates();
	virtual ~ReferenceCoordinates() { }
	void PreUserTrackingAction(const G4Track *track);
	void PostUserTrackingAction(const G4Track *track);
	void UserSteppingAction(const G4Step *step);
	void addStep(const G4ThreeVector &pos, G4double z, bool final=false);
};


std::vector<BLCoordinates::CLSegment> BLCoordinates::segmentCLVector;
unsigned BLCoordinates::currentCL = 0;
G4bool BLCoordinates::ring = false;
BLCoordinates::State BLCoordinates::state = UNUSED;
G4bool BLCoordinates::steppingVerbose = false;
G4bool BLCoordinates::useRC = false;
G4bool BLCoordinates::validRC = false;
std::vector<BLCoordinates::RCSegment> BLCoordinates::segmentRCVector;
G4double BLCoordinates::tolerance = 1.0e-4;
G4ThreeVector BLCoordinates::ReferenceCoordinates::unitZ(0.0,0.0,1.0);
G4double BLCoordinates::mostRecentReferenceZ = 0;

G4Allocator<BLCoordinates> BLCoordinates::allocator;

static G4RotationMatrix IdentityMatrix;


BLCoordinates::BLCoordinates() : G4VUserTrackInformation()
{
	// no state change, yet
	verify=VERIFY; 
	segmentCL=99999999; 
	segmentRC=99999999; 
	referenceCLsegment = -1;
	referenceTransformUsed = 0;
	global[0]=global[1]=global[2]=global[3]=0.0; 
	centerline[0]=centerline[1]=centerline[2]=centerline[3]=0.0; 
	reference[0]=reference[1]=reference[2]=reference[3]=0.0; 
	segmentCLRadiusFailure = -1;
}

BLCoordinates::BLCoordinates(const G4double _global[4]) : 
			G4VUserTrackInformation()
{
	state = TRACKING;
	verify = VERIFY;
	findCLSegment(_global);
}

G4bool BLCoordinates::isValid() const 
{
	return this != 0 && verify == VERIFY; 
}

// find lowest segmentCL containing _global[]; set segmentCL and centerline[]
G4bool BLCoordinates::findCLSegment(const G4double _global[4])
{
	for(segmentCL=0; segmentCL<segmentCLVector.size(); ++segmentCL) {
		segmentCLVector[segmentCL].transform.getLocal(centerline,_global);
//if(steppingVerbose) printf("FindCLSegment: seg=%d minZ=%.3f maxZ=%.3f Zcl=%.3f Rcl=%.3f rCut=%.3f\n",segmentCL,segmentCLVector[segmentCL].minZ,segmentCLVector[segmentCL].maxZ,centerline[2],sqrt(centerline[0]*centerline[0]+centerline[1]*centerline[1]),segmentCLVector[segmentCL].radiusCut);
		// note no check against minZ -- segmentCLs are ordered
		if(centerline[2] <= segmentCLVector[segmentCL].maxZ &&
		   (segmentCLVector[segmentCL].radiusCut <= 0.0 ||
		    (int)segmentCL == segmentCLRadiusFailure ||
		    sqrt(centerline[0]*centerline[0]+centerline[1]*centerline[1]) <=
	   				segmentCLVector[segmentCL].radiusCut)) {
			return true;
		}
	}
	--segmentCL;
	return false;
}

G4bool BLCoordinates::setGlobal(const G4double _global[4])
{
	global[0] = _global[0];
	global[1] = _global[1];
	global[2] = _global[2];
	global[3] = _global[3];

	if(segmentCL >= segmentCLVector.size()) {
//printf("setGlobal-0: calling findCLSegment\n");
		G4bool ret = findCLSegment(_global);
		if(ret && validRC) updateReferenceCoordinates();
//if(!ret) printf("setGlobal-0: findCLSegment FAILED\n");
		return ret;
	}

	segmentCLVector[segmentCL].transform.getLocal(centerline,_global);
//printf("setGlobal-1A: global=%.1f,%.1f,%.1f  cl=%.1f,%.1f,%.1f\n",_global[0],_global[1],_global[2],centerline[0],centerline[1],centerline[2]);
//printf("setGlobal-1: seg=%d minZ=%.3f maxZ=%.3f Zcl=%.3f Rcl=%.3f rCut=%.3f\n",segmentCL,segmentCLVector[segmentCL].minZ,segmentCLVector[segmentCL].maxZ,centerline[2],sqrt(centerline[0]*centerline[0]+centerline[1]*centerline[1]),segmentCLVector[segmentCL].radiusCut);
	while(centerline[2] > segmentCLVector[segmentCL].maxZ ||
	   (segmentCLVector[segmentCL].radiusCut > 0.0 &&
	    (int)segmentCL != segmentCLRadiusFailure &&
	    sqrt(centerline[0]*centerline[0]+centerline[1]*centerline[1]) >
	   				segmentCLVector[segmentCL].radiusCut)) {
		if(++segmentCL >= segmentCLVector.size()) {
			if(ring && findCLSegment(_global)) {
				segmentCLVector[segmentCL].transform.getLocal(centerline,
									_global);
				if(validRC) updateReferenceCoordinates();
				return true;
			}
			--segmentCL;
//printf("setGlobal-1a: FAILED\n");
			return false;
		}
		segmentCLVector[segmentCL].transform.getLocal(centerline,_global);
//printf("setGlobal-2: seg=%d minZ=%.3f maxZ=%.3f Zcl=%.3f Rcl=%.3f rCut=%.3f\n",segmentCL,segmentCLVector[segmentCL].minZ,segmentCLVector[segmentCL].maxZ,centerline[2],sqrt(centerline[0]*centerline[0]+centerline[1]*centerline[1]),segmentCLVector[segmentCL].radiusCut);
/***
		// note no check against minZ -- segmentCLs are ordered
		if(centerline[2] > segmentCLVector[segmentCL].maxZ ||
		   (segmentCLVector[segmentCL].radiusCut > 0.0 &&
		    segmentCL != segmentCLRadiusFailure &&
		    sqrt(centerline[0]*centerline[0]+centerline[1]*centerline[1]) >
					segmentCLVector[segmentCL].radiusCut)) {
			if(ring && findCLSegment(_global)) {
				segmentCLVector[segmentCL].transform.getLocal(centerline,
									_global);
				if(validRC) updateReferenceCoordinates();
				return true;
			}
//printf("setGlobal-2a: FAILED\n");
			return false;
		}
***/
	}

	if(validRC) updateReferenceCoordinates();
	return true;
}

G4bool BLCoordinates::setGlobal(const G4ThreeVector &pos, G4double time) 
{
	global[0]=pos[0],global[1]=pos[1],global[2]=pos[2];
	global[3]=time;
	return setGlobal(global);
}

void BLCoordinates::update(G4Track *track)
{
	// ensure there is a BLCoordiantes linked into the track
	BLCoordinates *c = (BLCoordinates *)track->GetUserInformation();
	if(!c || !c->isValid()) {
		c = new BLCoordinates();
		track->SetUserInformation(c);
	}

	// now update the object and apply radiusCut
	if(!c->setGlobal(track->GetPosition(),track->GetGlobalTime())) {
		track->SetTrackStatus(fStopAndKill);
		if(steppingVerbose)
			printf("Track Centerline radius exceeds radiusCut\n");
		c->segmentCLRadiusFailure = c->segmentCL;
	}
}

void BLCoordinates::getCoords(BLCoordinateType type, G4double local[4])
{
//printf("AAA\n"); fflush(stdout);
	prevType = type;

	switch(type) {
	case BLCOORD_GLOBAL:
//printf("BBB\n"); fflush(stdout);
		local[0] = global[0]; local[1] = global[1];
		local[2] = global[2]; local[3] = global[3];
//printf("getCoords(GLOBAL): %.1f,%.1f,%.1f\n",local[0],local[1],local[2]);
		return;
	case BLCOORD_CENTERLINE:
//printf("CCC\n"); fflush(stdout);
		BLAssert(segmentCL < segmentCLVector.size());
		local[0] = centerline[0]; local[1] = centerline[1];
		local[2] = centerline[2]; local[3] = centerline[3];
//printf("getCoords(CL): %.1f,%.1f,%.1f\n",local[0],local[1],local[2]);
		return;
	case BLCOORD_REFERENCE:
//printf("DDD1 segmentRCVector.size=%d segmentRC=%d\n",segmentRCVector.size(),segmentRC); fflush(stdout);
		BLAssert(validRC);
		// interpolate
		if(reference[2] < segmentRCVector[segmentRC].interpolatePrev &&
							   segmentRC > 0) {
//printf("DDD2\n"); fflush(stdout);
			referenceTransformUsed = 
				&segmentRCVector[segmentRC-1].interpolate;
			referenceTransformUsed->getLocal(reference,global);
//printf("DDD3\n"); fflush(stdout);
		} else if(reference[2] > 
				segmentRCVector[segmentRC].interpolateThis) {
//printf("DDD4\n"); fflush(stdout);
			referenceTransformUsed = 
				&segmentRCVector[segmentRC].interpolate;
			referenceTransformUsed->getLocal(reference,global);
//printf("DDD5\n"); fflush(stdout);
		}
//printf("DDD6\n"); fflush(stdout);
{ G4ThreeVector axis = referenceTransformUsed->getRotation().getAxis();
G4double angle = referenceTransformUsed->getRotation().getDelta();
//printf("segmentRC=%d minZ=%.4g maxZ=%.4g global=%.1f,%.1f,%.1f ref=%.1f,%.1f,%.1f axis=%.1f,%.1f,%.1f angle(deg)=%.1f\n",segmentRC,segmentRCVector[segmentRC].minZ,segmentRCVector[segmentRC].maxZ,global[0],global[1],global[2],reference[0],reference[1],reference[2],axis[0],axis[1],axis[2],angle/deg);
}
//printf("DDD7\n"); fflush(stdout);
		local[0] = reference[0]; local[1] = reference[1];
		local[2] = reference[2]; local[3] = reference[3];
//printf("getCoords(REF): %.1f,%.1f,%.1f\n",local[0],local[1],local[2]);
		return;
	case BLCOORD_LOCAL:
		if(localTransform != 0) {
			localTransform->getLocal(local,global);
		} else {
			G4Exception("BLCoordinates","No local transform",
							FatalException, "");
		}
		return;
	}
	G4Exception("BLCoordinates","Invalid Coordinate Type",FatalException,
									"");
}

void BLCoordinates::getCoords(BLCoordinateType type, G4ThreeVector &local)
{
//printf("getCoords type=%d &local=%p\n",(int)type,&local); fflush(stdout);
	G4double tmp[4];
//printf("A\n"); fflush(stdout);
	getCoords(type,tmp);
//printf("B\n"); fflush(stdout);
	local[0] = tmp[0];
	local[1] = tmp[1];
	local[2] = tmp[2];
//printf("getCoords(): %.1f,%.1f,%.1f\n",local[0],local[1],local[2]); fflush(stdout);
}

G4RotationMatrix &BLCoordinates::getRotation(BLCoordinateType type)
{
	switch(type) {
	case BLCOORD_GLOBAL:
		return IdentityMatrix;
	case BLCOORD_CENTERLINE:
		BLAssert(segmentCL < segmentCLVector.size());
		return segmentCLVector[segmentCL].transform.getRotation();
	case BLCOORD_REFERENCE:
		BLAssert(validRC && referenceTransformUsed != 0);
		return referenceTransformUsed->getRotation();
	case BLCOORD_LOCAL:
		if(localTransform != 0)
			return localTransform->getRotation();
		G4Exception("BLCoordinates","No local transform",
							FatalException, "");
	}
	G4Exception("BLCoordinates","Invalid Coordinate Type",FatalException,
									"");
	return IdentityMatrix;
}

void BLCoordinates::init()
{
	if(state == TRACKING) {
		G4Exception("BLCoordinates","Invalid init during tracking",
							FatalException, "");
	}
	if(state == UNUSED) {
		state = CONSTRUCTING;
		segmentCLVector.push_back(CLSegment(-DBL_MAX, +DBL_MAX, 0.0,
						BLCoordinateTransform()));
		currentCL = 0;
		steppingVerbose = 
			BLManager::getObject()->getSteppingVerbose() != 0;
	}
}

void BLCoordinates::start(const G4double _global[4], G4double z,
			G4RotationMatrix *rotation, G4double radiusCut,
			G4bool _ring)
{
	if(state != UNUSED || currentCL != 0) {
		G4Exception("BLCoordinates","Invalid start",FatalException, "");
	}
	init();
	
	G4ThreeVector pos(0,0,-z);
	if(rotation)
		pos = *rotation * pos;
	pos[0] += _global[0];
	pos[1] += _global[1];
	pos[2] += _global[2];
	segmentCLVector[0].transform = (rotation ?
					BLCoordinateTransform(rotation,pos) : 
					BLCoordinateTransform(pos));
	segmentCLVector[0].radiusCut = radiusCut;
	ring = _ring;
	steppingVerbose = BLManager::getObject()->getSteppingVerbose() != 0;
}

void BLCoordinates::corner(G4double z, G4RotationMatrix *rotation,
							G4double radiusCut)
{
	init();
	
	G4double local[4], glob[4];
	local[0]=local[1]=local[3]=0.0;
	local[2] = z;
	segmentCLVector[currentCL].transform.getGlobal(local,glob);

	segmentCLVector[currentCL].maxZ = z;
	G4RotationMatrix *p=0;
	if(rotation)
		// if C=currentCL rot, and R=this rot, we want R referenced to
		// the CURRENT centerline coordinates, so we need C R C^-1.
		// But we compose this with C (on the right), and get C R.
		p = new G4RotationMatrix(*getCurrentRotation() * *rotation);
	else if(getCurrentTransform().isRotated())
		p = getCurrentRotation();

	G4ThreeVector pos(0,0,-z);
	if(p)
		pos = *p * pos;
	pos[0] += glob[0];
	pos[1] += glob[1];
	pos[2] += glob[2];
	segmentCLVector.push_back(CLSegment(z, +DBL_MAX, radiusCut,
					BLCoordinateTransform(p,pos)));
	currentCL = segmentCLVector.size() - 1;
}

BLCoordinateTransform &BLCoordinates::getCurrentTransform()
{
	init();
	return segmentCLVector[currentCL].transform;
}

G4RotationMatrix *BLCoordinates::getCurrentRotation()
{
	init();
	return &segmentCLVector[currentCL].transform.getRotationInverse();
}

void BLCoordinates::getCurrentGlobal(const G4double local[4], G4double _global[4])
{
	init();
	segmentCLVector[currentCL].transform.getGlobal(local,_global);
}

void BLCoordinates::getCurrentGlobal(const G4ThreeVector& local, G4ThreeVector& _global)
{
	init();
	segmentCLVector[currentCL].transform.getGlobal(local,_global);
}

void BLCoordinates::getGlobalAnywhere(const G4ThreeVector& local, 
						G4ThreeVector& _global)
{
	for(unsigned int seg=0; seg<segmentCLVector.size(); ++seg) {
		// note no check against minZ -- segmentCLs are ordered
		if(local[2] <= segmentCLVector[seg].maxZ &&
		   (segmentCLVector[seg].radiusCut <= 0.0 ||
		    sqrt(local[0]*local[0]+local[1]*local[1]) <=
	   				segmentCLVector[seg].radiusCut)) {
			segmentCLVector[seg].transform.getGlobal(local,_global);
			return;
		}
	}
}

G4RotationMatrix *BLCoordinates::getGlobalRotation(const G4ThreeVector& local)
{
	for(unsigned int seg=0; seg<segmentCLVector.size(); ++seg) {
		// note no check against minZ -- segmentCLs are ordered
		if(local[2] <= segmentCLVector[seg].maxZ &&
		   (segmentCLVector[seg].radiusCut <= 0.0 ||
		    sqrt(local[0]*local[0]+local[1]*local[1]) <=
	   				segmentCLVector[seg].radiusCut)) {
			return &segmentCLVector[seg].transform.getRotation();
		}
	}

	return 0;
}

G4double BLCoordinates::getRadiusCut(const G4ThreeVector& local)
{
	for(unsigned int seg=0; seg<segmentCLVector.size(); ++seg) {
		// note no check against minZ -- segmentCLs are ordered
		if(local[2] <= segmentCLVector[seg].maxZ &&
		   (segmentCLVector[seg].radiusCut <= 0.0 ||
		    sqrt(local[0]*local[0]+local[1]*local[1]) <=
	   				segmentCLVector[seg].radiusCut)) {
			return segmentCLVector[seg].radiusCut;
		}
	}

	return 0;
}

BLCoordinateType BLCoordinates::getCoordinateType(G4String s)
{
	char c = tolower(s[0u]);
	switch(c) {
	case 'g':	return BLCOORD_GLOBAL;
	case 'c':	return BLCOORD_CENTERLINE;
	case 'r':	useReferenceCoordinates();
			return BLCOORD_REFERENCE;
	case 'l':	return BLCOORD_LOCAL;
	default:	BLCommand::printError("Unknown coordinate type '%s',"
				" Centerline used.\n",s.c_str());
			return BLCOORD_CENTERLINE;
	}
}

void BLCoordinates::useReferenceCoordinates()
{

	if(!useRC) {
		useRC = true;
		new ReferenceCoordinates();
	}
}


void BLCoordinates::updateReferenceCoordinates()
{
	// if we are lost, use CL coord segment as a starting point
	if(segmentRC >= segmentRCVector.size()) {
		segmentRC = segmentCLVector[segmentCL].firstRCSegment;
		if(segmentRC > 0) --segmentRC; // try 1 previous
	}
	if(segmentRC >= segmentRCVector.size())
		segmentRC = 0;

	for( ; segmentRC<segmentRCVector.size(); ++segmentRC) {
		segmentRCVector[segmentRC].transform.getLocal(reference,global);
//printf("segmentRC=%d minZ=%.4g maxZ=%.4g global=%.1f,%.1f,%.1f ref=%.1f,%.1f,%.3f\n",segmentRC,segmentRCVector[segmentRC].minZ,segmentRCVector[segmentRC].maxZ,global[0],global[1],global[2],reference[0],reference[1],reference[2]);
		referenceTransformUsed = &segmentRCVector[segmentRC].transform;
		if(segmentRC != 0) // no correction (minZ = -DBL_MAX)
		   reference[2] += (segmentRCVector[segmentRC].correction - 1.0)
		   	       * (reference[2]-segmentRCVector[segmentRC].minZ);
//printf("         corrected reference[2]=%.1f\n",reference[2]);
		if(reference[2] > segmentRCVector[segmentRC].maxZ)
			continue;
		return;
	}
	G4Exception("BLCoordinates","Cannot determine reference coordinate segment",
						FatalException,"");
}

BLCoordinates::ReferenceCoordinates::ReferenceCoordinates() : 
				prevDir(), prevPos(), stepPos(), prevStepPos()
{
	manager = BLManager::getObject();
	stepZ = -DBL_MAX;
	prevZ = 0.0;

	manager->registerReferenceParticleStep(0,this);
	manager->registerTrackingAction(this);
}

void BLCoordinates::ReferenceCoordinates::PreUserTrackingAction(const G4Track *track)
{
	BLCoordinates *coord = (BLCoordinates *)track->GetUserInformation();
	if(!coord || !coord->isValid()) {
		G4Exception("BLCoordinates","Coordinate Object Missing in Reference Track",
							FatalException, "");
	}

	coord->segmentRC = 99999999;

	if(BLManager::getObject()->getState() == BEAM) {
		if(!validRC) {
			G4Exception("BLCoordinates","Reference Coordinates not available",
							FatalException, "");
		}
		return;
	}

	if(BLManager::getObject()->getState() != REFERENCE) return;

	// @@@ need to handle multiple reference particles !!!!

	prevPos = track->GetPosition();
	prevDir = track->GetMomentumDirection();

	// initial Z is same as Centerline Coord value.
	G4ThreeVector cl;
	coord->getCoords(BLCOORD_CENTERLINE,cl);
	prevZ = cl[2];
	coord->setMostRecentReferenceZ(prevZ);

	// fake initial step of 1 mm 
	// (ensures sanity if reference bends in first step)
	stepPos = prevPos - prevDir;
	prevStepPos = stepPos - prevDir;
	stepZ = prevZ - 1.0;
}

void BLCoordinates::ReferenceCoordinates::PostUserTrackingAction(const G4Track *track)
{
	if(BLManager::getObject()->getState() != REFERENCE) return;

	// add 1 mm to final step
	// (ensures sanity if reference bends in final step)
	G4ThreeVector pos = track->GetPosition()+track->GetMomentumDirection();
	G4double z = track->GetTrackLength() + 1.0;
	addStep(pos,z,true);

	validRC = true;

	// for debugging -- dump segmentRCVector
	printf("\n# segmentRCVector\n#N minZ maxZ globX globY globZ axisX axisY axisZ angleDeg corr\n");
	G4double local[4];
	local[0] = local[1] = local[2] = local[3] = 0.0;
	for(unsigned i=0; i<segmentRCVector.size(); ++i) {
		G4double global[4];
		local[2] = (i==0 ? 0.0 : segmentRCVector[i].minZ);
		segmentRCVector[i].transform.getGlobal(local,global);
		G4ThreeVector axis = 
			segmentRCVector[i].transform.getRotation().getAxis();
		printf("%d %.4g %.4g %.4g %.4g %.4g %.4f %.4f %.4f %.1f %.3f\n",i,
			segmentRCVector[i].minZ,segmentRCVector[i].maxZ,
			global[0],global[1],global[2],axis[0],axis[1],axis[2],
			segmentRCVector[i].transform.getRotation().getDelta()/deg,
			segmentRCVector[i].correction);
	}
	printf("\n");
	fflush(stdout);
}

void BLCoordinates::ReferenceCoordinates::UserSteppingAction(const G4Step *step)
{
	// handle all states
	if(BLManager::getObject()->getState() != REFERENCE) return;

	G4Track *track = step->GetTrack();
	G4ThreeVector dir = track->GetMomentumDirection();
	if((dir-prevDir).mag2() > tolerance*tolerance)
		addStep(prevPos, prevZ);

	prevDir = dir;
	prevPos = track->GetPosition();
	prevZ = track->GetTrackLength();
	BLCoordinates::setMostRecentReferenceZ(prevZ);
}

void BLCoordinates::ReferenceCoordinates::addStep(const G4ThreeVector &pos, 
						G4double z, bool final)
{
	G4ThreeVector delta = pos - stepPos;
	G4double magDelta = delta.mag();

	// protect against delta being too small (should never happen,
	// as momentum direction changed significantly or first/final step
	// is >= 1 mm).
	if(magDelta < tolerance) return;

	// logic for rings is in CL coordinates; when we are lost, use CL
	// segment as starting point
	if(segmentCLVector[currentCL].firstRCSegment > segmentRCVector.size())
		segmentCLVector[currentCL].firstRCSegment = 
							segmentRCVector.size();

	G4ThreeVector rot = unitZ.cross(delta.unit());
	G4double a = asin(rot.mag());
	rot = rot.unit();
	// handle all 4 quadrants
	G4double b = unitZ.dot(delta.unit());
	if(a >= 0.0) {
		if(b >= 0.0)
			a = a;
		else
			a = M_PI - a;
	} else {
		if(b >= 0.0)
			a = a;
		else
			a = -M_PI + a;
	}

	G4RotationMatrix *p = 0;
	if(a > tolerance)
		p = new G4RotationMatrix(rot,a);
	G4ThreeVector tmp(0,0,-stepZ);
	if(p)
		tmp = *p * tmp;
	tmp += stepPos;

	G4double minZ = (segmentRCVector.size()==0 ? -DBL_MAX : stepZ);
	G4double maxZ = (final ? DBL_MAX : z);
	G4double interpPrev = stepZ + (z-stepZ)/4.0;
	G4double interpThis = z - (z-stepZ)/4.0;
	G4double correction = (z-stepZ)/magDelta;
//printf("addStep minZ=%.4g maxZ=%.4g start x,y,z=%.3f,%.3f,%.3f  axis=%.3f,%.3f,%.3f aDeg=%.1f  corr=%.6f\n",minZ,maxZ,stepPos.x(),stepPos.y(),stepPos.z(),rot[0],rot[1],rot[2],a/deg,correction);
	segmentRCVector.push_back(RCSegment(minZ,maxZ,
		BLCoordinateTransform(p,tmp),correction,interpPrev,interpThis));

	// set previous interpolate
	if(segmentRCVector.size() >= 2) {
		delta = (pos-prevStepPos);
		rot = unitZ.cross(delta.unit());
		a = asin(rot.mag());
		rot = rot.unit();
		// handle all 4 quadrants
		b = unitZ.dot(delta.unit());
		if(a >= 0.0) {
			if(b >= 0.0)
				a = a;
			else
				a = M_PI - a;
		} else {
			if(b >= 0.0)
				a = a;
			else
				a = -M_PI + a;
		}
		p = 0;
		if(a > tolerance)
			p = new G4RotationMatrix(rot,a);
		tmp = G4ThreeVector(0,0,-stepZ);
		if(p)
			tmp = *p * tmp;
		tmp += stepPos;
//printf("addStep prev interpolate x,y,z=%.3f,%.3f,%.3f  axis=%.3f,%.3f,%.3f aDeg=%.1f\n",stepPos.x(),stepPos.y(),stepPos.z(),rot[0],rot[1],rot[2],a/deg);
		segmentRCVector[segmentRCVector.size()-2].interpolate
						= BLCoordinateTransform(p,tmp);
	}

	// update after this segment
	prevStepPos = stepPos;
	stepPos = pos;
	stepZ = z;
}
