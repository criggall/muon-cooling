//	BLCoordinates.hh
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

#ifndef BLCOORDINATES_HH
#define BLCOORDINATES_HH

#include <vector>
#include "G4Allocator.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VUserTrackInformation.hh"
#include "G4UserTrackingAction.hh"
#include "G4UserSteppingAction.hh"
#include "BLCoordinateTransform.hh"
#include "BLManager.hh"

/**	BLCoordinateType enumerates the coordinate systems.
 **/
enum BLCoordinateType {BLCOORD_GLOBAL, BLCOORD_CENTERLINE, BLCOORD_REFERENCE,
			BLCOORD_LOCAL};

/**	class BLCoordinates manages the centerline coordinates.
 *
 *	Whenever an element is placed into the world group, the current 
 *	centerline coordinates are normally used, rather than global
 *	coordinates (the local coordinates of the world volume). 
 *
 *	Centerline coordinates are X,Y,Z, with Z along the nominal beam 
 *	direction, X is beam left, and Y is beam up. The centerline itself
 *	is always X=0,Y=0, and consists of a sequence of segmentCLs each of
 *	which is straight, and the centerline is continuous at each
 *	corner (bend in the centerline). Except within bending magnets
 *	the reference particle is supposed to travel directly down the
 *	centerline. Z is continuous along the centerline, including at
 *	corners; the first segmentCL extends to Z=-infinity, and the last
 *	segmentCL extends to Z=+infinity.
 *
 *	Note that the coordinate transform global->centerline is NOT unique
 *	-- near corners both segmentCLs may refer to a given point. The
 *	transform centerline->global is unique for a given segmentCL. The
 *	centerline itself is continuous at corners, but tracks will in
 *	general have a discontinuity in their centerline coordinates when
 *	transitioning from one segmentCL to the next. A track will transition
 *	from one segmentCL to the next when either its Z coordinate exceeds
 *	that of the initial segmentCL's end (corner), or when its radius
 *	exceeds the radiusCut for the initial segmentCL.
 *
 *	Initially the centerline coordinates are identical to global 
 *	coordinates. The start() function will change the starting location 
 *	and orientation of the first segmentCL, and the corner() function ends 
 *	the current segmentCL and begins a new one. start() must not be called
 *	after any element has been placed in the first segmentCL, or after
 *	corner() has been called.
 *
 *	Instances of this class represent the current centerline coordinates
 *	of a track, and are linked into the track via 
 *	G4Track::SetUserInformation() (an instance is registered with
 *	BLManager to manage the instance linked into the track,
 *	and they can also be used for individual points, if necessary).
 *	Static members and functions are used to access and create the
 *	centerline segmentCLs.
 *
 *	Because of the need to transition tracks from one segmentCL to the
 *	next, corners with acute angles are prohibited.
 *
 *	NOTE: BLManager expects G4Track::GetUserInformation() to return a 
 *	pointer to a BLCoordinates object, or to an object of a class derived
 *	from it. This means that any G4VUserTrackInformation class used in
 *	a simulation using BLManager should be derived from this class, and not
 *	directly from G4VUserTrackInformation (but verify is used to detect
 *	if it is not a BLCoordinates object).
 *
 *	Whenever local coordinates are used, the global-to-local coordinate
 *	transform must already have been set; the last one set is used.
 *
 *	Note the BLManager::UserSteppingAction() manages the updates to the
 *	current track's BLCoordinates before doing anything else.
 *
 *	NOTE: "CL" in a vairable means it is for centerline coordinates;
 *	"RC" means it is for reference coordinates.
 **/
class BLCoordinates : public G4VUserTrackInformation {
	enum { VERIFY=0x13579246 };
	int verify;
	unsigned int segmentCL;
	unsigned int segmentRC;
	BLCoordinateType prevType;
	G4double global[4];
	G4double centerline[4];
	G4double reference[4];
	int segmentCLRadiusFailure;
	BLCoordinateTransform *localTransform;
	G4bool findCLSegment(const G4double _global[4]);
	static G4Allocator<BLCoordinates> allocator;
public:
	/// Default constructor -- for UserSteppingAction or without point.
	BLCoordinates();

	/// Constructor -- finds the earliest segmentCL containing the point
	/// (global coordinates).
	BLCoordinates(const G4double _global[4]);

	/// Destructor.
	virtual ~BLCoordinates() { }

	/// operator new - using G4Allocator.
	/// NOTE: any derived class MUST re-implement this.
	void *operator new(size_t) {
		return (void *)allocator.MallocSingle();
	}

	/// operator delete - using G4Allocator.
	/// NOTE: any derived class MUST re-implement this.
	void operator delete(void *p) {
		allocator.FreeSingle((BLCoordinates *)p);
	}

	/// isValid() verifies that this is a valid BLCoordinates instance.
	/// NOTE: this is not virtual, and cannot be made virtual!
	G4bool isValid() const;

	/// Returns the global-to-local coordinate transform.
	BLCoordinateTransform *getLocalTransform() const
		{ return localTransform; }

	/// Sets the global-to-local coordinate transform.
	/// The coordinate trransform MUST be persistent through all uses.
	void setLocalTransform(BLCoordinateTransform *t)
		{ localTransform = t; }

	/// setGlobal() updates the centerline coordinates for the new point[]
	/// of the track; uses the current segmentCL #, and will update it if
	/// the track crosses into a new segmentCL.
	/// returns false if the track fails radiusCut in its current and 
	/// next segmentCL (or all segmentCL if ring is true).
	/// Also updates reference coordinates, if validRC is true and the
	/// centerline radiusCut succeeded.
	G4bool setGlobal(const G4double _global[4]);
	G4bool setGlobal(const G4ThreeVector &pos, G4double time);

	/// getCenterlineCoords() returns the centerline coordinates.
	/// DEPRECATED: use getCoords() instead.
	void getCenterlineCoords(G4double local[4]) const 
		{ local[0]=centerline[0]; local[1]=centerline[1];
		  local[2]=centerline[2]; local[3]=centerline[3]; }
	void getCenterlineCoords(G4ThreeVector& local) const
		{ local[0]=centerline[0]; local[1]=centerline[1];
		  local[2]=centerline[2]; }

	/// getCoords() returns the appropriate coordinate values.
	void getCoords(BLCoordinateType type, G4double local[4]);
	void getCoords(BLCoordinateType type, G4ThreeVector &local);

	/// getCLZ() returns the centerline z coordinate.
	G4double getCLZ() { return centerline[2]; }

	/// getRotation() returns the appropriate coordinate rotation for
	/// the previous call to getCoords().
	/// Given a G4ThreeVector B in global coordinates, this transforms
	/// it to centerline components:
	///   coord->getCoords(BLCOORD_CENTERLINE,local);
	///   B = coord->getRotation() * B;
	G4RotationMatrix &getRotation() { return getRotation(prevType); }

	/// getRotation() returns the rotation for the specified coordinates.
	/// Given a G4ThreeVector B in global coordinates, this transforms
	/// it to centerline components:
	///   B = coord->getRotation(BLCOORD_CENTERLINE) * B;
	G4RotationMatrix &getRotation(BLCoordinateType type);

	/// getCenterlineRotation() returns the rotation matrix used for the
	/// last setGlobal() of this instance. This is the rotation of
	/// the centerline wrt global coordinate axes.
	/// DEPRECATED -- use getRotation().
	G4RotationMatrix getCenterlineRotation() {
		if(segmentCL >= segmentCLVector.size())
			return G4RotationMatrix();
		return segmentCLVector[segmentCL].transform.getRotationInverse();
	}

	/// getCurrentRadiusCut() returns the radiusCut for the current CL
	/// segment.
	static double getCurrentRadiusCut() {
		if(currentCL >= segmentCLVector.size())
			return DBL_MAX;
		else
			return segmentCLVector[currentCL].radiusCut;
	}

	/// getTransformRotation() returns the rotation matrix of the
	/// global -> centerline coordinate transformation; it is the
	/// inverse of getCenterlineRotation().
	/// Given a G4ThreeVector B in global coordinates, this transforms
	/// it to centerline components:
	///   B = getTransformRotation() * B;
	/// DEPRECATED -- use getRotation().
	G4RotationMatrix getTransformRotation() {
		if(segmentCL >= segmentCLVector.size())
			return G4RotationMatrix();
		return segmentCLVector[segmentCL].transform.getRotation();
	}

	/// getSegmentCL() re turns the CL segment # (for debugging).
	int getSegmentCL() { return segmentCL; }

	/// update()
	/// This function ensures that a UserTrackInformation object is
	/// linked into the track; if none is, it creates a new
	/// BLCoordinates instance and links it in. It also calls setGlobal()
	/// for the BLCoordinates linked into the track.
	static void update(G4Track *track);

	/// Print from G4VUserTrackInformation.
	virtual void Print() const { }

	// static variables and functions to define and manipulate segmentCLs.
	// NOTE: these routines must all be called before the first instance
	// is created. That is, the centerline must be completely laid out
	// before any tracks are created.
private:
	enum State { UNUSED, CONSTRUCTING, TRACKING };
	struct CLSegment {
		G4double minZ;
		G4double maxZ;
		G4double radiusCut;
		BLCoordinateTransform transform;
		unsigned int firstRCSegment;
		CLSegment(G4double _minZ, G4double _maxZ, G4double _radiusCut,
				const BLCoordinateTransform &_transform) :
							transform(_transform)
		{ minZ=_minZ; maxZ=_maxZ; radiusCut=_radiusCut; 
		  firstRCSegment=99999999; }
		friend class BLCoordinates;
	};
	static std::vector<CLSegment> segmentCLVector;
	static unsigned currentCL;
	static G4bool ring;
	static State state;
	static G4bool steppingVerbose;
	static void init();
public:
	/// start() sets the starting point and direction for the
	/// first segmentCL. Must not be called after corner() is
	/// called, or after any element has been placed using centerline
	/// coordinates. point[] is in global coordinates, z is the starting
	/// value of z in CENTERLINE coordinates, _ring indicates whether
	/// a ring is present (so jumps in segmentCL # are permitted), and
	/// rotation should be the rotation of the first segmentCL's
	/// centerline wrt the global coordinate axes (i.e. it is an
	/// active rotation of the centerline as an object wrt the global
	/// coordinate axes).
	static void start(const G4double point[4], G4double z,
			G4RotationMatrix *rotation=0, G4double radiusCut=0.0,
			G4bool _ring=false);

	/// corner() creates a corner in the centerline by ending the
	/// current segmentCL and beginning a new one. z is the position
	/// in CENTERLINE coordinates of the current segmentCL, and rotation
	/// is relative to the current CENTERLINE coordinates (i.e. as an
	/// active rotation of the centerline as an object wrt the CENTERLINE
	/// coordinate axes of the current segmentCL). z must be greater than
	/// the centerline z of the previous corner or start. rotation can be
	/// an identity matrix (e.g. to merely set a new radiusCut).
	static void corner(G4double z, G4RotationMatrix *rotation,
							G4double radiusCut);

	/// getCurrentTransform() returns the current BLCoordinateTransform
	/// for the centerline coordinates (i.e. the coordinate transform 
	/// global->centerline of the current segmentCL).
	static BLCoordinateTransform &getCurrentTransform();

	/// getCurrentRotation() returns the current G4RotationMatrix for
	/// the centerline coordinates (i.e. the active rotation of the
	/// current segmentCL, as an object rotated wrt global coordinate
	/// axes).
	/// Given a G4ThreeVector B in global coordinates, this transforms
	/// it to centerline components:
	///   B = *getCurrentRotation() * B;
	static G4RotationMatrix *getCurrentRotation();

	/// getCurrentGlobal() returns global coordinates from centerline coordinates,
	/// using the current segmentCL.
	static void getCurrentGlobal(const G4double local[4], G4double _global[4]);

	/// getCurrentGlobal() returns global coordinates from centerline coordinates,
	/// using the current segmentCL.
	static void getCurrentGlobal(const G4ThreeVector& local, G4ThreeVector& _global);

	/// getGlobalAnywhere() returns global coordinates from centerline
	/// coordinates, using the lowest segmentCL containing local.
	static void getGlobalAnywhere(const G4ThreeVector& local, 
						G4ThreeVector& _global);

	/// getGlobalRotation() returns the local->global rotation matrix
	/// coordinates, using the lowest segmentCL containing local.
	static G4RotationMatrix *getGlobalRotation(const G4ThreeVector& local);

	/// getRadiusCut() returns the radius cut for the lowest
	/// segmentCL containing local.
	static G4double getRadiusCut(const G4ThreeVector& local);

	/// getCoordinateType() converts a G4String name to a
	/// BLCoordinateType. Names are: global, centerline, reference, local;
	/// only the first letter is used, case insensitive.
	/// Also enables the generation of Reference Coordinates, if they
	/// are selected by s.
	static BLCoordinateType getCoordinateType(G4String s);

	// static variables and function to manipulate segmentRC-s.
	// NOTE: reference coordinates can easily get onto the wrong segment.
	// This is best avoided by restricting the beam to a small radius from
	// the reference track so that ambiguities are avoided. But that is
	// not always possible (e.g. in a helical cooling channel) -- in that
	// case it is best to start tracking in a sement where the reference
	// particle is close to the centerline.
	//
	// updateReferenceCoordinates() always sets segmentRC to the highest
	// segment N such that Z<=seg[N].maxZ; it could have Z<seg[N].minZ 
	// but will have Z>seg[N-1].maxZ (here Z is evaluated for each seg[]).
	// 
	// In the plane of two successive segments N and N+1, there is a
	// triangle with Z>seg[N].maxZ and Z<seg[N+1].minZ on the outside of
	// the bend; there is also a triangle valid for both segments. This
	// is handled by creating an interpolation transform between every
	// pair of segments; seg[N-1].interpolate is used when 
	// Z<seg[N].interpolatePrev, and seg[N].interpolate is used when
	// Z>seg[N].interpolateThis (1/4 and 3/4 along the segment in Z).
	// The interpolation covers the entire omitted triangle, and much
	// of the overlap triangle (especially near the reference track).
	//
	// updateReferenceCoordinates() must keep up as the Track is tracked,
	// but the interpolation is performed only when the reference
	// coordinates are needed in getCoords().
	//
	// correction accounts for curvature in the reference track, in which
	// case the reference path length is larger than the straight line 
	// Segment. The difference is pro-rated in Z within the Segment, but
	// is ignored in X and Y (in keeping with the straight-line
	// approximation). Without this, values of Z would not be continuous,
	// and worse, Segment errors would occur near the downstream edge.
private:
	struct RCSegment {
		G4double minZ;
		G4double maxZ;
		G4double correction;
		BLCoordinateTransform transform;
		G4double interpolatePrev;
		G4double interpolateThis;
		BLCoordinateTransform interpolate;
		RCSegment(G4double _minZ, G4double _maxZ,
			const BLCoordinateTransform &_transform, G4double corr,
			G4double _interpolatePrev, G4double _interpolateThis) 
				: transform(_transform), interpolate(_transform)
		{ minZ=_minZ; maxZ=_maxZ; correction=corr;
		  interpolatePrev=_interpolatePrev;
		  interpolateThis=_interpolateThis; }
		friend class BLCoordinates;
	};
	static G4bool useRC;
	static G4bool validRC;
	static std::vector<RCSegment> segmentRCVector;
	static int currentRC;
	static G4double tolerance;
	static G4double mostRecentReferenceZ;
	class ReferenceCoordinates;
	int referenceCLsegment; // current CL segment of reference coords.
	BLCoordinateTransform *referenceTransformUsed;
	/// updateReferenceCoordinates() will update the reference coordinates,
	/// given the centerline coordintes are already updated and valid.
	/// Called by setGlobal().
	void updateReferenceCoordinates();
public:
	/// useReferenceCoordinates() will arrange to fill segmentRCVector[]
	/// with reference coordinate segments as the reference particle is 
	/// tracked. Must be called before the reference particle is tracked.
	/// Unless this is called, reference coordinates will not be available.
	/// Multiple calls are the same as one.
	static void useReferenceCoordinates();

	/// setMostRecentReferenceZ() sets the most recent z value of the
	/// reference coords. VALID ONLY WHILE TRACKING REFERENCE PARTICLE!
	static void setMostRecentReferenceZ(G4double z)
		{ mostRecentReferenceZ = z; } 

	/// getMostRecentReferenceZ() returns the most recent z value of the
	/// reference coords. VALID ONLY WHILE TRACKING REFERENCE PARTICLE!
	static G4double getMostRecentReferenceZ()
		{ return mostRecentReferenceZ; } 
};

#endif // BLCOORDINATES_HH
