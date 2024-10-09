//	BLCoordinateTransform.hh
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

#ifndef BLCOORDINATETRANSFORM_HH
#define BLCOORDINATETRANSFORM_HH

#include <stdio.h>

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#define DUMPROTATION(rot,str)						\
{ G4ThreeVector x(1.0,0.0,0.0), z(0.0,0.0,1.0);				\
  G4ThreeVector rx = rot * x;						\
  G4ThreeVector rz = rot * z;						\
  printf("%s: %.3f,%.3f,%.3f / %.3f,%.3f,%.3f / %.3f,%.3f,%.3f\n", 	\
    str,rot.xx(),rot.xy(),rot.xz(),rot.yx(),rot.yy(),rot.yz(), 		\
    rot.zx(),rot.zy(),rot.zz());					\
  printf("        rot*x=%.3f,%.3f,%.3f   rot*z = %.3f,%.3f,%.3f\n", 	\
    rx[0],rx[1],rx[2],rz[0],rz[1],rz[2]);				\
}

/**	class BLCoordinateTransform defines a linear coordinate transform.
 *
 *	An object is first defined in its local coordinates (i.e. with the
 *	beam centerline along +Z). The object is rotated in place, and is then 
 *	translated to position. That's how the constructors work; the
 *	actual transformation goes the other way. And normally the transform
 *	of each group is applied to the transforms of its elements, so each
 *	placed element's transform goes global->local in one fell swoop.
 *
 *	This class is completely inline.
 **/
class BLCoordinateTransform {
	G4RotationMatrix rot_p2l;	/// rotation parent->local
	G4RotationMatrix rot_l2p;	/// rotation local->parent
	G4ThreeVector position;		/// parent position of local (0,0,0)
	G4bool rotated;		/// for efficiency in un-rotated xforms
	void checkRotation() {
		// two test vectors is enough
		G4ThreeVector x(1.0,0.0,0.0), y(0.0,1.0,0.0);
		x = rot_p2l * x;
		y = rot_p2l * y;
		if(fabs(y[0]) > 1e-8 || fabs(y[1]-1.0) > 1e-8 || 
		   fabs(y[2]) > 1e-8 || fabs(x[0]-1.0) > 1e-8 ||
		   fabs(x[1]) > 1e-8 || fabs(x[2]) > 1e-8)
			rotated = true;
		else
			rotated = false;
	}
public:
	/// Default constructor. Identity transformation.
	BLCoordinateTransform() : rot_p2l(), rot_l2p(), position(),
				  rotated(false) { }

	/// Constructor for a simple unrotated position.
	/// pos is the parent location of the object.
	BLCoordinateTransform(const G4ThreeVector& pos) : rot_p2l(), 
				rot_l2p(), position(pos), rotated(false) { }

	/// Constructor for a rotation and position
	BLCoordinateTransform(const G4RotationMatrix &rot, 
				const G4ThreeVector &pos) :
				rot_p2l(rot), rot_l2p(rot),
				position(pos), rotated(false) 
		{ rot_p2l.invert(); checkRotation(); }

	/// Constructor for a rotation and position
	BLCoordinateTransform(const G4RotationMatrix *rot, 
				const G4ThreeVector &pos) :
				rot_p2l(), rot_l2p(),
				position(pos), rotated(false) 
		{ if(rot) rot_p2l = rot_l2p = *rot; rot_p2l.invert(); 
		  checkRotation(); 
		  }

	/// getLocal() - get local coordinates from global coordinates.
	/// This assumes that the BLCoordinateTransform transform for each
	/// enclosing group has been applied (in descending order).
	/// (otherwise "global" is really just "parent")
	/// Both arguments can be the same 4-vector (Point); the time coord
	/// is unchanged.
	void getLocal(G4double local[4], const G4double global[4]) const {
		G4ThreeVector l, g;	//@@ can probably optimize this....
		g[0] = global[0];
		g[1] = global[1];
		g[2] = global[2];
		getLocal(l,g);
		local[0] = l[0];
		local[1] = l[1];
		local[2] = l[2];
		local[3] = global[3];
	}

	/// getLocal() - get local coordinates from global coordinates.
	/// This assumes that the BLCoordinateTransform transform for each
	/// enclosing group has been applied (in descending order).
	/// (otherwise "global" is really just "parent")
	/// The two arguments must NOT be the same.
	void getLocal(G4ThreeVector& local, const G4ThreeVector& global) const
	{
		G4ThreeVector temp;
		if(rotated)
			local = rot_p2l * (global - position);
		else
			local = global - position;
	}

	/// getGlobal() - get Global coordinates from local coordinates.
	/// This assumes that the BLCoordinateTransform transform for each
	/// enclosing group has been applied (in descending order).
	/// (otherwise "global" is really just "parent")
	/// Both arguments can be the same 4-vector (Point); the time coord
	/// is unchanged.
	void getGlobal(const G4double local[4], G4double global[4]) const {
		G4ThreeVector l, g;	//@@ can probably optimize this....
		l[0] = local[0];
		l[1] = local[1];
		l[2] = local[2];
		getGlobal(l,g);
		global[0] = g[0];
		global[1] = g[1];
		global[2] = g[2];
		global[3] = local[3];
	}

	/// getGlobal() - get global coordinates from local coordinates.
	/// This assumes that the BLCoordinateTransform transform for each
	/// enclosing group has been applied (in descending order).
	/// (otherwise "global" is really just "parent")
	/// The two arguments should NOT be the same.
	void getGlobal(const G4ThreeVector& local, G4ThreeVector& global) const
	{
		if(rotated)
			global = rot_l2p * local + position;
		else
			global = local + position;
	}

	/// apply() will apply a BLCoordinateTransform to this one.
	/// This function is used to apply the transform for a parent to
	/// objects placed into the parent group -- the objects' rotations
	/// and positions are specified wrt the parent in the parent's local
	/// coordinates. This routine should be called during construction,
	/// not during tracking (it inverts a matrix).
	void apply(const BLCoordinateTransform& parentTransform) {
		if(parentTransform.rotated) {
			position = parentTransform.rot_l2p * position +
						parentTransform.position;
			if(rotated)
				rot_p2l = rot_p2l * parentTransform.rot_p2l;
			else
				rot_p2l = parentTransform.rot_p2l;
			rot_l2p = rot_p2l;
			rot_l2p.invert();
		} else {
			position += parentTransform.position;
		}
		checkRotation();
	}

	/// getPosition() returns the position of the transform.
	G4ThreeVector& getPosition() { return position; }

	/// getRotation() returns the Rotation of the transform.
	G4RotationMatrix& getRotation() { return rot_p2l; }

	/// getRotationInverse() returns the inverse Rotation of the transform
	/// (this is the active rotation of the object).
	/// OK to use during tracking (the inverse is computed just once).
	G4RotationMatrix& getRotationInverse() { return rot_l2p; }

	/// isRotated() returns true iff this transformation includes rotation.
	G4bool isRotated() const { return rotated; }
};

#endif // BLCOORDINATETRANSFORM_HH
