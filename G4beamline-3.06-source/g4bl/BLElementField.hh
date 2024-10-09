//	BLElementField.hh
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

#ifndef BLELEMENTFIELD_HH
#define BLELEMENTFIELD_HH

#include "globals.hh"

const G4double LARGE=9.0e99;


/**	class BLElementField - interface for the EM field of one element
 *
 *	This is the interface class used by BLGlobalField to compute the field
 *	value at a given point[].
 *
 *	An element command that represents an element with an EM field will
 *	derive a class from this one and implement the computation for the 
 *	element. The construct() function of the element command class will add
 *	the derived object into BLGlobalField. Such derived classes are
 *	typically local to the file implementing the element command.
 **/
class BLElementField {
	BLElementField(const BLElementField&);	// undefined
	BLElementField& operator=(const BLElementField&);	// undefined
protected:
	G4double minX, minY, minZ, maxX, maxY,maxZ;
public:
	/// Constructor.
	BLElementField()
		{ minX=minY=minZ=-LARGE; maxX=maxY=maxZ=LARGE; }

	/// Destructor.
	virtual ~BLElementField() { }

	/// setGlobalPoint() ensures that the point is within the global
	/// bounding box of this ElementField. global coordinates.
	/// Normally called 8 times for the corners of the local bounding
	/// box, after a local->global coordinate transform.
	/// If never called, the global bounding box is infinite.
	/// BEWARE: if called only once, the bounding box is just a point.
	void setGlobalPoint(const G4double point[4]) {
		if(minX == -LARGE || minX > point[0]) minX = point[0];
		if(minY == -LARGE || minY > point[1]) minY = point[1];
		if(minZ == -LARGE || minZ > point[2]) minZ = point[2];
		if(maxX == LARGE || maxX < point[0]) maxX = point[0];
		if(maxY == LARGE || maxY < point[1]) maxY = point[1];
		if(maxZ == LARGE || maxZ < point[2]) maxZ = point[2];
	}

	/// isInBoundingBox() returns true if the point is within the
	/// global bounding box. global coordinates.
	bool isInBoundingBox(const G4double point[4]) const {
		if(point[2] < minZ || point[2] > maxZ) return false;
		if(point[0] < minX || point[0] > maxX) return false;
		if(point[1] < minY || point[1] > maxY) return false;
		return true;
	}

	// getBoundingBox() returns the min/max values.
	void getBoundingBox(G4double &xMin, G4double &xMax, G4double &yMin,
			G4double &yMax, G4double &zMin, G4double &zMax) const {
		xMin=minX; xMax=maxX;
		yMin=minY; yMax=maxY;
		zMin=minZ; zMax=maxZ;
	}

	// setBoundingBox() sets the bounding box.
	void setBoundingBox(G4double &xMin, G4double &xMax, G4double &yMin,
			G4double &yMax, G4double &zMin, G4double &zMax) {
		minX=xMin; maxX=xMax;
		minY=yMin; maxY=yMax;
		minZ=zMin; maxZ=zMax;
	}

	/// addFieldValue() will add the field value for this element to
	/// field[].
	/// Implementations must be sure to verify that point[] is within 
	/// the field region, and do nothing if not.
	/// point[] is in global coordinates and geant4 units; x,y,z,t.
	/// field[] is in geant4 units; Bx,By,Bz,Ex,Ey,Ez.
	/// For efficiency, the caller may (but need not) call
	/// isInBoundingBox(point), and only call this function if that
	/// returns true.
	virtual void addFieldValue(const G4double point[4], G4double field[6])
								const = 0;
};

#endif // BLELEMENTFIELD_HH
