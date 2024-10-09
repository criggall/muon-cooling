//	BLWindowShape.hh
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

#ifndef BLWINDOWSHAPE_HH
#define BLWINDOWSHAPE_HH

#include <vector>
#include "globals.hh"

/**	struct BLWindowShape describes the shape of a window.
 *	This is intended to be used to create a G4PolyCone for the window, 
 *	a G4PolyCone for its contents, and a G4Tubs for the flange and
 *	possibly the pipe to which the flange is bolted (bolts are not 
 *	modeled).
 *
 *	This class may be subclassed to implement other methods of
 *	constructing a window shape.
 *
 * 	The standard constructor reads data from filename.
 *	The file format is a series of lines:
 *		First character # means comment, * means printed comment.
 *		Comment and blank lines are ignored. Units are mm.
 *		The first line contains the 4 flange variables:
 *			innerR outerR insideZ outsideZ
 *		The remaining lines contain 3 values for a given radius:
 *			r z t
 *		The first line must have r=0.0 and have the largest z value
 *		and the largest z+t value. z is the inside of the window,
 *		z+t is the outside of the window. All values must be positive.
 *		Successive r values must increase by at least 0.010 mm.
 *		Successive z values and z+t values must decrease by at least
 *		0.010 mm.
 *		flangeInnerRadius should equal the last r value.
 *		Any z origin may be used (it will be subtracted away).
 *	This is intended to be easy to interface to window design
 *	spreadsheets (export a list of 3 columns delimited by spaces, then
 *	add the appropriate comments and flange values at the top).
 **/
struct BLWindowShape {
	G4double flangeInnerRadius;	// should be same as last r[] value
	G4double flangeOuterRadius;	// must be larger than flangeInnerRadius
	G4double flangeInsideZ;		// in same Z coords as the window
	G4double flangeOutsideZ;	// in same Z coords as the window
	std::vector<G4double> r;	// first must be 0.0
	std::vector<G4double> z;	// Z coord of inside, first must be 
					// largest, and they must decrease
					// monotonically by at least 0.010 mm
	std::vector<G4double> t;	// thickness along z
public:
	/// default constructor.
	BLWindowShape() : r(), z(), t() { flangeInnerRadius=0.0;
		flangeOuterRadius=0.0; flangeOutsideZ=0.0; flangeInsideZ=0.0; }

	/// Constructor. Reads directory/filename for the data.
	BLWindowShape(G4String filename);

	/// destructor.
	virtual ~BLWindowShape() { }
};

#endif // BLWINDOWSHAPE_HH
