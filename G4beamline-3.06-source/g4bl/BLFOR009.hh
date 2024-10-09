//	BLFOR009.hh
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
#ifndef BLFOR009_HH
#define BLFOR009_HH

#include <stdio.h>
#include <map>

#include "globals.hh"
#include "G4ThreeVector.hh"

#include "BLTrackFile.hh" // BLTrackFileStatus only

/* A Region is placed in Z at the first track written to it. All tracks
   aving Z within REGION_SIZE of that first track are put into the same Region.
   If z is monotonically increasing during tracking, this handles a ring
   correctly (for N orbits around the ring, a single virtualdetector will
   generate N regions).
 */
const double REGION_SIZE = 1.0*cm;
struct Region;	// in BLFOR009.cc
struct RegionCompare {
	bool operator()(const double a, const double b) const
		{ return a < b-REGION_SIZE; }
};
typedef std::map<double,Region*,RegionCompare> RegionMap;

/**	class BLFOR009 creates an ASCII file containing tracks in ICOOL 
 *	FOR009.DAT format.
 *
 *	NOTE: Because the FOR009.DAT format requires the particles be sorted
 *	by region (JSRG), all tracks are stored in memory until close() is
 *	called; the file is written at that time. The ICOOL region is
 *	determined from the Z position of the track, and is simply 
 *	sequentially assigned as tracks at different Z positions are written.
 *	The region is +-REGION_SIZE in Z from the first track of each region.
 *
 *	XP[], PP[], BFLD[], and EFLD[] should all be converted to
 *	Centerline coordinates before calling write().
 *
 *	This data format comes from ICOOL v 2.77, User's Guide section 4.2.3.
 *	The first line of the file is the title.
 *	The second and third lines are comments, supposedly units and column 
 *	labels.
 *	There follow the tracks, one per line, sorted by "region". The first
 *	track of each region should be the "reference" particle in ICOOL
 *	parlance; in g4beamline it is the reference particle -- if multiple
 *	reference particles intersect the region, the first will be 
 *	"reference", and the following ones will be considered "beam".
 *	The variables for each track are:
 *		IEVT		(I) event #
 *		IPNUM		(I) track # for this event
 *		IPTYP		(I) particle type: >0 for positive charge, 
 *					<0 for negative
 *					1=e, 2=mu, 3=pi, 4=K, 5=proton
 *		IPFLG		(I) "flag", always 0
 *		JSRG		(I) region number (see above)
 *		TP		(F) time (sec)
 *		XP[3]		(F) position (meters)
 *		PP[3]		(F) momentum (GeV/c)
 *		BFLD[3]		(F) Magnetic field (Tesla)
 *		EVTWT		(F) weight
 *		EFLD[3]		(F) Electric field (V/meter)
 *		SARC		(F) arclength (meter) -- set to 0.0
 *		POL[3]		(F) spin -- set to 0.0
 *
 *	Note that event 0 will be ignored, as that value of IEVT is
 *	reserved for the reference particle.
 **/
class BLFOR009 {
	G4String filename;
	BLTrackFileStatus status;
	FILE *file;
	G4String mode;
	RegionMap *map;
public:
	/// Constructor. Opens the file; mode="r" or "w".
	/// title is an identifying title in the first line of the file
	/// (ignored for mode=r).
	/// READ MODE IS UNIMPLEMENTED.
	BLFOR009(G4String _filename, G4String title, G4String mode);

	/// Destructor. Writes any data and closes the file.
	~BLFOR009();

	/// close() writes the data to the file and closes it.
	void close();

	/// write() will write one track to the file (mode must have been "w").
	/// Units are standard geant4 internal units.
	/// Note that weight is optional.
	BLTrackFileStatus write(G4ThreeVector& pos, G4double time,
			G4ThreeVector& momentum, G4ThreeVector& Bfield,
			G4ThreeVector& Efield, int PDGid, int eventId,
			int trackId, int parentId, G4double weight=1.0);

	/// read() will read one track from the file (mode must have been "r").
	/// Units are standard geant4 internal units.
	/// Use the other version of read() if you want to use the weight.
	/// UNIMPLEMENTED.
	BLTrackFileStatus read(G4ThreeVector& pos, G4double& time,
			G4ThreeVector& momentum, G4ThreeVector& Bfield,
			G4ThreeVector& Efield, int& PDGid, int& eventId,
			int& trackId, int& parentId);
	/// read() will read one track from the file (mode must have been "r").
	/// Units are standard geant4 internal units.
	/// UNIMPLEMENTED.
	BLTrackFileStatus read(G4ThreeVector& pos, G4double& time,
			G4ThreeVector& momentum, G4ThreeVector& Bfield,
                        G4ThreeVector& Efield, int& PDGid, int& eventId,
			int& trackId, int& parentId, G4double &weight);
};

#endif // BLFOR009_HH
