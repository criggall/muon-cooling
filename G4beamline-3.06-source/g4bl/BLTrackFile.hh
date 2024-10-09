//	BLTrackFile.hh
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

#ifndef BLTRACKFILE_HH
#define BLTRACKFILE_HH

#include <stdio.h>
#include "globals.hh"
#include "G4ThreeVector.hh"

enum BLTrackFileStatus { BLTF_OK=0, BLTF_ENDOFFILE=1, BLTF_ERROR=2 , 
			BLTF_DUMMY=3};

/**	class BLTrackFile creates or reads an ASCII file containing tracks.
 *
 *      The file format is ASCII:
 *      Lines begining with # are comments
 *      First line is a structured comment (if not present the input
 *      routine issues a warning):
 *              #BLTrackFile ... user comment...
 *      Second line is a comment giving the column names:
 *              #x y z Px Py Pz t PDGid EventID TrackID ParentID weight
 *      Third line is a comment giving the units:
 *              #mm mm mm MeV/c MeV/c MeV/c ns - - - - -
 *      Thereafter follow the tracks, one per line.
 *	If a comment line containing 'cm cm cm' is read, the units are
 *	changed from mm to cm.
 *
 *	This class now supports two formats: the original BLTrackFile
 *	format, and an extended BLTrackFile2 format which appends additional
 *	fields to the original format:
 *		Bx, By, Bz (Tesla)
 *		Ex, Ey, Ez (Megavolts/meter)
 *		ProperTime (ns)
 *		PathLength (mm)
 *		PolX, PolY, PolZ (polarization)
 *		InitX, InitY, InitZ (position when track was created, mm)
 *		InitT (time when track was created, ns)
 *		InitialKE (MeV when track was created)
 *
 *	Note that calls to read/read2 and write/write2 must match the value
 *	of version given to the constructor. A version=1 file can be read
 *	as version=2, but the extra fields will be read as zero. A version=2
 *	file can be read as version=1, but the extra fields will be ignored.
 *
 *	While the input routine can handle initial spaces in the first
 *	column, it is STRONGLY suggested you not put any there (so cut/grep
 *	will work). Any fixed or floating-point format will do; PDGid,
 *      EventID, TrackID, and ParentID are integers (but are read as doubles,
 *	so .0 can be appended).
 *      Common PDGid-s:
 *              e-      11              e+      -11
 *              mu-     13              mu+     -13
 *              pi+     211             pi-     -211
 *              proton  2212            anti_proton -2212
 *              neutron 2112            anti_neutron -2112
 *              gamma   22
 *
 *	Removed features: the ability to omit Weight is removed.
 **/
class BLTrackFile {
	BLTrackFileStatus status;
	FILE *file;
	G4String mode;
	int version;
	double unit;
public:
	/// Default constructor. Reads return immediate EOF; writes do nothing.
	BLTrackFile();

	/// Constructor. Opens the file; mode="r" or "w".
	/// comment is an identifying comment in the first line of the file
	/// (ignored for mode=r).
	BLTrackFile(G4String filename, G4String comment, G4String _mode, 
								int _version=1);

	/// Destructor. Closes the file.
	~BLTrackFile();

	/// write() will write one track to the file (mode must have been "w").
	/// Units are standard geant4 internal units.
	/// This is BLTrackFile version 1.
	BLTrackFileStatus write(const G4ThreeVector& pos, G4double time,
			const G4ThreeVector& momentum, int PDGid, int eventId,
			int trackId, int parentId, G4double weight);

	/// write2() will write one track to the file (mode must have been "w").
	/// Units are standard geant4 internal units.
	/// This is BLTrackFile version 2.
	BLTrackFileStatus write2(const G4ThreeVector& pos, G4double time,
			const G4ThreeVector& momentum, int PDGid, int eventId,
			int trackId, int parentId, G4double weight,
			const G4ThreeVector& bfield,
			const G4ThreeVector& Efield,
			G4double properTime, G4double pathLength,
			const G4ThreeVector& polarization,
			const G4ThreeVector& initialPos,
			G4double initialT, G4double initialKE);

	/// read() will read one track from the file (mode must have been "r").
	/// Units are standard geant4 internal units.
	/// This is BLTrackFile version 1.
	BLTrackFileStatus read(G4ThreeVector& pos, G4double& time,
			G4ThreeVector& momentum, int& PDGid, int& eventId,
			int& trackId, int& parentId, G4double &weight);

	/// read2() will read one track from the file (mode must have been "r").
	/// Units are standard geant4 internal units.
	/// This is BLTrackFile version 2.
	BLTrackFileStatus read2(G4ThreeVector& pos, G4double& time,
			G4ThreeVector& momentum, int& PDGid, int& eventId,
			int& trackId, int& parentId, G4double &weight,
			G4ThreeVector& bfield, G4ThreeVector& Efield,
			G4double& properTime, G4double& pathLength,
			G4ThreeVector& polarization, G4ThreeVector& initialPos,
			G4double& initialT, G4double& initialKE);

	/// flush() will flush buffers to the file.
	void flush() { fflush(file); }
};

#endif // BLTRACKFILE_HH
