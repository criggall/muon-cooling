//	BLTrackNTuple.hh
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

#ifndef BLTRACKNTUPLE_HH
#define BLTRACKNTUPLE_HH

#include "globals.hh"
#include "G4Step.hh"

#include "BLNTuple.hh"
#include "BLEvaluator.hh"
#include "BLCoordinates.hh"


/**	class BLTrackNTuple is an interface class to an NTuple for tracks.
 *
 *	This class implements two new formats that add additional fields
 *	to the NTuple: asciiEnhanced, rootEnhanced.
 **/
class BLTrackNTuple {
	BLNTuple *ntuple;
	G4String require;
	BLEvaluator *eval;
	int noSingles;
	int nFields;
	BLCoordinateType coordinateType;
protected:
	/// Constructor.
	BLTrackNTuple();
public:
	/// Destructor.
	virtual ~BLTrackNTuple() { 
		if(ntuple) delete ntuple;
		if(eval) delete eval;
	}

	/// close() will close the ntuple
	void close() { if(ntuple) ntuple->close(); ntuple = 0; }

	/// create() will create a new NTuple for writing.
	/// Note that the field names are handled internally, and that the
	/// format can change the fields used.
	/// format must be a known format of NTuple ("" is default).
	/// "category/name" is the name of the NTuple (the "/" may be omitted 
	/// if category is null, depending on format).
	/// filename is a hint, and is not used by all formats.
	/// coordinateType is the type of coordiantes to use.
	/// require is an expression using track variables that must be nonzero
	/// to enter a Step into the NTuple.
	/// If noSingles is nonzero, the NTuple itself is not filled, but
	/// callbacks are made (to fill NTuples created by the ntuple command).
	static BLTrackNTuple *create(G4String format, G4String category,
			G4String name, G4String filename,
			BLCoordinateType coordinateType=BLCOORD_CENTERLINE, 
			G4String _require="", int _noSingles=0);

	/// appendTrack() appends the data for this Track to the
	/// NTuple. Uses arguments to create():
	/// If the require expression is nonempty and evaluates to zero,
	/// no entry is made to the NTuple or callbacks.
	/// If noSingles is nonzero, the NTuple itself is not filled, but
	/// callbacks are made (to fill NTuples created by the ntuple command).
	/// If coordinateType==BLCOORD_LOCAL (in create()), the local
	/// coordinates must be put into the track before calling this routine.
	/// returns true, unless the require test fails.
	bool appendTrack(const G4Track *track);

	/// appendTrack() appends the data for this Track, with updated
	/// values of the other arguments
	/// returns true, unless the require test fails.
	bool appendTrack(const G4Track *track, double time, 
				G4ThreeVector position, G4ThreeVector momentum,
				double properTime, double trackLength);

	/// getName() returns the name of the NTuple (not including category)
	G4String getName() const { return ntuple->getName(); }

	/// getFields() returns a colon-separated list of the fields.
	G4String getFields() { return ntuple->getFields(); }

	/// annotate() will add an annotation line to any ASCII format.
	/// Normally line should begin with a "#", and have no "\r" or "\n".
	/// line == "" will output an empty line to ASCII files.
	/// Ignored for NTuples that don't support annotations.
	virtual void annotate(G4String line) { ntuple->annotate(line); }

	/// needsReference() returns true if this NTuple needs the Reference
	/// particle.
	virtual bool needsReference() { return ntuple->needsReference(); }

	/// getFormatList() returns a list of known formats.
	/// This class implements asciiEnhanced and rootEnhanced.
	static G4String getFormatList();
};

#endif // BLTRACKNTUPLE_HH
