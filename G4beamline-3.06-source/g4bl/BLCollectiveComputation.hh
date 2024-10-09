//	BLCollectiveComputation.hh
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

#ifndef BLCOLLECTIVECOMPUTATION_H
#define BLCOLLECTIVECOMPUTATION_H

#include "G4Event.hh"
#include "G4Track.hh"

/**	struct BLTrackData -- data related to a track (for collective tracking)
 **/
struct BLTrackData {
	G4Event *event;
	G4Track *track;
public:
	/// Constructor.
	BLTrackData(G4Event *e, G4Track *t) { event = e;  track=t; }

	/// destructor. DOES NOT DELETE ANY OTHER OBJECTS.
	~BLTrackData() { }
};

/**	class BLCollectiveComputation -- base class for a collective computation
 *
 *	Normally a command class that implements a collective computation 
 *	will be derived from this class as well as BLCommand.
 *	Its command() function should register itself with
 *	    BLRunManager::getObject()->registerCollectiveComputation(this);
 *	and then call
 *	    BLRunManager::getObject()->setCollectiveMode(true);
 *
 *	If steps in time are also desired, it should call
 *	    BLRunManager::getObject()->setDeltaT(0.1*ns);
 *	(with desired delta-T). Note that steps in time are not exact, so
 *	the computation should probably check that the variation in time
 *	among the tracks is small enough (in practive, 1 ps accuracy is
 *	normal).
 *
 *	Note that setDeltaT() can be called at any time, and will affect the
 *	next step. It is probably desirable to adjust deltaT depending on the
 *	importance of the collective computation (when the bunch is diffuse,
 *	bigger deltaT won't affect the accuracy but will reduce the CPU time).
 *	Caveat: if steps in time are to be used, setDeltaT() must be called
 *	before tracking begins, and should never be <= 0.0 (values <= 0.0
 *	means each track takes individual steps in space, which is probably
 *	not useful).
 *
 *	Note that runmanager->rejectCollectiveStep() can be called during
 *	tracking or in collectiveStep() to reject the current step and 
 *	re-process it. This only makes sense if setDeltaT() has been
 *	called with a smaller value of deltaT. Attempts to reject the very
 *	first step are silently ignored (that step is used to process all
 *	tracks to the same time).
 *
 *	The basic sequence of collective tracking, including callbacks from
 *	BLRunManger to each registered collective computation is:
 *
 *	(generate beam tracks from beam commands, filling trackVector)
 *	stepTime = (largest time from beam tracks)
 *	beginCollectiveTracking()
 *	first=true
 *	loop over steps {
 *		rejected=false
 *		nActive=0
 *		loop over tracks in trackVector {
 *			if(rejected && !first) break;
 *			beginTrack()
 *			(discard remaining secondaries, if any)
 *			process track to stepTime or killed
 *			if(track is active) ++nActive
 *			if(keepSecondaries)
 *				append secondaries to trackVector
 *			else
 *				discard secondaries
 *		}
 *		if(!rejected || first) collectiveStep()
 *		if(rejected && !first) {
 *			restore saved tracks and stepTime
 *		} else {
 *			save tracks and stepTime
 *			stepTime += deltaT
 *		}
 *		first=false
 *	} while(nActive > 0)
 *	endCollectiveTracking()
 *
 *	Note: if multiple BLCollectiveComputation instances are registered
 *	with BLRUnManager, the callbacks to them will occur in the order
 *	in which they were registered (normally the order in which their
 *	commands appear in the input file).
 *
 **/
class BLCollectiveComputation {
public:
	BLCollectiveComputation() { };
	virtual ~BLCollectiveComputation() { }

	/// beginCollectiveTracking() is called after the track vector is 
	/// constructed but before trackng begins.
	/// Can be used to allocate internal arrays or create NTuple-s.
	/// Note that additional tracks can be added to v as tracking proceeds;
	/// no tracks are ever removed (but they will change status).
	/// The vector v will not change during the run.
	virtual void beginCollectiveTracking(std::vector<BLTrackData>& v) { }

	/// beginTrack() is called just before the track v[index] is
	/// tracked for the next step. Note that the user class can be
	/// derived from G4UserTrackingAction and registered to be called
	/// at the start and end of tracking for each track. The difference
	/// is this routine relates the track to the entry in v[], and
	/// also: BLRunManager::processOneTrack() can be called here, but
	/// not in PreUserTrackingAction().
	virtual void beginTrack(std::vector<BLTrackData>& v, int index) { }

	/// collectiveStep() is called after each collective step
	/// to perform the computation.
	/// It can store data for later use, and can modify the tracks.
	/// Note that collectiveStep() is not called until after the
	/// first step is taken, so you might want to call your
	/// collectiveStep() function at the end of beginCollectiveTracking()
	/// to store the state of tracks at the earliest time.
	/// Note: collectiveStep() can call runmanager->rejectCollectiveStep(),
	/// in which case the tracking since the previous return from
	/// collectiveStep() will be discarded and re-done; for this to make
	/// sense, runmanager->setDeltaT() should also be called with a
	/// smaller value -- this is one way to adapt deltaT dynamically.
	virtual void collectiveStep(std::vector<BLTrackData>& v) = 0;

	/// endCollectiveTracking() is called after tracking is complete.
	/// Can be used to delete internal arrays and print a summary.
	virtual void endCollectiveTracking(std::vector<BLTrackData>& v) { }
};

#endif // BLCOLLECTIVECOMPUTATION_H
