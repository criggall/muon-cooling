//	BLRunManager.hh
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

#ifndef BLRUNMANAGER_H
#define BLRUNMANAGER_H

#include <vector>
#include <map>
#include <setjmp.h>

#include "G4RunManager.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4StackManager.hh"
#include "G4SteppingManager.hh"
#include "G4Navigator.hh"
#include "G4PropagatorInField.hh"

#include "BLAssert.hh"
#include "BLCollectiveComputation.hh"
#include "BLCallback.hh"


/**	class BLRunManager -- enhanced version of G4RunManager.
 *
 *	There are three basic modes of BLRunManager: normal, collective,
 *	and enhanced.
 *
 *	Normal mode:
 *	This class inherents from G4RunManager, and can thus be used
 *	just like it, for standard control of the Geant4 simulation
 *	(via function BeamOn()).
 *
 *	Collective mode:
 *
 *	In collective mode a vector<G4TrackData*> is generated from the
 *	input file's beam commands, and then all tracks are stepped for
 *	one time step after which the BLCollectiveComputation::collectiveStep()
 *	is called; this is continued until all tracks are killed. Any
 *	secondaries generated can either be kept or discarded (see
 *	setCollectiveMode() below). The TrackID-s of primary tracks are
 *	preserved, and TrackIDs for secondary tracks are assigned sequentially
 *	from 1001, ignoring EventID.
 *
 *	Collective mode is selected by calling setCollectiveMode(), and
 *	then behaving externally just like normal mode (i.e. just call
 *	BeamOn()).
 *
 *	Enhanced mode:
 *	This class also implements functions that permit enhanced control
 *	of the simulation -- in particular there is a processOneTrack()
 *	function that simply processes a single track.
 *
 *	NOTE: the word "track" is ambiguous, and is commonly used as both
 *	a noun and a verb. Here it is a noun only (referring to a G4Track),
 *	and the verb is "process". Processing a track consists of stepping
 *	it through the simulation until it is killed or suspended, while
 *	holding all secondaries it generates. The secondaries are themselves
 *	tracks, and must be dealt with before another call to processOneTrack()
 *	(see below).
 *
 *	Enhanced run control: All of the usual G4beamline initialization
 *	must occur, so the following code is best put into a BLCallback
 *	function registered to replace the main loop -- that means it gets
 *	called after the post-Reference callbacks are called, and that after
 *	all replace-main-loop callbacks are called the program exits (thus
 *	avoiding the main loop).
 *
 *		BLRunManager *runmgr = BLRunManager::getObject();
 *		runmgr->beginRun(int runid);
 *		--- loop over events ---
 *			runmgr->beginEvent(int evid);
 *			--- loop over tracks ---
 *				--- get a track
 *				runmgr->processOneTrack(track);
 *				--- handle secondaries (see below)
 *			runmgr->endEvent();
 *		runmgr->endRun();
 *		// write out NTuples, etc. (if appropriate)
 *
 *	There can be a loop over runs, if necessary.
 *
 *	Note that these enhanced run control routines do not create or handle
 *	trajectories, and ignore them if present. Hence no visualization
 *	of tracks is possible.
 *
 *	Note that secondaries are tracks, but they have not yet been assigned
 *	either a trackid or a parentid. Those should be assigned before the
 *	track is deferred, as all deferred tracks should have those id-s.
 *
 *	Handling secondaries -- do one of:
 *		runmgr->discardAllSecondaries();
 *	or
 *		runmgr->processAllSecondariesAndDeferredTracks();
 *	or
 *		runmgr->deferAllSecondaries();
 *		// followed sometime later by:
 *		runmgr->processAllDeferredTracksAndTheirSecondaries();
 *	or
 *		while((track=runmgr->popOneSecondary()) != 0) {
 *			// NOTE: you cannot call processOneTrack() here!
 *			// -- that would discard  any remaining secondaries.
 *			if(--- condition ---)
 *				runmgr->deferOneTrack(track);
 *			else
 *				delete track;
 *		}
 *		// followed sometime later by:
 *		runmgr->processAllDeferredTracksAndTheirSecondaries();
 **/
class BLRunManager : public G4RunManager, public BLCallback {
	static BLRunManager *singleton;
	G4TrackingManager *trackManager;
	G4StackManager *stackManager;
	jmp_buf jmpBuf;
	bool validJmpBuf;
	std::vector<BLTrackData> trackVector;
	std::vector<BLCollectiveComputation *> computeVector;
	bool collectiveMode;
	bool keepSecondaries;
	G4double stepTime; // collective mode only
	G4double deltaT;       // collective mode only
	G4int currentTrackIndex; // collective mode only
	bool rejected;		// collective mode only
	int nextSecondaryTrackID; // collective mode only
	void appendTrack(G4Track *track); // collective mode only
public:
	/// returns the (singleton) BLRunManager.
	/// Note that G4RunManager::GetRunManager() returns a pointr to
	/// the base class, not this class.
	static BLRunManager *getObject();

	/// Constructor
	BLRunManager();

	/// Destructor
	virtual ~BLRunManager();

	/// BeamOn() will track beam events, in normal or collective mode.
	/// See setCollectiveMode() below.
	virtual void BeamOn(int nEvents, const char *macroFile=0,
							G4int n_select=-1);

	/// DoEventLoop() is REPLACED, so abondonCurrentEvent() works.
	void DoEventLoop(G4int n_event,const char* macroFile,G4int n_select);

	/// abandonCurrentEvent() will abandon the current event by performing
	/// a longjmp into DoEventLoop(). Not guaranteed to work, but probably
	/// does, even if called from BLAlarm::signal().
	/// NOTE: If the current event can be abandoned, this function never
	/// returns; if the event cannot be abandoned, it returns immediately.
	void abandonCurrentEvent();

	/// callback() from BLCallback.
	virtual void callback(int type);

	// Collective run control:
	// NONE of these routines are called in normal or enhanced mode.

	/// beamOnCollective() processes an entire run in collective mode.
	/// In collective mode, a vector of BLTrackData is constructed, and
	/// all tracks in it are stepped one time step at a time. This is 
	/// essentially an inversion of the usual loop over steps inside
	/// loops over tracks and events.
	/// All secondaries are immediately added to the end of trackVector,
	/// if keepSecondaries is true.
	/// The Tune and Reference particles are always tracked in normal mode.
	void beamOnCollective(int nEvents);

	/// registerCollectiveComputation() will register a computation to
	/// be called after every time step in collective mode. Multiple
	/// computations are called in the order they are registered.
	/// Must be called before beamOnCollective() or steps may be missed
	/// (usually called in the command() function of a BLCommand, so the
	/// order of commands in the input file is preserved).
	void registerCollectiveComputation(BLCollectiveComputation *compute)
		{ computeVector.push_back(compute); }

	/// getTrackVector() returns the trackVector.
	/// Note: the track vector and its tracks can be modified inside the
	/// collective computation (takes effect immediately).
	/// As BLRunManager is a singleton, the trackVector never changes.
	/// No track is ever removed from the track vector, but they can
	/// have their status set to fStopAndKill which effectively does that.
	std::vector<BLTrackData> &getTrackVector() { return trackVector; }

	/// getCurrentTrackIndex() returns the index in the track vector
	/// of the track currently being tracked. Returns -1 if not in
	/// collective mode, or if tracking is not in progress.
	int getCurrentTrackIndex() { return currentTrackIndex; }

	/// getCollectiveMode() returns true if in collective mode.
	bool getCollectiveMode() { return collectiveMode; }

	/// setCollectiveMode() arranges for calls to BeamOn() to call
	/// beamOnCollective() rather than the usual G4RunManager function.
	/// This is for BEAM only; TUNE and REFERENCE particles are always
	/// tracked normally via G4RunManager::BeamOn().
	void setCollectiveMode(bool flag=true) { collectiveMode=flag; }

	/// getStepTime() will get the ste time.
	/// NOTE: inside collectiveStep() this function returns the time
	/// of the current step; outside collectiveStep() it returns the
	/// next step time (i.e. the time at which tracks should be suspended).
	G4double getStepTime() { return stepTime; }

	/// setStepTime() will set the time for the next step in time.
	/// If called from within collectiveStep(), it may confuse other
	/// classes collectiveStep()-s. If called from beginCollectiveTracking()
	/// it sets the time of the first step.
	void setStepTime(G4double t) { stepTime = t; }

	/// getDeltaT() will get the time interval between time steps
	G4double getDeltaT() { return deltaT; }

	/// setDeltaT() will set the time interval between time steps.
	/// All tracks will be stepped to the same value of global time.
	/// Initiall value is 1.0*ns. deltaT <= 0.0 is an error.
	/// Can be called at any time; if called within collectiveStep()
	/// it will affect the immediately following time step.
	void setDeltaT(G4double dt) { BLAssert(dt>0.0); deltaT = dt; 
	}

	/// rejectCollectiveStep() will cause the run manager to abandon
	/// tracking (if in progress) and to discard all of the current
	/// tracks, replacing the vector with a copy saved immediately
	/// after the previous call to collectiveStep() returned.
	/// Thus the current collective step will be repeated as if
	/// no tracking had been performed since the previous return from
	/// collectiveStep(). For this to make sense, setDeltaT() 
	/// should have been called with a smaller value
	/// -- this is intended for a collective algorithm to dynamically
	/// adapt the value of deltaT.
	/// Note: the first step is to bring all beam tracks to a common
	/// value of global time, and rejecting this step is ignored.
	void rejectCollectiveStep() { rejected = true; }



	// Enhanced run control:
	// NONE of these routines are called in normal run-control mode.
	// NONE of these routines are called in collective run-control mode.
	// In most cases these routines should be called from a callback()
	// function registered to replace the main loop (exception: calling
	// processOneTrack() as "pre-tracking" for a collective computation).

	/// beginRun() begins a run in enhanced run-control mode.
	void beginRun(int runid=0);

	/// endRun() ends a run in enhanced run-control mode.
	void endRun();

	/// beginEvent() begins an event in enhanced run-control mode.
	void beginEvent(int evid=0);

	/// endEvent() ends an event in enhanced run-control mode.
	void endEvent();

	/// getNextBeamEventAndTrack() will get the next beam event and track.
	/// returns true if one is returned, false if none are left.
	/// Neithe event actions nor track actions are performed.
	/// Both *pevent and *ptrack must be deleted.
	bool getNextBeamEventAndTrack(G4Event **pevent, G4Track **ptrack);

	/// processOneTrack() tracks a single track until it is suspended or 
	/// killed.
	/// Note that before tracking begins it deletes any held secondaries,
	/// so if secondaries are not to be discarded they must be deferred
	/// or processed before the next call.
	/// You may want to call BLManager::setState(SPECIAL) to prevent
	/// this track from being entered into NTuples for virtualdetectors
	/// (etc.) encountered during processing.
	/// This function MUST NOT BE CALLED if any track is being processed
	/// via any means (e.g. the normal BeamOn() is executing).
	void processOneTrack(G4Track *track);

	/// discardAllSecondaries() will discard all secondaries.
 	void discardAllSecondaries();

	/// deferAllSecondaries() will defer all secondaries until
	/// deferred tracks are processed. Returns # tracks deferred.
	int deferAllSecondaries(int secondaryid=10000, int parentid=-1);

	/// deferOneTrack() will defer a single track until deferred tracks
	/// are processed.
	void deferOneTrack(G4Track *track);

	/// processAllSecondariesAndDeferredTracks() will defer all secondaries
	/// and then process all deferred tracks and their secondaries.
	/// returns # tracks processed.
	int processAllSecondariesAndDeferredTracks(int secondaryid=10000,
							int parentid=-1) {
		int first = secondaryid;
		secondaryid += deferAllSecondaries(secondaryid,parentid);
		secondaryid +=
		       processAllDeferredTracksAndTheirSecondaries(secondaryid);
		return secondaryid-first;
	}

	/// processAllDeferredTracksAndTheirSecondaries() will process all
	/// deferred tracks, including processing their secondaries.
	/// Returns # tracks processed.
	int processAllDeferredTracksAndTheirSecondaries(int trackid=10000);

	/// popOneSecondary() will return a pointer to one secondary, removing
	/// it from the list of secondaries; returns NULL if no more.
	/// Intended to be called in a loop immediately after a call to
	/// processOneTrack(). The returned pointer should eventually be
	/// deleted. Order is LIFO.
	G4Track *popOneSecondary();

	/// popOneDeferredTrack() will return a pointer to one deferred
	/// track, removing it from the list of deferred tracks; returns
	/// NULL if no more. The returned pointer should eventually be
	/// deleted. Note that a loop calling this function can defer other
	/// tracks or secondaries, which simply extends the loop.
	/// Order is LIFO.
	G4Track *popOneDeferredTrack();

	// protected functions made public (for testing and special purposes):
	void RunInitialization() { G4RunManager::RunInitialization(); }
	void RunTermination() { G4RunManager::RunTermination(); }
	G4EventManager *getEventManager() { return eventManager; }
	G4TrackingManager *getTrackingManager() { return trackManager; }
	G4StackManager *getStackManager() { return stackManager; }
	void setCurrentEvent(G4Event *ev) { currentEvent = ev; }
	G4bool getRunAborted() { return runAborted; }
};

#endif // BLRUNMANAGER_H
