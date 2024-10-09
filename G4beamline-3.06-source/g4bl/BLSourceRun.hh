//	BLSourceRun.hh
#ifndef BLSOURCERUN_HH
#define BLSOURCERUN_HH

#include "G4TrackingManager.hh"
#include "G4VUserDetectorConstruction.hh"

/**	class BLSourceRun is a base class for a particle source which needs a
 *	special Geant4 run to generate events. For instance, rdecaysource.
 *
 *	This class implements all of the BLManager callbacks with do-nothing
 *	routines. In operation during a source run, the BLManager calls ONLY
 *	callbacks of this class, so whichever routines are implemented in a
 *	derived class will be called.
 *
 *	By calling BLManager::registerSourceRun(this), a derived class can
 *	arrange the BLManager to perform a Geant4 run in which only the 
 *	callbacks from the class itself are called. This happens after the
 *	geometry and physics are setup, and after the Tune and Reference
 *	particles are tracked (if any), but before the "post-Reference"
 *	(pre-beam tracking) callbacks are called. This includes the Beam
 *	callbacks to generate events. The result is that the derived class
 *	is in complete control of the source run -- essentially the only
 *	things setup by G4beamline are the geometry and the physics processes.
 *
 *	Note the derived-class must normally register itself as a beam; its
 *	nextBeamEvent() must must read-back the events it generated in the
 *	source run. An alternate is to issue a beam command that reads a
 *	file written during the source run.
 *
 *	Note that UserZSteppingAction() cannot be used (because the BLManager
 *	UserSteppingAction does not get that far in SOURCE mode).
 *
 *	Technical note: this class does not inherit from all the BLManager
 *	callback classes, it simply implements their virtual functions.
 *	As this file must be #included by BLManager.hh, this approach avoids
 *	problems with recursive includes.
 **/
class BLSourceRun {
protected:
	G4TrackingManager *fpTrackingManager;
public:
	/// Constructor
	BLSourceRun() { fpTrackingManager = 0; }

	/// Derived classes prepare for their SOURCE run.
	/// Called just before the beamOn() call for the SOURCE run.
	virtual void prepare() { }

	/// Derived classes shut down.
	/// Called just after the SOURCE run has ended.
	virtual void shutDown() { }

	/// Callback functions for RunAction
	virtual void BeginOfRunAction(const G4Run *run) { }
	virtual void EndOfRunAction(const G4Run *run) { }

	/// Callback functions for EventAction
	virtual void BeginOfEventAction(const G4Event* event) { }
	virtual void EndOfEventAction(const G4Event* event) { }

	/// Callback functions for TrackingAction
	void SetTrackingManagerPointer(G4TrackingManager *p) 
		{fpTrackingManager=p;}
	virtual void PreUserTrackingAction(const G4Track *track) { }
	virtual void PostUserTrackingAction(const G4Track *track) { }

	/// Callback functions for SteppingAction
	virtual void UserSteppingAction(const G4Step *step) { }

	/// Callback functions for StackingAction
	virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track*) 
		{ return fUrgent; } // keeps all tracks by default
	virtual void NewStage() { }
	virtual void PrepareNewEvent() { }

	/// callbacks for generating the beam during the source run.
	virtual bool nextSourceEvent(G4Event *event) { return false; }
};

#endif // BLSOURCERUN_HH
