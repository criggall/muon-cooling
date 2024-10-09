/*	BLRunManager.cc - based on G4RunManager.cc of Geant4 9.0. */

#include "G4StateManager.hh"
#include "G4TransportationManager.hh"
#include "G4ProcessTable.hh"
#include "G4VProcess.hh"
#include "G4UImanager.hh"
#include "G4Timer.hh"

#include "BLAssert.hh"
#include "BLRunManager.hh"
#include "BLManager.hh"
#include "BLMPI.hh"

/**	class TStepLimiter limits steps based on global time.
 *
 *	The steps in time are not exact, as Geant4 inherently limits steps
 *	in space, not time. This class assumes the velocity of the track
 *	does not change during the step. This is approximate, but for a
 *	50 MeV/c mu+ in a 20 MegaVolt/meter electric field, it is accurate
 *	to 4 ps for steps of 1 ns; for a 200 MeV/c mu+ it is accurate to 0.5
 *	ps for steps of 1 ns. For better timing accuracy, use smaller steps.
 *	Obviously the accuracy will be better for tracks with beta near 1.
 *	Due to the algorithm used, these timing errors do NOT accumulate.
 **/
class TStepLimiter : public G4VProcess {
	static G4TrackStatus prevStatus;
	class ParticleChange : public G4VParticleChange {
	public:
		ParticleChange() : G4VParticleChange() { }
	};
	ParticleChange change;
	BLRunManager *runManager;
public:
	TStepLimiter() : G4VProcess("TStepLimiter",fUserDefined),
		change() { runManager = BLRunManager::getObject(); }

	virtual ~TStepLimiter() { }

	G4VProcess *clone() { return new TStepLimiter(); }

	/// getPrevStatus() returns the previous status of the most recently
	/// suspended track.
	static G4TrackStatus getPrevStatus() { return prevStatus; }

	// G4VProcess functions from G4VProcess
	virtual G4bool IsApplicable(const G4ParticleDefinition&)
		{ return true; }

	virtual G4double AlongStepGetPhysicalInteractionLength(
			const G4Track& track, G4double  previousStepSize,
			G4double  currentMinimumStep, G4double& proposedSafety,
			G4GPILSelection* selection)
		{ return -1.0; }

	virtual G4VParticleChange* AlongStepDoIt(const G4Track& track,
			const G4Step& stepData)
		{ change.Initialize(track);
		  return &change;
		}

	virtual G4double AtRestGetPhysicalInteractionLength(
			const G4Track& track, G4ForceCondition* condition)
		{ *condition = NotForced;
		  G4double remaining = runManager->getStepTime() -
		  					track.GetGlobalTime();
		  if(remaining <= 0.0) remaining = DBL_MIN;
		  return remaining;
		}

	virtual G4double PostStepGetPhysicalInteractionLength(
			const G4Track& track, G4double   previousStepSize,
			G4ForceCondition* condition)
		{ *condition = NotForced;
		  G4double remaining = runManager->getStepTime() -
		  					track.GetGlobalTime();
		  if(remaining <= 0.0) remaining = DBL_MIN;
		  return remaining * track.GetVelocity();
		}

	virtual G4VParticleChange *PostStepDoIt(const G4Track& track,
							const G4Step& aStep)
		{ prevStatus = fAlive;
		  change.Initialize(track);
		  change.ProposeTrackStatus(fSuspend);
		  return &change;
		}

	virtual G4VParticleChange* AtRestDoIt(const G4Track& track,
							const G4Step&  aStep)
		{ prevStatus = fStopButAlive;
		  change.Initialize(track);
		  change.ProposeTrackStatus(fSuspend);
		  return &change;
		}
};
G4TrackStatus TStepLimiter::prevStatus = fPostponeToNextEvent;


BLRunManager *BLRunManager::singleton = 0;

BLRunManager *BLRunManager::getObject() 
{
	if(!singleton) {
		new BLRunManager();
		// startup sequencing prevents registering callback, so it
		// is hard-coded in BLManager::handleCallbacks(). This is due
		// to the BLManager constructor calling this function.
	}

	return singleton; 
}

BLRunManager::BLRunManager() : G4RunManager(), trackVector(), computeVector()
{
	BLAssert(singleton==0);
	singleton = this;

	trackManager = 0;
	stackManager = 0;
	validJmpBuf = false;
	collectiveMode = false;
	keepSecondaries = true;
	stepTime = 0.0;
	deltaT = 1.0*ns;
	currentTrackIndex = -1;
	rejected = false;
	nextSecondaryTrackID = -9999;
}

BLRunManager::~BLRunManager()
{
	trackVector.clear();
	computeVector.clear();
}

void BLRunManager::BeamOn(int nEvents, const char *macroFile, G4int n_select) 
{
	if(collectiveMode && BLManager::getObject()->getState() == BEAM)
		beamOnCollective(nEvents);
	else
		G4RunManager::BeamOn(nEvents,macroFile,n_select);
}

void BLRunManager::DoEventLoop(G4int n_event,const char* macroFile,G4int n_select)
{
	if(verboseLevel>0) { timer->Start(); }

	G4String msg;
	if(macroFile!=0) { 
		if(n_select<0) n_select = n_event;
		msg = "/control/execute ";
		msg += macroFile;
	} else {
		n_select = -1;
	}

	int local_n_select = n_select; // argument can be clobbered by longjmp

	// Event loop
	volatile G4int i_event = 0;
	if(setjmp(jmpBuf) != 0) ++i_event;
	validJmpBuf = true;
	while(i_event<n_event) {
		currentEvent = GenerateEvent(i_event);
		eventManager->ProcessOneEvent(currentEvent);
		AnalyzeEvent(currentEvent);
		UpdateScoring();
		if(i_event<local_n_select)
			G4UImanager::GetUIpointer()->ApplyCommand(msg);
		StackPreviousEvent(currentEvent);
		currentEvent = 0;
		if(runAborted) break;
		++i_event;
	}
	validJmpBuf = false;

	if(verboseLevel>0) {
		timer->Stop();
		G4cout << "Run terminated." << G4endl;
		G4cout << "Run Summary" << G4endl;
		if(runAborted) {
			G4cout << "  Run Aborted after " << i_event + 1 <<
						" events processed." << G4endl;
		} else {
			G4cout << "  Number of events processed : " << n_event
								<< G4endl; 
		}
		G4cout << "  "  << *timer << G4endl;
	}
}

void BLRunManager::abandonCurrentEvent()
{
	if(validJmpBuf) {
		G4Exception("BLRunManager","Event Time Limit Alarm",
							EventMustBeAborted,"");
		G4StateManager::GetStateManager()->
						SetNewState(G4State_GeomClosed);
		longjmp(jmpBuf,1);
	}
}

void BLRunManager::callback(int type)
{
	// Only perform actions before Tune particle
	if(type != 0) return;

	// in collectiveMode, add a TStepLimiter to every particle
	if(collectiveMode) {
		if(BLMPI::isMPI())
			G4Exception("BLRunManager",
				"Collective mode is incompatible with MPI",
				FatalException,"");
		printf("BLRunManager: adding TStepLimiter processes to all particles\n");
	    	G4ParticleTable::G4PTblDicIterator *myParticleIterator
			= G4ParticleTable::GetParticleTable()->GetIterator();
	    	myParticleIterator->reset();
	    	while((*myParticleIterator)()) {
			G4ParticleDefinition *pd = 
					myParticleIterator->value();
			if(pd->IsShortLived()) continue;
			G4ProcessManager *pmgr = pd->GetProcessManager();
			if(!pmgr) continue;
			TStepLimiter *tsl = new TStepLimiter();
			pmgr->AddProcess(tsl,ordDefault,-1,ordDefault);
	    	}
		stepTime = DBL_MAX;
	}
}

void BLRunManager::beamOnCollective(int nEvents)
{
	BLManager *manager = BLManager::getObject();
	G4TransportationManager *transportManager =
			G4TransportationManager::GetTransportationManager();
	G4StateManager* stateManager = G4StateManager::GetStateManager();
	numberOfEventToBeProcessed = nEvents;

	if(!ConfirmBeamOnCondition()) {
		G4Exception("BLRunManager","Cannot run beam",FatalException,"");
	}
	printf("================== In Collective Mode ==================\n");
	ConstructScoringWorlds();
	RunInitialization();
	G4ApplicationState currentState = stateManager->GetCurrentState();
	if(currentState != G4State_GeomClosed) {
		G4Exception("BLRunManager","Geometry not closed",FatalException,
									"");
	}

	manager->setState(BEAM);
	stateManager->SetNewState(G4State_EventProc);
	trackManager = 0; // unused in collective mode
	stackManager = eventManager->GetStackManager();
	runAborted = false;

	// the per-event time limit makes no sense in collective mode
	manager->setEventTimeLimit(-1); // infinite

	// event loop to create trackVector
	G4PrimaryTransformer *transformer = new G4PrimaryTransformer();
	G4Navigator* navigator = transportManager->GetNavigatorForTracking();
	stepTime = -DBL_MAX;
	int nev;
	for(nev=0; nev<nEvents; ++nev) {
		currentEvent = GenerateEvent(nev);
		if(!currentEvent || runAborted) break;
		stackManager->PrepareNewEvent();
		manager->BeginOfEventAction(currentEvent);
		G4TrackVector *tv=transformer->GimmePrimaries(currentEvent,0);
		int evId=currentEvent->GetEventID();
		if(evId <= 0) {
			evId = nev;
			currentEvent->SetEventID(evId);
		}
		nextSecondaryTrackID = manager->getNextSecondaryTrackID();
		for(unsigned j=0; j<tv->size(); ++j) {
			G4Track *track = (*tv)[j];
			BLAssert(j==0); // cannot handle multiple primaries 
			track->SetTrackID(manager->getPrimaryTrackID());
			appendTrack(track);
			G4double t = track->GetGlobalTime();
			if(t > stepTime) stepTime = t;
		}
		if(tv->size() == 0) delete currentEvent;
		tv->clear();
		if(runAborted) break;
		transportManager->SetNavigatorForTracking(navigator);
		manager->setEventID(evId+1);
	}

	printf("=========== Collective: %d Events, %ld Tracks ==============\n",
				nev,(long)trackVector.size());

	for(unsigned i=0; i<computeVector.size(); ++i)
		computeVector[i]->beginCollectiveTracking(trackVector);

	// Step loop
	runAborted = false;
	int nActive;
	bool first=true;
	G4double saveStepTime = stepTime;
	std::vector<G4Track*> saveVector;
	do {
		if(runAborted) goto RunAborted;
		nActive = 0;
		rejected = false;
		// note: secondaries can be appended to trackVector during loop
		for(unsigned i=0; i<trackVector.size(); ++i) {
			if(rejected && !first) break;
			BLTrackData *td = &trackVector[i];
			G4Track *track = td->track;
			G4TrackStatus trackStatus = track->GetTrackStatus();
//printf("Collective Track loop: ev=%d trk=%d status=%d fAlive=%d fSuspend=%d\n",td->event->GetEventID(),track->GetTrackID(),(int)trackStatus, fAlive,fSuspend);
			if(trackStatus != fAlive && 
						trackStatus != fStopButAlive)
				continue;
			++nActive;
			currentEvent = td->event;
			manager->setEventID(td->event->GetEventID());
			currentTrackIndex = i;
			for(unsigned j=0; j<computeVector.size(); ++j)
				computeVector[j]->beginTrack(trackVector,i);
			if(runAborted) goto RunAborted;
			processOneTrack(track);
			currentTrackIndex = -1;
			if(td->event->IsAborted())
				track->SetTrackStatus(fKillTrackAndSecondaries);
			if(runAborted) goto RunAborted;
			trackStatus = track->GetTrackStatus();
			if(trackStatus == fSuspend)
			   	track->SetTrackStatus(TStepLimiter::getPrevStatus());
			G4TrackVector *tv = trackManager->GimmeSecondaries();
			if(keepSecondaries && trackStatus != fKillTrackAndSecondaries) {
				for(unsigned j=0; j<tv->size(); ++j)
					appendTrack((*tv)[j]);
				tv->clear();
			} else {
				for(unsigned j=0; j<tv->size(); ++j)
					delete (*tv)[j];
				tv->clear();
			}
		}
		if(runAborted) goto RunAborted;
		// don't call collectiveStep() if rejected during tracking
		if(!rejected || first) {
			for(unsigned i=0; i<computeVector.size(); ++i)
				computeVector[i]->collectiveStep(trackVector);
		}
		if(rejected && !first) {
			// discard current tracks, and restore saved ones
			unsigned n=saveVector.size();
			// first, discard any added secondaries
			while(trackVector.size() > n) {
				delete trackVector.back().track;
				delete trackVector.back().event;
				trackVector.pop_back();
			}
			BLAssert(trackVector.size() == saveVector.size());
			for(unsigned i=0; i<n; ++i) {
				delete trackVector[i].track;
				trackVector[i].track = saveVector[i];
			}
			stepTime = saveStepTime;
		} else {
			// delete saved tracks, and save current ones
			for(unsigned i=0; i<saveVector.size(); ++i)
				delete saveVector[i];
			saveVector.clear();
			for(unsigned i=0; i<trackVector.size(); ++i)
				saveVector.push_back(
					new G4Track(*trackVector[i].track));
			saveStepTime = stepTime;
			stepTime += deltaT;
		}
		first = false;
	} while(nActive > 0);

RunAborted:
	for(unsigned i=0; i<computeVector.size(); ++i)
		computeVector[i]->endCollectiveTracking(trackVector);

	RunTermination();
	manager->setState(IDLE);
}

void BLRunManager::appendTrack(G4Track *track)
{
	G4int evId = currentEvent->GetEventID();
	if(track->GetTrackID() <= 0) {
		track->SetTrackID(nextSecondaryTrackID++);
	}

	G4Event *tmp = currentEvent;
	currentEvent = new G4Event(currentEvent->GetEventID());

	trackVector.push_back(BLTrackData(currentEvent,track));

	currentEvent = tmp;
}

void BLRunManager::beginRun(int runid)
{
	BLManager::getObject()->setState(BEAM);
	ConstructScoringWorlds();
	RunInitialization();
	SetRunIDCounter(runid);
	runAborted = false;
}

void BLRunManager::endRun()
{
	if(BLManager::getObject()->getState() == IDLE) return;
	RunTermination();
	BLManager::getObject()->setState(IDLE);
}

void BLRunManager::beginEvent(int evid)
{
	if(currentEvent) delete currentEvent;
	currentEvent = new G4Event(evid);
	trackManager = eventManager->GetTrackingManager();
	stackManager = eventManager->GetStackManager();
	BLManager::getObject()->BeginOfEventAction(currentEvent);
}

void BLRunManager::endEvent()
{
	BLManager::getObject()->EndOfEventAction(currentEvent);
	if(currentEvent) delete currentEvent;
	currentEvent = 0;
	trackManager = 0;
	stackManager = 0;
}

bool BLRunManager::getNextBeamEventAndTrack(G4Event **pevent, G4Track **ptrack)
{
	static bool more=true;
	static int nev=0;
	static BLManager *manager=0;
	static G4PrimaryTransformer *transformer=0;
	if(!manager) {
		manager = BLManager::getObject();
		transformer = new G4PrimaryTransformer();
	}

	if(!more) return false;

	manager->setState(BEAM);
	currentEvent = GenerateEvent(nev);
	if(!currentEvent || runAborted) {
		more = false;
		return false;
	}
	*pevent = currentEvent;
	G4TrackVector *tv=transformer->GimmePrimaries(currentEvent,0);
	int evId=currentEvent->GetEventID();
	if(evId <= 0) {
		evId = nev;
		currentEvent->SetEventID(evId);
	}
	for(unsigned j=0; j<tv->size(); ++j) {
		G4Track *track = (*tv)[j];
		BLAssert(j==0); // cannot handle multiple primaries 
		track->SetTrackID(manager->getPrimaryTrackID());
		track->SetParentID(manager->getPrimaryParentID());
		*ptrack = track;
	}
	tv->clear();
	if(runAborted) {
		more = false;
		return false;
	}
	++nev;
	manager->setState(SPECIAL);
	return true;
}

void BLRunManager::processOneTrack(G4Track *track)
{
	if(!trackManager || !stackManager) {
		trackManager = eventManager->GetTrackingManager();
		stackManager = eventManager->GetStackManager();
	}
	BLAssert(track != 0);
	trackManager->ProcessOneTrack(track);
}

void BLRunManager::discardAllSecondaries()
{
	G4TrackVector *secondaries = trackManager->GimmeSecondaries();
	for(unsigned j=0; j<secondaries->size(); ++j)
		delete (*secondaries)[j];
	secondaries->clear();
}

int BLRunManager::deferAllSecondaries(int secondaryid, int parentid)
{
	int n=0;
	G4TrackVector *secondaries = trackManager->GimmeSecondaries();
	for(unsigned j=0; j<secondaries->size(); ++j) {
		G4Track *track = (*secondaries)[j];
		track->SetTrackID(secondaryid++);
		track->SetParentID(parentid); 
		deferOneTrack(track);
		++n;
	}
	secondaries->clear();

	return n;
}

void BLRunManager::deferOneTrack(G4Track *track)
{
	stackManager->PushOneTrack(track);
}

int BLRunManager::processAllDeferredTracksAndTheirSecondaries(int secondaryid)
{
	int first=secondaryid;
	for(;;) {
		G4VTrajectory *trajectory=0;
		G4Track *track = stackManager->PopNextTrack(&trajectory);
		if(!track) break;
		processOneTrack(track);
		int parentid = track->GetTrackID();
		delete track;
		secondaryid += deferAllSecondaries(secondaryid,parentid);
	}

	return secondaryid-first;
}

G4Track *BLRunManager::popOneSecondary()
{
	G4TrackVector *secondaries = trackManager->GimmeSecondaries();
	int i = secondaries->size();
	if(i == 0) return 0;
	G4Track *track = secondaries->back();
	secondaries->pop_back();
	return track;
}

G4Track *BLRunManager::popOneDeferredTrack()
{
	G4VTrajectory *trajectory=0;
	return stackManager->PopNextTrack(&trajectory);
}
