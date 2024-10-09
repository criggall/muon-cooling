//	BLManager.cc
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

#ifdef G4BL_GSL
#include <gsl/gsl_errno.h>
#endif

#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"
#include "G4GeometryManager.hh"
#include "G4StateManager.hh"
#include "G4UserTrackingAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserStackingAction.hh"
#include "G4VProcess.hh"
#include "G4VParticleChange.hh"
#include "G4ProcessTable.hh"
#include "G4ParticleTable.hh"
#include "G4HadronicProcessStore.hh"

#include "BLAssert.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLMPI.hh"
#include "BLRunManager.hh"
#include "BLBeam.hh"
#include "BLNTuple.hh"
#include "BLPhysics.hh"
#include "BLGroup.hh"
#include "BLGlobalField.hh"
#include "BLCommand.hh"
#include "BLCoordinates.hh"
#include "BLAlarm.hh"
#include "BLSignal.hh"
#include "BLTime.hh"
#include "BLZStep.hh"
#include "BLTrackInfo.hh"
#ifdef G4BL_GUI
#include "BLQt.h"
#endif

#ifdef G4BL_VISUAL
#include "BLVisManager.hh"
#endif

#include "mysnprintf.hh"

extern void g4bl_exit(int); // in g4beamline.cc

// Param definitions are all moved to the BLManager constructor, because
// other initializers use BLmanager.

/** class ZStepLimiter limits the ZStep
 **/
class ZStepLimiter : public G4VProcess {
	static G4double maxStep;
	class ParticleChange : public G4VParticleChange {
	public:
		ParticleChange() : G4VParticleChange() { }
	};
	ParticleChange change;
public:
	ZStepLimiter() : G4VProcess("ZStepLimiter",fUserDefined), change() {
	    G4ParticleTable::G4PTblDicIterator *myParticleIterator =
			G4ParticleTable::GetParticleTable()->GetIterator();
	    myParticleIterator->reset();
	    while((*myParticleIterator)()) {
		G4ParticleDefinition *pd = myParticleIterator->value();
		G4ProcessManager *pmgr = pd->GetProcessManager();
		if(!pmgr) {
			printf("ZStepLimiter: particle '%s' has no ProcessManager!\n",
				pd->GetParticleName().c_str());
			continue;
		}
		pmgr->AddProcess(this,-1,-1,ordDefault);
	    }
	}

	G4VProcess *clone() { return new ZStepLimiter(*this); }

	static void setMaxStep(G4double ms) { maxStep = ms; }

	virtual G4double PostStepGetPhysicalInteractionLength(
			const G4Track& track, G4double   previousStepSize, 
			G4ForceCondition* condition)
		{ *condition = NotForced;
		  double limit = maxStep-0.010*mm;
		  if(limit < 0.0) limit = DBL_MIN;
		  return limit;
		}
	G4VParticleChange* PostStepDoIt( const G4Track& track, 
      			const G4Step&  stepData)
		{ change.Initialize(track);
		  return &change;
		}
	virtual G4double AlongStepGetPhysicalInteractionLength(
			const G4Track& track, G4double  previousStepSize,
			G4double  currentMinimumStep, G4double& proposedSafety,
			G4GPILSelection* selection)
		{ return -1.0; }
	virtual G4VParticleChange* AlongStepDoIt( const G4Track& track,
			const G4Step& stepData)
		{ change.Initialize(track);
		  return &change;
		}
	virtual G4double AtRestGetPhysicalInteractionLength( 
			const G4Track& track, G4ForceCondition* condition)
		{ return -1.0; }
	virtual G4VParticleChange* AtRestDoIt( const G4Track& track,
			const G4Step& stepData)
		{ change.Initialize(track);
		  return &change;
		}
};
G4double ZStepLimiter::maxStep = DBL_MAX;


BLManager *BLManager::blManager = 0;
bool BLManager::initialized=false;

/*	User Action classes
 *	These derive from the G4User* classes, and simply call the
 *	corresponding function of BLManager.
 *
 *	These exist merely to give the Geant4 code something to delete
 *	at the end of the program; if BLManager were derived from all
 *	of these it would be deleted multiple times by Geant4 code.
 */
class BLManager_UserRunAction : public G4UserRunAction {
	BLManager *mgr;
public:
	BLManager_UserRunAction(BLManager *p) { mgr = p; }
	void BeginOfRunAction(const G4Run *run) {
		mgr->BeginOfRunAction(run);
	}
	void EndOfRunAction(const G4Run *run) {
		mgr->EndOfRunAction(run);
	}
};
class BLManager_UserEventAction : public G4UserEventAction {
	BLManager *mgr;
public:
	BLManager_UserEventAction(BLManager *p) { mgr = p; }
	void BeginOfEventAction(const G4Event* event) {
		mgr->BeginOfEventAction(event);
	}
	void EndOfEventAction(const G4Event* event) {
		mgr->EndOfEventAction(event);
	}
};
class BLManager_UserTrackingAction : public G4UserTrackingAction {
	BLManager *mgr;
public:
	BLManager_UserTrackingAction(BLManager *p) { mgr = p; }
	void PreUserTrackingAction(const G4Track *track) {
		mgr->setTrackingManager(fpTrackingManager);
		mgr->PreUserTrackingAction(track);
	}
	void PostUserTrackingAction(const G4Track *track) {
		mgr->PostUserTrackingAction(track);
	}
};
class BLManager_UserSteppingAction : public G4UserSteppingAction {
	BLManager *mgr;
public:
	BLManager_UserSteppingAction(BLManager *p) { mgr = p; }
	void UserSteppingAction(const G4Step *step) {
		mgr->setSteppingManager(fpSteppingManager);
		mgr->UserSteppingAction(step);
	}
};
class BLManager_PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
	BLManager *mgr;
public:
	BLManager_PrimaryGeneratorAction(BLManager *p) { mgr = p; }
	void GeneratePrimaries(G4Event *event) {
		mgr->GeneratePrimaries(event);
	}
};
class BLManager_UserStackingAction : public G4UserStackingAction {
	BLManager *mgr;
public:
	BLManager_UserStackingAction(BLManager *p) { mgr = p; }
	G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track *track) {
		return mgr->ClassifyNewTrack(track);
	}
	void NewStage() { mgr->NewStage(); }
	void PrepareNewEvent() { mgr->PrepareNewEvent(); }
};

BLManager::BLManager() : G4VUserDetectorConstruction(), 
		G4VExceptionHandler(), keepEventList(), skipEventList(),
		beamVector(), referenceVector(),
		runActionVector(), beamRunActionVector(),
		eventActionVector(), beamEventActionVector(),
		trackingActionVector(), preReferenceCallbackVector(),
		postReferenceCallbackVector(), postTrackingCallbackVector(),
		replaceMainLoopCallbackVector(), visualizationCallbackVector(),
		physicsCallbackVector(),
		allStepVector(), allStepMap(), tpStepMap(), rpStepMap(),
		tpStepVector(), rpStepVector(),
		beamStepMap(), beamStepVector(), verboseFormat(),
		tuneZStep(), referenceZStep(), beamZStep(), 
		stackingActionVector(), trackIDMap(), userCodeVector(),
		exceptionCount(), sourceRun(0)
{
	// Parameter definitions (variables are all unused):
	BLSetParam unused_1("histoFile","g4beamline",
				"Default (Root) NTuple output filename");
	BLSetParam unused_2("histoUpdate","0",
				"Output update interval (events)");
	BLSetParam unused_3("viewer","none",
				"Visualization driver selected (default=none)");
	BLSetParam unused_4("eventTimeLimit","30",
				"CPU Time Limit (sec)");
	BLSetParam unused_5("steppingVerbose","0",
				"Set nonzero to print each step");
	BLSetParam unused_6("steppingFormat", "N GLOBAL CL KE STEP VOL PROCESS",
				"Format for printing steps");
	BLSetParam unused_7("zTolerance","2.0",
				"Tolerance for Z steps (mm)");
	BLSetParam unused_8("wallClockLimit","-1",
			"Limit on wall clock time in seconds; -1 is infinite");

	if(blManager)
		G4Exception("BLManager","Object Already Exists",FatalException,
									"");
	blManager = this;

	state = IDLE;
	steppingVerbose = 0;

	// ensure ZStep vectors have bookends
	tuneZStep.push_back(ZStep(-DBL_MAX,0));
	tuneZStep.push_back(ZStep(DBL_MAX,0));
	referenceZStep.push_back(ZStep(-DBL_MAX,0));
	referenceZStep.push_back(ZStep(DBL_MAX,0));
	beamZStep.push_back(ZStep(-DBL_MAX,0));
	beamZStep.push_back(ZStep(DBL_MAX,0));
	currentZStep = &beamZStep;
	indexZStep = 1;
	zTolerance = Param.getDouble("zTolerance");
	prevZ = -DBL_MAX;

	// default verbose format
	verboseFormat.push_back(NSTEP);
	verboseFormat.push_back(GLOBAL);
	verboseFormat.push_back(CL);
	verboseFormat.push_back(KE);
	verboseFormat.push_back(STEP);
	verboseFormat.push_back(VOL);
	verboseFormat.push_back(PROCESS);
	verboseFormat.push_back(NEWLINE);

	runManager = 0;
	state = IDLE;
	beamIndex = 0;
	histoUpdate = 0;
	startProgram = BLTime::time();
	startRun = 0;
	startEvent = 0;
	eventTimeLimit = 30;
	physics = 0;
	worldPhysicalVolume = 0;
	eventID = 1;
	trackID=-9999;
	currentTrack = 0;
	allExceptions = false;
	fatalExceptions = 0;
	eventsAborted=0;
	stuckTracks = 0;
	warnings=0;
	prevEventID = 0;
	eventsProcessed = 0;
	endRun = false;
	seedMethod = EVENT_NUMBER;
	nStuckSteps = 0;
	primaryTrackID = -1;
	primaryParentID = -1;
	nextSecondaryTrackID = 1001;

	BLSignal::init();
	BLAlarm::init();

	G4StateManager * stateManager = G4StateManager::GetStateManager();
	stateManager->SetExceptionHandler(this);
#ifdef G4BL_GSL
	gsl_set_error_handler(gsl_error_handler);
#endif
}

void BLManager::delayedConstruction()
{
	runManager = new BLRunManager();
	G4StateManager * stateManager = G4StateManager::GetStateManager();
	stateManager->SetExceptionHandler(this);
}

BLManager::~BLManager()
{
	// cannot delete runManager, as it deletes the detector, 
	// steppingaction, eventaction, ..., and that becomes multiple
	// deletes of this object. But we are exiting, so we don't really
	// need to delete anything.
	//delete runManager;

	blManager = 0;
}

void BLManager::initialize()
{
	if(initialized) {
		G4Exception("BLManager","Already Initialized",FatalException,
									"");
	}
	initialized = true;

	G4UImanager* UI = G4UImanager::GetUIpointer();

	if(!physics) {
		G4Exception("BLManager","No physics registered",JustWarning,
					"Using default physics list.");
		G4String cmd("physics default");
		BLCommand::doCommand(cmd);
		BLAssert(physics != 0);
	}

	// get various parameter values (may have changed since constructor)
	histoUpdate = Param.getInt("histoUpdate");
	eventTimeLimit = Param.getInt("eventTimeLimit");
	zTolerance = Param.getDouble("zTolerance");
	steppingVerbose = Param.getInt("steppingVerbose");
	wallClockLimit = Param.getInt("wallClockLimit");

	setSteppingFormat();

	// set initialization and user classes
	runManager->SetUserInitialization((G4VUserDetectorConstruction*)this);
	/// physics->getPhysicsList() is set in registerPhysicsList()
	runManager->SetUserAction(new BLManager_UserSteppingAction(this));
	runManager->SetUserAction(new BLManager_UserTrackingAction(this));
	runManager->SetUserAction(new BLManager_UserRunAction(this));
	runManager->SetUserAction(new BLManager_UserEventAction(this));
	runManager->SetUserAction(new BLManager_PrimaryGeneratorAction(this));
	runManager->SetUserAction(new BLManager_UserStackingAction(this));

	// apply initial range cuts
	UI->ApplyCommand("/range/cutG 2 mm");
	UI->ApplyCommand("/range/cutE 2 mm");

	// initialize the RunManager, construct the apparatus, etc.
	runManager->Initialize();

	// set initial verbosities
	UI->ApplyCommand("/control/verbose 0");
	UI->ApplyCommand("/run/verbose 0");
	UI->ApplyCommand("/event/verbose 0");
	UI->ApplyCommand("/tracking/verbose 0");
	UI->ApplyCommand("/hits/verbose 0");
	UI->ApplyCommand("/material/verbose 0");
	UI->ApplyCommand("/process/setVerbose 0 all");
	UI->ApplyCommand("/process/verbose 0");
	UI->ApplyCommand("/process/eLoss/verbose 0");
	UI->ApplyCommand("/process/had/verbose 0");
	UI->ApplyCommand("/process/em/verbose 0");
	G4HadronicProcessStore::Instance()->SetVerbose(0);

	// initialize the beam
	for(beamIndex=0; beamIndex<beamVector.size(); ++beamIndex)
		beamVector[beamIndex]->init();
	for(unsigned i=0; i<referenceVector.size(); ++i)
		referenceVector[i]->init();
	beamIndex = 0;

	// create the ZStepLimiter
	new ZStepLimiter();
}

BLManager *BLManager::getObject()
{
	if(!blManager) new BLManager();
	return blManager;
}

void BLManager::insertZStep(std::vector<ZStep>& vector, G4double z,
							ZSteppingAction *action)
{
	// the bookends in vector (see BLManager()) ensure this works
	std::vector<ZStep>::iterator i;
	for(i=vector.begin(); i<vector.end(); ++i) {
		if(z < i->z) {	// preserve order for ==
			ZStep zs(z,action);
			vector.insert(i,zs);
			break;
		}
	}
}

void BLManager::registerZStep(G4double z, ZSteppingAction *sa, G4int when)
{
	if(when & 1) insertZStep(tuneZStep,z,sa);
	if(when & 2) insertZStep(referenceZStep,z,sa);
	if(when & 4) insertZStep(beamZStep,z,sa);
}

void BLManager::trackTuneAndReferenceParticles()
{
	if(referenceVector.size() == 0)
		return;

	if(!physics) {
		G4Exception("BLManager","No physics registered",FatalException,
									"");
	}

	// Tune and Reference particles cannot use collective mode
	bool collectiveMode = runManager->getCollectiveMode();
	runManager->setCollectiveMode(false);

	printf("================= Prepare Tune Particle(s) ===========\n");
	physics->setDoStochastics(FORCE_OFF,0);
	runManager->Initialize();

	printf("================= Begin Tune Particle(s) =============\n");
	state = TUNE;
	setEventID(-2);
	beamIndex = 0;
	runManager->BeamOn(referenceVector.size());
	state = IDLE;

	// now track center particle
	printf("================== Begin Reference Particle(s) ===============\n");
	state = REFERENCE;
	setEventID(-1);
	beamIndex = 0;
	runManager->BeamOn(referenceVector.size());
	state = IDLE;
	beamIndex = 0;

	physics->setDoStochastics(NORMAL,0);
	runManager->setCollectiveMode(collectiveMode);
}

void BLManager::handleSourceRun()
{
	for(int i=0; i<sourceRunVector.size(); ++i) {
		// setting sourceRun makes all BLManager callbacks just call it
		sourceRun = sourceRunVector[i];
		printf("================== Source Run ==================\n");
		state = SOURCE;
		sourceRun->prepare();
		runManager->BeamOn(0x7FFFFFFF);
		sourceRun->shutDown();
		state = IDLE;
		sourceRun = 0;
	}
}

void BLManager::trackBeam()
{
	if(beamVector.size() == 0) {
		G4Exception("BLManager","No beam registered",JustWarning,
									"");
		return;
	}
	if(!physics) {
		G4Exception("BLManager","No physics registered",FatalException,
									"");
	}

	printf("================== Prepare Tracking Beam ==================\n");
	physics->setDoStochastics(NORMAL);
	runManager->Initialize();

	printf("================== Begin Tracking Beam ===============\n");
	state = BEAM;
	setEventID(1);
	beamIndex = 0;
	runManager->BeamOn(0x7FFFFFFF);
	for(beamIndex=0; beamIndex<beamVector.size(); ++beamIndex)
		beamVector[beamIndex]->summary();
	beamIndex = 0;
	state = IDLE;
}

void BLManager::displayVisual()
{
#ifdef G4BL_VISUAL
	printf("================== Prepare Visualization ==================\n");
	state = VISUAL;
	physics->setDoStochastics(NORMAL);
	runManager->Initialize();
	setEventID(1);
	beamIndex = 0;

	// handle viewer=name and viewer=name,100,5 (5 images with 100 events)
	G4String viewer = Param.getString("viewer");
	std::vector<G4String> parts = BLCommand::splitString(viewer,",",true);
	int evPerImage = (parts.size() > 1 ? atoi(parts[1].c_str()) : 1);
	int nImages = (parts.size() > 2 ? atoi(parts[2].c_str()) : 1);
	viewer = parts[0];
	if(viewer == "best") nImages = 999999;
#ifdef G4BL_GUI
	if(viewer == "best") viewer = "Qt";
	if(viewer == "Qt") {
		// instantiate G4UIExecutive for Qt
		BLQt::init1();
		// re-arange Qt UI window for G4beamline (not Geant4 commands)
		BLQt::init2();
	}
#endif
	BLVisManager *visManager = new BLVisManager(viewer);
	visManager->init();

	handleCallbacks(4);

	visManager->generateImages(evPerImage, nImages);

	delete visManager;
	state = IDLE;
	beamIndex = 0;
#else
	G4Exception("BLManager","G4beamline was built without visualization",
							FatalException,"");
#endif
}

void BLManager::displayGeometry(G4VPhysicalVolume *phys, int level)
{
	if(phys == 0) phys = worldPhysicalVolume;
	for(int i=0; i<level; ++i) printf("  ");
	G4ThreeVector pos = phys->GetObjectTranslation();
	printf("%s x=%.3f y=%.3f z=%.3f\n",phys->GetName().c_str(),
			pos[0],pos[1],pos[2]);
	G4RotationMatrix rot = phys->GetObjectRotationValue();
	if(!rot.isIdentity()) {
		for(int i=0; i<level; ++i) printf("  ");
		BLCommand::dumpRotation(&rot,"");
	}
	G4LogicalVolume *log = phys->GetLogicalVolume();
	int n = log->GetNoDaughters();
	for(int i=0; i<n; ++i) {
		G4VPhysicalVolume *p = log->GetDaughter(i);
		if(p)
			displayGeometry(p,level+1);
	}
}

G4VPhysicalVolume *BLManager::Construct()
{
	// ensure the global field is initialized
	(void)BLGlobalField::getObject();

	// construct the world
	worldPhysicalVolume =  BLGroup::constructWorld();

	// make it invisible
	worldPhysicalVolume->GetLogicalVolume()->SetVisAttributes(
					BLCommand::getVisAttrib("invisible"));

	return worldPhysicalVolume;
}

void BLManager::handleCallbacks(int type)
{
	// startup sequencing prevents RunManager from registering
	// so just do it
	BLRunManager::getObject()->callback(type);

	if(type == 0) {
		for(unsigned i=0; i<preReferenceCallbackVector.size(); ++i)
			preReferenceCallbackVector[i]->callback(type);
	} else if(type == 1) {
		for(unsigned i=0; i<postReferenceCallbackVector.size(); ++i)
			postReferenceCallbackVector[i]->callback(type);
	} else if(type == 2) {
		for(unsigned i=0; i<postTrackingCallbackVector.size(); ++i)
			postTrackingCallbackVector[i]->callback(type);
		exceptionSummary();
	} else if(type == 3) {
		for(unsigned i=0; i<replaceMainLoopCallbackVector.size(); ++i)
			replaceMainLoopCallbackVector[i]->callback(type);
		if(replaceMainLoopCallbackVector.size() > 0) {
			BLAlarm::clear();
			// we have replaced the main loop, so closeup and exit
			BLNTuple::summary();
			BLNTuple::closeAll();
			// handle post-tracking callbacks
			handleCallbacks(2);
			blManager = 0;
			g4bl_exit(0);
		}
	} else if(type == 4) {
		for(unsigned i=0; i<visualizationCallbackVector.size(); ++i)
			visualizationCallbackVector[i]->callback(type);
	} else if(type == -1) {
		for(unsigned i=0; i<physicsCallbackVector.size(); ++i)
			physicsCallbackVector[i]->callback(type);
	}
}

void BLManager::BeginOfRunAction(const G4Run *run)
{
	if(sourceRun) {
		sourceRun->BeginOfRunAction(run);
		return;
	}

	startEvent = startRun = BLTime::time();
	eventsProcessed = 0;
	endRun = false;

	for(unsigned int i=0; i<runActionVector.size(); ++i) {
		runActionVector[i]->BeginOfRunAction(run);
	}

	if(state == BEAM) {
		for(unsigned int i=0; i<beamRunActionVector.size(); ++i)
			beamRunActionVector[i]->BeginOfRunAction(run);
	}
}

void BLManager::EndOfRunAction(const G4Run *run)
{
	if(sourceRun) {
		sourceRun->EndOfRunAction(run);
		return;
	}

	if(state == BEAM) {
		for(unsigned int i=0; i<beamRunActionVector.size(); ++i)
			beamRunActionVector[i]->EndOfRunAction(run);
	}

	for(unsigned int i=0; i<runActionVector.size(); ++i)
		runActionVector[i]->EndOfRunAction(run);

	if(state == TUNE || state == REFERENCE)
		++eventsProcessed;

	if(state != VISUAL)
		printf("Run complete  %d Events  %ld seconds\n",
			eventsProcessed,(long)(BLTime::time()-startRun));
}

void BLManager::BeginOfEventAction(const G4Event* event)
{
	if(sourceRun) {
		sourceRun->BeginOfEventAction(event);
		return;
	}

if(endRun || event->IsAborted()) return;

	if(state == BEAM) {

	    int evId = event->GetEventID();
	    if(steppingVerbose && !runManager->getCollectiveMode()) 
		printf("\n\n=================== Event %d ==================\n",
								evId);

	    // call all registered action-s
	    std::vector<BLManager::EventAction*>::iterator i;
	    for(i=beamEventActionVector.begin(); i<beamEventActionVector.end(); ++i) {
		(*i)->BeginOfEventAction(event);
	    }
	}

	// call all registered action-s
	std::vector<BLManager::EventAction*>::iterator i;
	for(i=eventActionVector.begin(); i<eventActionVector.end(); ++i) {
		(*i)->BeginOfEventAction(event);
	}

	startEvent = BLTime::time();
	// add 10 seconds so the test in UserSteppingAction() can try first
	// (it kills 1 event, BLAlarm kills the entire job).
	if(eventTimeLimit > 0) BLAlarm::set(eventTimeLimit+10);
}

void BLManager::EndOfEventAction(const G4Event* event)
{
	if(sourceRun) {
		sourceRun->EndOfEventAction(event);
		return;
	}

	BLAlarm::clear();

	if(endRun || event->IsAborted()) return;

	// call all registered action-s
	std::vector<BLManager::EventAction*>::iterator i;
	for(i=eventActionVector.begin(); i<eventActionVector.end(); ++i) {
		(*i)->EndOfEventAction(event);
	}

	if(state != BEAM) return;

	// call all registered action-s
	for(i=beamEventActionVector.begin(); i<beamEventActionVector.end(); ++i) {
		(*i)->EndOfEventAction(event);
	}
	BLNTuple::update();

	if(endRun) return; // (might have been set in a user action)

	incrEventsProcessed(event->GetEventID());

	// update histogram file, if appropriate
	if((histoUpdate > 0 && eventsProcessed%histoUpdate == 0 && 
	   eventsProcessed > 0) || BLSignal::sigusr1()) {
		G4Exception("BLManager","NTuples flushed",JustWarning,"");
		BLNTuple::flushAll();
	}
	if(BLSignal::sigusr2())
		G4Exception("BLManager","SIGUSR2 received - Exiting",
							FatalException,"");
}

void BLManager::PreUserTrackingAction(const G4Track *track)
{
	static bool first=true;
	if(first) {
		first=false;
		std::vector<BLManager::TrackingAction*>::iterator i;
		for(i=trackingActionVector.begin(); 
					 i<trackingActionVector.end(); ++i) {
			(*i)->SetTrackingManagerPointer(fpTrackingManager);
		}
	}

	nStuckSteps = 0;

	// link a BLTrackInfo into the track (called coord -- historical)
	BLTrackInfo *coord = (BLTrackInfo *)track->GetUserInformation();
	if(!coord || !coord->isValid()) {
		coord = new BLTrackInfo();
		((G4Track *)track)->SetUserInformation(coord);
	}
	coord->setGlobal(track->GetPosition(),track->GetGlobalTime());

	// set the external trackID and parentID
	trackID = coord->getExternalTrackID();
	if(trackID <= 0) {
		if(primaryTrackID >= 0) {
			trackID = primaryTrackID;
			coord->setExternalTrackID(trackID);
			coord->setExternalParentID(primaryParentID);
			primaryTrackID = primaryParentID = -1;
		} else {
			trackID = nextSecondaryTrackID++;
			coord->setExternalTrackID(trackID);
			coord->setExternalParentID(
					trackIDMap[track->GetParentID()]);
		}
	}
	trackIDMap[track->GetTrackID()] = trackID;
	currentTrack = track;

	// handle sourceRun
	if(sourceRun) {
		sourceRun->PreUserTrackingAction(track);
		return;
	}

	// print header, if appropriate
	if(steppingVerbose && !runManager->getCollectiveMode()) {
	    printf("=========== EventID %d TrackID %d %s    ParentID %d  KE=%.3fMeV",
			eventID,trackID,
			track->GetDefinition()->GetParticleName().c_str(),
			getExternalParentID(track),track->GetKineticEnergy());
	    if(getExternalParentID(track) == 0) {
	    	printf("   CreatorProcess=Beam");
	    } else {
	    	const G4VProcess *proc = track->GetCreatorProcess();
		if(proc != 0)
			printf("   CreatorProcess=%s",
					proc->GetProcessName().c_str());
	    }
	    printf(" ===========\n");
	}

	// set currentZStep
	if(state == TUNE)
		currentZStep = &tuneZStep;
	else if(state == REFERENCE)
		currentZStep = &referenceZStep;
	else
		currentZStep = &beamZStep;
	indexZStep = 1;
	prevZ = coord->getCLZ();
	while(prevZ > (*currentZStep)[indexZStep].z) {
		++indexZStep;
	}
	BLAssert(indexZStep < currentZStep->size());
	G4double dz = (*currentZStep)[indexZStep].z - prevZ;
	if(dz < zTolerance) dz = zTolerance;
	ZStepLimiter::setMaxStep(dz);

	// loop over all registered user tracking actions
	std::vector<BLManager::TrackingAction*>::iterator i;
	for(i=trackingActionVector.begin(); i<trackingActionVector.end(); ++i) {
		(*i)->PreUserTrackingAction(track);
	}
}

void BLManager::PostUserTrackingAction(const G4Track *track)
{
	if(sourceRun) {
		sourceRun->PostUserTrackingAction(track);
		return;
	}

	std::vector<BLManager::TrackingAction*>::iterator i;
	for(i=trackingActionVector.begin(); i<trackingActionVector.end(); ++i) {
		(*i)->PostUserTrackingAction(track);
	}

	// delete the BLcoordinates
	G4VUserTrackInformation *ti = track->GetUserInformation();
	if(ti && !runManager->getCollectiveMode()) {
		delete ti;
		((G4Track *)track)->SetUserInformation(0);
	}

	trackID = -9999;
	currentTrack = 0;
}

void BLManager::UserSteppingAction(const G4Step *step)
{
	if(sourceRun) {
		sourceRun->UserSteppingAction(step);
		return;
	}

	if(BLSignal::received()) {
		char tmp[64];
		snprintf(tmp,sizeof(tmp),"Signal %d received",
							BLSignal::value());
		G4Exception("BLManager",tmp,FatalException,"");
	}

	G4Track *track = step->GetTrack();
	G4StepPoint *prePoint = step->GetPreStepPoint();
	G4StepPoint *postPoint = step->GetPostStepPoint();
	G4VPhysicalVolume *preVol=0;
	G4VPhysicalVolume *postVol=0;
	if(prePoint) preVol = prePoint->GetPhysicalVolume();
	if(postPoint) postVol = postPoint->GetPhysicalVolume();

	// need to update BLCoordinates before steppingVerbose print.
	// So BLCoordinates is not registered with SteppingAction, it
	// is handled specially, right here.
	BLCoordinates::update(track);
	BLCoordinates *coord = (BLCoordinates *)track->GetUserInformation();
	if(!coord || !coord->isValid()) {
		G4Exception("BLManager","Coordinates Got Lost",FatalException,
									"");
	}

	// steppingVerbose print moved to the top of this routine
	if(steppingVerbose > 0) { 
		static bool first=true;
		int nstep = track->GetCurrentStepNumber();
		if(first || (nstep <= 1 && !runManager->getCollectiveMode()))
			steppingVerbosePrint(0,0,1);
		steppingVerbosePrint(step,track,0);
		first = false;
	}

	G4TrackStatus status = track->GetTrackStatus();
	if(status == fStopAndKill || status == fKillTrackAndSecondaries)
		goto quit;

	/* check for stuck track */
	if(step->GetStepLength() > 0.0001*mm ||
	   step->GetDeltaTime() > 0.0001*ns) {
		nStuckSteps = 0;
	} else if(status == fAlive && ++nStuckSteps >= 100) {
		if(track->GetVelocity()/c_light < 0.010) {
		    ++stuckTracks;
		    track->SetTrackStatus(fStopButAlive);
		    if(steppingVerbose > 0)
		    	printf("BLManager: Stuck Track -- stopped\n");
		} else {
		    G4Exception("BLManager Stepping Action",
		    "Stuck Track -- Killed",JustWarning,
		    "100 steps in a row, each less than 0.1 micron and 0.1 ps");
		    track->SetTrackStatus(fStopAndKill);
		}
	}

	// check event time limit
	if(eventTimeLimit > 0 && BLTime::time()-startEvent > eventTimeLimit) {
		G4Exception("BLManager","Event Time Limit",EventMustBeAborted,"");
		goto quit;
	}
	if(wallClockLimit>0 && BLTime::time()-startProgram>=wallClockLimit) {
		G4Exception("BLManager","Wall Clock Limit - Exiting",
						FatalException,"");
		goto quit;
	}

	// call ZStep actions
	// NOTE: indexZstep points to the next entry in the +Z direction.
	// Note also the bookends constrain indexZStep to always be valid.
	if(currentZStep->size() > 2) { // i.e. not just 2 bookends
		G4double thisZ = coord->getCLZ();
		if(thisZ == prevZ) goto noZstep;
		G4double minZ = (thisZ<prevZ ? thisZ : prevZ);
		G4double maxZ = (thisZ>prevZ ? thisZ : prevZ);
		// find first entry that might be spanned by the step
		while(indexZStep > 0 && minZ <= (*currentZStep)[indexZStep-1].z)
			--indexZStep;
		while(minZ > (*currentZStep)[indexZStep].z)
			++indexZStep;
		int indexPrev = indexZStep-1;
		// loop over all entries actually spanned by this step
		for( ; indexZStep<currentZStep->size()-1; ++indexZStep) {
			G4double z=(*currentZStep)[indexZStep].z;
			if(maxZ <= z) break;
			if((*currentZStep)[indexZStep].action == 0) continue;
			// save current values from track
			G4ThreeVector pos = track->GetPosition();
			G4ThreeVector mom = track->GetMomentum();
			G4double time = track->GetGlobalTime();
			G4double ke = track->GetKineticEnergy();
			G4double properTime = track->GetProperTime();
			G4double trackLength = track->GetTrackLength();
			G4double localTime = track->GetLocalTime();
			G4ThreeVector pol = track->GetPolarization();
			// interpolate linearly to z, global coords
			G4ThreeVector deltaPos=step->GetDeltaPosition();
			G4ThreeVector deltaMom=(postPoint->GetMomentum() -
						prePoint->GetMomentum());
			G4double deltaTime=step->GetDeltaTime();
			G4double deltaE=(postPoint->GetKineticEnergy() -
						prePoint->GetKineticEnergy());
			G4double deltaProperTime=(postPoint->GetProperTime() -
						prePoint->GetProperTime());
			G4double deltatrackLength=step->GetStepLength();
			G4ThreeVector deltaPol=(postPoint->GetPolarization() -
						prePoint->GetPolarization());
			G4double f=(z-thisZ)/(prevZ-thisZ);
			track->SetPosition(pos-f*deltaPos);
			track->SetMomentumDirection((mom-f*deltaMom).unit());
			track->SetGlobalTime(time-f*deltaTime);
			track->SetKineticEnergy(ke-f*deltaE);
			track->SetProperTime(properTime-f*deltaProperTime);
			track->AddTrackLength(-f*deltatrackLength);
			track->SetLocalTime(localTime-f*deltaTime);
			track->SetPolarization(pol-f*deltaPol);
			BLCoordinates::update(track);
			(*currentZStep)[indexZStep].action->
						UserZSteppingAction(track);
			// restore current values to track
			track->SetPosition(pos);
			track->SetMomentumDirection(mom.unit());
			track->SetGlobalTime(time);
			track->SetKineticEnergy(ke);
			track->SetProperTime(properTime);
			track->AddTrackLength(f*deltatrackLength);
			track->SetLocalTime(localTime);
			track->SetPolarization(pol);
			BLCoordinates::update(track);
		}
		double dz = (thisZ>prevZ ? (*currentZStep)[indexZStep].z-thisZ : 
					   thisZ-(*currentZStep)[indexPrev].z);
		dz = (dz>zTolerance*2.0 ? dz : zTolerance*2.0);
		ZStepLimiter::setMaxStep(dz);
		prevZ = thisZ;
	}
noZstep:

	// call all-state stepping actions before per-state actions.
	{ std::vector<BLManager::SteppingAction*>::iterator i;
	  for(i=allStepVector.begin(); i!=allStepVector.end(); ++i) {
		  (*i)->UserSteppingAction(step);
		  if(track->GetTrackStatus() != fAlive)
			  goto quit;
	  }
	  if(preVol && allStepMap.count(preVol) > 0) {
		allStepMap[preVol]->UserSteppingAction(step);
	  	if(track->GetTrackStatus() != fAlive)
			goto quit;
	  }
	  if(postVol  && preVol != postVol && allStepMap.count(postVol) > 0) {
		allStepMap[postVol]->UserSteppingAction(step);
	  	if(track->GetTrackStatus() != fAlive)
			goto quit;
	  }
	}

	// call Tune Particle stepping actions
	if(state == TUNE) {
		std::vector<BLManager::SteppingAction*>::iterator i;
		for(i=tpStepVector.begin(); i!=tpStepVector.end(); ++i) {
			(*i)->UserSteppingAction(step);
			if(track->GetTrackStatus() != fAlive)
				goto quit;
		}
		if(preVol && tpStepMap.count(preVol) > 0) {
			tpStepMap[preVol]->UserSteppingAction(step);
			if(track->GetTrackStatus() != fAlive)
				goto quit;
		}
		if(postVol && preVol != postVol && 
						tpStepMap.count(postVol) > 0) {
			tpStepMap[postVol]->UserSteppingAction(step);
			if(track->GetTrackStatus() != fAlive)
				goto quit;
		}
	}

	// call Reference Particle stepping actions
	if(state == REFERENCE) {
		std::vector<BLManager::SteppingAction*>::iterator i;
		for(i=rpStepVector.begin(); i!=rpStepVector.end(); ++i) {
			(*i)->UserSteppingAction(step);
			if(track->GetTrackStatus() != fAlive)
				goto quit;
		}
		if(preVol && rpStepMap.count(preVol) > 0) {
			rpStepMap[preVol]->UserSteppingAction(step);
			if(track->GetTrackStatus() != fAlive)
				goto quit;
		}
		if(postVol && preVol != postVol && 
						rpStepMap.count(postVol) > 0) {
			rpStepMap[postVol]->UserSteppingAction(step);
			if(track->GetTrackStatus() != fAlive)
				goto quit;
		}
	}

	// call beam stepping actions
	if(state == BEAM) {
		std::vector<BLManager::SteppingAction*>::iterator i;
		for(i=beamStepVector.begin(); i!=beamStepVector.end(); ++i) {
			(*i)->UserSteppingAction(step);
			if(track->GetTrackStatus() != fAlive)
				goto quit;
		}
		if(preVol && beamStepMap.count(preVol) > 0) {
			beamStepMap[preVol]->UserSteppingAction(step);
			if(track->GetTrackStatus() != fAlive)
				goto quit;
		}
		if(postVol  && preVol != postVol && 
					     beamStepMap.count(postVol) > 0) {
			beamStepMap[postVol]->UserSteppingAction(step);
			if(track->GetTrackStatus() != fAlive)
				goto quit;
		}
	}

quit:	;
}

void BLManager::GeneratePrimaries(G4Event *event)
{
	if(sourceRun) {
		eventID = prevEventID + 1;
		event->SetEventID(eventID);
		if(!sourceRun->nextSourceEvent(event)) goto end_run;
		eventID = event->GetEventID();
		prevEventID = eventID; // (this increments it, too)
		return;
	}

	if(beamVector.size() == 0) {
		G4Exception("BLManager","No beam registered",FatalException,
									"");
	}

	switch(state) {
	case IDLE:
	case SPECIAL:
		G4Exception("BLManager","Erroneous call to GeneratePrimaries",
							FatalException,"");
	case TUNE:
		setEventID(-2);
		event->SetEventID(-2);
		while(beamIndex < referenceVector.size()) {
			if(referenceVector[beamIndex++]->
					generateReferenceParticle(event))
			   	return;
		}
		goto end_run;
	case REFERENCE:
		setEventID(-1);
		event->SetEventID(-1);
		while(beamIndex < referenceVector.size()) {
			if(referenceVector[beamIndex++]->
					generateReferenceParticle(event))
			   	return;
		}
		goto end_run;
	case VISUAL:
	case BEAM:
		// Checking skipEvent() has been moved into nextBeamEvent()
		// default event # (nextBeamEvent() can change it)
		eventID = prevEventID + 1;
		event->SetEventID(eventID);
		// generate the event, indexing through beamVector[]
		if(beamIndex >= beamVector.size()) {
			goto end_run;
		}
		while(!beamVector[beamIndex]->nextBeamEvent(event)){
			if(++beamIndex >= beamVector.size())
				goto end_run;
		}
		// update eventID (may have changed)
		eventID = event->GetEventID();
		prevEventID = eventID; // (this increments it, too)
		break;
	case SOURCE:		// never get here (handled above)
		break;
	}
	return;
end_run:
	// RunManager cannot abort the event from inside
	// UserGeneratePrimaries(), so we do a soft abort
	// to the RunManager, and abort the event ourself.
	// The result is the same as a hard abort.
	runManager->AbortRun(true);
	event->SetEventAborted();
	beamIndex = 0;
	endRun = true;
}

G4ClassificationOfNewTrack BLManager::ClassifyNewTrack(const G4Track *track)
{
	if(sourceRun) {
		return sourceRun->ClassifyNewTrack(track);
	}

	std::vector<BLManager::StackingAction*>::iterator i;
	for(i=stackingActionVector.begin(); i<stackingActionVector.end(); ++i) {
		G4ClassificationOfNewTrack c = (*i)->ClassifyNewTrack(track);
		if(c == fKill) return c;
	}
	return fUrgent;
}

void BLManager::NewStage()
{
	if(sourceRun) {
		sourceRun->NewStage();
		return;
	}

	std::vector<BLManager::StackingAction*>::iterator i;
	for(i=stackingActionVector.begin(); i<stackingActionVector.end(); ++i) {
		(*i)->NewStage();
	}
}

void BLManager::PrepareNewEvent()
{
	if(sourceRun) {
		sourceRun->PrepareNewEvent();
		return;
	}

	std::vector<BLManager::StackingAction*>::iterator i;
	for(i=stackingActionVector.begin(); i<stackingActionVector.end(); ++i) {
		(*i)->PrepareNewEvent();
	}
}

int BLManager::getExternalTrackID(const G4Track *track)
{
	if(runManager->getCollectiveMode()) return track->GetTrackID();
	BLTrackInfo *p = (BLTrackInfo *)track->GetUserInformation();
	if(p && p->isValid()) 
		return p->getExternalTrackID();
	return trackIDMap[track->GetTrackID()];
}

int BLManager::getExternalParentID(const G4Track *track)
{
	if(runManager->getCollectiveMode()) return track->GetParentID();
	BLTrackInfo *p = (BLTrackInfo *)track->GetUserInformation();
	if(p && p->isValid()) 
		return p->getExternalParentID();
	return trackIDMap[track->GetParentID()];
}

void BLManager::setExternalTrackID(G4Track *track, int trackID, int parentID)
{
	if(trackID < 0) trackID = nextSecondaryTrackID++;

	BLTrackInfo *ti = (BLTrackInfo *)track->GetUserInformation();
	if(!ti || !ti->isValid()) {
		ti = new BLTrackInfo();
		track->SetUserInformation(ti);
	}
	ti->setExternalTrackID(trackID);
	ti->setExternalParentID(parentID);
}

void BLManager::incrEventsProcessed(int eventID)
{
	// print event number, if appropriate
	++eventsProcessed;
	if(!runManager->getCollectiveMode()) {
		if(eventsProcessed <= 10 ||
		   (eventsProcessed < 100 && eventsProcessed%10 == 0) ||
		   (eventsProcessed < 1000 && eventsProcessed%100 == 0) ||
		   eventsProcessed%1000 == 0) {
			printf("Event %d Completed",eventID);
			int t = BLTime::time() - startRun;
			if(t <= 0) t = 1;
			printf("  %d events  realTime=%d sec  %.1f ev/sec",
				eventsProcessed,t,(double)eventsProcessed/t);
			printf("\n");
			fflush(stdout);
		}
	}
}

G4bool BLManager::Notify(const char* originOfException,
	const char* exceptionCode, G4ExceptionSeverity severity,
	const char* description)
{
	fflush(stdout);

	G4bool abortProgram = false;
	const char *p="UNKNOWN";
	switch(severity) {
	case FatalException:
     		p = "Fatal Exception";		
     		abortProgram = true;
		++fatalExceptions;
		break;
	case FatalErrorInArgument:
     		p = "Fatal Error in Argument";	
     		abortProgram = true;
		++fatalExceptions;
		break;
	case RunMustBeAborted:
     		p = "Run Must Be Aborted";
		runManager->AbortRun(false);
     		abortProgram = false;
		++fatalExceptions;
		break;
	case EventMustBeAborted:
		BLAlarm::clear();
     		p = "Event must be Aborted";
		runManager->AbortEvent();
     		abortProgram = false;
		++eventsAborted;
		break;
	default:
     		p = "Warning";
     		abortProgram = false;
		++warnings;
		break;
	}

	const char *pname="";
	G4double kineticEnergy=0.0;
	int evid=0, trkid=0;
	if(G4StateManager::GetStateManager()->GetCurrentState() ==
				G4State_EventProc && currentTrack != 0) {
		pname = currentTrack->GetDefinition()->
						GetParticleName().c_str();
		kineticEnergy = currentTrack->GetKineticEnergy();
		evid = eventID;
		trkid = trackID;
	}
	if(BLMPI::isMPI() && !BLMPI::isRank0())
		BLMPI::sendExceptionMessage(originOfException,exceptionCode,p,
		      description,evid,trkid,pname,kineticEnergy,abortProgram);
	printException(originOfException,exceptionCode,p,
		      description,evid,trkid,pname,kineticEnergy,abortProgram);

	if(abortProgram) {
		BLAlarm::clear();
		fprintf(stdout,"g4beamline: attempting to close up after fatal G4Exception...\n");
		int value=99;
		if(strstr(exceptionCode,"Exiting") != 0)
			value = 0;
		if(BLMPI::isMPI())
			BLMPI::closeupAndExit(value);
		BLNTuple::summary();
		BLNTuple::closeAll();
		handleCallbacks(2);
		g4bl_exit(value);
	}

	return abortProgram;
}

int BLManager::tallyException(const char *code)
{
	return ++exceptionCount[code];
}

void BLManager::printException(const char *origin, const char *code,
		const char *severity, const char *description,
		int evid, int trkid, const char *particleName,
		G4double kineticEnergy, bool abortProgram, bool printOnce)
{
	// thin out the printing of many similar exceptions
	int n = tallyException(code);
	if(printOnce && n>1) return;
	int interval = 1;
	if(n >= 1000)
		interval = 1000;
	else if(n >= 100)
		interval = 100;
	else if(n >= 10)
		interval = 10;
	if(allExceptions) interval = 1;
	if(n%interval != 0)
		return;

	fprintf(stdout,"**************************************************************************\n");
	fprintf(stdout,"*** G4Exception: %s\n",code);
	fprintf(stdout,"***    severity: %s\n",severity);
	fprintf(stdout,"***   issued by: %s\n",origin);
	if(strlen(description) > 0)
		fprintf(stdout,"*** description: %s\n",description);
	if(strstr(description,"Missing mandatory data")) {
		fprintf(stdout,
		    "***              Run 'g4bldata --install' to fix this.\n");
	} else if(particleName[0] != '\0') {
		fprintf(stdout,"***     EventID: %d     TrackID: %d   %s  KE=%.4g MeV\n",
				evid,trkid,particleName,kineticEnergy);
	}
	if(interval > 1)
		fprintf(stdout,"***    printing: every %d-th occurrence\n",
			interval);
	fprintf(stdout,"**************************************************************************\n");

	fflush(stdout);
}

void BLManager::exceptionSummary()
{
	printf("\nExceptions: %d Fatal, %d Events Aborted, "
			"%d Stuck Tracks (stopped), %d Warnings\n",
			fatalExceptions,eventsAborted,stuckTracks,warnings);
	std::map<G4String,int>::iterator it;
	for(it=exceptionCount.begin(); it!=exceptionCount.end(); ++it) {
		G4String except = it->first;
		int count = it->second;
		printf(" %6d times: %s\n",count,except.c_str());
	}
}

void BLManager::setSteppingFormat()
{
	G4String fmt = Param.getString("steppingFormat");
	BLArgumentVector args;
	BLArgumentMap namedArgs;
	
	// convert "," to " "
	G4String s;
	for(unsigned i=0; i<fmt.size(); ++i)
		s += (fmt[i]==',' ? ' ' : fmt[i]);

	BLCommand::parseArgs(s,args,namedArgs);
	appendVerboseFormat("");
	for(unsigned i=0; i<args.size(); ++i)
		appendVerboseFormat(args[i]);
	appendVerboseFormat("NEWLINE");
}

G4String BLManager::getFormatHelp()
{
	G4String s = 
		"        EXT       toggle extended precision (3 more digits)\n"
		"        TAG       print a '>' (useful to grep output)\n"
		"        N         step number\n"
		"        NSTEP     Synonym of N\n"
		"        GLOBAL    X,Y,Z,T in global coords\n"
		"        XYZT      Synonym of GLOBAL\n"
		"        CL        X,Y,Z,dxdz,dydz in CL coords\n"
		"        CLX       extended precision CL\n"
		"        KE        kinetic energy\n"
		"        STEP      step length\n"
		"        STEPLEN   Synonym of STEP\n"
		"        VOL       volume name\n"
		"        VOLNAME   Synonym of VOL\n"
		"        PROCESS   process name\n"
		"        B         magnetic field\n"
		"        E         electric field\n"
		"        P         3-momentum\n"
		"        MAT       material name\n"
		"        ID        event ID, track ID, parent ID\n"
		"        PART      particle name\n"
		"        SEG       centerline coord segment number\n"
		"        WT        weight\n"
		"        POLAR     polarization\n"
		"        NL        <newline>\n"
		"        NEWLINE   Synonym of NL\n"
		"        \\n        Synonym of NL\n";

	return s;
}

void BLManager::appendVerboseFormat(G4String fmt)
{
	G4String f;
	for(const char *p=fmt.c_str(); *p; ++p)
		f += toupper(*p);
	if(f == "") verboseFormat.clear();
	else if(f == "TAG") verboseFormat.push_back(TAG);
	else if(f == "EXT") verboseFormat.push_back(EXT);
	else if(f == "N") verboseFormat.push_back(NSTEP);
	else if(f == "NSTEP") verboseFormat.push_back(NSTEP);
	else if(f == "GLOBAL") verboseFormat.push_back(GLOBAL);
	else if(f == "XYZT") verboseFormat.push_back(GLOBAL);
	else if(f == "CL") verboseFormat.push_back(CL);
	else if(f == "CLX") verboseFormat.push_back(CLX);
	else if(f == "KE") verboseFormat.push_back(KE);
	else if(f == "STEP") verboseFormat.push_back(STEP);
	else if(f == "STEPLEN") verboseFormat.push_back(STEP);
	else if(f == "VOL") verboseFormat.push_back(VOL);
	else if(f == "VOLNAME") verboseFormat.push_back(VOL);
	else if(f == "PROCESS") verboseFormat.push_back(PROCESS);
	else if(f == "B") verboseFormat.push_back(B);
	else if(f == "E") verboseFormat.push_back(E);
	else if(f == "P") verboseFormat.push_back(P);
	else if(f == "MAT") verboseFormat.push_back(MAT);
	else if(f == "ID") verboseFormat.push_back(ID);
	else if(f == "PART") verboseFormat.push_back(PART);
	else if(f == "SEG") verboseFormat.push_back(SEG);
	else if(f == "WT") verboseFormat.push_back(WT);
	else if(f == "POLAR") verboseFormat.push_back(POLAR);
	else if(f == "\n") verboseFormat.push_back(NEWLINE);
	else if(f == "NL") verboseFormat.push_back(NEWLINE);
	else if(f == "NEWLINE") verboseFormat.push_back(NEWLINE);
	else printf("BLManager: invalid verbose format '%s'\n",
		fmt.c_str());
}

void BLManager::steppingVerbosePrint(const G4Step *step, const G4Track *track, 
							int header)
{
	G4StepPoint *prePoint=0;
	G4VPhysicalVolume *preVol=0;
	G4StepPoint *postPoint=0;
	G4double time=0;
	G4ThreeVector position;
	G4Material *mat=0;
	int eventNum = -999;
	int nstep=0;
	const char *procName="??";
	BLCoordinates *coord = (track!=0 ? (BLCoordinates *)track->GetUserInformation() : 0);
	if(coord && !coord->isValid()) coord = 0;

	if(!header) {
		prePoint = step->GetPreStepPoint();
		preVol = prePoint->GetPhysicalVolume();
		postPoint = step->GetPostStepPoint();
		time = track->GetGlobalTime();
		position = track->GetPosition();
		mat = preVol->GetLogicalVolume()->GetMaterial();
		nstep = track->GetCurrentStepNumber();
		procName = (postPoint->GetProcessDefinedStep() ?
		  postPoint->GetProcessDefinedStep()->GetProcessName().c_str() :
		  "UserLimit"); 
		if(track->GetTrackStatus() == fStopButAlive) {
			static char name[128];
			sprintf(name,"Stopped:%s",procName);
			procName = name;
		}
		if(track->GetTrackStatus() == fStopAndKill) {
			static char name[128];
			sprintf(name,"Killed:%s",procName);
			procName = name;
		}
		const G4Event* event = runManager->GetCurrentEvent();
		eventNum = event->GetEventID();
	}

	bool extended=false;
	for(unsigned i=0; i<verboseFormat.size(); ++i) {
		switch(verboseFormat[i]) {
		case EXT:
			extended = !extended;
			break;
		case TAG:
			if(header) printf(" ");
			else printf(">");
			break;
		case NSTEP:
			if(header) printf(" Step");
			else printf(" %4d",nstep);
			break;
		case GLOBAL:
		    if(!extended) {
			if(header) printf("   X(mm)   Y(mm)   Z(mm)    T(ns)");
			else printf(" %7.1f %7.1f %7.1f %8.2f",
				position[0],position[1],position[2],
				time);
		    } else {
			if(header) printf("      X(mm)      Y(mm)      Z(mm)      T(ns) ");
			else printf(" %10.4f %10.4f %10.4f %11.5f",
				position[0],position[1],position[2],
				time);
		    }
			break;
		case CL:
		    if(extended) goto clx;
		    if(header) {
			printf("   CL: X        Y        Z     dxdz    dydz");
		    } else {
			G4ThreeVector dir = track->GetMomentumDirection();
			// transform to centerline coordinates
			G4ThreeVector clpos;
			coord->getCoords(BLCOORD_CENTERLINE,clpos);
			dir = coord->getRotation() * dir;
			printf(" %8.1f %8.1f %8.1f %7.4f %7.4f",
					clpos[0],clpos[1],clpos[2],
					dir[0]/dir[2],dir[1]/dir[2]);
		    }
			break;
		case CLX:
clx:		    if(header) {
			printf("   CL: X           Y           Z        dxdz       dydz   ");
		    } else {
			G4ThreeVector dir = track->GetMomentumDirection();
			// transform to centerline coordinates
			G4ThreeVector clpos;
			coord->getCoords(BLCOORD_CENTERLINE,clpos);
			dir = coord->getRotation() * dir;
			printf(" %11.4f %11.4f %11.4f %10.7f %10.7f",
					clpos[0],clpos[1],clpos[2],
					dir[0]/dir[2],dir[1]/dir[2]);
		    }
			break;
		case KE:
		    if(!extended) {
			if(header) printf("  KE(MeV)");
			else printf(" %8.1f",track->GetKineticEnergy());
		    } else {
			if(header) printf("    KE(MeV) ");
			else printf(" %11.4f",track->GetKineticEnergy());
		    }
			break;
		case STEP:
		    if(!extended) {
			if(header) printf("  StepLen");
			else printf(" %8.2f",step->GetStepLength());
		    } else {
			if(header) printf("    StepLen ");
			else printf(" %11.5f",step->GetStepLength());
		    }
			break;
		case VOL:
			if(header) printf(" This Volume     ");
			else printf(" %-16s",preVol->GetName().c_str());
			break;
		case PROCESS:
			if(header) printf(" Process        ");
			else printf(" %-16s",procName);
			break;
		case B:
		    if(header) {
		    	printf("       Bx,By,Bz (Tesla)     ");
		    } else {
			G4double point[4], field[6];
			point[0]=position[0],point[1]=position[1],point[2]=position[2];
			point[3]=time;
			BLGlobalField::getObject()->GetFieldValue(point,field);
			G4ThreeVector B(field[0],field[1],field[2]);
			B = coord->getRotation(BLCOORD_CENTERLINE) * B;
			printf(" %8.4f %8.4f %8.4f",
				B[0]/tesla,B[1]/tesla,B[2]/tesla);
		    }
			break;
		case E:
		    if(header) {
		    	printf("      Ex,Ey,Ez (MV/meter)  ");
		    } else {
			G4double point[4], field[6];
			point[0]=position[0],point[1]=position[1],point[2]=position[2];
			point[3]=time;
			BLGlobalField::getObject()->GetFieldValue(point,field);
			G4ThreeVector E(field[3],field[4],field[5]);
			E = coord->getRotation(BLCOORD_CENTERLINE) * E;
			printf(" %8.4f %8.4f %8.4f",E[0]/(megavolt/meter),
				E[1]/(megavolt/meter),E[2]/(megavolt/meter));
		    }
			break;
		case P:
		  if(!extended) {
		    if(header) {
		    	printf("      Px,Py,Pz (MeV/c)     ");
		    } else {
			G4ThreeVector P = track->GetMomentum();
			P = coord->getRotation(BLCOORD_CENTERLINE) * P;
			printf(" %8.1f %8.1f %8.1f",P[0]/(MeV),
				P[1]/(MeV),P[2]/(MeV));
		    }
		  } else {
		    if(header) {
		    	printf("          Px,Py,Pz (MeV/c)          ");
		    } else {
			G4ThreeVector P = track->GetMomentum();
			P = coord->getRotation(BLCOORD_CENTERLINE) * P;
			printf(" %11.4f %11.4f %11.4f",P[0]/(MeV),
				P[1]/(MeV),P[2]/(MeV));
		    }
		  }
			break;
		case MAT:
		    if(header) {
		    	printf(" Material ");
		    } else {
			printf(" %-9s",mat->GetName().c_str());
		    }
			break;
		case ID:
		    if(header) {
		    	printf(" Event# Trk Prnt");
		    } else {
			printf(" %6d %3d %4d",eventNum,
				BLManager::getObject()->
				    getExternalTrackID(track),
				BLManager::getObject()->
				    getExternalParentID(track));
		    }
			break;
		case PART:
		    if(header) {
		    	printf(" Particle");
		    } else {
			printf(" %8s",track->GetDefinition()->GetParticleName().c_str());
		    }
			break;
		case SEG:
		    if(header) {
		    	printf(" Seg");
		    } else {
			printf(" %3d",coord->getSegmentCL());
		    }
			break;
		case WT:
			if(header) {
			    printf("  Wt  ");
			} else {
			    printf(" %.3f",track->GetWeight());
			}
			break;
		case POLAR:
			if(header) {
				printf("    Polarization   ");
			} else {
				G4ThreeVector polar=track->GetPolarization();
				printf(" %5.3f,%5.3f,%5.3f ",
						polar[0],polar[1],polar[2]);
			}
			break;
		case NEWLINE:
			printf("\n");
			break;
		}
	}

	fflush(stdout); // help debugging location of crashes
}

#ifdef G4BL_GSL
void BLManager::gsl_error_handler(const char *reason, const char *file,
					int lineno, int gsl_errno)
{
	// ignore underflows
	if(gsl_errno == GSL_EUNDRFLW) return;

	char tmp[256];
	snprintf(tmp,sizeof(tmp)-1,"file=%s line=%d gsl_errno=%d",file,lineno,
								gsl_errno);
	G4Exception("BLManager",reason,EventMustBeAborted,tmp);
}
#endif
