//	BLCMDparticlecolor.cc
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
//

#include <map>
#include <vector>

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4VisAttributes.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"

#include "BLAssert.hh"
#include "BLCommand.hh"
#include "BLCommandAlias.hh"
#include "BLManager.hh"
#include "BLParam.hh"

/**	class BLCMDparticlecolor -- implement the particlecolor command
 **/
class BLCMDparticlecolor : public BLCommand, public BLManager::TrackingAction,
						public BLManager::EventAction {
	static std::map<const G4ParticleDefinition*,const G4VisAttributes*> pd2va;
	const G4VisAttributes *plus;
	const G4VisAttributes *minus;
	const G4VisAttributes *neutral;
	const G4VisAttributes *reference;
public:
	/// Constructor.
	BLCMDparticlecolor();

	/// Destructor.
	~BLCMDparticlecolor() { }

	/// commandName() returns "particlecolor"
	G4String commandName()  { return "particlecolor"; }

	/// command() executes the command associated with this element.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command. 
	void defineNamedArgs() { }

	/// getVisAttributes()  will return the G4VisAttributes for a
	/// particle.
	const G4VisAttributes *getVisAttributes(const G4ParticleDefinition *pd) const;

	/// PreUserTrackingAction() from BLManager::TrackingAction
	void PreUserTrackingAction(const G4Track* track);

	/// PostUserTrackingAction from BLManager::TrackingAction
	void PostUserTrackingAction(const G4Track* track);

	/// BeginOfEventAction from BLManager::EventAction
	void BeginOfEventAction(const G4Event* event);

	/// EndOfEventAction from BLManager::EventAction
	void EndOfEventAction(const G4Event* event);
};
std::map<const G4ParticleDefinition*,const G4VisAttributes*> BLCMDparticlecolor::pd2va;

BLCMDparticlecolor defaultParticleColor;

BLCommandAlias aliasParticleColor("trackcolor",defaultParticleColor);


/// drawTrajectory() is a static function to draw a trajectory
/// (making it static is the simplest way to share it between classes
/// BLTrajectory and BLSaveTrajectory)
static void drawTrajectory(const G4VTrajectory *traj, const G4VisAttributes &va, G4int i_mode=0);


/**	class BLTrajectory implements particle-color drawing.
 *
 * 	Note this class is derived from G4Trajectory (not G4VTrajectory),
 * 	so it is efficient. But G4Trajectory cannot be saved across runs,
 * 	so there is also a BLSaveTrajectory class.
 *
 *	Note we need to supply operators new and delete, because G4Trajectory
 *	has versions of them specific to sizeof(G4Trajectory).
 **/
class BLTrajectory : public G4Trajectory {
	G4VisAttributes visAttrib;
public:
	/// Constructors
	BLTrajectory() : G4Trajectory() , visAttrib()
		{ } 
	BLTrajectory(const G4Track *track) : G4Trajectory(track) , visAttrib()
		{ const G4VisAttributes *va = defaultParticleColor.getVisAttributes(track->GetDefinition()); 
		  if(va) visAttrib = *va;
		}
	BLTrajectory(BLTrajectory &t) : G4Trajectory(t) , visAttrib(t.visAttrib)
		{ }

	/// Destructor
	~BLTrajectory() { }

	/// getVisAttributes() returns the vis attributes for the trajectory.
	G4VisAttributes& getVisAttributes() { return visAttrib; }

	/// DrawTrajectory will draw the trajectory using the selected color.
	void DrawTrajectory(G4int i_mode=0) const 
		{ drawTrajectory(this,visAttrib,i_mode); }

	/// Operators
	inline void* operator new(size_t ignored) {
		void* aTrajectory;
		aTrajectory = (void*)allocator.MallocSingle();
		return aTrajectory;
	}
	inline void  operator delete(void *p) {
		allocator.FreeSingle((BLTrajectory*)p);
	}

	static G4Allocator<BLTrajectory> allocator;
};
G4Allocator<BLTrajectory> BLTrajectory::allocator;

/**	class BLSaveTrajectory is used to save a trajectory from the Reference
 *	particle run to the visual run. Neither G4Trajectory nor BLTrajectory
 *	can be saved, due to their new/delete functions. The constructor here
 *	does a deep copy of the trajectory, including each point.
 *
 *	Much of this source is copied from geant4.6.2.p02.
 *
 *	We don't need operators new and delete, because this class derives
 *	from G4Vtrajectory (not G4Trajectory).
 **/
class BLSaveTrajectory : public G4VTrajectory {
	TrajectoryPointContainer* positionRecord;
	G4int                     fTrackID;
	G4int                     fParentID;
	G4int                     PDGEncoding;
	G4double                  PDGCharge;
	G4String                  ParticleName;
	G4ThreeVector             initialMomentum;
	G4VisAttributes visAttrib;
public:
	BLSaveTrajectory(G4VTrajectory *traj);
	~BLSaveTrajectory();
	G4bool operator == (const G4VTrajectory& right) const { ::abort(); }
	G4int GetTrackID() const { return fTrackID; }
	G4int GetParentID() const { return fParentID; }
	G4String GetParticleName() const { return ParticleName; }
	G4double GetCharge() const { return PDGCharge; }
	G4int GetPDGEncoding() const { return PDGEncoding; }
	G4ThreeVector GetInitialMomentum() const { return initialMomentum; }
	int GetPointEntries() const { return positionRecord->size(); }
	G4VTrajectoryPoint* GetPoint(G4int i) const 
		{ return (*positionRecord)[i]; }
	void DrawTrajectory(G4int i_mode=0) const
		{ drawTrajectory(this,visAttrib,i_mode); }
	void AppendStep(const G4Step* aStep) { ::abort(); }
	void MergeTrajectory(G4VTrajectory* secondTrajectory) { ::abort(); }
};

/// vector of Reference particle trajectories (to be added to the first event
/// during VISUAL tracking)
std::vector<BLSaveTrajectory*> ReferenceTrajectories;


BLCMDparticlecolor::BLCMDparticlecolor() : BLCommand(), BLManager::TrackingAction(), BLManager::EventAction()
{
	registerCommand(BLCMDTYPE_AUX);
	setSynopsis("Set the colors for particle tracks.");
	setDescription("Arguments are of the form 'name=1,1,0', where name is\n"
			"the standard name of a particle, and 1,1,0 is the\n"
			"R,G,B value desired for its color ('' for invisible)\n"
			"The special names plus, minus, and neutral will set\n"
			"colors for unnamed particles of each charge.\n"
			"The name reference will apply to the reference\n"
			"track (defaults to invisible).");
	
	G4String red("1,0,0"), green("0,1,0"), blue("0,0,1"), white("1,1,1");
	plus = getVisAttrib(blue);
	minus = getVisAttrib(red);
	neutral = getVisAttrib(green);
	reference = getVisAttrib(white);
}

int BLCMDparticlecolor::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	int retval = 0;

	if(argv.size() > 0) {
		printError("particlecolor: unnamed arguments are invalid.");
		retval = 1;
	}

	BLArgumentMap::iterator i;
	for(i=namedArgs.begin(); i!=namedArgs.end(); ++i) {
		G4String name = i->first;
		G4String value = i->second;
		const G4VisAttributes *va = getVisAttrib(value);
		if(name == "plus") {
			plus = va;
			continue;
		} else if(name == "minus") {
			minus = va;
			continue;
		} else if(name == "neutral") {
			neutral = va;
			continue;
		} else if(name == "reference") {
			reference = va;
			continue;
		}
		G4ParticleDefinition *pd = 
			G4ParticleTable::GetParticleTable()->FindParticle(name);
		if(!pd) {
			printError("particlecolor: particle '%s' not found",
								name.c_str());
			retval = 1;
			continue;
		}
		pd2va[pd] = va;
	}

	BLManager::getObject()->registerTrackingAction(this);
	BLManager::getObject()->registerEventAction(this,false);

	print("",namedArgs);

	return retval;
}

const G4VisAttributes *BLCMDparticlecolor::getVisAttributes(const G4ParticleDefinition *pd) const
{
	BLManagerState state = BLManager::getObject()->getState();
	if(state == REFERENCE)
		return reference;
	if(state != VISUAL)
		return 0;

	if(pd2va.count(pd) > 0)
		return pd2va[pd];
	G4double charge = pd->GetPDGCharge();
	if(charge > 0.0)
		return plus;
	else if(charge < 0.0)
		return minus;
	return neutral;
}

void BLCMDparticlecolor::PreUserTrackingAction(const G4Track* track)
{
	if(Param.getString("viewer") == "none")
		return;
	if(fpTrackingManager)
		fpTrackingManager->SetTrajectory(new BLTrajectory(track));
	BLManagerState state = BLManager::getObject()->getState();
	if(state == REFERENCE)
		fpTrackingManager->SetStoreTrajectory(true);
}

void BLCMDparticlecolor::PostUserTrackingAction(const G4Track* track)
{
	if(!fpTrackingManager) return;
	// get Reference particle trajectories
	BLManagerState state = BLManager::getObject()->getState();
	if(state == REFERENCE) {
		G4VTrajectory *t = fpTrackingManager->GimmeTrajectory();
		// this is a deep copy, so let the tracking manager delete t
		if(t) ReferenceTrajectories.push_back(new BLSaveTrajectory(t));
	}
}

void BLCMDparticlecolor::BeginOfEventAction(const G4Event* event)
{
}

void BLCMDparticlecolor::EndOfEventAction(const G4Event* event)
{
	BLManagerState state = BLManager::getObject()->getState();
	if(state == VISUAL && ReferenceTrajectories.size() > 0) {
		// add ReferenceTrajectories to this event
		G4TrajectoryContainer *tc = event->GetTrajectoryContainer();
		for(unsigned i=0; i<ReferenceTrajectories.size(); ++i) {
			if(tc)
				tc->push_back(ReferenceTrajectories[i]);
		}
		// clear ReferenceTrajectories so only one copy is drawn, and
		// because the tracking manager is going to delete them.
		ReferenceTrajectories.clear();
	}
}



BLSaveTrajectory::BLSaveTrajectory(G4VTrajectory *traj) : G4VTrajectory(),
								visAttrib()
{
	positionRecord = new TrajectoryPointContainer();
	fTrackID = traj->GetTrackID();
	fParentID = traj->GetParentID();
	PDGEncoding = traj->GetPDGEncoding();
	PDGCharge = traj->GetCharge();
	ParticleName = traj->GetParticleName();
	initialMomentum = traj->GetInitialMomentum();
	G4ParticleDefinition *pd = 
		G4ParticleTable::GetParticleTable()->FindParticle(PDGEncoding);
	const G4VisAttributes *va = defaultParticleColor.getVisAttributes(pd);
	if(va) visAttrib = *va;

	unsigned n = traj->GetPointEntries();
	for(unsigned i=0; i<n; ++i) {
		G4TrajectoryPoint *p = 
			dynamic_cast<G4TrajectoryPoint*>(traj->GetPoint(i));
		BLAssert(p != 0);
		G4TrajectoryPoint *q = new G4TrajectoryPoint(*p);
		positionRecord->push_back(q);
	}
}
BLSaveTrajectory::~BLSaveTrajectory()
{
	unsigned n = positionRecord->size();
	for(unsigned i=0; i<n; ++i)
		delete GetPoint(i);
	positionRecord->clear();
	delete positionRecord;
}




// copied from geant4 5.0 source/tracking/src/G4VTrajectory.cc
// and the colors were then changed.
static void drawTrajectory(const G4VTrajectory *traj,
					const G4VisAttributes &va, G4int i_mode)
{
  // If i_mode>=0, draws a trajectory as a polyline (blue for
  // positive, red for negative, green for neutral) and, if i_mode!=0,
  // adds markers - yellow circles for step points and magenta squares
  // for auxiliary points, if any - whose screen size in pixels is
  // given by abs(i_mode)/1000.  E.g: i_mode = 5000 gives easily
  // visible markers.

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager) return;

  const G4double markerSize = abs(i_mode)/1000;
  G4bool lineRequired (i_mode >= 0);
  G4bool markersRequired (markerSize > 0.);

  G4Polyline trajectoryLine;
  G4Polymarker stepPoints;
  G4Polymarker auxiliaryPoints;

  for (G4int i = 0; i < traj->GetPointEntries() ; i++) {
    G4VTrajectoryPoint* aTrajectoryPoint = traj->GetPoint(i);
    const std::vector<G4ThreeVector>* auxiliaries
      = aTrajectoryPoint->GetAuxiliaryPoints();
    if (auxiliaries) {
      for (size_t iAux = 0; iAux < auxiliaries->size(); ++iAux) {
	const G4ThreeVector pos((*auxiliaries)[iAux]);
	if (lineRequired) {
	  trajectoryLine.push_back(pos);
	}
	if (markersRequired) {
	  auxiliaryPoints.push_back(pos);
	}
      }
    }
    const G4ThreeVector pos(aTrajectoryPoint->GetPosition());
    if (lineRequired) {
      trajectoryLine.push_back(pos);
    }
    if (markersRequired) {
      stepPoints.push_back(pos);
    }
  }

  if (lineRequired) {
    if(!va.IsVisible()) return;
    trajectoryLine.SetVisAttributes(va);
    pVVisManager->Draw(trajectoryLine);
  }
  if (markersRequired) {
    auxiliaryPoints.SetMarkerType(G4Polymarker::squares);
    auxiliaryPoints.SetScreenSize(markerSize);
    auxiliaryPoints.SetFillStyle(G4VMarker::filled);
    G4VisAttributes auxiliaryPointsAttribs(G4Colour(0.,1.,1.));  // Magenta
    auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);
    pVVisManager->Draw(auxiliaryPoints);

    stepPoints.SetMarkerType(G4Polymarker::circles);
    stepPoints.SetScreenSize(markerSize);
    stepPoints.SetFillStyle(G4VMarker::filled);
    G4VisAttributes stepPointsAttribs(G4Colour(1.,1.,0.));  // Yellow.
    stepPoints.SetVisAttributes(&stepPointsAttribs);
    pVVisManager->Draw(stepPoints);
  }
}
