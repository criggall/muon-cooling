//	BLCMDemfactor.cc
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

#include <set>

#include "G4WrapperProcess.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleChangeForMSC.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4VMultipleScattering.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4TransportationManager.hh"

#include "BLManager.hh"
#include "BLCommand.hh"
#include "BLCallback.hh"
#include "BLAssert.hh"

#include "fnmatch.h"

/*	Reverse-engineering muMsc and muIoni (background for the overall design)
 *
 *	For all charged particles, the first three entries in the process list
 *	are fixed:
 *		[0] Transportation
 *		[1] multiple scattering
 *		[2] ionization energy loss (incl generate delta-rays)
 *
 *	These routines change the following quantities (listed in order called):
 *	Routine				momDir	Pos	KE
 *	--------------------		-----	---	--
 *	Transp::AlongStepDoIt		Y	Y	Y
 *	muMsc::AlongStepDoIt		N/Y	N/?	N
 *	muIoni::AlongStepDoIt		N	N	Y
 *	Transp::PostStepDoit		N	N	N
 *	muMsc::PostStepDoit		Y	Y	N
 *	muIoni::PostStepDoit		Y	N	Y
 *	(None of the GPIL routines change anything.)
 *	(Two answers are for Geant4 releases 9.5/9.6.p02)
 *
 *	Note that inside AlongStepDoit(track,step), the track is from the
 *	start of the step, while the "current" track is in 
 *	step->GetPostStepPoint()->GetTrack().
 *
 *	Inside PostStepDoit(track,step) however, the track is updated
 *	after each processes' call, as is step->GetPostStepPoint()->GetTrack().
 *
 *	Handle Pos and momentumDirection in WrapMsc AlongStepDoit() and
 *	PostStepDoIt(), and WrapEloss PostStepDoit. Handle KE in WrapEloss
 *	AlongStepDoit and PostStepDoit.
 *
 *	Further complication: muMsc sets condition=Forced, but muIoni sets
 *	condition=NotForced -- so muIoni::PostStepDoit is not always called.
 *
 *	We cannot access the copy constructor of G4particleChange, so we
 *	cannot copy BLCMDemfactor; just use defaultBLCMDemfactor.
 */

static void dumpProcesses(const G4ParticleDefinition *pd, const char *tag)
{
	printf("\nparticle %s processes at %s\n",pd->GetParticleName().c_str(),
								tag);
	G4ProcessManager *pm = pd->GetProcessManager();
	G4ProcessVector *pv = pm->GetProcessList();
	for(int j=0; j<pv->size(); ++j) {
		G4VProcess *ppp = (*pv)[j];
		printf("j=%d %s\n",j,ppp->GetProcessName().c_str());
	}
}

/**	class WrapMsc wraps all multiple-scattering processes.
 **/
class WrapMsc : public G4WrapperProcess {
	double msc_ratio;
	double eLoss_ratio;
	double deltaRay_ratio;
	G4SafetyHelper *safetyHelper;
public:
	WrapMsc() : G4WrapperProcess("emfactor_")
		{ msc_ratio = eLoss_ratio = deltaRay_ratio = 1.0; 
		  safetyHelper = 0; }

	void setRatio(G4double msc, G4double eLoss, G4double deltaRay) 
		{ msc_ratio = msc; eLoss_ratio = eLoss; 
		  deltaRay_ratio = deltaRay;  }

	virtual void RegisterProcess(G4VProcess *p)
		{ if(p) theProcessName = "emfactor_";
		  if(p) G4WrapperProcess::RegisterProcess(p);
		  else pRegProcess = 0; }
	virtual void StartTracking(G4Track *track) 
		{ pRegProcess->StartTracking(track); } 
	virtual void EndTracking() 
		{ pRegProcess->EndTracking(); }
	virtual G4bool IsApplicable(const G4ParticleDefinition &pd) 
		{ BLAssert(false); return false; }
	virtual void BuildPhysicsTable(const G4ParticleDefinition &pd) 
		{ BLAssert(false); }
	virtual void PreparePhysicsTable(const G4ParticleDefinition &pd) 
		{ BLAssert(false);}
	virtual G4bool StorePhysicsTable(const G4ParticleDefinition* ,
	                                 const G4String& directory,
					 G4bool          ascii = false) 
		{ BLAssert(false); return false; }
	virtual G4bool RetrievePhysicsTable( const G4ParticleDefinition* ,
	                                     const G4String& directory,
					     G4bool          ascii = false) 
		{ BLAssert(false); return false; }
	virtual void ResetNumberOfInteractionLengthLeft() 
		{ pRegProcess->ResetNumberOfInteractionLengthLeft(); }

	virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double  previousStepSize,
                             G4double  currentMinimumStep,
                             G4double& proposedSafety,
                             G4GPILSelection* selection);

	virtual G4VParticleChange* AlongStepDoIt(
                             const G4Track& track,
			     const G4Step& stepData);

	virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition);

	virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& track,
                             const G4Step&  stepData);

	virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition);

	virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& track,
                             const G4Step& stepData);
};


/**	class WrapEloss wraps all Eloss processes
 **/
class WrapEloss : public G4WrapperProcess {
	double msc_ratio;
	double eLoss_ratio;
	double deltaRay_ratio;
	G4ParticleChange myParticleChange;
	G4SafetyHelper *safetyHelper;
public:
	WrapEloss() : G4WrapperProcess("emfactor_Eloss"), myParticleChange()
		{ msc_ratio = eLoss_ratio = deltaRay_ratio = 1.0; 
		  safetyHelper = 0; }

	void setRatio(G4double msc, G4double eLoss, G4double deltaRay) 
		{ msc_ratio = msc; eLoss_ratio = eLoss; 
		  deltaRay_ratio = deltaRay;  }

	virtual void RegisterProcess(G4VProcess *p)
		{ if(p) theProcessName = "emfactor_";
		  if(p) G4WrapperProcess::RegisterProcess(p);
		  else pRegProcess = 0; }
	virtual void StartTracking(G4Track *track) 
		{ pRegProcess->StartTracking(track); } 
	virtual void EndTracking() 
		{ pRegProcess->EndTracking(); }
	virtual G4bool IsApplicable(const G4ParticleDefinition &pd) 
		{ return pRegProcess->IsApplicable(pd); }
	virtual void BuildPhysicsTable(const G4ParticleDefinition &pd) 
		{ BLAssert(false); }
	virtual void PreparePhysicsTable(const G4ParticleDefinition &pd) 
		{ BLAssert(false);}
	virtual G4bool StorePhysicsTable(const G4ParticleDefinition* ,
	                                 const G4String& directory,
					 G4bool          ascii = false) 
		{ BLAssert(false); return false; }
	virtual G4bool RetrievePhysicsTable( const G4ParticleDefinition* ,
	                                     const G4String& directory,
					     G4bool          ascii = false) 
		{ BLAssert(false); return false; }
	virtual void      ResetNumberOfInteractionLengthLeft() 
		{ pRegProcess->ResetNumberOfInteractionLengthLeft(); }

	virtual G4double AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double  previousStepSize,
                             G4double  currentMinimumStep,
                             G4double& proposedSafety,
                             G4GPILSelection* selection);

	virtual G4VParticleChange* AlongStepDoIt(
                             const G4Track& track,
			     const G4Step& stepData);

	virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition);

	virtual G4VParticleChange* PostStepDoIt(
                             const G4Track& track,
                             const G4Step&  stepData);

	virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition);

	virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& track,
                             const G4Step& stepData);
};

/**	class BLCMDemfactor will re-weight multiple scattering and eLoss.
 *
 *	This command permits the user to artificially decrease (or increase)
 *	the amount of multiple scattering and/or the energy loss fluctuations
 *	in materials (the mean energy loss is not affected).
 *
 *	The command works by wrapping the multiple-scattering and energy-loss
 *	processes for the selected particles. This wrapping must be deferred
 *	until PreUserTrackingAction(), so the regular processes get properly
 *	initialized. They must then be un-wrapped in PostUserTrackingAction().
 **/
class BLCMDemfactor : public BLCommand, public BLCallback,
					public BLManager::TrackingAction {
	G4String particle;
	G4double msc;
	G4double eLoss;
	G4double deltaRay;
	WrapMsc wrapMsc;
	WrapEloss wrapEloss;
	std::set<const G4ParticleDefinition*> particleDefs;
	std::set<const G4ParticleDefinition*> warnings;
	bool first;
	// internal utility functions
	bool isMatchingParticle(G4ParticleDefinition *pd)
		{ return matchList(pd->GetParticleName(),particle); }
	void replaceProcess(G4ProcessManager *pmgr, const G4VProcess *orgProc,
						const G4VProcess *newProc);
public:
	/// Constructor.
	BLCMDemfactor();

	/// commandName() returns "emfactor".
	virtual G4String commandName() { return "emfactor"; }
	
	/// command() implements the emfactor command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// callback() from BLCallback.
	void callback(int type);

	// from BLManager::TrackingAction
	virtual void PreUserTrackingAction(const G4Track *track);
	virtual void PostUserTrackingAction(const G4Track *track);
};
BLCMDemfactor defaultBLCMDemfactor;

BLCMDemfactor::BLCMDemfactor() : BLCommand(), BLCallback()
{
	registerCommand(BLCMDTYPE_PHYSICS);
	setSynopsis("Multiply multiple scattering and energy loss by factors.");
	setDescription("The multiple scattering and ionization energy loss "
	"processes are not really independent; the msc factor applies to "
	"scattering in both, and the eLoss factor applies to the fluctuations "
	"of energy loss in both. The 'mean' energy loss is not affected; "
	"note that in some EM models, including the one used for muons, the "
	"fluctuations only increase the loss from the 'mean', so the average "
	"energy loss is affected, but only slightly.\n\n"
	"The deltaRay parameter is the probability of keeping any given "
	"delta-ray. The 3-momentum of each discarded delta-ray is added back "
	"into the particle's 3-momentum, and its K.E. is recomputed.\n\n"
	"Note that with msc=0 tracks can still be scattered if deltaRay>0.\n\n"
	"This command applies ONLY to multiple scattering and ionization "
	"energy loss. Other processes still apply, which can scatter and cause "
	"a reduction of the track's energy. For muons you may wish to disable "
	"these processes (in the physics command): muBrems, muPairProd, "
	"CoulombScat, Decay; similarly for other particles.\n\n"
	"Problems may occur if any factor is > 1.0; in particular, the "
	"position may be wrong - it might even cross a volume boundary "
	"illegally, causing major tracking errors (code is optimized for "
	"0.0 <= msc,eLoss,DeltaRay <= 1.0)."
	);

	particle = "";
	msc = 1.0;
	eLoss = 1.0;
	deltaRay = 1.0;
	first = true;
}

int BLCMDemfactor::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(!first) G4Exception("emfactor", "Only one instance permitted",
							FatalException, "");
	first = false;

	int retval = handleNamedArgs(namedArgs);

	wrapMsc.setRatio(msc,eLoss,deltaRay);
	wrapEloss.setRatio(msc,eLoss,deltaRay);

	// too early to enumerate particles or processes here
	BLManager::getObject()->registerCallback(this,0);

	// do the wrapping/unwrapping in TrackingAction
	BLManager::getObject()->registerTrackingAction(this);

	print("");

	return retval;
}

void BLCMDemfactor::defineNamedArgs()
{
	argString(particle,"particle","Comma-separated list of particle "
					"names or patterns (default=none).");
	argDouble(msc,"msc","Ratio for multiple scattering (1.0).");
	argDouble(eLoss,"eLoss","Ratio for ionization energy loss fluctuations (1.0).");
	argDouble(deltaRay,"deltaRay","Probability to keep delta-rays (1.0).");
	argDouble(eLoss,"Eloss","Synonym for eLoss.");
	argDouble(eLoss,"eloss","Synonym for eLoss.");
}

void BLCMDemfactor::callback(int type)
{
	if(type == 0) {
	    // fill particleDefs (avoid string comparisons during tracking)
	    G4ParticleTable::G4PTblDicIterator *myParticleIterator =
			G4ParticleTable::GetParticleTable()->GetIterator();
	    myParticleIterator->reset();
	    while((*myParticleIterator)()) {
		G4ParticleDefinition *pd = myParticleIterator->value();
		if(!isMatchingParticle(pd) || pd->GetPDGCharge() == 0.0 ||
		   pd->IsShortLived())
			continue;
		particleDefs.insert(pd);
	    }
	}

	return;
}

void BLCMDemfactor::PreUserTrackingAction(const G4Track *track)
{
	BLAssert(wrapMsc.GetRegisteredProcess() == 0);
	BLAssert(wrapEloss.GetRegisteredProcess() == 0);

	// check if we apply to this track
	const G4ParticleDefinition *pd = track->GetParticleDefinition();
	if(particleDefs.count(pd) == 0) {
		return;
	}
dumpProcesses(pd,"PreUserTrackingAction start");
	// wrap msc and Eloss processes
	G4ProcessManager *pm = pd->GetProcessManager();
	G4ProcessVector *pv = pm->GetProcessList();
	int nMsc=0, nEloss=0;
	for(int j=1; j<=2; ++j) { // just *Msc and *Ioni
G4VProcess *ppp = (*pv)[j];
printf("PreUserTrackingAction j=%d %s\n",j,ppp->GetProcessName().c_str());
		G4VMultipleScattering *p =
				dynamic_cast<G4VMultipleScattering*>((*pv)[j]);
		if(p != 0) {
printf("wrap msc\n");
			wrapMsc.RegisterProcess(p);
			replaceProcess(pm,p,&wrapMsc);
			++nMsc;
		}
		G4VEnergyLossProcess *q =
				dynamic_cast<G4VEnergyLossProcess*>((*pv)[j]);
		if(q != 0) {
printf("wrap eLoss\n");
			wrapEloss.RegisterProcess(q);
			replaceProcess(pm,q,&wrapEloss);
			++nEloss;
		}
	}
	if(nMsc != 1)
		G4Exception("emfactor", "Invalid # msc processes",
				FatalException, pd->GetParticleName());
	if(nEloss != 1)
		G4Exception("emfactor", "Invalid # Eloss processes",
				FatalException, pd->GetParticleName());
	if(warnings.count(pd) == 0) {
		warnings.insert(pd);
		G4Exception("emfactor", "Wrapped msc process",
				JustWarning,pd->GetParticleName());
		G4Exception("emfactor", "Wrapped Eloss process",
				JustWarning,pd->GetParticleName());
	}
dumpProcesses(pd,"PreUserTrackingAction end");
}

void BLCMDemfactor::PostUserTrackingAction(const G4Track *track)
{
	const G4ParticleDefinition *pd = track->GetParticleDefinition();
	G4ProcessManager *pm = pd->GetProcessManager();
dumpProcesses(pd,"PostUserTrackingAction start");
	const G4VProcess *p = wrapMsc.GetRegisteredProcess();
	if(p != 0) {
printf("unwrap msc\n");
		//#G4ProcessManager *pm = 
		//#	(G4ProcessManager*)wrapMsc.GetProcessManager();
		replaceProcess(pm,&wrapMsc,p);
		wrapMsc.RegisterProcess(0);
	}

	p = wrapEloss.GetRegisteredProcess();
	if(p != 0) {
printf("unwrap eLoss\n");
		//#G4ProcessManager *pm = 
		//#	(G4ProcessManager*)wrapEloss.GetProcessManager();
		replaceProcess(pm,&wrapEloss,p);
		wrapEloss.RegisterProcess(0);
	}
dumpProcesses(pd,"PostUserTrackingAction end");
}

void BLCMDemfactor::replaceProcess(G4ProcessManager *pmgr,
		const G4VProcess *orgProc, const G4VProcess *newProc)
{
printf("replaceProcess pmgr=%p orgProc=%p newProc=%p\n",pmgr,orgProc,newProc);
	// loop over all process vectors in the process manager
	for(int i=0; i<=6; ++i) {
		G4ProcessVector *pv = 0;
		switch(i) {
		case 0:	continue; // not present in AtRestGPIL
		case 1:	continue; // not present in AtRestDoIt
		case 2:	pv = pmgr->GetProcessVector(idxAlongStep,typeGPIL);
			break;
		case 3:	pv = pmgr->GetProcessVector(idxAlongStep,typeDoIt);
			break;
		case 4:	pv = pmgr->GetProcessVector(idxPostStep,typeGPIL);
			break;
		case 5:	pv = pmgr->GetProcessVector(idxPostStep,typeDoIt);
			break;
		case 6:	pv = pmgr->GetProcessList();
			break;
		}
		int k=0;
		for(int j=0; j<pv->size(); ++j) {
			if((*pv)[j] == orgProc) {
				(*pv)[j] = (G4VProcess *)newProc;
				++k;
printf("replaceProcess i=%d j=%d k=%d\n",i,j,k);
			}
		}
		BLAssert(k <= 1);
	}
}

G4double WrapMsc::AlongStepGetPhysicalInteractionLength(
			const G4Track& track, G4double  previousStepSize,
			G4double  currentMinimumStep, G4double& proposedSafety,
			G4GPILSelection* selection)
{
	G4double value = pRegProcess->AlongStepGetPhysicalInteractionLength(
			track,previousStepSize,currentMinimumStep,
			proposedSafety,selection);
	return value;
}

G4VParticleChange* WrapMsc::AlongStepDoIt(const G4Track& track,
						const G4Step& stepData)
{
	G4ThreeVector initialPos = stepData.GetPostStepPoint()->GetPosition();
	G4ThreeVector initialMD =
			stepData.GetPostStepPoint()->GetMomentumDirection();

	G4ParticleChangeForMSC *pc = dynamic_cast<G4ParticleChangeForMSC*>(
			pRegProcess->AlongStepDoIt(track,stepData));

	BLAssert(pc != 0);
	BLAssert(pc->GetNumberOfSecondaries() == 0);
G4ThreeVector procPos = *pc->GetPosition();
G4ThreeVector procMD = *pc->GetMomentumDirection();

	// get new momentum direction
	G4ThreeVector orgMD = *pc->GetMomentumDirection();
	BLAssert(fabs(orgMD.mag2()-1.0) < 1.0e12);
	G4ThreeVector newMD = orgMD;
printf("WrapMsc::AlongStepDoIt: orgMD: %.6f,%.6f,%.6f initialMD:%.6f,%.6f,%.6f\n",orgMD[0],orgMD[1],orgMD[2],initialMD[0],initialMD[1],initialMD[2]);
	double angle = initialMD.dot(orgMD); // can sometimes be 1.000000000001
	angle = (fabs(angle) < 1.0 ? acos(angle) : 0.0);
printf("WrapMsc::AlongStepDoIt: angle=%.6f\n",angle);
	if(fabs(angle) > 1E-6) {
		G4ThreeVector tmp = initialMD;
		G4ThreeVector axis = tmp.cross(orgMD).unit(); // modifies tmp
		tmp = initialMD;
		newMD = tmp.rotate(angle*msc_ratio,axis); // modifies tmp
printf("WrapMsc::AlongStepDoIt: orgMD: %.6f,%.6f,%.6f newMD:%.6f,%.6f,%.6f\n",orgMD[0],orgMD[1],orgMD[2],newMD[0],newMD[1],newMD[2]);
	}
	BLAssert(fabs(newMD.mag2()-1.0) < 1.0e12);
	pc->ProposeMomentumDirection(newMD);

	return pc;
}

G4double WrapMsc::PostStepGetPhysicalInteractionLength(const G4Track& track,
		G4double previousStepSize, G4ForceCondition* condition)
{
	G4double value = pRegProcess->PostStepGetPhysicalInteractionLength(
					track,previousStepSize,condition);
	return value;
}

G4VParticleChange* WrapMsc::PostStepDoIt(const G4Track& track,
						const G4Step&  stepData)
{
	if(safetyHelper == 0) {
	    safetyHelper = G4TransportationManager::GetTransportationManager()
							->GetSafetyHelper();
	    safetyHelper->InitialiseHelper();
	}

	G4ThreeVector initialPos = stepData.GetPostStepPoint()->GetPosition();
	G4ThreeVector initialMD =
			stepData.GetPostStepPoint()->GetMomentumDirection();

	G4ParticleChangeForMSC *pc = dynamic_cast<G4ParticleChangeForMSC*>(
			pRegProcess->PostStepDoIt(track,stepData));

	BLAssert(pc != 0);
	BLAssert(pc->GetNumberOfSecondaries() == 0);
G4ThreeVector procPos = *pc->GetPosition();
G4ThreeVector procMD = *pc->GetMomentumDirection();

	// get new momentum direction
	G4ThreeVector orgMD = *pc->GetMomentumDirection();
	BLAssert(fabs(orgMD.mag2()-1.0) < 1.0e12);
	G4ThreeVector newMD = orgMD;
printf("WrapMsc::PostStepDoIt: orgMD: %.6f,%.6f,%.6f initialMD:%.6f,%.6f,%.6f\n",orgMD[0],orgMD[1],orgMD[2],initialMD[0],initialMD[1],initialMD[2]);
	double angle = initialMD.dot(orgMD); // can sometimes be 1.000000000001
	angle = (fabs(angle) < 1.0 ? acos(angle) : 0.0);
printf("WrapMsc::PostStepDoIt: angle=%.6f\n",angle);
	if(fabs(angle) > 1E-6) {
		G4ThreeVector tmp = initialMD;
		G4ThreeVector axis = tmp.cross(orgMD).unit(); // modifies tmp
		tmp = initialMD;
		newMD = tmp.rotate(angle*msc_ratio,axis); // modifies tmp
printf("WrapMsc::PostStepDoIt: orgMD: %.6f,%.6f,%.6f newMD:%.6f,%.6f,%.6f\n",orgMD[0],orgMD[1],orgMD[2],newMD[0],newMD[1],newMD[2]);
	}
	BLAssert(fabs(newMD.mag2()-1.0) < 1.0e12);
	pc->ProposeMomentumDirection(newMD);

	// get new position, staying within the current PhysicalVolume
	// NOTE: in principle this can have an error in the scaled position;
	// in practice, the Geant4 msc model has already made this error:
	// Near a geometrical boundary, the displacement is limited by the
	// safety (which can be 0.0 right on a volume boundary).
	// As the msc process has already limited the jump from un-scattered
	// to scattered position, it is guaranteed that initialPos is within
	// the safety of orgPos. Note the msc process limits the step length
	// to keep this error reasonable.
	// Using a smaller minStep will reduce the error.
	// This code is modeled after G4VMscModel::ComputeDisplacement().
	G4double safety = stepData.GetPostStepPoint()->GetSafety();
	G4ThreeVector orgPos = *pc->GetPosition();
	G4ThreeVector deltaPos = (initialPos - orgPos) * (1.0 - msc_ratio);
	if(deltaPos.mag() > safety)
		safety = safetyHelper->ComputeSafety(orgPos);
	if(deltaPos.mag() > safety) {
		deltaPos = deltaPos.unit() * safety * 0.999;
	}
	G4ThreeVector newPos = orgPos + deltaPos;
	safetyHelper->ReLocateWithinVolume(newPos);
	pc->ProposePosition(newPos);

G4ThreeVector finalPos = *pc->GetPosition();
G4ThreeVector finalMD = *pc->GetMomentumDirection();
printf("WrapMsc::PostStepDoIt: procPos: %.3f,%.3f,%.3f finalPos:%.3f,%.3f,%.3f\n",procPos[0],procPos[1],procPos[2],finalPos[0],finalPos[1],finalPos[2]);
printf("WrapMsc::PostStepDoIt: procMD: %.3f,%.3f,%.3f finalMD:%.3f,%.3f,%.3f\n",procMD[0],procMD[1],procMD[2],finalMD[0],finalMD[1],finalMD[2]);

	return pc;
}

G4double WrapMsc::AtRestGetPhysicalInteractionLength(
		const G4Track& track, G4ForceCondition* condition)
{
	G4double value = pRegProcess->
			AtRestGetPhysicalInteractionLength(track,condition);

	return value;
}

G4VParticleChange* WrapMsc::AtRestDoIt(const G4Track& track,
                             			const G4Step& stepData)
{
	G4VParticleChange *pc = pRegProcess->AtRestDoIt(track,stepData);
	return pc;
}

G4double WrapEloss::AlongStepGetPhysicalInteractionLength(
			const G4Track& track, G4double  previousStepSize,
			G4double  currentMinimumStep, G4double& proposedSafety,
			G4GPILSelection* selection)
{
	G4double value = pRegProcess->AlongStepGetPhysicalInteractionLength(
			track,previousStepSize,currentMinimumStep,
			proposedSafety,selection);
	return value;
}

G4VParticleChange* WrapEloss::AlongStepDoIt(const G4Track& track,
						const G4Step& stepData)
{
	G4ParticleChangeForLoss *pc=0;

	// We cannot get inside the eLoss process to isolate the fluctuation,
	// so call it twice, once without fluctuations, and once with
	// fluctuations. In principle AlongStepDoIt could change things in
	// the original process and screw this up; in practice this works...
	// Moreover, only KE is affected by AlongStepDoIt.
	G4VEnergyLossProcess *proc = 
			dynamic_cast<G4VEnergyLossProcess*>(pRegProcess);
	BLAssert(proc != 0);

	if(fabs(eLoss_ratio-1.0) < 1E-12) {
		// shortcut if 1.0 -- makes testing msc easier
		pc = dynamic_cast<G4ParticleChangeForLoss*>(
				   	proc->AlongStepDoIt(track,stepData));
	} else {
		proc->SetLossFluctuations(false);
		pc = dynamic_cast<G4ParticleChangeForLoss*>(
				   	proc->AlongStepDoIt(track,stepData));
		BLAssert(pc != 0);
		G4double mean = pc->GetProposedKineticEnergy();
		proc->SetLossFluctuations(true);
		pc = dynamic_cast<G4ParticleChangeForLoss*>(
				   	proc->AlongStepDoIt(track,stepData));
		BLAssert(pc != 0);
		G4double ke = pc->GetProposedKineticEnergy();
		ke = mean + (ke - mean) * eLoss_ratio;
		pc->SetProposedKineticEnergy(ke);
	}

	return pc;
}

G4double WrapEloss::PostStepGetPhysicalInteractionLength(const G4Track& track,
		G4double previousStepSize, G4ForceCondition* condition)
{
	G4double value = pRegProcess->PostStepGetPhysicalInteractionLength(
					track,previousStepSize,condition);

	return value;
}

G4VParticleChange* WrapEloss::PostStepDoIt(const G4Track& track,
						const G4Step&  stepData)
{
	// shortcut makes testing msc easier
	if(fabs(deltaRay_ratio-1.0) < 1E-12) {
		return pRegProcess->PostStepDoIt(track,stepData);
	}

	myParticleChange.Initialize(track);

	if(safetyHelper == 0) {
	    safetyHelper = G4TransportationManager::GetTransportationManager()
							->GetSafetyHelper();
	    safetyHelper->InitialiseHelper();
	}

	G4ThreeVector initialPos = stepData.GetPostStepPoint()->GetPosition();

	if(stepData.GetPostStepPoint()->GetProcessDefinedStep() == this) {
		G4ParticleChangeForLoss *pc =
				dynamic_cast<G4ParticleChangeForLoss*>(
				pRegProcess->PostStepDoIt(track,stepData));
		BLAssert(pc != 0);
		myParticleChange.ProposeEnergy(pc->GetProposedKineticEnergy());
		myParticleChange.ProposeMomentumDirection(
					pc->GetProposedMomentumDirection());
		myParticleChange.ProposePolarization(
					pc->GetProposedPolarization());
		// copy secondaries, applying deltaRay_ratio, and adding
		// the momentum of omitted secondaries back into deltaMom
		G4int n=pc->GetNumberOfSecondaries();
		G4ThreeVector deltaMom(0.0,0.0,0.0);
		myParticleChange.SetNumberOfSecondaries(n);
		for(G4int i=0; i<n; ++i) {
			const G4Track *t = pc->GetSecondary(i);
			if(G4UniformRand() < deltaRay_ratio)
				myParticleChange.AddSecondary(new G4Track(*t));
			else
				deltaMom += t->GetMomentum();
		}
		// add deltaMom back into the track (myParticleChange)
		if(deltaMom.mag2() != 0.0) {
			G4double ke = myParticleChange.GetEnergy();
			G4double mass = myParticleChange.GetMass();
			G4ThreeVector md =
				*myParticleChange.GetMomentumDirection();
			G4ThreeVector ptot = md * sqrt(ke*ke+2*ke*mass);
			ptot += deltaMom;
			ke = sqrt(ptot.mag2()+mass*mass) - mass;
			myParticleChange.ProposeEnergy(ke);
			myParticleChange.ProposeMomentumDirection(ptot.unit());
		}
		// remove secondaries from original particleChange
		pc->SetNumberOfSecondaries(n); // yes, it deletes them
	}

	// Eloss fluctuations have already been handled in AlongStepDoIt

	// get new position, staying within the current PhysicalVolume
	// NOTE: in principle this can have an error in the scaled position;
	// in practice, the Geant4 msc model has already made this error:
	// Near a geometrical boundary, the displacement is limited by the
	// safety (which can be 0.0 right on a volume boundary).
	// As the msc process has already limited the jump from un-scattered
	// to scattered position, it is guaranteed that initialPos is within
	// the safety of orgPos. Note the msc process limits the step length
	// to keep this error reasonable.
	// Using a smaller minStep will reduce the error.
	// This code is modeled after G4VMscModel::ComputeDisplacement().
	G4double safety = stepData.GetPostStepPoint()->GetSafety();
	G4ThreeVector orgPos = *myParticleChange.GetPosition();
	G4ThreeVector deltaPos = (initialPos - orgPos) * (1.0 - msc_ratio);
	if(deltaPos.mag() > safety)
		safety = safetyHelper->ComputeSafety(orgPos);
	if(deltaPos.mag() > safety) {
		deltaPos = deltaPos.unit() * safety * 0.999;
	}
	G4ThreeVector newPos = orgPos + deltaPos;
	safetyHelper->ReLocateWithinVolume(newPos);
	myParticleChange.ProposePosition(newPos);

	return &myParticleChange;
}

G4double WrapEloss::AtRestGetPhysicalInteractionLength(
		const G4Track& track, G4ForceCondition* condition)
{
	G4double value = pRegProcess->
			AtRestGetPhysicalInteractionLength(track,condition);

	return value;
}

G4VParticleChange* WrapEloss::AtRestDoIt(const G4Track& track,
                             			const G4Step& stepData)
{
	G4VParticleChange *pc = pRegProcess->AtRestDoIt(track,stepData);
	return pc;
}

