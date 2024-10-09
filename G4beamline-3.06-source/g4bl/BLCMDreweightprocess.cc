//	BLCMDreweightprocess.cc
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

#include <vector>
#include <map>
#include <stdio.h>

#include "G4WrapperProcess.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"
#include "G4Track.hh"

#include "BLManager.hh"
#include "BLCommand.hh"
#include "BLCallback.hh"
#include "BLAssert.hh"
#include "BLParam.hh"

#include "fnmatch.h"

#define PREFIX "Reweight_"

static int DebugFlag=0;

/**	class BLCMDreweightprocess will re-weight physics processes.
 *
 *	This command permits the user to artificially increase (or decrease)
 *	the cross-section of one or more physics processes applied to one or
 *	more particles. This can be used to reduce the variance of a 
 *	simulation intended to measure a rare process. The track weights are
 *	adjusted accordingly.
 *
 *	For example, in a simulation with 1E6 beam particles, essentially
 *	none will get through an absorber that is 10 interaction lengths
 *	long. To determine how many make it through for a beam of 1E12
 *	particles, the interaction cross-section must be reduced.
 **/
class BLCMDreweightprocess : public BLCommand, public BLCallback {
	G4String particle;
	G4String process;
	G4double ratio;
	// internal utility functions
	bool isMatchingParticle(G4ParticleDefinition *pd)
	    { return particle.size()==0 || 
	    			matchList(pd->GetParticleName(),particle); }
	bool isMatchingProcess(const G4VProcess *proc)
	    { return matchList(proc->GetProcessName(),process); }
	bool wrapProcess(G4ProcessManager *pmgr, G4VProcess *proc);
public:
	/// Constructor.
	BLCMDreweightprocess();

	/// commandName() returns "reweightprocess".
	virtual G4String commandName() { return "reweightprocess"; }
	
	/// command() implements the reweightprocess command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// callback() from BLCallback.
	void callback(int type);
};
BLCMDreweightprocess defaultBLCMDreweightprocess;

/**	class ReweightProcess does the re-weighting.
 **/
class ReweightProcess : public G4WrapperProcess {
	double ratio;
	G4ProcessManager *pmgr;
public:
	ReweightProcess(G4ProcessManager *_pmgr, G4VProcess *proc, 
								double _ratio);

	virtual void StartTracking(G4Track *track) {
		currentInteractionLength = -1.0;
		G4VProcess::ResetNumberOfInteractionLengthLeft();
		pRegProcess->StartTracking(track);
	}

	virtual void EndTracking() {
		G4VProcess::ClearNumberOfInteractionLengthLeft();
		currentInteractionLength = -1.0;
		pRegProcess->EndTracking();
	}
	virtual G4bool IsApplicable(const G4ParticleDefinition &pd) { return pRegProcess->IsApplicable(pd); }
	virtual void BuildPhysicsTable(const G4ParticleDefinition &pd) { pRegProcess->BuildPhysicsTable(pd); }
	virtual void PreparePhysicsTable(const G4ParticleDefinition &pd) { pRegProcess->PreparePhysicsTable(pd);}
	virtual G4bool StorePhysicsTable(const G4ParticleDefinition* ,
                                         const G4String& directory,
                                         G4bool          ascii = false)
		{ BLAssert(false); return false; }
	virtual G4bool RetrievePhysicsTable( const G4ParticleDefinition* ,
                                             const G4String& directory,
                                             G4bool          ascii = false)
		{ BLAssert(false); return false; }
	virtual void SetProcessManager(const G4ProcessManager *p) { pRegProcess->SetProcessManager(p); }
	virtual  const G4ProcessManager* GetProcessManager() { return pRegProcess->GetProcessManager(); }
	virtual void      ResetNumberOfInteractionLengthLeft() { pRegProcess->ResetNumberOfInteractionLengthLeft(); }


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



BLCMDreweightprocess::BLCMDreweightprocess() : BLCommand(), BLCallback()
{
	registerCommand(BLCMDTYPE_PHYSICS);
	setSynopsis("modify the cross-section of a physics process.");
	setDescription("This command will modify the cross-section of a "
		"physics process, modifying the track weights so that "
		"weighted histograms give the same statistical result as if "
		"the command were not used. Used properly, this can greatly "
		"reduce the variance of the result.\n\n"
		"Care should be taken to ensure that all regions that should "
		"be sampled actually are sampled. For instance, if ratio>1 "
		"the interaction length will be reduced, and deep inside an "
		"absorber there may be no sampling because no simulated tracks "
		"ever get there, even though real tracks will. If ratio<1 then "
		"the upstream regions of absorbers will be under sampled; this "
		"is usually OK, as the desired result is to increase the "
		"sampling deep inside the absorber.\n\n"
		"For ratio>>1 this command can be used to examine rare "
		"processes. This can induce multiple rare interactions in an "
		"event when normally none would be expected; the weights will "
		"still correctly correspond to the real interaction, even "
		"though the event topologies don't.\n\n"
		"This command cannot reweight any continuous process (e.g. "
		"multiple scattering, ionization energy loss, etc.) -- it "
		"will issue a fatal exception if applied to such a process.\n\n"
		"This command should not be applied to other processes that "
		"re-weight tracks (e.g. the neutrino command). "
		"Indeed it probably won't give the correct weights if any "
		"such process applies to the particle (except for itself, it "
		"cannot determine the unmodified interaction length of such "
		"procsses, which is needed to compute the weight).\n\n"
		"The re-weighting applies to both PostStep and AtRest "
		"processes, but some AtRest processes do not work with this "
		"re-weighting; for instance, Decay works properly in PostStep "
		"for a moving particle, but not for a stopped one AtRest. "
		"This is related to the Geant4 limitation that exactly one "
		"process be active AtRest, and exactly one step be taken.\n\n"
		"Be sure to test your use of this commnd for a simple physical "
		"situation before believing its results."
	);

	particle = "";
	process = "";
	ratio = 1.0;
}

int BLCMDreweightprocess::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	BLCMDreweightprocess *p = 
			new BLCMDreweightprocess(defaultBLCMDreweightprocess);

	int retval = p->handleNamedArgs(namedArgs);

	// too early to modify processes here
	BLManager::getObject()->registerCallback(p,0);

	p->print("");

	return retval;
}

void BLCMDreweightprocess::defineNamedArgs()
{
	argString(particle,"particle","Comma-separated list of particle "
						"patterns ('' => all).");
	argString(process,"process","Comma-separated list of process patterns.");
	argDouble(ratio,"ratio","Ratio of artificial to real cross-section.");
}

void BLCMDreweightprocess::callback(int type)
{
	BLManager *blmanager = BLManager::getObject();
	if(type == 0) {
		// wrap matching processes for matching particles
		G4ParticleTable::G4PTblDicIterator *myParticleIterator =
			G4ParticleTable::GetParticleTable()->GetIterator();
		myParticleIterator->reset();
		while((*myParticleIterator)()) {
			G4ParticleDefinition *pd = myParticleIterator->value();
			if(!isMatchingParticle(pd)) continue;
			G4ProcessManager *pmgr = pd->GetProcessManager();
			BLAssert(pmgr != 0);
			int n=pmgr->GetProcessListLength();
			G4ProcessVector *pv=pmgr->GetProcessList();
			for(int i=0; i<n; ++i) {
			    G4VProcess *p=(*pv)[i];
			    if(!isMatchingProcess(p)) continue;
			    char tmp[256];
			    sprintf(tmp,"process=%s particle=%s",
				    		p->GetProcessName().c_str(),
						pd->GetParticleName().c_str());
			    if(wrapProcess(pmgr,p))  {
				bool flag=blmanager->showAllExceptions(true);
				G4Exception("reweightprocess",
				        	"Process re-weighted",
						JustWarning, tmp);
				blmanager->showAllExceptions(flag);
			    } else {
				G4Exception("reweightprocess",
				        	"Process cannot be re-weighted",
						FatalException, tmp);
			    }
			}
		}
	}
}


bool BLCMDreweightprocess::wrapProcess(G4ProcessManager *pmgr,G4VProcess *proc)
{
	// cannot re-weight any continuous processes (don't know how to do it)
	if(pmgr->GetAlongStepIndex(proc) >= 0) return false;

	ReweightProcess *rp = new ReweightProcess(pmgr,proc,ratio);

	// fix the positions of all processes (prevent them from moving)
	for(int i=idxAtRest; i<=idxPostStep; ++i) {
		G4ProcessVector *pv = pmgr->GetProcessVector(
					(G4ProcessVectorDoItIndex)i,typeDoIt);
		for(int k=1; k<pv->size(); ++k) { // 0 is Transportation
			pmgr->SetProcessOrdering((*pv)[k],
					(G4ProcessVectorDoItIndex)i,k);
		}
	}

	int orderAtRest = pmgr->GetProcessOrdering(proc,idxAtRest);
	int orderAlongStep = ordLast;
	int orderPostStep = pmgr->GetProcessOrdering(proc,idxPostStep);

	pmgr->RemoveProcess(proc);

	pmgr->AddProcess(rp,orderAtRest,orderAlongStep,orderPostStep);

	return true;
}

ReweightProcess::ReweightProcess(G4ProcessManager *_pmgr,G4VProcess *proc,
				double _ratio) : G4WrapperProcess(PREFIX)
{
	ratio = _ratio;
	pmgr = _pmgr;
	RegisterProcess(proc);
}

G4double ReweightProcess::AlongStepGetPhysicalInteractionLength(
			const G4Track& track, G4double  previousStepSize,
			G4double  currentMinimumStep, G4double& proposedSafety,
			G4GPILSelection* selection)
{
	// original process has no AlongStepGetPhysicalInteractionLength()
	// our AlongStepDoIt() always gets called.
	return DBL_MAX;
}

G4VParticleChange* ReweightProcess::AlongStepDoIt(const G4Track& track,
						const G4Step& stepData)
{
	G4VParticleChange *pc = &aParticleChange;
	pc->Initialize(track);

	// Re-weight the track because it survived the artificial interaction
	// length. Do it here, because all GPIL routines have been called,
	// but no PostStepDoIt() or AtRestDoIt() have been called. So all
	// processes have computed their currentInteractionLength, and this
	// re-weighting will apply to all discrete processes.
	// Note that only the first ReweightProcess in the list will do the
	// re-weight; others return with no change.
	G4double effectiveIL=0.0, realIL=0.0;
	int n=pmgr->GetProcessListLength();
	G4ProcessVector *pv=pmgr->GetProcessList();
	int nReWeight = 0;
	for(int i=0; i<n; ++i) {
		G4VProcess *p=(*pv)[i];
		G4double currentIL=p->GetCurrentInteractionLength();
		if(currentIL <= 0.0) continue;
		effectiveIL += 1.0/currentIL;
		ReweightProcess *rw = dynamic_cast<ReweightProcess*>(p);
		if(rw != 0) {
		    if(nReWeight++ != 0 && rw == this) return pc;
		    currentIL = rw->pRegProcess->GetCurrentInteractionLength();
		    BLAssert(currentIL > DBL_MIN);
		}
		realIL += 1.0/currentIL;
	}
	effectiveIL = 1.0/effectiveIL;
	realIL = 1.0/realIL;
	G4double stepLength = stepData.GetStepLength();
	G4double factor = exp(stepLength/effectiveIL - stepLength/realIL);
	pc->ProposeParentWeight(track.GetWeight()*factor);

	return pc;
}

G4double ReweightProcess::PostStepGetPhysicalInteractionLength(const G4Track& track,
		G4double previousStepSize, G4ForceCondition* condition)
{
	if(previousStepSize > 0.0) {
		SubtractNumberOfInteractionLengthLeft(previousStepSize);
		if(theNumberOfInteractionLengthLeft < 0.0) {
			theNumberOfInteractionLengthLeft=perMillion;
		}
	}

	pRegProcess->PostStepGetPhysicalInteractionLength(track,
						previousStepSize, condition);
	currentInteractionLength = pRegProcess->GetCurrentInteractionLength();
	if(currentInteractionLength < 0.0)
		G4Exception("reweightprocess","Process cannot be re-weighted",
			FatalException, pRegProcess->GetProcessName().c_str());
	currentInteractionLength /= ratio;

	G4double value = theNumberOfInteractionLengthLeft *
						currentInteractionLength;
	return value;
}

G4VParticleChange* ReweightProcess::PostStepDoIt(const G4Track& track,
						const G4Step&  stepData)
{
	G4VParticleChange *pc = pRegProcess->PostStepDoIt(track,stepData);

	// re-weight track
	pc->ProposeParentWeight(track.GetWeight()/ratio);
	// (Geant4 bug requires the following update to the weight, sometimes.)
	stepData.GetPostStepPoint()->SetWeight(track.GetWeight()/ratio);

	// re-weight secondaries
	int n=pc->GetNumberOfSecondaries();
	for(int i=0; i<n; ++i) {
		G4Track *t=pc->GetSecondary(i);
		t->SetWeight(t->GetWeight()/ratio);
	}

	ClearNumberOfInteractionLengthLeft();
	return pc;
}

G4double ReweightProcess::AtRestGetPhysicalInteractionLength(
		const G4Track& track, G4ForceCondition* condition)
{
	G4double value = pRegProcess->
			AtRestGetPhysicalInteractionLength(track,condition);

	value /= ratio;

	return value;
}

G4VParticleChange* ReweightProcess::AtRestDoIt(const G4Track& track,
                             			const G4Step& stepData)
{
	G4VParticleChange *pc = pRegProcess->AtRestDoIt(track,stepData);

	// re-weight track 
	pc->ProposeParentWeight(track.GetWeight()/ratio);
	// (Geant4 bug requires the following update to the weight, sometimes.)
	stepData.GetPostStepPoint()->SetWeight(track.GetWeight()/ratio);

	// re-weight secondaries
	int n=pc->GetNumberOfSecondaries();
	for(int i=0; i<n; ++i) {
		G4Track *t=pc->GetSecondary(i);
		t->SetWeight(t->GetWeight()/ratio);
	}

	ClearNumberOfInteractionLengthLeft();
	return pc;
}

