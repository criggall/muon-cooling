//	BLCMDphysics.cc
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

#include <vector>
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "G4MuIonisation.hh"
#include "G4ProcessTable.hh"
#include "G4ParticleTable.hh"
#include "G4StepLimiter.hh"
#include "G4PhysListFactory.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayTable.hh"
#include "G4VDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"
#include "G4SynchrotronRadiation.hh"
#include "G4GammaConversionToMuons.hh"
#include "G4Transportation.hh"
#include "G4PionDecayMakeSpin.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleChange.hh"
// for muon spin decay from examples/extended/field/field05
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"
#include "G4DecayWithSpin.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "BLAssert.hh"
#include "BLPhysics.hh"
#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLGlobalField.hh"

const G4String Indent1(IndentDesc);
const G4String Indent2(IndentDesc+"     ");
const G4String Indent3(IndentDesc+"         ");

/**	class BLCMDphysics implements a command to select from the geant4
 *	physics lists of physics processes.
 **/
class BLCMDphysics : public BLPhysics, public BLCommand, public BLCallback,
			public BLManager::TrackingAction {
	G4VUserPhysicsList *physicsList;
	G4int doStochastics;
	G4int fluct;
	G4String disable;
	// G4double minRangeCut; is protected in BLPhysics
	G4int list;
	G4int synchrotronRadiation;
	G4int synchrotronRadiationMuon;
	G4int gammaToMuPair;
	G4int spinTracking;
	std::vector<G4String> disableList;
	G4PhysListFactory *physListFactory;
	BLManager *manager;
public:
	/// Default Constructor.
	BLCMDphysics();

	/// Destructor.
	~BLCMDphysics() { }
	
	/// setDoStochastics() sets whether or not stochastic processes
	/// are to be enabled.
	void setDoStochastics(BLSetValue value, G4int warn=1);

	/// getPhysicsList() returns the G4PhysicsList.
	G4VUserPhysicsList *getPhysicsList() { return physicsList; }

	/// commandName() returns "physics".
	G4String commandName() { return "physics"; }
	
	/// command() implements the command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments.
	void defineNamedArgs();

	/// help() provides special help text for the physics command.
	void help(bool detailed);

	virtual bool isSpinTrackingEnabled() { return spinTracking != 0; }

	/// callback() is from BLCallback; it performs modifictions to
	/// particle physics lists.
	void callback(int type);

	/// PreUserTrackingAction() and PostUserTrackingAction() from
	/// BLManager::TrackingAction.
	void PreUserTrackingAction(const G4Track *track);
	void PostUserTrackingAction(const G4Track *track) { }

	/// isLocalPhysicsList() returns true if thsi is a local physics list.
	bool isLocalPhysicsList(G4String name);

	/// getLocalPhysicsList() returns a local physics list.
	G4VUserPhysicsList *getLocalPhysicsList(G4String name);

	/// listPhysicsLists() will list all known physics lists, with synopsis.
	void listPhysicsLists();

	/// printCase() prints the description of a physics list.
	void printCase(G4String name, G4String description);
	
	/// getSynopsis() returns the synopsis for a physics list, by name.
	G4String getSynopsis(G4String listName);
};

BLCMDphysics defaultPhysicsCmd;

/**	class TrackTimeLimit imposes maxTime as a limit on the global time
 *	of a track. Applies to both PostStep and AtRest.
 *
 *	This process is added to every particle; because it is not a stochastic
 *	process, this ensures there is always a valid process for At Rest.
 **/
class TrackTimeLimiter : public G4VProcess {
	G4double maxTime;
	G4double limit;
	class ParticleChange : public G4ParticleChange {
	public:
		ParticleChange() : G4ParticleChange() { }
	};
	ParticleChange change;
public:
	TrackTimeLimiter(double _maxTime) : 
			G4VProcess("TrackTimeLimiter",fUserDefined), change() {
		maxTime = _maxTime;
		limit = 0.0;
	}

	G4VProcess *clone() { return new TrackTimeLimiter(*this); }

	virtual G4double PostStepGetPhysicalInteractionLength(
			const G4Track& track, G4double   previousStepSize, 
			G4ForceCondition* condition)
		{ *condition = NotForced;
		  limit = maxTime - track.GetGlobalTime();
		  limit *= track.GetVelocity();
		  if(limit < 0.0) limit = DBL_MIN;
		  return limit;
		}
	G4VParticleChange* PostStepDoIt( const G4Track& track, 
      			const G4Step&  stepData)
		{ change.Initialize(track);
		  change.ProposeTrackStatus(fStopAndKill);
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
		{ *condition = NotForced;
		  limit = maxTime - track.GetGlobalTime();
		  if(limit < 0.0) limit = DBL_MIN;
		  return limit;
		}
	virtual G4VParticleChange* AtRestDoIt( const G4Track& track,
			const G4Step& stepData)
		{ change.Initialize(track);
		  change.ProposeGlobalTime(track.GetGlobalTime()+limit);
		  change.ProposeTrackStatus(fStopAndKill);
		  return &change;
		}
};

BLCMDphysics::BLCMDphysics() : BLPhysics(), BLCommand(), disableList()
{
	registerCommand(BLCMDTYPE_PHYSICS);
	setSynopsis("Defines the physics processes and controls them.");
	setDescription("Exactly one physics command must be present.\n"
		"This command implements the geant4 physics lists\n"
		"of physics processes.\n"
		"The command is 'physics QGSP' for the QGSP list, and\n"
		"similarly for the other lists. With no argument it prints\n"
		"the available physics lists. 'default' will select a standard "
		"list, currently FTFP_BERT.\n\n"
		"Note that stochastic processes are always disabled while\n"
		"tracking the tune and reference particles.\n"
		"The only non-stochastic processes are Transportation and\n"
		"ionization energy loss (with fluctuations disabled).\n\n"
		"For muon beam studies you may want to disable the Decay "
		"process.\n\n"
		"spinTracking=1 will enable tracking the spins of e+, e-, mu+, "
		"and mu- (only), including pion decays to polarized muons, "
		"polarized muon decays, and the force from the particle's "
		"magnetic moment. The rare pion decay to e nu is included, "
		"but the e are unpolarized. The muon decays give the correct "
		"distribution for the electron, but only approximate "
		"distributions for the neutrinos, and all are unpolarized.\n\n"
		"NOTE: this command defines the particles used throughout "
		"the simulation, so this command must come before others "
		"that use particle names.\n\n"
		"NOTE: the rare decay mode for pi+/pi- to e nu is always "
		"added, with branching ratio 1.230E-4.\n\n"
		"The default all-around physics list for HEP is called "
		"'FTFP_BERT'.\n\n"
		"NOTE: At present there is a bug in Geant4 tracking and the "
		"magnetic moment is not used. For >~ 0.1 GeV muons and "
		"practical fields this is a very small error.\n\n"
		"NOTE: the RadioactiveDecay process is added to all ions.");

	physicsList = 0;
	doStochastics = 1;
	fluct = 1;
	list = 0;
	synchrotronRadiation = 0;
	synchrotronRadiationMuon = 0;
	gammaToMuPair = 0;
	spinTracking = 0;
	disable = "";
	minRangeCut = 1.0 * mm;
	physListFactory = 0;
	manager = BLManager::getObject();
}

int BLCMDphysics::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(physicsList != 0 && argv.size() > 0) {
		printError("physics: Multiple physics commands not allowed!");
		return -1;
	}

	if(argv.size() == 0)
		return -1;

	// handle default and obsolete physics lists
	if(argv[0] == "default") {
		argv[0] = "FTFP_BERT";
		printf("default physics list is FTFP_BERT\n");
	} else if(argv[0] == "QGSC") {
		argv[0] = "QGSC_BERT";
		printf("Physics list QGSC is obsolete, replaced by QGSC_BERT\n");
	} else if(argv[0] == "QGSP") {
		argv[0] = "FTFP_BERT";
		printf("Physics list QGSP is obsolete, replaced by FTFP_BERT\n");
	}

	if(!physListFactory) physListFactory = new G4PhysListFactory();
	if(isLocalPhysicsList(argv[0]))
		physicsList = getLocalPhysicsList(argv[0]);
	if(!physicsList) {
		if(!physListFactory->IsReferencePhysList(argv[0])) {
			printError("physics: Unknown physics list '%s'",
							argv[0].c_str());
			printError("The choices are:");
			listPhysicsLists();
			return -1;
		}
		physicsList = physListFactory->GetReferencePhysList(argv[0]);
	}

	// Always have ions decay radioactively
	G4VModularPhysicsList *pl = 
			dynamic_cast<G4VModularPhysicsList*>(physicsList);
	BLAssert(pl != 0);
	pl->RegisterPhysics(new G4RadioactiveDecayPhysics());

	handleNamedArgs(namedArgs);

	if(spinTracking != 0) {
		BLGlobalField::setSpinTracking(true);
		manager->registerTrackingAction(this);
	}

	disableList = splitString(disable,',',true);

	manager->registerPhysics(this);

	// now is too early to manage the physics processes; do it in a callback
	manager->registerCallback(this,0);

	print(argv[0]);

	return 0;
}

void BLCMDphysics::defineNamedArgs()
{
	argString(disable,"disable","Comma-separated list of processes to disable (e.g. 'Decay,msc').");
	argString(disable,"inactivate","Synonym of disable.");
	argString(disable,"deactivate","Synonym of disable.");
	argInt(doStochastics,"doStochastics","Set to zero to disable all stochastic processes.");
	argInt(fluct,"fluct","Set to zero to disable ionization energy loss fluctuations.");
	argDouble(minRangeCut,"minRangeCut","Minimum range cut for particle production (1 mm)");
	argInt(list,"list","Nonzero to list the processes (later in output).");
	argInt(gammaToMuPair,"gammaToMuPair","Nonzero to add gamma->mu+mu- (0).");
	argInt(spinTracking,"spinTracking","Nonzero to track particle spins (0).");
	argInt(synchrotronRadiation,"synchrotronRadiation","Nonzero to add synchrotron radiation to the physics list for e- and e+.");
	argInt(synchrotronRadiationMuon,"synchrotronRadiationMuon","Nonzero to add synchrotron radiation to the physics list for mu- and mu+ NOTE: This is experimental, and may not work correctly.");
	argDouble(maxTime,"maxTime","Maximum time for tracking (1000000ns).",1.0*ns);
}

void BLCMDphysics::help(bool detailed)
{
	BLCommand::help(detailed);
	if(detailed) {
		printf("\n             ----- PHYSICS LISTS -----\n\n");
		G4String s;
s="Further guidance in selecting physics lists is available at:\n"
"http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsListGuide/html/index.html\n\n";
		printf(wrapWords(s,Indent1,Indent1).c_str());
s="The default all-around physics list for HEP is 'FTFP_BERT'.\n\n";
		printf(wrapWords(s,Indent1,Indent1).c_str());
s="Each basic physics list uses the standard Geant4 electromagnetic physics Opt0. "
"This can be changed by adding one of these suffixes (the best but slowest is _EMZ):\n";
		printf(wrapWords(s,Indent1,Indent1).c_str());
s="_EMV EM Opt1 less precise but faster.\n";
		printf(wrapWords(s,Indent1,Indent2).c_str());
s="_EMX EM Opt2 less precise but faster.\n";
		printf(wrapWords(s,Indent1,Indent2).c_str());
s="_EMY EM Opt3 it uses a set of EM processes with accurate simulation of gamma and charged particle transport. Only the Urban multiple scattering model is used for all charged particles and all energies. Slower.\n";
		printf(wrapWords(s,Indent1,Indent2).c_str());
s="_EMZ EM Opt4 the best set of electromagnetic physics models selected from the low energy and standard packages. With its concentration on the best possible physics, electromagnetic option 4 is slower than the standard EM package.\n";
		printf(wrapWords(s,Indent1,Indent2).c_str());
s="_LIV is made on top of electromagnetic option 3 by substitution of standard models for gamma and electrons from Livermore set of models.\n";
		printf(wrapWords(s,Indent1,Indent2).c_str());
s="_PEN is made on top of electromagnetic option 3 by substitution of standard models for gamma, electrons and positrons from Penelope-2008 set of models.\n";
		printf(wrapWords(s,Indent1,Indent2).c_str());
s="_GS is made on top of the default electromagnetic configurations by substitution of the Urban multiple scattering model for electrons and positrons by the Goudsmit-Saunderson model.\n";
		printf(wrapWords(s,Indent1,Indent2).c_str());
s="_LE is made on top of the default electromagnetic configurations by substitution of the Urban multiple scattering model for electrons and positrons by the LowEWentzelVI model. Also, using 5D gamma conversion model and Lindhard-Sorensent model for ion ionisation.\n";
		printf(wrapWords(s,Indent1,Indent2).c_str());
s="_WVI is made on top of the default electromagnetic configurations by substitution of the Urban multiple scattering model for electrons and positrons by the WentzelVI model and ATIMA ion ionisation model.\n";
		printf(wrapWords(s,Indent1,Indent2).c_str());
s="_SS is made on top default electromagnetic configurations by substitution of all multiple scattering models by single scattering models.\n\n";
		printf(wrapWords(s,Indent1,Indent2).c_str());
s="The basic physics lists specify the hadronic modeling (can add EM suffix):\n";
		printf(wrapWords(s,Indent1,Indent1).c_str());
s="             FTFP_BERT The current Geant4 default, good for most HEP users.\n" ;
		printf(s);
s="             QBBC      Recommended for medical and space physics simulations.\n";
		printf(s);
s="             QGSP_BERT The former Geant4 default, good for most HEP users.\n";
		printf(s);
s="             QGSP_BIC  Another list for HEP.\n";
		printf(s);
s="Shielding Recommended for simulation of deep shielding. Neutrons of 20 MeV and lower use the High Precision neutron models and cross sections to describe elastic and inelastic scattering, capture and fission. The G4NDL database is required for this physics list.\n";
		printf(wrapWords(s,Indent1,Indent3).c_str());
		printf("\n             These physics lists are available:\n");
		listPhysicsLists();
		printf("\n             NOTE: all _HP lists require the G4NDL dataset.\n");
	}
}

void BLCMDphysics::setDoStochastics(BLSetValue value, G4int warn)
{
	if(value == NORMAL) {
		if(doStochastics)
			value = FORCE_ON;
		else
			value = FORCE_OFF;
	}

	G4UImanager* UI = G4UImanager::GetUIpointer();
	G4String cmd;
	char rangeCmd[64];
	switch(value) {
	case FORCE_OFF:
		// turn off delta rays (and other particle production) by
		// setting production cut very high in energy by setting
		// the required range large.
		UI->ApplyCommand("/run/setCut 1 km");
		// Turn off enegry-loss straggling
		UI->ApplyCommand("/process/eLoss/fluct false");
		cmd = "/process/inactivate ";
		switch(warn) {
		case 0:	break;
		case 1:
		case 2:
			G4Exception("physics",
				"All stochastic processes disabled",
				JustWarning,"");
			break;
		}
		stochasticsEnabled = false;
		break;
	case FORCE_ON:
		// turn on delta rays (and other particle production)
		sprintf(rangeCmd,"/run/setCut %.6g mm",minRangeCut/mm);
		UI->ApplyCommand(rangeCmd);
		// Turn on enegry-loss straggling
		UI->ApplyCommand("/process/eLoss/fluct true");
		cmd = "/process/activate ";
		if(warn != 0) printf("Stochastic processes are enabled.\n");
		stochasticsEnabled = true;
		break;
	case NORMAL:
		return;
	}

	if(fluct == 0) UI->ApplyCommand("/process/eLoss/fluct false");

	// run through list of processes, applying cmd to stochastic ones
	G4ProcessTable *pt = G4ProcessTable::GetProcessTable();
	G4ProcessTable::G4ProcNameVector *nv = pt->GetNameList();
	if(list) printf("List of physics processes:\n");
	for(unsigned int i=0; i<nv->size(); ++i) {
		G4String name = (*nv)[i];
		if(list) printf("\t%s\n",name.c_str());
		if(isStochasticProcess(name))
			UI->ApplyCommand(cmd+name);
	}
	if(list) printf("\n");
	list = 0;	// print list only once

	// now disable all proceses in the disableList
	cmd = "/process/inactivate ";
	for(unsigned i=0; i<disableList.size(); ++i) {
		switch(warn) {
		case 0:	break;
		case 1:
		case 2:	G4Exception("physics","Physics process disabled",
					JustWarning,disableList[i].c_str());
			break;
		}
		UI->ApplyCommand(cmd+disableList[i]);
	}
}

void BLCMDphysics::callback(int type)
{
	G4ProcessTable* processTbl = G4ProcessTable::GetProcessTable();
	G4ParticleTable *particleTbl = G4ParticleTable::GetParticleTable();
	G4ParticleTable::G4PTblDicIterator *iter = particleTbl->GetIterator();

	// add G4StepLimiter and TrackTimeLImiter to all particles
	iter->reset();
	while( (*iter)() ) {
	    G4ParticleDefinition* particle = iter->value();
	    G4ProcessManager* pmgr = particle->GetProcessManager();
	    if(!pmgr) {
		printf("BLCMDphysics: particle '%s' has no ProcessManager!\n",
			particle->GetParticleName().c_str());
		continue;
	    }
	    pmgr->AddProcess(new G4StepLimiter(),-1,-1,ordDefault);
	    pmgr->AddProcess(new TrackTimeLimiter(maxTime),ordDefault,-1,
	    							ordDefault);
	}

	// add pi+/pi- decay to e nu
	double br=1.230e-4;
	{ G4ParticleDefinition *pi=particleTbl->FindParticle("pi+");
	  G4DecayTable* dt = pi->GetDecayTable();
	  BLAssert(dt->entries() == 1); /// ensure only mu-nu is present
	  G4VDecayChannel *munu = dt->GetDecayChannel(0);
	  munu->SetBR(1.0-br);  
	  dt->Insert(new G4PhaseSpaceDecayChannel("pi+",br,2,"e+","nu_e"));
	}
	{ G4ParticleDefinition *pi=particleTbl->FindParticle("pi-");
	  G4DecayTable* dt = pi->GetDecayTable();
	  BLAssert(dt->entries() == 1); /// ensure only mu-nu is present
	  G4VDecayChannel *munu = dt->GetDecayChannel(0);
	  munu->SetBR(1.0-br);  
	  dt->Insert(new G4PhaseSpaceDecayChannel("pi-",br,2,"e-","anti_nu_e"));
	}

	if(synchrotronRadiation) {
		particleTbl->FindParticle("e-")->GetProcessManager()->
			AddDiscreteProcess(new G4SynchrotronRadiation);
		particleTbl->FindParticle("e+")->GetProcessManager()->
			AddDiscreteProcess(new G4SynchrotronRadiation);
	}

	if(synchrotronRadiationMuon) {
		particleTbl->FindParticle("mu-")->GetProcessManager()->
			AddDiscreteProcess(new G4SynchrotronRadiation);
		particleTbl->FindParticle("mu+")->GetProcessManager()->
			AddDiscreteProcess(new G4SynchrotronRadiation);
	}

	if(gammaToMuPair) {
		particleTbl->FindParticle("gamma")->GetProcessManager()->
			AddDiscreteProcess(new G4GammaConversionToMuons);
	}

	if(spinTracking) {
		// enable tracking with magnetic moment
		iter->reset();
		while( (*iter)() ) {
			G4ParticleDefinition* particle = iter->value();
			assert(particle != 0);
			if(particle->GetPDGMagneticMoment() == 0.0) continue;
			G4ProcessManager* pmgr=particle->GetProcessManager();
			assert(pmgr != 0);
			G4ProcessVector *pv = pmgr->GetProcessList();
			assert(pv != 0);
			for(int i=0; i<pv->size(); ++i) {
				G4VProcess *p = (*pv)[i];
				assert(p != 0);
				if(p->GetProcessType() != fTransportation) 
					continue;
				G4Transportation *t = 
					dynamic_cast<G4Transportation*>(p);
				if(!t) continue;
				t->EnableUseMagneticMoment(true);
			}
		}

		// handling pion and muon decays with spin is modeled after
		// examples/extended/field/field05
		G4ProcessManager* pmgr;

		// pi+ decay with spin
		G4VProcess *newDecay =
				 new G4PionDecayMakeSpin("PionDecayMakeSpin");
		G4VProcess *decay =
		     processTbl->FindProcess("Decay",G4PionPlus::PionPlus());
		pmgr = G4PionPlus::PionPlus()->GetProcessManager();
		assert(pmgr != 0);
		if(decay) pmgr->RemoveProcess(decay);
		pmgr->AddProcess(newDecay);
		pmgr ->SetProcessOrdering(newDecay, idxPostStep);
		pmgr ->SetProcessOrdering(newDecay, idxAtRest);

		// pi- decay with spin
		decay =
		    processTbl->FindProcess("Decay",G4PionMinus::PionMinus());
		pmgr = G4PionMinus::PionMinus()->GetProcessManager();
		assert(pmgr != 0);
		if(decay) pmgr->RemoveProcess(decay);
		pmgr->AddProcess(newDecay);
		pmgr ->SetProcessOrdering(newDecay, idxPostStep);
		pmgr ->SetProcessOrdering(newDecay, idxAtRest);

		// muon decay modes with spin
		G4DecayTable* dt = new G4DecayTable();
		dt->Insert(new G4MuonDecayChannelWithSpin("mu+",0.986));
		dt->Insert(new G4MuonRadiativeDecayChannelWithSpin("mu+",.014));
		G4MuonPlus::MuonPlusDefinition()->SetDecayTable(dt);
		dt = new G4DecayTable();
		dt->Insert(new G4MuonDecayChannelWithSpin("mu-",0.986));
		dt->Insert(new G4MuonRadiativeDecayChannelWithSpin("mu-",.014));
		G4MuonMinus::MuonMinusDefinition()->SetDecayTable(dt);

		// replace the Decay process for mu+
		G4DecayWithSpin* decayWithSpin = new G4DecayWithSpin();
		decay =
		     processTbl->FindProcess("Decay",G4MuonPlus::MuonPlus());
		pmgr = G4MuonPlus::MuonPlus()->GetProcessManager();
        	assert(pmgr != 0);
		if(decay) pmgr->RemoveProcess(decay);
		pmgr->AddProcess(decayWithSpin);
		pmgr ->SetProcessOrdering(decayWithSpin, idxPostStep);
		pmgr ->SetProcessOrdering(decayWithSpin, idxAtRest);

		// replace the Decay process for mu-
		decay =
		     processTbl->FindProcess("Decay",G4MuonMinus::MuonMinus());
		pmgr = G4MuonMinus::MuonMinus()->GetProcessManager();
		assert(pmgr != 0);
		if(decay) pmgr->RemoveProcess(decay);
		pmgr->AddProcess(decayWithSpin);
		pmgr ->SetProcessOrdering(decayWithSpin, idxPostStep);
		pmgr ->SetProcessOrdering(decayWithSpin, idxAtRest);
	}
}

void BLCMDphysics::PreUserTrackingAction(const G4Track *track)
{
	if(spinTracking) {
		double anomaly=0.0;
		switch(track->GetParticleDefinition()->GetPDGEncoding()) {
		case -11: case 11:	// e+ e-
			anomaly = 0.0011596521811;
			break;
		case -13: case 13:	// mu+ mu-
			anomaly = 0.0011659208;
			break;
		}
		BLGlobalField::setSpinAnomaly(anomaly);
		if(manager->getSteppingVerbose() != 0 && anomaly != 0.0)
			printf("physics: setSpinAnomaly(%.14f)\n",anomaly);
	}
}

bool BLCMDphysics::isLocalPhysicsList(G4String name)
{
	return false;
}

G4VUserPhysicsList *BLCMDphysics::getLocalPhysicsList(G4String name)
{
	return 0;
}


void BLCMDphysics::listPhysicsLists()
{
	if(!physListFactory) physListFactory = new G4PhysListFactory();
	const std::vector<G4String> v = physListFactory->AvailablePhysLists();
	for(unsigned i=0; i<v.size(); ++i) {
		G4String n = v[i];
		G4String s = getSynopsis(n);
		printCase(n,s);
	}
}

void BLCMDphysics::printCase(G4String name, G4String description)
{
	G4String indent1("             "), indent("       ");
	indent += indent1;
	indent1 += name;
	do { indent1 += " "; } while(indent1.size() < indent.size());
	printf("%s",wrapWords(description,indent1,indent).c_str());
}

G4String BLCMDphysics::getSynopsis(G4String listName)
{
#define CASE(NAME,SYNOPSIS)				\
	if(listName == #NAME) return SYNOPSIS;

	CASE(CHIPS,"No synopsis available.");
	CASE(FTFP_BERT,"For calorimetry, colliders, and general HEP. Based on the Fritiof description of string excitation and fragmentation. Uses Geant4 Bertini cascade for primary protons, neutrons, pions and Kaons below ~10GeV.");
	CASE(FTFP_BERT_HP,"identical to FTFP_BERT except that neutrons of 20 MeV and lower use the High Precision neutron models and cross sections to describe elastic and inelastic scattering, capture and fission. The G4NDL database is required.");
	CASE(FTFP_BERT_ATL,"identical to FTFP_BERT except that transition between the FTF model and the Bertini cascade is changed for Atlas.");
	CASE(FTFP_BERT_TRV,"identical to FTFP_BERT except that several electromagnetic and hadronic models are substituted by their alternative.");
	CASE(FTFP_INCLXX,"identical to FTFP_BERT except that the Bertini cascade is substituted by the INCL++ cascade.");
	CASE(FTFP_INCLXX_HP,"identical to FTFP_INCLXX except that neutrons of 20 MeV and lower use the High Precision neutron models and cross sections to describe elastic and inelastic scattering, capture and fission. The G4NDL database is required.");
	CASE(FTFQGSP_BERT,"identical to FTFP_BERT except that QGS string model is used.");
	CASE(FTF_BIC,"identical to FTFP_BERT except that the Binary cascade is used as by the FTF model instead of internal FTF cascade code for rescattering of secondary particles within nucleus.");
	CASE(QBBC,"Similar to FTFP_BERT but optimized for medical and space physics. Recomended for applications where accurate simulation for low-energy transport of protons and neutrons is needed. It usually produces the best agreement in the energy range below 1 GeV for thin target experiments.");
	CASE(QGSP_BERT,"Uses the QGS quark-gluon model with the Bertini cascade model.");
	CASE(QGSP_BERT_HP,"identical to QGSP_BERT except that neutrons of 20 MeV and lower use the High Precision neutron models and cross sections to describe elastic and inelastic scattering, capture and fission. The G4NDL database is required.");
	CASE(QGSP_FTFP_BERT,"identical to QGSP_BERT except that neutrons of 20 MeV and lower use the High Precision neutron models and cross sections to describe elastic and inelastic scattering, capture and fission. The G4NDL database is required.");
	CASE(QGSP_INCLXX,"identical to QGSP_BERT except that the Bertini cascade is substituted by the INCL++ cascade.");
	CASE(QGSP_INCLXX_HP,"identical to QGSP_BERT_HP except that the Bertini cascade is substituted by the INCL++ cascade. The G4NDL database is required.");
	CASE(QGSP_BIC,"Similar to QGSP_BERT, except the BIC binary cascade model is used.");
	CASE(QGSP_BIC_HP,"identical to QGSP_BIC except that neutrons of 20 MeV and lower use the High Precision neutron models and cross sections to describe elastic and inelastic scattering, capture and fission. The G4NDL database is required.");
	CASE(QGSP_BIC_AllHP,"identical to QGSP_BIC except that neutrons, protons, and light ions of 20 MeV and lower use the High Precision models and cross sections to describe elastic and inelastic scattering, capture and fission. The G4NDL database is required for this physics list.");
	CASE(Shielding,"Recommended for simulation of deep shielding. Neutrons of 20 MeV and lower use the High Precision neutron models and cross sections to describe elastic and inelastic scattering, capture and fission. The G4NDL database is required for this physics list.");
	CASE(ShieldingM,"different transition from the Bertini model to the FTF model: from 9.5 to 9.9 GeV.");
	CASE(ShieldingLEND,"different configurations of low-energy neutron transport with LEND option.");

	// These are from previous version (earlier CASEs take precedence):
	CASE(LBE,"For low background experiments (e.g. underground)");
	CASE(LHEP,"For calorimetry -- is the fastest, when it comes to CPU. It uses the LEP and HEP parametrized models for inelastic scattering. The modeling parametrizes the final states individual inelastic reactions, so you will not see resonances, and the detailed secondary angular distributions for O(100MeV) reactions may not be described perfectly. The average quantities will be well described.");
	CASE(LHEP_BERT,"For 'typical' HEP collider detectors and low-energy penetration shielding (< 5 GeV) -- is as LHEP for energies above 3GeV. Below 3 GeV the Bertini cascade is used for nucleon and pion induced reactions. The price to pay is reduced CPU performance. It does not include gamma and electro nuclear physics.");
	CASE(QBBC,"No synopsis available.");
	CASE(QGSC,"For calorimetry and high energy physics trackers -- is as QGSP for the initial reaction, but uses chiral invariant phase-space decay (multi-quasmon fragmentation) to model the behavior of the system's fragmentation.");
	CASE(QGSC_BERT,"For calorimetry and high energy physics trackers -- is as QGSP for the initial reaction, but uses chiral invariant phase-space decay (multi-quasmon fragmentation) to model the behavior of the system's fragmentation. Uses Geant4 Bertini cascade for nucleon and pion induced reactions.");
	CASE(QGSP,"For calorimetry and high energy physics trackers and high-energy and medium-energy production targets -- uses theory driven modeling for the reactions of energetic pions, kaons, and nucleons. It employs quark gluon string model for the 'punch-through' interactions of the projectile with a nucleus, the string excitation cross-sections being calculated in quasi-eikonal approximation. A pre-equilibrium decay model with an extensive evaporation phase to model the behavior of the nucleus 'after the punch'. It uses current best pion cross-section.");
	CASE(QGSP_BERT,"Like QGSP, but using Geant4 Bertini cascade for primary protons, neutrons, pions and Kaons below ~10GeV. In comparison to experimental data we find improved agreement to data compared to QGSP which uses the low energy parameterised (LEP) model for all particles at these energies. The Bertini model produces more secondary neutrons and protons than the LEP model, yielding a better agreement to experimental data.");
	CASE(QGSP_BERT_HP,"Like QGSP_BERT but with _HP modeling for neutrons. The G4NDL database is required.");
	CASE(QGSP_BIC,"Like QGSP, but using Geant4 Binary cascade for primary protons and neutrons with energies below ~10GeV, thus replacing the use of the LEP model for protons and neutrons In comparison to teh LEP model, Binary cascade better describes production of secondary particles produced in interactions of protons and neutrons with nuclei.");
	CASE(QGSP_BIC_HP,"Like QGSP_BIC but with _HP modeling for neutrons. The G4NDL database is required.");

	return "No synopsis available.";
}
