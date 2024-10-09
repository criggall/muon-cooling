//	BLCMDmuonium.cc
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

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayTable.hh"
#include "G4Decay.hh"
#include "G4DecayProducts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4MuonDecayChannel.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4TransportationManager.hh"
#include "G4Navigator.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"

#include "BLAssert.hh"
#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLCallback.hh"
#include "BLParam.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

extern void g4bl_exit(int);

const double IonizationEnergy = 13.6 * eV;
const double StripMFP = 6.55E-7 * (gram/cm2); // mean free path for stripping
// [1] Prog. Theor. Exp. Phys. 2013, 103C01 (10 pages) DOI: 10.1093/ptep/ptt080
//   Bakule et al, Measurement of muonium emission from silica aerogel
// This assumes the MFP in g/cm^2 for silica aerogel applies to all materials.



/**	class Muonium - particle definition for muonium.
 **/
class Muonium : public G4ParticleDefinition {
	static Muonium *object;
	Muonium(const G4String&     aName,        G4double        mass,
		G4double            width,        G4double        charge,
		G4int               iSpin,        G4int           iParity,
		G4int               iConjugation, G4int           iIsospin,
		G4int               iIsospin3,    G4int           gParity,
		const G4String&     pType,        G4int           lepton,
		G4int               baryon,       G4int           encoding,
		G4bool              stable,       G4double        lifetime,
		G4DecayTable        *decaytable );
	static Muonium *create();
public:
	static Muonium *getObject() {
		if(!object) object = create();
		return object;
	}
};
Muonium *Muonium::object = 0;

/**	class DecayingMuPlus -- particle definition for a shortlived muon
 **/
class DecayingMuPlus : public G4ParticleDefinition
{
 private:
   static DecayingMuPlus* object;
   DecayingMuPlus(){}
   ~DecayingMuPlus(){}

 public:
   static DecayingMuPlus* Definition();
   static DecayingMuPlus* MuonPlusDefinition() { return Definition(); }
   static DecayingMuPlus* MuonPlus() { return Definition(); }
};
DecayingMuPlus *DecayingMuPlus::object = 0;


/** class MuoniumDecay - physics process to decay a Muonium atom.
 *
 **/
class MuoniumDecay : public G4VProcess {
	class ParticleChange : public G4ParticleChange {
	public:
		ParticleChange() : G4ParticleChange() { }
	};
	class BLCMDmuonium *cmd;
	ParticleChange change;
	double maxProperTime;
public:
	MuoniumDecay(class BLCMDmuonium *p);
	void StartTracking(G4Track *track);
	virtual G4double PostStepGetPhysicalInteractionLength(
			const G4Track& track, G4double   previousStepSize,
			G4ForceCondition* condition) ;
	virtual G4VParticleChange* PostStepDoIt( const G4Track& track,
			const G4Step&  stepData);
	virtual G4double AlongStepGetPhysicalInteractionLength(
			const G4Track& track, G4double  previousStepSize,
			G4double  currentMinimumStep, G4double& proposedSafety,
			G4GPILSelection* selection) {
		return -1.0; 
	}
	virtual G4VParticleChange* AlongStepDoIt( const G4Track& track,
			const G4Step& stepData) {
		change.Initialize(track);
		return &change;
	}
	virtual G4double AtRestGetPhysicalInteractionLength(
			const G4Track& track, G4ForceCondition* condition);
	virtual G4VParticleChange* AtRestDoIt( const G4Track& track,
			const G4Step& stepData);
};

/** class MuoniumSurface - physics process to strip, adsorb, or reflect a
 *  Muonium atom.
 **/
class MuoniumSurface : public G4VProcess {
	enum Mode { TRANSPORT, STRIP, ADSORB, REFLECT };
	class ParticleChange : public G4ParticleChange {
	public:
		ParticleChange() : G4ParticleChange() { }
	};
	class BLCMDmuonium *cmd;
	ParticleChange change;
	Mode mode;	// set in PostStepGPIL
	double maxPath;	// set in PostStepGPIL
public:
	MuoniumSurface(class BLCMDmuonium *p);
	void StartTracking(G4Track *track);
	virtual G4double PostStepGetPhysicalInteractionLength(
			const G4Track& track, G4double   previousStepSize,
			G4ForceCondition* condition) ;
	virtual G4VParticleChange* PostStepDoIt( const G4Track& track,
			const G4Step&  stepData);
	virtual G4double AlongStepGetPhysicalInteractionLength(
			const G4Track& track, G4double  previousStepSize,
			G4double  currentMinimumStep, G4double& proposedSafety,
			G4GPILSelection* selection) {
		return -1.0; 
	}
	virtual G4VParticleChange* AlongStepDoIt( const G4Track& track,
			const G4Step& stepData) {
		change.Initialize(track);
		return &change;
	}
	virtual G4double AtRestGetPhysicalInteractionLength(
			const G4Track& track, G4ForceCondition* condition) {
		return -1.0;
	}
	virtual G4VParticleChange* AtRestDoIt( const G4Track& track,
			const G4Step& stepData) {
		change.Initialize(track);
		return &change;
	}
};


/**	class BLCMDmuonium -- implement the muonium command
 **/
class BLCMDmuonium : public BLCommand , public BLCallback {
	G4double E1;
	G4double E2;
	G4double E3;
	G4double E4;
	G4double stripMFP;
	G4double minDensity;
	G4int verbose;
	friend class MuoniumDecay;
	friend class MuoniumSurface;
public:
	/// Constructor.
	BLCMDmuonium();

	/// Destructor.
	~BLCMDmuonium() { }

	/// commandName() returns "muonium"
	G4String commandName()  { return "muonium"; }

	/// command() executes the command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command. 
	void defineNamedArgs();

	/// callack() from BLCallback.
	void callback(int type);

};
BLCMDmuonium defaultMuonium;

BLCMDmuonium::BLCMDmuonium() : BLCommand(), BLCallback()
{
	registerCommand(BLCMDTYPE_PHYSICS);
	setSynopsis("Define the muonium (Mu) particle and processes.");
	setDescription("Muonium is an atom, and there are lots of detailed "
	"atomic physics processes that are NOT modeled here. This is good "
	"enough to model muonium transport in vacuum, and adsorbing, "
	"reflecting, or stripping at a surface; nothing else is modeled, and "
	"these are only approximated and parameterized. The ionization energy "
	"of Mu is known to be within 0.5% of Hydrogen (13.6 eV).\n\n"
	"MuoniumDecay: the Mu lifetime is the same as for mu+; it is split "
	"into an e- and a decaying_mu+. In the Mu rest frame they have equal "
	"and opposite momenta in a random direction; they share 13.6 eV of KE, "
	"with the ratio from a classical orbit. The decaying_mu+ is a mu+ "
	"which decays immediately.\n\n"
	"MuoniumSurface: handles stripping, reflecting, and adsorbing at a "
	"surface; parameters are stripMFP, minDensity, E1, E2, E3, and E4. "
	"The mean-free-path for stripping is ASSUMED to be independent of "
	"KE and material; the default value is from [1] (silica aerogel). "
	"If the material being entered has a density < minDensity, the Mu is "
	"transported WITHOUT ENERGY LOSS. Otherwise its KE determines which "
	"process applies: if KE<E1, it is adsorbed or reflected; if KE>E2 "
	"it is stripped; between E1 and E2 the probabilities are linear. "
	"If it is not stripped, then if KE<E3 it is adsorbed; if KE>E4 it is "
	"reflected, and between E3 and E4 the probabilities are linear.\n\n"
	"For stripping and decay the Mu polarization is transferred to the "
	"mu+.\n\n"
	"NOTE: This command must come before the 'physics' command.\n\n"
	"[1] Prog. Theor. Exp. Phys. 2013, 103C01 (10 pages) "
	"DOI:10.1093/ptep/ptt080 "
	" Bakule et al, Measurement of muonium emission from silica aerogel\n\n"
	);

	E1 = IonizationEnergy - 2.0*eV;
	E2 = IonizationEnergy * 2.0;
	E3 = IonizationEnergy / 10.0;
	E4 = E1;
	stripMFP = StripMFP;
	minDensity = 0.01 * (gram/cm3);
	verbose = 0;
}

int BLCMDmuonium::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	static bool first=true;
	if(!first) {
		printf("muonium: Command cannot be used more than once.\n");
		return -1;
	}
	first = false;

	// check that no physics command has been executed
	if(BLManager::getObject()->getPhysics() != 0)
		G4Exception("muonium","Used after 'physics' command",
			FatalException,"The 'muonium' command must come BEFORE "
					"the 'physics' command.");

	int retval = 0;

	retval = handleNamedArgs(namedArgs);

	print("");

	// too early to define particles and physics processes
	BLManager::getObject()->registerCallback(this,-1);
	BLManager::getObject()->registerCallback(this,0);

	return retval;
}

void BLCMDmuonium::defineNamedArgs()
{
	argDouble(E1,"E1","The KE below which all Mu adsorb or reflect (11.6 eV).",eV);
	argDouble(E2,"E2","The KE above which all Mu strip (27.2 eV).",eV);
	argDouble(E3,"E3","The KE below which all un-stripped Mu adsorb (1.36 eV).",eV);
	argDouble(E4,"E4","The KE above which all un-stripped Mu reflect (11.6 eV).",eV);
	argDouble(stripMFP,"stripMFP","The mean free path for stripping (gram/cm2).", gram/cm2);
	argDouble(minDensity,"minDensity","The minimum density below which Mu do not interact in any way (gram/cm3).",gram/cm3);
	argInt(verbose,"verbose","Set nonzero for verbose prints (0).");
}

void BLCMDmuonium::callback(int type)
{
fprintf(stderr,"BLCMDmuonium::callback(%d)\n",type);

	if(type == -1) {
		G4ParticleTable *table = G4ParticleTable::GetParticleTable();
		// define particle: decaying_mu+
		G4ParticleDefinition *pd1 = DecayingMuPlus::Definition();
		G4ParticleDefinition *pd2 = table->FindParticle("decaying_mu+");
		G4ParticleDefinition *pd3 = table->FindParticle(-23456);
		BLAssert(pd1 != 0 && pd1 == pd2 && pd1 == pd3);
		// define particle: Mu
		pd1 = Muonium::getObject();
		pd2 = table->FindParticle("Mu");
		pd3 = table->FindParticle(12345);
		BLAssert(pd1 != 0 && pd1 == pd2 && pd1 == pd3);
		pd1->SetProcessManager(new G4ProcessManager(pd1));
		if(verbose >= 1) printf("muonium: particles decaying_mu+ and "
								"Mu created\n");
	} else if(type == 0) {
		// in Muonium,  replace Decay process with MuoniumDecay
		G4ParticleTable *table = G4ParticleTable::GetParticleTable();
		G4ParticleDefinition *pd = table->FindParticle("Mu");
		G4ProcessManager *pmgr = pd->GetProcessManager();
		G4ProcessVector *pv = pmgr->GetProcessList();
		int n = pv->size();
		for(int i=0; i<n; ++i) {
			G4VProcess *process = (*pv)[i];
			if(process->GetProcessName().find("Decay") != 
							G4String::npos) {
				pmgr->RemoveProcess(process);
				break;
			}
		}
		MuoniumDecay *md = new MuoniumDecay(this);
		pmgr->AddProcess(md);
		pmgr->SetProcessOrdering(md,idxPostStep);
		pmgr->SetProcessOrdering(md,idxAtRest);
		// add MuoniumSurface
		MuoniumSurface *ms = new MuoniumSurface(this);
		pmgr->AddProcess(ms);
		pmgr->SetProcessOrdering(ms,idxPostStep);
		// decaying_mu+ has correct processes already
		if(verbose >= 1) printf("muonium: MuoniumDecay and "
					"MuoniumSurface added to Mu.\n");
	}
}

Muonium::Muonium(
	    const G4String&     aName,        G4double            mass,
	    G4double            width,        G4double            charge,
	    G4int               iSpin,        G4int               iParity,
	    G4int               iConjugation, G4int               iIsospin,
	    G4int               iIsospin3,    G4int               gParity,
	    const G4String&     pType,        G4int               lepton,
	    G4int               baryon,       G4int               encoding,
	    G4bool              stable,       G4double            lifetime,
	    G4DecayTable        *decaytable )
	    : G4ParticleDefinition( aName, mass, width, charge, iSpin, iParity,
		    iConjugation, iIsospin, iIsospin3, gParity, pType,
		    lepton, baryon, encoding, stable, lifetime, decaytable )
{
}

Muonium *Muonium::create()
{
	G4ParticleTable *table = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition *muon = table->FindParticle("mu+");
	G4ParticleDefinition *electron = table->FindParticle("e-");
	BLAssert(muon != 0 && electron != 0);

	double mass = muon->GetPDGMass() + electron->GetPDGMass() -
							IonizationEnergy;
	double lifetime = muon->GetPDGLifeTime();

//    Arguments for constructor are as follows
//               name             mass          width         charge
//             2*spin           parity  C-conjugation        
//          2*Isospin       2*Isospin3       G-parity        
//               type    lepton number  baryon number   PDG encoding 
//             stable         lifetime    decay table        
	return new Muonium("Mu",  mass,           0.0,          0.0,
	            0,               0,             0,
		    0,               0,             0,
	       "atom",               0,             0,          12345,
	        false,        lifetime,             0);
}

MuoniumDecay::MuoniumDecay(class BLCMDmuonium *p) : 
			G4VProcess("MuoniumDecay",fUserDefined), change()
{
	cmd = p;
	maxProperTime = 0.0;
}

void MuoniumDecay::StartTracking(G4Track *track)
{
	double lifetime = Muonium::getObject()->GetPDGLifeTime();
	maxProperTime = -lifetime * std::log(G4UniformRand());
}

G4double MuoniumDecay::PostStepGetPhysicalInteractionLength(
		const G4Track& track, G4double   previousStepSize,
		G4ForceCondition* condition) 
{
	double properTimeLeft = maxProperTime - track.GetProperTime();
	if(properTimeLeft < 0.0) properTimeLeft = DBL_MIN;
	double betaGamma = track.GetMomentum().mag() / 
					Muonium::getObject()->GetPDGMass();
	if(cmd->verbose >= 1) printf("MuoniumDecay::PostStepGPIL=%.5g\n",
					properTimeLeft * c_light * betaGamma);
	return properTimeLeft * c_light * betaGamma;
}

G4VParticleChange* MuoniumDecay::PostStepDoIt( const G4Track& track,
		const G4Step&  stepData)
{
	static G4ParticleTable *table = G4ParticleTable::GetParticleTable();
	static G4ParticleDefinition *muon = table->FindParticle("decaying_mu+");
	static G4ParticleDefinition *electron = table->FindParticle("e-");

	change.Initialize(track);
	change.ProposeTrackStatus(fStopAndKill);

	// random direction in rest frame 
	double cosTheta = -1.0 + 2.0*G4UniformRand();
	double sinTheta = sqrt(1.0 - cosTheta*cosTheta);
	double phi = 2.0*M_PI*G4UniformRand();
	G4ThreeVector dir = G4ThreeVector(sinTheta*sin(phi),
					  sinTheta*cos(phi),
					  cosTheta);

	// non-rel decay, f = KE(e-) / ( KE(e-)+KE(mu+) )
	double f = electron->GetPDGMass() / muon->GetPDGMass();
	f = 1.0 / (1.0 + f*f);

	// create two DynamicParticles, boost them to lab
	G4DecayProducts prod(*track.GetDynamicParticle());
	G4DynamicParticle *dp = new G4DynamicParticle(electron,dir,
						f*IonizationEnergy);
	prod.PushProducts(dp);
	dp = new G4DynamicParticle(muon,-dp->GetMomentum());
	prod.PushProducts(dp);
	prod.Boost(track.GetTotalEnergy(),track.GetMomentumDirection());

	// add the secondaries, transferring the Mu polarization to the mu+
	G4Track *mu = new G4Track(prod.PopProducts(),
				track.GetGlobalTime(), track.GetPosition());
	mu->SetPolarization(track.GetPolarization());
	G4Track *e = new G4Track(prod.PopProducts(),
				track.GetGlobalTime(), track.GetPosition());
	change.AddSecondary(mu);
	change.AddSecondary(e);

	if(cmd->verbose >= 1) printf("MuoniumDecay: KE=%.5g eV   t=%.3f ns\n",
			track.GetKineticEnergy()/eV,track.GetGlobalTime()/ns);

	return &change;
}

G4double MuoniumDecay::AtRestGetPhysicalInteractionLength(
			const G4Track& track, G4ForceCondition* condition) 
{
	double properTimeLeft = maxProperTime - track.GetProperTime();
	if(properTimeLeft < 0.0) properTimeLeft = DBL_MIN;
	if(cmd->verbose >= 1) printf("MuoniumDecay::AtRestGPIL=%.5g\n",
								properTimeLeft);
	return properTimeLeft;
}

G4VParticleChange* MuoniumDecay::AtRestDoIt( const G4Track& track,
			const G4Step& stepData) 
{
	return PostStepDoIt(track,stepData);
}

MuoniumSurface::MuoniumSurface(class BLCMDmuonium *p) : 
			G4VProcess("MuoniumSurface",fUserDefined), change()
{
	cmd = p;
	mode = TRANSPORT;
	maxPath = DBL_MAX;
}

void MuoniumSurface::StartTracking(G4Track *track)
{
	G4ForceCondition condition;
	maxPath = PostStepGetPhysicalInteractionLength(*track,0.0,&condition);
}

G4double MuoniumSurface::PostStepGetPhysicalInteractionLength(
		const G4Track& track, G4double   previousStepSize,
		G4ForceCondition* condition) 
{
	double density = (track.GetNextMaterial()!= 0 ?
			  track.GetNextMaterial()->GetDensity() :
			  1E-12*(gram/cm3));
	double strip_mfp = cmd->stripMFP/density;
	if(cmd->verbose>=1) printf("MuoniumSurface::PostStepGPIL KE=%.5g eV\n",
						track.GetKineticEnergy()/eV);
	if(density < cmd->minDensity) {
		//@ APPROXIMATION: no model for energy loss, so just TRANSPORT
		mode = TRANSPORT;
		if(cmd->verbose >= 1) printf("MuoniumSurface::PostStepGPIL: "
					"transport in low-density material\n");
		return DBL_MAX;
	}

	double KE = track.GetKineticEnergy();
	if(KE <= cmd->E1) {
		mode = ADSORB; // changed below
	} else if(KE >= cmd->E2) {
		mode = STRIP;
	} else {
		if(G4UniformRand() < (KE - cmd->E1)/(cmd->E2 - cmd->E1))
			mode = ADSORB;
		else
			mode = STRIP;
	}

	if(mode == ADSORB) {
		if(KE <= cmd->E3) {
			mode = ADSORB;
		} else if(KE >= cmd->E4) {
			mode = REFLECT;
		} else {
		if(G4UniformRand() < (KE - cmd->E3)/(cmd->E4 - cmd->E3))
			mode = ADSORB;
		else
			mode = REFLECT;
		}
	}

	switch(mode) {
	case TRANSPORT:
		if(cmd->verbose >= 1) printf("MuoniumSurface::PostStepGPIL: "
								"transport\n");
		return DBL_MAX;
	case STRIP:
		if(cmd->verbose >= 1) printf("MuoniumSurface::PostStepGPIL: "
					"strip MFP=%.5g mm\n",strip_mfp/mm);
		return strip_mfp;
	case ADSORB:
		if(cmd->verbose >= 1) printf("MuoniumSurface::PostStepGPIL: "
								"adsorb\n");
		return DBL_MIN;
	case REFLECT:
		if(cmd->verbose >= 1) printf("MuoniumSurface::PostStepGPIL: "
								"reflect\n");
		return DBL_MIN;
	}
	// never get here but compiler complains
	return DBL_MAX;
}

G4VParticleChange* MuoniumSurface::PostStepDoIt( const G4Track& track,
		const G4Step&  stepData)
{
	static G4ParticleTable *table = G4ParticleTable::GetParticleTable();
	static G4ParticleDefinition *muon = table->FindParticle("mu+");
	static G4ParticleDefinition *electron = table->FindParticle("e-");

	change.Initialize(track);

	switch(mode) {
	case TRANSPORT:
		break;
	case STRIP:
		{ change.ProposeTrackStatus(fStopAndKill);
		  G4ThreeVector dir = track.GetMomentumDirection();
		  // electron velocity = 0
		  G4DynamicParticle *dp = new G4DynamicParticle(electron,dir,
									0.0);
		  G4Track *e = new G4Track(dp, track.GetGlobalTime(),
		  					track.GetPosition());
		  // muon velocity = Mu velocity
		  double KE = track.GetKineticEnergy() * muon->GetPDGMass() /
		  			Muonium::getObject()->GetPDGMass();
		  dp = new G4DynamicParticle(muon,dir,KE);
		  G4Track *mu = new G4Track(dp, track.GetGlobalTime(),
		  					track.GetPosition());
		  mu->SetPolarization(track.GetPolarization());
		  change.AddSecondary(mu);
		  change.AddSecondary(e);
		  if(cmd->verbose >= 1) printf("MuoniumSurface::PostStepDoIt: "
		  						"strip\n");
		}
		break;
	case ADSORB:
		change.ProposeVelocity(0.0);
		change.ProposeTrackStatus(fStopButAlive);
		  if(cmd->verbose >= 1) printf("MuoniumSurface::PostStepDoIt: "
		  						"adsorb\n");
		break;
	case REFLECT:
		{ G4bool valid;
		  // get normal (from G4OpBoundaryProcess.cc)
		  G4int hNavId = G4ParallelWorldProcess::GetHypNavigatorID();
		  std::vector<G4Navigator*>::iterator iNav =
		  	G4TransportationManager::GetTransportationManager()->
						GetActiveNavigatorsIterator();
		  G4ThreeVector normal = (iNav[hNavId])->
		  		GetGlobalExitNormal(track.GetPosition(),&valid);
		  if(valid && normal.mag2() > 1E-3) {
			// reflect - whether normal is in or out does not matter
			normal /= normal.mag();
		  	G4ThreeVector dir = track.GetMomentumDirection();
			dir -= 2.0*dir.dot(normal)*normal;
			change.ProposeMomentumDirection(dir);
		  	if(cmd->verbose >= 1) printf("MuoniumSurface::PostStepDoIt: "
		  						"reflect\n");
		  } else {
		  	// no valid normal -- treat as ADSORB
			G4Exception("MuoniumSurface","No valid normal",
								JustWarning,"");
			change.ProposeVelocity(0.0);
			change.ProposeTrackStatus(fStopButAlive);
		  	if(cmd->verbose >= 1) printf("MuoniumSurface::PostStepDoIt: "
		  					"reflect => adsorb\n");
		  }
		}
		break;
	}

	return &change;
}

DecayingMuPlus* DecayingMuPlus::Definition()
{
  if (object !=0) return object;
  const G4String name = "decaying_mu+";
  // search in particle table]
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance ==0) {
  // create particle
  //
  //    Arguments for constructor are as follows
  //               name             mass          width         charge
  //             2*spin           parity  C-conjugation
  //          2*Isospin       2*Isospin3       G-parity
  //               type    lepton number  baryon number   PDG encoding
  //             stable         lifetime    decay table
  //             shortlived      subType    anti_encoding
  anInstance = new G4ParticleDefinition(
                 name, 0.1056583715*GeV, 2.99598e-16*MeV,  +1.*eplus, 
		    1,               0,                0,          
		    0,               0,                0,             
	     "lepton",              -1,                0,        -23456,
		false,      DBL_MIN,             NULL,
             true,           "mu"
              );
   // Bohr Magnetron
   G4double muB =  0.5*eplus*hbar_Planck/(anInstance->GetPDGMass()/c_squared) ;
   
   anInstance->SetPDGMagneticMoment( muB * 1.0011659209);

  //create Decay Table 
  G4DecayTable* table = new G4DecayTable();
  // create a decay channel
  G4VDecayChannel* mode = new G4MuonDecayChannel("mu+",1.00);
  table->Insert(mode);
  anInstance->SetDecayTable(table);
  }
  object = reinterpret_cast<DecayingMuPlus*>(anInstance);
  return object;
}
