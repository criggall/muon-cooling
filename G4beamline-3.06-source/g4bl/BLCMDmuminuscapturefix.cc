//	BLCMDmuminuscapturefix.cc
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

#include <stdio.h>

#include "G4VProcess.hh"
#include "G4UImanager.hh"
#include "G4ProcessVector.hh"
#include "G4ProcessManager.hh"
#include "G4Neutron.hh"
#include "G4MuonMinusCapture.hh"
#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

/// Default values are for Al, based on agreement with MARS.
const double DEFAULT_NEUTRONMEANNUMBER=2.5;
const double DEFAULT_NEUTRONMEANENERGY=10.0*MeV;

/**	class MuMinusCaptureFix will fixup the G4MuonMinusCapture
 *	process.
 *
 *	The constructor replaces muMinusCaptureAtRest with this class, so all
 *	that is needed is to create an instance of this class with the desired
 *	arguments. That must occur AFTER the physics processes are completely
 *	set up but before any tracking.
 *
 *	This class adds additional neutrons to mu- capture. The neutrons are
 *	added with a Poisson distribution with mean=neutronMeanNumber,
 *	and have an exponential distribution in kinetic energy:
 *		(1/neutronMeanEnergy)*exp(-KE/neutronMeanEnergy)
 *	As the muonic atom cascades to its ground state it forgets the
 *	incident mu- direction, so neutrons are generated isotropically
 *	in the lab.
 *
 *	The extra neutrons are added only to those captures that are hadronic
 *	(i.e. not decay in orbit). The value of neutronMeanNumber should
 *	reflect this. The default values correspond to Al.
 *
 *	Author: Tom Roberts
 *	Date: September 8, 2010
 *	This class is completely independent of G4beamline and can be copied
 *	into other code (keep the above include-s and definitions).
 **/
class MuMinusCaptureFix : public G4MuonMinusCapture {
	double neutronMeanNumber;
	double neutronMeanEnergy;
	G4ParticleChange change;
public:
	/// Constructor.
	MuMinusCaptureFix(double _neutronMeanNumber=DEFAULT_NEUTRONMEANNUMBER,
			  double _neutronMeanEnergy=DEFAULT_NEUTRONMEANENERGY);

	/// AtRestDoit() from G4MuonMinusCapture.
	virtual G4VParticleChange* AtRestDoIt(const G4Track&, const G4Step&);

	/// createExtraNeutron() creates an extra neutron with the given KE 
	/// distribution, isotropically at position.
	G4Track *createExtraNeutron(const G4ThreeVector &position, double time);
};

MuMinusCaptureFix::MuMinusCaptureFix(double _neutronMeanNumber, double _neutronMeanEnergy)
			: G4MuonMinusCapture(), change()
{
	neutronMeanNumber = _neutronMeanNumber;
	neutronMeanEnergy = _neutronMeanEnergy;

	// replace G4MuonMinusCapture with this MuMinusCaptureFix.
	G4ParticleDefinition *mu = G4MuonMinus::MuonMinus();
	G4ProcessManager *pm = mu->GetProcessManager();
again:	G4ProcessVector *pv = pm->GetProcessList();
	int n = pv->size();
	for(int i=0; i<n; ++i) {
		G4VProcess *process = (*pv)[i];
		if(process->GetProcessName() == "muMinusCaptureAtRest") {
			pm->RemoveProcess(process);
			goto again;
		}
	}
	pm->AddProcess(this,ordDefault,-1,-1);
	theProcessName = "MuMinusCaptureFix";
}


G4VParticleChange* MuMinusCaptureFix::AtRestDoIt(const G4Track& track, 
							const G4Step& step)
{
	int nExtraNeutrons = CLHEP::RandPoissonQ::shoot(neutronMeanNumber);

	// localtime will become the earliest time of any secondary hadron
	// from the original process; extra neutrons are not added to
	// decay in orbit captures (events with no secondary hadrons).
	// Note the capture gives prompt gammas which must not be considered.
	double localtime=DBL_MAX;
	bool haveHadron=false;

	// call original process's AtRestDoit()
	G4VParticleChange *p=G4MuonMinusCapture::AtRestDoIt(track,step);
	int n0=p->GetNumberOfSecondaries();
	G4ThreeVector position=track.GetPosition();

	// setup a ParticleChange with room for the extra neutrons (this is a
	// copy of all modifications to aParticleChange in the original code)
	change.Initialize(track);
	change.SetNumberOfSecondaries(n0+nExtraNeutrons);
	for(int i=0; i<n0; ++i) {
		G4Track *t=p->GetSecondary(i);
		change.AddSecondary(new G4Track(*t)); // avoid double deletes
		G4String type=t->GetDefinition()->GetParticleType();
		if(type == "baryon" || type == "meson" || type == "nucleus") {
			haveHadron = true;
			double time=t->GetGlobalTime();
			if(time < localtime) localtime = time;
		}
	}
	change.ProposeLocalEnergyDeposit(0.0);
	change.ProposeTrackStatus(fStopAndKill);
	p->Clear();

	// store extra neutrons
	if(haveHadron) {
	    for(int i=0; i<nExtraNeutrons; ++i) {
		G4Track *neutron = createExtraNeutron(position,localtime);
		change.AddSecondary(neutron);
	    }
	}

	return &change;
}

G4Track *MuMinusCaptureFix::createExtraNeutron(const G4ThreeVector &position,
								double time)
{
	G4ParticleDefinition *neutron=G4Neutron::Neutron();
	double KE=CLHEP::RandExponential::shoot(neutronMeanEnergy);
	double costheta=2.0*(G4UniformRand()-0.5);
	double sintheta=sqrt(1.0-costheta*costheta);
	double phi=2.0*M_PI*G4UniformRand();
	G4ThreeVector direction(sintheta*cos(phi),sintheta*sin(phi),costheta);
	G4DynamicParticle *dyn = new G4DynamicParticle(neutron,direction,KE);
	G4Track *track = new G4Track(dyn,time,position);
	return track;
}





// Below is the code for the Gbeamline command to apply MuMinusCaptureFix
// in G4beamline.

#include "BLManager.hh"
#include "BLCommand.hh"
#include "BLCallback.hh"

/**	class BLCMDmuminuscapturefix fixes up the G4MuonMinusCapture
 *	process.
 *
 **/
class BLCMDmuminuscapturefix : public BLCommand, public BLCallback {
	double neutronMeanNumber;
	double neutronMeanEnergy;
	static bool first;
public:
	/// Constructor.
	BLCMDmuminuscapturefix();

	/// commandName() returns "muminuscapturefix".
	virtual G4String commandName() { return "muminuscapturefix"; }
	
	/// command() implements the muminuscapturefix command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// callback() from BLCallback.
	void callback(int type);

};
BLCMDmuminuscapturefix defaultMuminuscapturefix;
bool BLCMDmuminuscapturefix::first=true;

BLCMDmuminuscapturefix::BLCMDmuminuscapturefix()
{
	registerCommand(BLCMDTYPE_PHYSICS);
	setSynopsis("Fixes up the G4MuonMinusCapture process.");
	setDescription(" This class adds extra neutrons to mu- capture. "
		"The neutrons are added with a Poisson distribution having "
		"a mean of neutronMeanNumber, and with an exponential "
		"distribution in kinetic energy:\n"
		"    (1/neutronMeanEnergy)*exp(-KE/neutronMeanEnergy)\n"
		"As the muonic atom cascades to its ground state it forgets "
		"the incident mu- direction, so the extra neutrons are "
		"generated isotropically in the lab.\n\n"
		"The extra neutrons are added only to those captures that are "
		"hadronic (i.e. not decay in orbit). The value of "
		"neutronMeanNumber should reflect this.\n\n"
		"The default values correspond to Aluminum.\n\n"
		"NOTE: the Geant4 processes have changed since this command "
		"was designed and validated, and they believe this command is "
		"no longer necessary; the user must ensure correct physics is "
		"implemented.");

	// values for neutrons from mu- stopping in Al
	// this is the "high-energy tail" that is missing from the G4 process.
	neutronMeanNumber = DEFAULT_NEUTRONMEANNUMBER;
	neutronMeanEnergy = DEFAULT_NEUTRONMEANENERGY;
}

int BLCMDmuminuscapturefix::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(!first) {
		printError("muminuscapturefix: Multiple commands are ignored");
		return -1;
	}
	first = false;

	int retval = handleNamedArgs(namedArgs);

	BLManager::getObject()->registerCallback(this,0);

	print("");

	return retval;
}

void BLCMDmuminuscapturefix::defineNamedArgs()
{
	argDouble(neutronMeanNumber,"neutronMeanNumber",
		"Mean mumber of extra neutrons per nuclear capture (2.5).");
	argDouble(neutronMeanEnergy,"neutronMeanEnergy",
		"Mean energy of neutron spectrum (MeV)",MeV);
}

void BLCMDmuminuscapturefix::callback(int type)
{
	(void)new MuMinusCaptureFix(neutronMeanNumber,neutronMeanEnergy);
}
