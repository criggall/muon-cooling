//	BLCMDsetdecay.cc
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
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VDecayChannel.hh"
#include "G4PhaseSpaceDecayChannel.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "BLAssert.hh"
#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLCallback.hh"
#include "BLParam.hh"

extern void g4bl_exit(int);

/**	class BLCMDsetdecay -- implement the setdecay command
 **/
class BLCMDsetdecay : public BLCommand {
public:
	/// Constructor.
	BLCMDsetdecay();

	/// Destructor.
	~BLCMDsetdecay() { }

	/// commandName() returns "setdecay"
	G4String commandName()  { return "setdecay"; }

	/// command() executes the command associated with this element.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command. 
	/// Because arg names can be decay modes, this is handled in command().
	void defineNamedArgs() { }
};
BLCMDsetdecay defaultSetDecay;

class SetDecayInstance : public BLCallback {
	G4String particleName;
	double lifetime;
	std::vector<G4String> channel;
	std::vector<double> br;
	friend class BLCMDsetdecay;
public:
	SetDecayInstance(G4String _particleName) : BLCallback(), channel(), br()
		{ particleName=_particleName; lifetime=-1.0; }

	/// callack() from BLCallback.
	void callback(int type);

	// isMatch() returns true if name matches chan's daughters
	bool isMatch(G4String name, G4VDecayChannel *chan);

	/// decayChannelName() returns the name of the channel.
	G4String decayChannelName(const G4VDecayChannel *channel);

	/// fatalError() arranges to abort the simulation after all 
	/// callbacks have run.
	void fatalError();
};


BLCMDsetdecay::BLCMDsetdecay() : BLCommand()
{
	registerCommand(BLCMDTYPE_PHYSICS);
	setSynopsis("Set lifetime, decay channels, and branching ratios for a particle's decay.");
	setDescription("The particle is specified by name as the first "
		"positional argument.\n\n"
		"The lifetime of the particle can be set, unless it is "
		"a short-lived particle (for which lifetime is fixed at 0 -- "
		"these are particles like quarks, Zs, and Ws). "
		"Units are ns.\n\n"
		"Decay channels are specified 'daughter1,daughter2=BR', where "
		"the daughter names are separated by commas, and the branching "
		"ratio is a value between 0 and 1 (inclusive); the order of "
		"daughters does not matter. The sum of all BRs must be 1.0. "
		"It is best to use existing channels for the particle, because "
		"the code for the decay distribution is retained; new decay "
		"channels are given a default phase-space distribution, "
		"which is probably valid only for a 2-body decay of a spin 0 "
		"particle. New channels are limited to 4 daughters. "
		"Note that all desired decay channels must be listed.\n\n"
		"Example to force fast decay (0.1 ns) of pi+ to a positron:\n"
		"    setdecay pi+ lifetime=0.1 e+,nu_e=1.0"
		);
}

int BLCMDsetdecay::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	int retval = 0;

	if(argv.size() < 1) {
		printError("setdecay: no decaying particle name.");
		retval = 1;
	}
	if(argv.size() > 1) {
		printError("setdecay: too many positional arguments.");
		retval = 1;
	}

	SetDecayInstance *instance = new SetDecayInstance(argv[0]);

	BLArgumentMap::iterator i;
	BLEvaluator eval;
	double totalBR=0.0;
	for(i=namedArgs.begin(); i!=namedArgs.end(); ++i) {
		G4String name = i->first;
		G4String value = i->second;
		if(name == "lifetime") {
			instance->lifetime = eval.evaluate(value);
			if(eval.status() != HepTool::Evaluator::OK)
			    printError("setdecay: invalid value for lifetime");
			continue;
		}
		double v = eval.evaluate(value);
		if(eval.status() != HepTool::Evaluator::OK) {
			printError("setdecay: invalid BR for %s",name.c_str());
			continue;
		}
		instance->channel.push_back(name);
		instance->br.push_back(v);
		totalBR += v;
	}
	if(instance->channel.size() > 0 && fabs(totalBR-1.0) > 1.0e-6)
		printError("setdecay: branching ratios do not total to 1");

	print(argv[0],namedArgs);

	BLManager::getObject()->registerCallback(instance,0);

	return retval;
}

void SetDecayInstance::callback(int type)
{
	BLAssert(type==0);

	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition *pd = particleTable->FindParticle(particleName);
	
	if(!pd) {
		BLCommand::printError("setdecay: invalid particle '%s'\n",
						particleName.c_str());
		fatalError();
		return;
	}

	char tmp[128];
	sprintf(tmp,"particle %s",particleName.c_str());
	G4Exception("setdecay","Decay has been modified",JustWarning, tmp);
	printf("setdecay %s",particleName.c_str());
	if(lifetime >= 0.0) {
		if(pd->IsShortLived()) {
			printf(" (lifetime ignored for short-lived particles)");
		} else {
			pd->SetPDGLifeTime(lifetime); 
			printf(" lifetime=%.4g ns",lifetime/ns);
		}
	}
	printf("\n");

	// handle stable particles, adding Decay process if necessary
	if(pd->GetPDGStable())
		pd->SetPDGStable(false);
	if(pd->GetDecayTable() == 0)
		pd->SetDecayTable(new G4DecayTable());
	G4ProcessManager* pmanager = pd->GetProcessManager(); 
	G4ProcessVector *pv = pmanager->GetProcessList();
	bool hasDecay=false;
	for(int i=0; i<pv->entries(); ++i) {
		G4VProcess *p=(*pv)[i];
		if(p->GetProcessName().find("Decay") != G4String::npos) {
			hasDecay = true;
			break;
		}
	}
	if(!hasDecay) {
		if(pd->GetPDGLifeTime() < 0.0) {
			printf("  *** Must specify lifetime!\n");
			fatalError();
			return;
		}
		printf("  Decay process added to previously-stable particle\n");
        	G4VProcess* aDecay = new G4Decay();
        	pmanager->AddProcess(aDecay);
        	pmanager->SetProcessOrdering(aDecay, idxPostStep, 5);
        	pmanager->SetProcessOrdering(aDecay, idxAtRest);
	}

	// no channels specified means keep the existing ones.
	if(channel.size() == 0) {
		printf("\n");
		return;
	}

	// Create new decay table (based on the original decay table).
	// Channels found in orgDecayTable are re-used with updated branching 
	// ratio; channels not found are created ("generic" phase space decay).
	// These loops handle multiple matching channels in orgDecayTable.
	G4DecayTable *orgDecayTable = pd->GetDecayTable();
	G4DecayTable *newDecayTable = new G4DecayTable();
	int nOrgChannels = orgDecayTable->entries();

	// save orgBR because we will be changing the BR
	std::vector<double> orgBRvect;
	for(int i=0; i<nOrgChannels; ++i)
		orgBRvect.push_back(orgDecayTable->GetDecayChannel(i)->GetBR());

	for(unsigned i=0; i<channel.size(); ++i) {
		G4String channelName=channel[i];
		// find total BR for all org channels that match channelName
		double totalBR=0.0;
		for(int j=0; j<nOrgChannels; ++j) {
			G4VDecayChannel *p=orgDecayTable->GetDecayChannel(j);
			if(isMatch(channelName,p))
				totalBR += p->GetBR();
		}
		if(totalBR > 0.0) {
			double ratio = br[i]/totalBR;
			// update channels' BR and insert into newDecayTable
			for(int j=0; j<nOrgChannels; ++j) {
				G4VDecayChannel *p=orgDecayTable->
							GetDecayChannel(j);
				if(isMatch(channelName,p)) {
					p->SetBR(p->GetBR()*ratio);
					newDecayTable->Insert(p);
				}
			}
			continue;
		}
		// Create the channel (default phase space decay -- valid for
		// 2-particle decays, but probably incorrect for 3 or more
		// daughters or non-zero spin). Caveat utilitor.
		std::vector<G4String> daughter = 
					BLCommand::splitString(channelName,",");
		int n=daughter.size();
		if(n > 4) {
			n = 4;
			BLCommand::printError("setdecay: too many daughters - "
						"simulation aborted.");
			fatalError();
		}
		while(daughter.size() < 4) daughter.push_back("");
		newDecayTable->Insert(new G4PhaseSpaceDecayChannel(particleName,
		    br[i],n,daughter[0],daughter[1],daughter[2],daughter[3]));

	}
	pd->SetDecayTable(newDecayTable);

	// print new and org decay tables
	printf("  Original Decay Channels                 New Decay Channels\n");
	for(int i=0; ; ++i) {
		double orgBR=-1.0, newBR=-1.0;
		G4String orgName="", newName="";
		char orgBRstring[32]="", newBRstring[32]="";
		if(i < orgDecayTable->entries()) {
			orgBR = orgBRvect[i];
			orgName = 
			    decayChannelName(orgDecayTable->GetDecayChannel(i));
			sprintf(orgBRstring,"%.6f",orgBR);
		}
		if(i < newDecayTable->entries()) {
			newBR = newDecayTable->GetDecayChannel(i)->GetBR();
			newName = 
			    decayChannelName(newDecayTable->GetDecayChannel(i));
			sprintf(newBRstring,"%.6f",newBR);
		}
		if(orgBR < 0.0 && newBR < 0.0) break;
		printf(" %9.9s %-26.26s    %9.9s %-26.26s\n",
		    orgBRstring,orgName.c_str(),newBRstring,newName.c_str());
	}
	printf("\n");

	// deliberate memory leak: don't delete orgDecayChannel (it might delete
	// the G4DecayChannel-s it uses, which we re-used).
}

bool SetDecayInstance::isMatch(G4String name, G4VDecayChannel *chan)
{
	std::vector<G4String> list = BLCommand::splitString(name,",");
	int nDaughters=chan->GetNumberOfDaughters();
	if(list.size() != (unsigned)nDaughters) return false;
	for(unsigned i=0; i<list.size(); ++i) {
		bool ok=false;
		for(int j=0; j<nDaughters; ++j) {
			if(list[i] == chan->GetDaughterName(j)) {
				ok = true;
				break;
			}
		}
		if(!ok) return false;
	}
	return true;
}

G4String SetDecayInstance::decayChannelName(const G4VDecayChannel *channel)
{
	int n=channel->GetNumberOfDaughters();
	G4String name=channel->GetDaughterName(0);
	for(int i=1; i<n; ++i) {
		name += ",";
		name += channel->GetDaughterName(i);
	}
	return name;
}

void SetDecayInstance::fatalError()
{
	class Abort : public BLCallback {
	public:
		Abort() : BLCallback() { }
		void callback(int type) {
			BLCommand::printError("setdecay: FATAL ERROR -- "
						"simulation aborted.\n");
			g4bl_exit(99);
		}
	};
	BLManager::getObject()->registerCallback(new Abort(),0);
}
