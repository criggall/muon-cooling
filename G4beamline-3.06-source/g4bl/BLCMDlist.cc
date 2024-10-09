//	BLCMDlist.cc
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

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "BLCommand.hh"
#include "BLManager.hh"

extern void g4bl_exit(int);

/**	class BLCMDlist implements the list command.
 *
 **/
class BLCMDlist : public BLCommand {
	G4String particle;
public:
	BLCMDlist();

	G4String commandName() { return "list"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	void listParticles();

	void listParticleDetails();
};

BLCMDlist defineList;

BLCMDlist::BLCMDlist()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("provides interactive list of interesting internal tables.");
	setDescription("list with no arguments lists all lists except processes.\n"
		"'list name' lists that specific one.\n"
		"'list -exit name(s)' will exit after listing.\n"
		"List names are:\n"
		"    commands    all commands\n"
		"    materials   currently known materials\n"
		"    physics     all physics lists\n"
		"    particles   currently known particles\n"
		"    processes   currently known physics processes ***\n"
		"NOTE: the particles and processes lists are "
		"not populated until the physics list is selected (via "
		"the physics command). Different physics lists use different "
		"processes and particles.\n\n"
		"***NOTE: listing processes will prevent any simulating, "
		"as will a non-empty particle list.");
}


int BLCMDlist::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	bool exit = (argv.size() > 0 && argv[0] == "-exit");
	if(exit) argv.erase(argv.begin());

	handleNamedArgs(namedArgs);

	if((argv.size() == 0 || argv[0] == "*" || argv[0] == "all") &&
						particle.size() == 0) {
		argv.clear();
		argv.push_back("commands");
		argv.push_back("materials");
		argv.push_back("physics");
		argv.push_back("particles");
		//argv.push_back("processes");
	}

	for(unsigned i=0; i<argv.size(); ++i) {
		if(argv[i].compareTo("commands",G4String::ignoreCase) == 0) {
			printf("\nCommands:\n---------\n");
			G4String s("help");
			doCommand(s);
		} else if(argv[i].compareTo("particles",G4String::ignoreCase) == 0) {
			printf("\nParticles:\n----------\n");
			listParticles();
		} else if(argv[i].compareTo("physics",G4String::ignoreCase) == 0) {
			printf("\nPhysics lists:\n--------------\n");
			G4String s("physics");
			doCommand(s);
		} else if(argv[i].compareTo("materials",G4String::ignoreCase) == 0) {
			G4String s("material");
			doCommand(s);
		} else if(argv[i].compareTo("processes",G4String::ignoreCase) == 0) {
			if(BLManager::getObject()->getPhysics() == 0) {
				printf("cannot list processes -- no physics registered.\n");
				return 0;
			}
			if(!BLManager::isInitialized()) {
				BLManager::getObject()->initialize();
				BLManager::getObject()->handleCallbacks(0);
			}
			printf("\nlist: Note that simulation will fail due to multiple initializations.\n");
			printf("Physics processes:\n------------------\n");
			G4String s("/process/list");
			doCommand(s);
		}
	}

	if(particle.size() > 0) {
		printf("\nlist: Note that simulation will fail due to multiple initializations.\n");
		printf("Particle details:\n----------------\n");
		if(!BLManager::isInitialized()) {
			BLManager::getObject()->initialize();
			BLManager::getObject()->handleCallbacks(0);
		}
		listParticleDetails();
	}

	if(exit) {
		fflush(stdout);
		g4bl_exit(0);
	}

	particle = ""; // in case another list command is issued

	return 0;
}

void BLCMDlist::defineNamedArgs()
{
	argString(particle,"particle","Comma-separated list of particles for "
				"which details will be printed");
}

void BLCMDlist::listParticles()
{
	if(BLManager::getObject()->getPhysics() == 0) {
		printf("cannot list particles -- no physics registered.\n");
		return;
	}

	G4ParticleTable::G4PTblDicIterator *i =
			G4ParticleTable::GetParticleTable()->GetIterator();
	i->reset();
	int n=0;
	while((*i)()) {
		G4ParticleDefinition *pd = i->value();
		printf("%15s PDGid=%d\n", pd->GetParticleName().c_str(),
			 pd->GetPDGEncoding());
		++n;
	}
	if(n == 0)
		printf("No particles listed -- select physics list first\n");
}

void BLCMDlist::listParticleDetails()
{
	G4ParticleTable *table = G4ParticleTable::GetParticleTable();
	std::vector<G4String> v = splitString(particle,',');
	for(unsigned i=0; i<v.size(); ++i) {
		G4ParticleDefinition *pd = table->FindParticle(v[i]);
		if(!pd) {
			printError("list: cannot find particle '%s'\n",
								v[i].c_str());
			continue;
		}
		G4ProcessManager *pmgr = pd->GetProcessManager();
		if(!pmgr) continue;
		G4String line(v[i] + " ");
		while(line.size() < 10) line += " ";
		printf("%sPDGid=%d  mass=%.3f MeV/c^2\n",line.c_str(),
			pd->GetPDGEncoding(),pd->GetPDGMass());
		G4ProcessVector *pv = pmgr->GetProcessList();
		line = "Processes: ";
		for(int i=0; i<pv->size(); ++i) {
			G4String d("");
			if(!pmgr->GetProcessActivation(i)) d = "(disabled)";
			line += (*pv)[i]->GetProcessName() + d + " ";
		}
		printf("%s",wrapWords(line,"          ","                     ").c_str());
	}


/*** 
	// add "bookends" to simplify the comparison
	G4String comma(",");
	particle = comma + particle + comma;

	G4ParticleTable::G4PTblDicIterator *myParticleIterator
			= G4ParticleTable::GetParticleTable()->GetIterator();
	myParticleIterator->reset();
	while((*myParticleIterator)()) {
		G4ParticleDefinition *pd = myParticleIterator->value();
		G4ProcessManager *pmgr = pd->GetProcessManager();
		if(!pmgr) continue;
		G4String name = pd->GetParticleName();
		if(particle.find(comma+name+comma) == particle.npos) continue;
		G4String line(name + " ");
		while(line.size() < 10) line += " ";
		printf("%sPDGid=%d  mass=%.3f MeV/c^2\n",line.c_str(),
			pd->GetPDGEncoding(),pd->GetPDGMass());
		G4ProcessVector *pv = pmgr->GetProcessList();
		line = "Processes: ";
		for(int i=0; i<pv->size(); ++i) {
			if(!pmgr->GetProcessActivation(i)) continue;
			line += (*pv)[i]->GetProcessName() + " ";
		}
		printf("%s",wrapWords(line,"          ","                     ").c_str());
	}
***/
}

