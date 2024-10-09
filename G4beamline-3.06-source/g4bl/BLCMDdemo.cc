//	BLCMDdemo.cc
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

#include "G4NistManager.hh"

#include "BLCommand.hh"
#include "BLTime.hh"

class BLCMDdemo : public BLCommand {
	G4String s1;
	G4String s2;
	G4double d1;
	G4double d2;
public:
	BLCMDdemo();

	G4String commandName() { return "demo"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	void display() { printf("s1='%s' s2='%s' d1=%g d2=%g\n",
				s1.c_str(),s2.c_str(),d1,d2); }
};

BLCMDdemo defaultDemoCommand;	// registers the command, and holds
					// default values of the arguments.

BLCMDdemo::BLCMDdemo()
{
	registerCommand(BLCMDTYPE_OTHER);
	setSynopsis("demo command.");
	setDescription("This demo command takes both positional and named args,\n"
		"and is the prototype class for all commands.\n"
		"All argument values are merely displayed.\n"
		"'demo default name=value...' sets default values.");
}

int BLCMDdemo::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	// handle setting default argument values
	if(argv.size() > 0 && argv[0] == "default") {
		handleNamedArgs(namedArgs);
		return 0;
	}

	// test the time() function (problem on Windows?)
	if(argv.size() > 0 && argv[0] == "time") {
		time_t start = BLTime::time();
		for(;;) {
			BLTime::sleepms(200);
			printf("time=%ld\r",BLTime::time()-start);
			fflush(stdout);
		}
	}

	BLCMDdemo *p = new BLCMDdemo(defaultDemoCommand);
	p->handleNamedArgs(namedArgs);
	for(unsigned int i=0; i<argv.size(); ++i)
		printf("arg%d='%s' ",i,argv[i].c_str());
	p->display();
	delete p;

	G4NistManager *m = G4NistManager::Instance();
	const std::vector<G4String> *v=&m->GetNistMaterialNames();
	for(unsigned i=0; i<v->size(); ++i)
		printf("%s\n",(*v)[i].c_str());

	m->FindOrBuildMaterial("G4_lH2");
	m->FindOrBuildMaterial("G4_Be");
	m->FindOrBuildMaterial("G4_Al");
	m->FindOrBuildMaterial("G4_Fe");
	m->FindOrBuildMaterial("G4_Cu");
	m->FindOrBuildMaterial("G4_W");
	m->FindOrBuildMaterial("G4_Pb");
	m->FindOrBuildMaterial("G4_AIR");
	m->FindOrBuildMaterial("G4_Galactic");
	m->FindOrBuildMaterial("G4_WATER");

	return 0;
}

void BLCMDdemo::defineNamedArgs()
{
	argString(s1,"s1","a demo string argument.");
	argString(s2,"s2","a demo string argument.");
	argDouble(d1,"d1","a demo double argument.");
	argDouble(d2,"d2","a demo double argument.");
}
