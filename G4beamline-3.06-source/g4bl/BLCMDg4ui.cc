//	BLCMDg4ui.cc
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

#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UIterminal.hh"

#include "BLManager.hh"
#include "BLParam.hh"

const int N_TYPE=5;	// # callback types

/**	class BLCMDg4ui arranges to present the user with the standard
 *	Geant4 user interface session.
 **/
class BLCMDg4ui : public BLCommand, public BLCallback {
	bool alreadyRegistered;
	std::vector<G4String> cmds[N_TYPE];
	int when;
public:
	/// Constructor.
	BLCMDg4ui();

	/// commandName() returns "g4ui".
	virtual G4String commandName() { return "g4ui"; }
	
	/// command() implements the g4ui command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// callback() from BLCallback.
	void callback(int type);
};

BLCMDg4ui defaultG4ui;

BLCMDg4ui::BLCMDg4ui()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("Accesses the Geant4 user interface");
	setDescription("Each positional argument is executed as a Geant4 UI\n"
		"command, according to the when parameter. No positional\n"
		"arguments means open a UI session on stdout/stdin.\n"
		"For a given value of when, the UI commands from all g4ui "
		"commands are executed in order.");

	alreadyRegistered = false;
	when = 0;
}

int BLCMDg4ui::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	when = 0;
	int retval = handleNamedArgs(namedArgs);

	if(when < 0 || when == 3 || when >= N_TYPE) {
		printError("g4ui: invalid when value -- command ignored");
		return -1;
	}

	if(argv.size() == 0) {
		cmds[when].push_back("--interactive--");
	} else {
		for(unsigned i=0; i<argv.size(); ++i)
			cmds[when].push_back(argv[i]);
	}

	if(!alreadyRegistered) {
		BLManager::getObject()->registerCallback(this, 0);
		BLManager::getObject()->registerCallback(this, 1);
		BLManager::getObject()->registerCallback(this, 2);
		BLManager::getObject()->registerCallback(this, 4);
		alreadyRegistered = true;
	}

	print("");

	return retval;
}

void BLCMDg4ui::defineNamedArgs()
{
	argInt(when,"when","0=before reference, 1=before beam, 2=after beam, "
				"3=cannot be used, 4=visualization.");
}

void BLCMDg4ui::callback(int _type)
{
	if(_type < 0 || _type >= N_TYPE || cmds[_type].size() == 0) return;

	static const char *where[N_TYPE] = {"before reference",
						"before beam", "after beam"};
	printf("\n");
	printf("Geant4 User Interface %s\n", where[_type]);
	fflush(stdout);

	G4UImanager *UI = G4UImanager::GetUIpointer();
	for(unsigned i=0; i<cmds[_type].size(); ++i) {
		if(cmds[_type][i] == "--interactive--") {
			printf("Type exit<CR> to exit the geant4 UI and resume G4beamline\n");
			fflush(stdout);
			G4UIterminal *terminal = new G4UIterminal(0,false);
			terminal->SessionStart();
			delete terminal;
		} else {
			printf("    %s\n",cmds[_type][i].c_str());
			fflush(stdout);
			int code = UI->ApplyCommand(cmds[_type][i]);
			if(code > 0)
				printf("Command Failed code=%d\n",code);
		}
	}
}
