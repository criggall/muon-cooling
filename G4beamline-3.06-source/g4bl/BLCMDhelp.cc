//	BLCMDhelp.cc
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

#include "BLCommand.hh"
#include "BLCommandAlias.hh"
#include "BLManager.hh"
#include "BLParam.hh"
#include "BLParam.hh"

extern void g4bl_exit(int);

/**	class BLCMDhelp implements the help command.
 *
 **/
class BLCMDhelp : public BLCommand {
public:
	BLCMDhelp();

	G4String commandName() { return "help"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	// no arguments.
};

BLCMDhelp defineHelp;

BLCommandAlias aliasHelp("man",defineHelp);	// Alias "man"

BLCMDhelp::BLCMDhelp()
{
	registerCommand(BLCMDTYPE_HELP);
	setSynopsis("provides interactive help.");
	setDescription("help with no arguments lists all commands.\n"
		"'help command' gives more detailed help on that command.\n"
		"'help *' gives detailed help for all commands.\n"
		"If the first argument is -exit\n"
		"the program will exit after printing the help.");
}

int BLCMDhelp::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	bool exit = (argv.size() > 0 && argv[0] == "-exit");
	if(exit) argv.erase(argv.begin());

	if(argv.size() == 0) {
		for(int j=0; j<=(int)BLCMDTYPE_OTHER; ++j) {
			switch((BLCmdType)j) {
			case BLCMDTYPE_HELP:
			    printf("\nThe help command:\n");
			    break;
			case BLCMDTYPE_CONTROL:
			    printf("\nProgram control commands:\n");
			    break;
			case BLCMDTYPE_LAYOUT:
			    printf("\nCenterline layout commands:\n");
			    break;
			case BLCMDTYPE_BEAM:
			    printf("\nBeam definition commands:\n");
			    break;
			case BLCMDTYPE_AUX:
			    printf("\nAuxiliary definition commands:\n");
			    break;
			case BLCMDTYPE_CUTS:
			    printf("\nTrack and Event cuts:\n");
			    break;
			case BLCMDTYPE_ELEMENT:
			    printf("\nBeamline element definition commands:\n");
			    break;
			case BLCMDTYPE_PLACE:
			    printf("\nThe place command:\n");
			    break;
			case BLCMDTYPE_DATA:
			    printf("\nData output commands:\n");
			    break;
			case BLCMDTYPE_PHYSICS:
			    printf("\nPhysics commands:\n");
			    break;
			case BLCMDTYPE_OTHER:
			    printf("\nOther commands:\n");
			    break;
			}
			std::map<G4String,BLCommand*>::iterator i;
			for(i=(*BLCommand::mapCommand).begin();
				    i != (*BLCommand::mapCommand).end(); ++i) {
				if((BLCmdType)j == i->second->getCmdType()) {
					printf("    ");
					i->second->help(false);
				}
			}
		}
		printf("\nProgram control Parameters:\n%s",
						Param.getHelpText().c_str());
		printf("steppingFormat is a space- or comma-separated list of items:\n%s\n\n",
			BLManager::getObject()->getFormatHelp().c_str());
	} else if(argv[0] == "*") {
		std::map<G4String,BLCommand*>::iterator i;
		for(i=(*BLCommand::mapCommand).begin();
				i != (*BLCommand::mapCommand).end(); ++i) {
			i->second->help(true);
			printf("\n");
		}
	} else {
		for(unsigned i=0; i<argv.size(); ++i) {
			if((*BLCommand::mapCommand).count(argv[i]) == 0) {
				printf("Command '%s' not found.\n",
							argv[i].c_str());
			} else {
				(*BLCommand::mapCommand)[argv[i]]->help(true);
			}
		}
	}

	if(exit) {
		fflush(stdout);
		g4bl_exit(0);
	}

	return 0;
}
