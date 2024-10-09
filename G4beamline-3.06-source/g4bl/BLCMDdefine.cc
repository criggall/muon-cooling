//	BLCMDdefine.cc
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

#include "BLCommand.hh"
#include "BLParam.hh"

/**	class BLCMDdefine implements the "define" command to define a macro
 *	(an argument-expanded set of commands).
 **/
class BLCMDdefine : public BLCommand {
public:
	BLCMDdefine();

	G4String commandName() { return "define"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);
};

BLCMDdefine defaultDefine;	// registers the command

/**	class Macro implements a macro command defined by the define command.
 **/
class Macro : public BLCommand {
	G4String name;
	G4String body;
	static int nExpansions;
public:
	Macro(G4String _name, G4String _body) : name(_name), body(_body) 
		{ deleteCommand(name);
		  registerCommand(BLCMDTYPE_OTHER);
		  setSynopsis("Defined macro.");
		}

	G4String commandName() { return name; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs() { }
};
int Macro::nExpansions = 0;


BLCMDdefine::BLCMDdefine() : BLCommand()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("defines a macro (argument-expanded set of commands).");
	setDescription("The first argument is the macro name, additional arguments\n"
	  "become lines in the body of the expanded macro. The macro name\n"
	  "becomes a command with up to 9 positional arguments. When the\n"
	  "command is issued, the body is expanded and executed, with\n"
	  "these substitutions:\n"
 	  "  $0     MacroName\n"
	  "  $1-$9  Positional arguments of the command\n"
	  "  $#     # macro expansions (for generating unique names)\n"
	  "NOTE: $paramname is expanded in the define command, but $$paramname "
	  "is expanded when the macro is invoked.");
}

int BLCMDdefine::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() < 1) {
		printError("define: need command name being defined");
		return -1;
	}

	// the arguments are the body
	G4String body;
	for(unsigned i=1; i<argv.size(); ++i) {
		if(argv[i].size() > 0) {body += argv[i]; body += "\n";}
	}

	new Macro(argv[0],body);

	// print the definition
	G4String indent1(commandName());
	do { indent1 += " "; } while(indent1.size() < IndentDesc.size());
	printf("%s%s:\n",indent1.c_str(),argv[0].c_str());
	for(unsigned i=1; i<argv.size(); ++i) {
		if(argv[i].size() > 0) 
			printf("%s%s\n",IndentDesc.c_str(),argv[i].c_str());
	}

	return 0;
}

int Macro::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	++nExpansions;

	if(namedArgs.size() > 0) {
		printError("macro: named arguments are ignored");
	}

	// print the macro command, with args
	G4String line(commandName());
	do { line += " "; } while(line.size() < IndentArg.size());
	for(unsigned i=0; i<argv.size(); ++i)
		{line += argv[i]; line += " ";}
	printf("%s\n",line.c_str());

	// ensure argv has 10 entries. Changes argv in caller, but nobody cares.
	while(argv.size() < 10)
		argv.push_back("");

	// replace: $0 with name, $1 - $9 with argv[0-8], $# with nExpansions
	G4String tmp = Param.expand(body);
	int place=0;
	for(;;) {
		G4String::size_type i = tmp.find('$',place);
		if(i == tmp.npos) break;
		if(tmp.c_str()[i+1] == '#') {
			char tmp2[16];
			sprintf(tmp2,"%d",nExpansions);
			tmp.replace(i,2,tmp2);
		} else if(isdigit(tmp.c_str()[i+1])) {
			int j = tmp.c_str()[i+1] - '0';
			if(j == 0)
				tmp.replace(i,2,name);
			else
				tmp.replace(i,2,argv[j-1]);
		} else {
			place = i+1;
		}
	}

	// execute the macro body (now in tmp)
	G4String::size_type i = 0;
	while(i < tmp.size()) {
		G4String::size_type j = tmp.find('\n',i);
		if(j == tmp.npos) break;
		G4String tmp3(tmp.substr(i,j-i));
		BLCommand::doCommand(tmp3);
		i = j + 1;
	}

	return 0;
}
