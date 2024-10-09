//	BLCMDparam.cc
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

#include <math.h>

#include "globals.hh"
#include "BLParam.hh"
#include "BLEvaluator.hh"

#ifdef WIN32
#define snprintf _snprintf
#endif
#ifdef __CYGWIN__
#include "mysnprintf.hh"
#endif


/**	class BLCMDparam implements the param command.
 **/
class BLCMDparam : public BLCommand {
	static std::map<G4String,G4String> *paramMap;
	void init() {if (!paramMap) paramMap = new std::map<G4String,G4String>;}
public:
	/// Constructor. registers the command and provides help text.
	BLCMDparam();

	/// commandName() returns "param";
	virtual G4String commandName() { return "param"; }

	/// command() implements the param command. It does not use the
	/// standard handleNamedArgs(), but accepts any named argument
	/// and defines a parameter of that name. Unnamed args are ignored
	/// with a warning. If no args are present, calls printParam().
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// printParam() will display all parameters on stdout.
	void printParam();

	virtual void help(bool detailed) {
		BLCommand::help(detailed);
		if(detailed) {
		    printf("\nProgram control Parameters:\n%s",
						Param.getHelpText().c_str());
		    printf("steppingFormat is a space- or comma-separated list of items:\n%s\n\n",
			BLManager::getObject()->getFormatHelp().c_str());
		}
	}
};

BLCMDparam defaultParam;

BLCMDparam::BLCMDparam()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("Defines parameter values.");
	setDescription("Parameters are named strings used to parameterize "
	    "the input file; a few are used for program control.\n"
	    "Parameters are set by the param commend, and on the command "
	    "line (all program arguments after the first are interpreted "
	    "as name=value).\n\n"
	    "'param name=value ...' defines parameters.\n"
	    "If no arguments are present, all parameters are displayed. "
	    "If the first argument is '-unset', this command will not\n"
	    "overwrite parameters that are already set (e.g. from the\n"
	    "command line - this permits the input.file to set defaults "
	    "that can be overridden on the command line).\n\n"
	    "Parameters are expanded only in the arguments of commands:\n"
	    "'cmd argname=$paramname [...]'; real-valued expressions for "
	    "arguments can use $paramname as a value, as long as paramname "
	    "contains a valid real expression.\n\n"
	    "Parameters are most useful for setting global things like viewer\n"
	    "and histoFile, or as parameters used in the input.file.\n\n"
	    "When a parameter is used, if it has not been defined, it will "
	    "be defined from the environment if possible; if it is not defined "
	    "in the environment then this generates an error message.\n\n"
	    "The values of parameters are strings, but if the value of a "
	    "parameter is set to a valid real expression including at least "
	    "one operator {+-*/^<>=()!~&%|?}, the parameter value will "
	    "be set to the numerical value of the expression to 8 significant "
	    "digits.\n\n"
	    "NOTE: pre-defined Program control parameters (listed below) are "
	    "defined before the command-line and are not affected by -unset.");
}

int BLCMDparam::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	bool unset = false;

	if(argv.size() > 0 && argv[0] == "-unset") {
		unset = true;
	}
	if(argv.size() > (unset ? 1 : 0)) {
		G4Exception("param command","un-named arguments ignored",
					JustWarning, "");
	}

	if(namedArgs.size() > 0) {
		BLArgumentMap::iterator i;
		for(i=namedArgs.begin(); i!=namedArgs.end(); ++i) {
			bool defined = Param.isDefined(i->first);
			if(!unset || !defined) {
				// check if value is a valid expression
				{ BLEvaluator e;
				  G4double v = e.evaluate(i->second);
				  // e has many symbols that conflict with
				  // simple symbols (e.g. "C" for Coulomb), so
				  // require an operator before evaluating as
				  // an expression
				  if(e.status() == HepTool::Evaluator::OK &&
				     i->second.find_first_of("+-*/^<>=()!~&%|?")
							!= i->second.npos) {
					char tmp[64];
					snprintf(tmp,sizeof(tmp),"%.8g",v);
					Param.setParam(i->first,tmp);
				  } else {
					Param.setParam(i->first,i->second);
				  }
				}
			}
			// print param definition
			G4String indent1(commandName());
			do { indent1 += " "; }
				while(indent1.size() < IndentDesc.size());
			indent1 += i->first;
			do { indent1 += " "; }
				while(indent1.size() < IndentArg.size());
			printf("%s%s%s\n",indent1.c_str(),
				Param.getString(i->first).c_str(),
				(unset&&defined ? "  (already defined)" : ""));
		}
	} else {
		printParam();
	}

	return 0;
}

void BLCMDparam::printParam()
{
	Param.printParam();
}
