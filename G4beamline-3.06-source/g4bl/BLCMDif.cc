//	BLCMDif.cc
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

#include <stdlib.h>
#include "BLCommand.hh"
#include "BLParam.hh"
#include "BLEvaluator.hh"

/**	class BLCMDif implements the if command for conditionals
 *
 *	Implements both a single-command if that can execute command(s),
 *	and if / elseif / else /endif
 **/
class BLCMDif : public BLCommand {
public:
	BLCMDif();

	G4String commandName() { return "if"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);
};

BLCMDif defineIf;

BLCMDif::BLCMDif()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("Conditional execution of command(s), and if/elseif/else/endif.");
	setDescription("Single-line format:\n"
		"    if $i==1 CMD1 CMD2 ...\n\n"
		"If the expression is true (nonzero), the commands are "
		"executed. The commands usually need to be quoted.\n\n"
		"Multi-line format:\n"
		"    if $i==1\n"
		"        CMD1 ...\n"
		"    elseif $i==2\n"
		"        CMD2 ...\n"
		"    else\n"
		"        CMD3 ...\n"
		"    endif\n\n"
		"The commands are executed in the usual way; elseif and else "
		"are optional, but endif is mandatory. Any number of "
		"elseif-s and commannds can be used. Do-s, "
		"for-s, and multi-line if-s can be nested in any manner.\n\n"
		"If there are spaces in the expression, it must be quoted.");
}

int BLCMDif::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() < 1 || namedArgs.size() != 0) {
		printError("Invalid if command");
		return -1;
	}

	BLEvaluator eval;
	int retval=0;

	if(argv.size() > 1) {
		// simple (old) if
		double val = eval.evaluate(argv[0]);
		if(eval.status() != HepTool::Evaluator::OK) {
			printError("if: Invalid expression '%s'\n",
							argv[0].c_str());
			return -1;
		}
		if(fabs(val) > 1.0e-12) {
			for(unsigned i=1; i<argv.size(); ++i) {
				if(doCommand(argv[i]) != 0)
					retval = -1;
			}
		}
	} else {
		// complex (multi-line) if
		double val = eval.evaluate(argv[0]);
		if(eval.status() != HepTool::Evaluator::OK) {
			printError("if: Invalid expression '%s'\n",
							argv[0].c_str());
			val = 0.0;
		}
		bool skipping = (fabs(val) < 1.0e-12);
		bool haveExecuted = !skipping;
		int level = 1;
		for(;;) {
			G4String *line = BLCommand::getNextCommand();
			if(!line) {
				printError("Unexpected EOF in if statement.");
				break;
			}
			if((*line) == "endif") {
				if(--level == 0) break;
				skipping = true;
				continue;
			} else if((*line) == "else") {
				if(level == 1) {
					if(!haveExecuted) {
						skipping = false;
						haveExecuted = true;
					} else {
						skipping = true;
					}
				}
				continue;
			} else if(line->find("elseif") == 0 &&
					line->find_first_of(" \t") == 6) {
				if(haveExecuted) {
					skipping = true;
					continue;
				}
				// overwrites argv and namedArgs (nobody cares)
				argv.clear();
				namedArgs.clear();
				if(parseArgs(line->substr(7),argv,namedArgs) < 0
				   || argv.size() != 1 || namedArgs.size() != 0) {
					printError("elseif: invalid elseif.");
					skipping = true;
					continue;
				}
				G4String t=Param.expand(argv[0]);
				double val = eval.evaluate(t.c_str());
				if(eval.status() != HepTool::Evaluator::OK) {
					printError("elseif: Invalid expression '%s'\n",
						line->substr(7).c_str());
					val = 0.0;
				}
				if(fabs(val) > 1.0e-12) {
					skipping = false;
					haveExecuted = true;
				}
				continue;
			} else if(line->find("if") == 0 &&
					line->find_first_of(" \t") == 2) {
				if(skipping) {
					// overwrites argv and namedArgs (nobody cares)
					argv.clear();
					namedArgs.clear();
					if(parseArgs(line->substr(3),argv,namedArgs) < 0) {
						printError("if: invalid if.");
						continue;
					}
					if(argv.size() == 1) ++level;
				}
			}
			if(!skipping)
				BLCommand::doCommand(*line);
		}
	}
	
	return retval;
}
