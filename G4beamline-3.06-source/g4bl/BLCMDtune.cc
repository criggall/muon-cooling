//	BLCMDtune.cc
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

#include <set>
#include "CLHEP/Evaluator/Evaluator.h"

#include "G4Track.hh"

#include "BLParam.hh"
#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLEvaluator.hh"
#include "BLTune.hh"
#include "BLCoordinates.hh"

enum State { INIT, TUNING, DONE, FIXED };

/**	class BLCMDtune - tunes a variable used as argument(s) to other elements
 **/
class BLCMDtune : public BLCommand, public BLManager::ZSteppingAction,
					public BLManager::RunAction {
	G4String name;
	G4double z0;
	G4double z1;
	G4double initial;
	G4double initialStep;
	G4String start;
	G4String expr;
	G4double tolerance;
	G4int maxIter;
	G4Track saveTrack;
	G4double minStep;
	BLEvaluator eval;
	State state;
	G4int count;
	G4double thisX, prevX, thisY, prevY;
public:
	/// Default constructor.
	BLCMDtune();

	/// Destructor.
	virtual ~BLCMDtune();

	/// Copy constructor.
	BLCMDtune(const BLCMDtune& r);

	/// commandName() returns "tune".
	G4String commandName() { return "tune"; }

	/// command() implements the tune command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// from BLManager::ZSteppingAction.
	void UserZSteppingAction(const G4Track *track);

	/// from BLManager::RunAction
	void BeginOfRunAction(const G4Run* run) { }
	void EndOfRunAction(const G4Run* run);
};
BLCMDtune defaultTune;


BLCMDtune::BLCMDtune() : BLCommand(), BLManager::ZSteppingAction(), 
			BLManager::RunAction(), saveTrack(), eval()
{
	registerCommand(BLCMDTYPE_CONTROL);
	setSynopsis("Tune a variable used as argument to other elements.");
	setDescription("This command samples the Tune particle track at z0, "
		"samples it again at z1, and varies its tune variable in "
		"order to bring the expression to zero. The tune variable "
		"is the first positional argument, and "
		"can be used in the argument expression(s) for tunable "
		"arguments located after z0. Due to the simple "
		"solver used, there should be an approximately linear "
		"dependence between the tune variable and the expression. "
		"This is suitable for tuning the By field of a genericbend "
		"or the maxGradient of a pillbox. "
		"At each solver step the saved Tune particle is re-started from "
		"z0, and when it reaches z1 the next step in the solver is "
		"taken.\n\n"
		"Note that multiple tune commands can be used together, as long "
		"as their z0-z1 regions are properly nested by at least 0.5 mm "
		"(i.e. each pair of regions must either not overlap at all, or "
		"one must be wholly contained in the other). "
		"This command can be complicated to use; see the User's Guide "
		"for more description and examples.\n\n"
		"This command places itself into the geometry.");
	
	name = "";
	z0 = -DBL_MAX;
	z1 = -DBL_MAX;
	initial = 0.0;
	initialStep = 0.1;
	start = "1";
	expr = "";
	tolerance = 0.0001;
	maxIter = 10;
	minStep = 0.01*mm;
	state = INIT;
	count = 0;
	thisX = prevX = thisY = prevY = 0.0;
}

BLCMDtune::~BLCMDtune()
{
}

BLCMDtune::BLCMDtune(const BLCMDtune& r) : BLCommand(r), 
			BLManager::ZSteppingAction(r),
			BLManager::RunAction(r), saveTrack(), eval()
{
	name = r.name;
	z0 = r.z0;
	z1 = r.z1;
	initial = r.initial;
	initialStep = r.initialStep;
	start = r.start;
	expr = r.expr;
	tolerance = r.tolerance;
	maxIter = r.maxIter;
	minStep = r.minStep;
	state = INIT;
	count = 0;
	thisX = prevX = thisY = prevY = 0.0;
}

int BLCMDtune::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("tune: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultTune.handleNamedArgs(namedArgs);
	}

	if(BLTune::isSet(argv[0])) {
		printError("tune: name '%s' is already in use!",name.c_str());
		return -1;
	}

	BLCMDtune *t = new BLCMDtune(defaultTune);
	t->minStep = Param.getDouble("minStep");
	int retval = t->handleNamedArgs(namedArgs);

	t->name = argv[0];
	BLTune::set(t->name,t->initial);
	t->eval.setVariable(t->name,t->initial);
	t->thisX = t->initial;
	if(fabs(t->initialStep) < 1.0e-12) {
		t->state = FIXED;
	} else {
		t->state = INIT;
		BLManager::getObject()->registerZStep(t->z0,t,1);
		BLManager::getObject()->registerZStep(t->z1,t,1);
		BLManager::getObject()->registerRunAction(t,false);
	}

	t->print(t->name);

	return retval;
}

void BLCMDtune::defineNamedArgs()
{
	argDouble(z0,"z0","The starting z position in CL coordinates.");
	argDouble(z1,"z1","The ending z position in CL coordinates.");
	argDouble(initial,"initial","Initial value of the variable 'name'");
	argDouble(initialStep,"initialStep","Initial step (0 to disable tuning)");
	argDouble(initialStep,"step","Synonym for initialStep");
	argString(start,"start","An expression that must be nonzero to start tuning (default=1)");
	argString(expr,"expr","The expression to tune to zero");
	argDouble(tolerance,"tolerance","The tolerance for expr to be zero");
	argInt(maxIter,"maxIter","The maximum number of iterations (10).");
}

void BLCMDtune::UserZSteppingAction(const G4Track *track)
{
	if(state == FIXED) return;

	BLCoordinates *coord =  
			(BLCoordinates *)track->GetUserInformation();
	G4ThreeVector clpos;
	coord->getCoords(BLCOORD_CENTERLINE,clpos);

	if(fabs(clpos[2]-z0) < minStep && state != TUNING) {
		// we just entered the tuning region -- save the track
		eval.setTrackVariables(track,BLCOORD_CENTERLINE,"0");
		G4double t = eval.evaluate(start);
		if(eval.status() != HepTool::Evaluator::OK) {
			char tmp[128];
			sprintf(tmp,"%s (variables are x0,y0,...,x1,y2,...)",
								start.c_str());
			G4Exception("tune","Invalid Start Expr",FatalException,
				tmp);
		}
		if(t == 0.0) goto quit;
		saveTrack.CopyTrackInfo(*track);
		saveTrack.SetUserInformation(0);
		state = TUNING;
		count = 0;
		printf("tune '%s' begun  %s=%.4f\n",name.c_str(),name.c_str(),
								thisX);
	} else if(fabs(clpos[2]-z1) < minStep && state == TUNING) {
		// we just reached the end of the tuning region
		if(++count > maxIter)
			G4Exception("tune","Iteration Limit",FatalException,
								name.c_str());
		eval.setTrackVariables(track,BLCOORD_CENTERLINE,"1");
		thisY = eval.evaluate(expr);
		if(eval.status() != HepTool::Evaluator::OK) {
			char tmp[128];
			sprintf(tmp,"%s (variables are x0,y0,...,x1,y2,...)",
								start.c_str());
			G4Exception("tune","Invalid Tune Expr",FatalException,
				tmp);
		}
		if(fabs(thisY) < tolerance) {
			state = DONE;
			printf("tune '%s' complete in %d steps  expr=%.6f  %s=%.4f\n",
				name.c_str(),count,thisY,name.c_str(),thisX);
			goto quit;
		}
		G4double newX=thisX+initialStep;
		if(count > 1) {
			G4double dydx = (prevY-thisY)/(prevX-thisX);
			newX = thisX - thisY/dydx;
		}
		printf("tune '%s' step %d  %s=%.4f  expr=%.6f  new %s=%.4f\n",
				name.c_str(),count,name.c_str(),thisX,thisY,
				name.c_str(),newX);
		prevX = thisX;
		prevY = thisY;
		thisX = newX;
		BLTune::set(name,thisX);
		eval.setVariable(name,thisX);
		// restore saveTrack and kill the current track
		BLManager::getObject()->getSteppingManager()->GetfSecondary()->
				push_back(new G4Track(saveTrack));
		((G4Track *)track)->SetTrackStatus(fStopAndKill);
	}
quit:	;
}

void BLCMDtune::EndOfRunAction(const G4Run* run)
{
	if(state != DONE && state != FIXED) {
		G4Exception("tune","Failed to Converge",FatalException,
							name.c_str());
	}
}
