//	g4bl.cc - Main program for g4beamline
// DEBUG_INPUT is for testing, used only if no command-line parameters 
// (e.g. in debugger)
//#define DEBUG_INPUT "C:/Users/tjrob/test.g4bl"
//#define G4BL_DIR "C:/Users/tjrob/G4beamline-3.03-debug"
//#define VIEWER_ARG "viewer=best"

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

#ifdef WIN32
#include <windows.h>
#undef min
#undef max
#include <tchar.h> 
#include <strsafe.h>
#include <process.h>
/***
#include "crtdbg.h"
***/
#else
#include <unistd.h>
#endif

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "G4NistManager.hh"
#include "G4GeometryManager.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4VisExecutive.hh"

#include "Util.hh"
#include "BLMPI.hh"
#include "BLNTuple.hh"
#include "BLParam.hh"
#include "BLCommand.hh"
#include "BLGroup.hh"
#include "BLManager.hh"
#include "BLTime.hh"
#include "BLSignal.hh"
#include "BLWriteAsciiFile.hh"
#include "Geant4Data.hh"

// nested macros are needed to stringize the value of a preprocessor symbol
#define QUOTE(A) QUOTE_(A)
#define QUOTE_(A) #A

// Global Variables
G4String g4bl_dir = "";
G4String input_file = "";

BLSetParam maxStep("maxStep","100.0","Maximum physics step size (mm)");
BLSetParam worldMaterial("worldMaterial","Vacuum","Material of the World volume");

void startupPrint()
{
    printf("*************************************************************\n");
    printf(" g4beamline version: %s                        (%s)\n",
					    QUOTE(G4BLVERSION),__DATE__);
    printf("                      Copyright : Tom Roberts, Muons, Inc.\n");
    printf("                        License : Gnu Public License\n");
    printf("                            WWW : http://g4beamline.muonsinc.com");
    fflush(stdout); // for checking the program started up
}

void g4bl_exit(int value)
{
	static bool first=true;
	if(!first) {
	    fprintf(stderr,"g4beamline: g4bl_exit re-entered, exiting fast\n");
	    exit(90);
	}
	first = false;

	if(BLMPI::isMPI())
		BLMPI::closeupAndExit(value);

	// Instead of deleting things, just open the geometry to avoid 
	// warnings when exit() is called.
	G4GeometryManager::GetInstance()->OpenGeometry();

	BLWriteAsciiFile::closeAll();

	if(value == 0)
		printf("g4beamline: simulation complete\n");
	else
		printf("g4beamline: simulation aborted\n");
	fflush(stdout);

	BLSignal::writeStackTrace(2); // does nothing if no stack trace
	exit(value);
}

#ifdef STUB
extern "C" void abort()
{
	BLSignal::generateStackTrace();
	G4Exception("abort()","abort() called",FatalException,"");
	_exit(99);
}
#endif // STUB

void printEnvironment()
{
	extern char **environ;
	int i = 1;
	char *s = *environ;

	printf("\nProcess Environment:\n");
	for (; s; i++) {
		printf("%s\n", s);
		s = *(environ+i);
	}
	printf("\n\n");
}

int main(int _argc, char *_argv[])
{
	// get program arguments
	std::vector<G4String> argV;
	for(int i=0; i<_argc; ++i) {
		argV.push_back(_argv[i]);
	}
	assert(argV.size() >= 1);
	if(argV.size() > 1) input_file = argV[1];

#ifdef DEBUG_INPUT
	// Set the input file to DEBUG_INPUT, but only if no args are given
	if(argV.size() < 2) {
		argV.push_back(DEBUG_INPUT);
		printf("DEBUG_INPUT: InputFile set to '" DEBUG_INPUT "'\n");
#ifdef VIEWER_ARG
		argV.push_back(VIEWER_ARG);
		printf("DEBUG_INPUT: %s\n",VIEWER_ARG);
#endif
#ifdef G4BL_DIR
		g4bl_dir = G4BL_DIR;
		printf("DEBUG_INPUT: g4bl_dir=%s\n",G4BL_DIR);
#endif
	}
#endif

	// get g4bl_dir, the install path. This executable is in bin under it.
	if(g4bl_dir == "") g4bl_dir = getG4blDir();
	printf("G4BL_DIR=%s\n",g4bl_dir.c_str());

	// find the Geant4 data sets, put them into my environment
	if(!Geant4Data::setup()) Geant4Data::runG4bldata();
	Geant4Data::printEnv();

	// initialize and set up MPI, if appropriate.
	// Sets stderr and stdout to /dev/null except for rank 0 or non-MPI.
	BLMPI::init(argV);

	// print my process ID (for g4blgui to properly terminate me)
#ifdef WIN32
	printf("G4beamline Process ID %d\n\n",_getpid());
#else
	printf("G4beamline Process ID %d\n\n",getpid());
#endif

    // Identify ourself (styled after the Geant4 identification banner)
    startupPrint();
    if(argV.size() < 2) {
      printf("\nUSAGE: g4beamline inputfile [viewer=TYPE] [param=value...]\n");
      printf("       viewer=best   Display detector and events visually\n");
      printf("                     (can be any valid viewer name).\n");
      printf("       viewer=none   (default) Track beam.\n");
      printf("       Any other parameters can be given.\n");
      ::exit(1);
    }

	// set PRNG seed. Note the seed is probably re-set
	// at the start of each event, according to the randomseed command
	// (event number if no such command).
	CLHEP::HepRandom::setTheSeed(0x7FFFFFFEL);

	// Ensure the BLManager has been created, and do its delayed constuction
	BLManager *mgr = BLManager::getObject();
	mgr->delayedConstruction();

	BLMPI::scan();

	// Handle arguments 2 thru argc-1 as arguments to the param command
	// They should be of the form name=value or name="a value with spaces"
	// Arguments that consist of a valid command will be kept until all
	// parameters are defined, and then will be executed in order.
	// Remember the shell will eat " and '
	G4String cmd("param");
	std::vector<G4String> argCommands;
	for(unsigned i=2; i<argV.size(); ++i) {
		cmd += " ";
		G4String s(argV[i]);
		G4String::size_type p = s.find("=");
		if(p == s.npos || s.find(" ") < p) {
			p = s.find(" ");
			G4String c = s.substr(0,p);
			if(BLCommand::isValidCommand(c))
			    argCommands.push_back(s);
			else
			    printf("*** g4beamline: invalid arg '%s' ignored\n",
					s.c_str());
			continue;
		}
		cmd += s.substr(0,p+1);
		if(s.find("'") == s.npos) {
			cmd += "'";
			cmd += s.substr(p+1);
			cmd += "'";
		} else if(s.find('"') == s.npos) {
			cmd += '"';
			cmd += s.substr(p+1);
			cmd += '"';
		} else {
			cmd += s;
		}
	}
	if(cmd.find("=") != cmd.npos)
		BLCommand::doCommand(cmd);
	for(unsigned i=0; i<argCommands.size(); ++i)
		BLCommand::doCommand(argCommands[i]);

	// arrange to perform a geometry test
	cmd = "geometry";
	BLCommand::doCommand(cmd);

	BLMPI::scan();

	// don't write output files in visual mode
	G4String viewer = Param.getString("viewer");
	if(viewer != "none")
		BLNTuple::disableOutputFiles();

	// read input file and end the World volume
	if(BLCommand::readFile(input_file) != 0) {
		char tmp[64];
		sprintf(tmp,"There were %d errors in the input",
				BLCommand::getErrorCount());
		G4Exception("main","Input Errors",FatalException,tmp);
	}
	BLGroup::getWorld()->end();

	// quit if world is empty
	if(BLGroup::getWorld()->getNChildren() == 0) {
		G4Exception("main","Empty World",FatalException,
			"No elements have been placed into the World");
	}

	if(viewer != "none" && BLMPI::isMPI())
		G4Exception("main","Viewer and MPI",FatalException,
				"Visualization cannot be used in MPI mode.");

	BLMPI::scan();

	// verify that viewer did not change from command-line (because we
	// disabled NTuple output files if viewer != none).
	if(viewer != Param.getString("viewer"))
		G4Exception("main","Viewer changed",FatalException,
		  "The viewer parameter may only be set on the command-line.");

	Param.printParam();

	// use the worldMaterial.
	BLGroup::getWorld()->setMaterial(Param.getString("worldMaterial"));

	// initialize BLManager (construct the geometry, etc.)
	mgr->initialize();

	BLMPI::scan();

	// handle pre-reference callbacks
	mgr->handleCallbacks(0);

	// Track tune and reference particle(s), if present
	if(mgr->nReference() > 0)
		mgr->trackTuneAndReferenceParticles();

	BLMPI::scan();

	// Perform any source runs.
	mgr->handleSourceRun();

	// handle post-reference callbacks
	mgr->handleCallbacks(1);

	// handle MPI (never returns if rank 0 in MPI mode).
	BLMPI::main();

	// handle replace main loop callbacks
	mgr->handleCallbacks(3);

	if(viewer != "none") {
		if(!BLElement::allOK())
			G4Exception("main","Element Error",JustWarning,
			    "Not all BLElement-s are OK -- continuing anyway");
		mgr->displayVisual();
	} else {
		if(!BLElement::allOK()) {
			G4Exception("main","Element Error",FatalException,
				"Not all BLElement-s are OK");
		}
		mgr->trackBeam();
		if(BLMPI::isMPI())
			BLMPI::closeupAndExit(0);
		BLNTuple::summary();
		BLNTuple::closeAll();
	}

	// handle post-tracking callbacks
	mgr->handleCallbacks(2);

	delete mgr;

	g4bl_exit(0);
}

