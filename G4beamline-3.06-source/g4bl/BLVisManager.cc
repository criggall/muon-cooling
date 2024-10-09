//	BLVisManager.cc
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

#ifdef G4BL_VISUAL

#include <stdlib.h>
#include <stdio.h>

#ifdef G4BL_GUI
#include "BLQt.h"
#endif

#include "G4UIterminal.hh"
#include "G4VUserVisAction.hh"

#include "G4Text.hh"

#include "BLVisManager.hh"
#include "BLManager.hh"
#include "BLCommand.hh"
#include "BLAssert.hh"
#include "BLMarkers.hh"

extern G4String g4bl_dir; // in g4bl.cc

BLVisManager *BLVisManager::manager = 0;

class DrawMarkers : public G4VUserVisAction {
	virtual void Draw() {
		BLMarkers::displayAll();
	}
};

BLVisManager::BLVisManager(G4String _viewer) :
	G4VisExecutive(), BLManager::RunAction(), initCommands(), beginRunCommands(),
	endRunCommands()
{
	viewer = _viewer; 
	UI = G4UImanager::GetUIpointer();
	manager = this;

	// read VISUAL_DEF_FILENAME, and find the section for our viewer
	readSection(viewer);

#ifndef WIN32
	// Some OpenGL implementations need this
	setenv("LIBGL_ALLOW_SOFTWARE","1",0);
#endif // !WIN32
}

void BLVisManager::init()
{
	Initialize();	// calls RegisterGraphicsSystems()

	printf("\nSelected visualization viewer: %s\n",viewer.c_str());

	doInitCommands();

	// register us for UserRunAction
	BLManager::getObject()->registerRunAction(this,false);
}

void BLVisManager::doInitCommands()
{
	for(unsigned int i=0; i<initCommands.size(); ++i) {
	    printf("BLVisManager init: %s\n",initCommands[i].c_str());
	    UI->ApplyCommand(initCommands[i]);
	}
	fputs("\n",stdout);
}

void BLVisManager::generateImages(int evPerImage, int nImages)
{
	// arrange to draw all BLMarkers
	SetUserAction(new DrawMarkers,G4VisExtent());
	UI->ApplyCommand("/vis/scene/add/userAction SetUserAction");

#ifdef G4BL_GUI
	// For the Qt viewer, just run the Qt UI and let the user start runs
	if(viewer == "Qt") {
		BLQt::setEvPerImage(evPerImage);
		BLQt::SessionStart();
		return;
	}
#endif

	if(nImages > 1) {
		G4UImanager* UI = G4UImanager::GetUIpointer();
		char cmd[64];
		sprintf(cmd,"/run/beamOn %d",evPerImage);
		for(int i=0; i<nImages; ++i)
			UI->ApplyCommand(cmd);
	} else if(BLManager::getObject()->getBeamVector()->size() > 0) {
		printf("To display a run with 10 events in the viewer, type '/run/beamOn 10<cr>'\n");
		G4UIterminal *terminal = new G4UIterminal(0,false);
		terminal->SessionStart();
		delete terminal;
	} else { // no beam, just show the world
		G4UImanager* UI = G4UImanager::GetUIpointer();
		UI->ApplyCommand("/vis/viewer/update");
	}
}

void BLVisManager::readSection(G4String section)
{
	FILE *in = fopen(VISUAL_DEF_FILENAME,"r");
	if(!in) {
		G4String name = g4bl_dir;
		name += "/share/g4beamline/" VISUAL_DEF_FILENAME;
		in = fopen(name.c_str(),"r");
		if(!in) {
			G4Exception("BLVisManager","Cannot read " 
					VISUAL_DEF_FILENAME,FatalException, "");
		}
	}
	G4String seek = "["; {seek += section; seek += "]";}
	char line[256];
	for(;;) {
		if(!fgets(line,sizeof(line),in)) {
			G4Exception("BLVisManager","Viewer not defined" 
					,FatalException, viewer.c_str());
		}
		if(line[0] == '#' || line[0] == '\n') continue;
		if(strncmp(line,seek.c_str(),seek.size()) == 0) break;
	}

	// Now read our commands into the command vectors
	while(fgets(line,sizeof(line),in)) {
		if(line[0] == '#' || line[0] == '\n') continue;
		if(line[0] == '[') break;
		char *p = strchr(line,'\n'); // remove newline
		if(p) *p = '\0';
		if(strncmp(line,"init:",5) == 0) {
			initCommands.push_back(G4String(line+5));
		} else if(strncmp(line,"include:",8) == 0) {
			readSection(line+8);
		} else if(strncmp(line,"help:",5) == 0) {
			puts(line+5);
#ifdef USE_SHARED_OBJECTS
		} else if(strncmp(line,"load:",5) == 0) {
			if(BLLoad::load(line+5))
				printf("BLVisManager: loaded '%s'\n",line+5);
			else
				BLCommand::printError("BLVisManager: Cannot "
					"load '%s'",line+5);
#endif // USE_SHARED_OBJECTS
		} else if(strncmp(line,"beginRun:",9) == 0) {
			beginRunCommands.push_back(G4String(line+9));
		} else if(strncmp(line,"endRun:",7) == 0) {
			endRunCommands.push_back(G4String(line+7));
		} else {
			printf("BLVisManager::init invalid command for viewer"
				"'%s':\n\t%s\n",viewer.c_str(),line);
		}
	}

	fclose(in);
}

void BLVisManager::BeginOfRunAction(const G4Run *run)
{
	for(unsigned int i=0; i<beginRunCommands.size(); ++i)
		UI->ApplyCommand(beginRunCommands[i]);
}

void BLVisManager::EndOfRunAction(const G4Run *run)
{
	for(unsigned int i=0; i<endRunCommands.size(); ++i)
		UI->ApplyCommand(endRunCommands[i]);
}

#else // G4BL_VISUAL
int dummyBLVisManager = 0;
#endif // G4BL_VISUAL
