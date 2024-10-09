//	BLVisManager.hh
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

#ifndef BLVISMANAGER_H
#define BLVISMANAGER_H

#include <vector>

#include "G4VisExecutive.hh"
#include "G4UImanager.hh"
#include "G4Run.hh"

#include "BLManager.hh"

#define VISUAL_DEF_FILENAME "viewer.def"

/**	class BLVisManager manages the visualization drivers.
 **/
class BLVisManager: public G4VisExecutive, public BLManager::RunAction {
	G4String viewer;
	G4UImanager* UI;
	std::vector<G4String> initCommands;
	std::vector<G4String> beginRunCommands;
	std::vector<G4String> endRunCommands;
	static BLVisManager *manager;
public:
	/// Default Constructor.
	BLVisManager() : G4VisExecutive(), BLManager::RunAction() {
		viewer = ""; UI = 0;
	}

	/// Constructor given the name of a viewer.
	/// Reads VISUAL_DEF_FILENAME and selects _viewer.
	BLVisManager(G4String _viewer);

	/// get pointer to the singleton object
	static BLVisManager *getObject() { return manager; }

	/// init() will Initialize() the graphics system, and execute the 
	/// "init:" commands of the selected viewer.
	void init();

	/// doInitCommands() will perform the init commands for the ciewer.
	void doInitCommands();

	/// generateImages() will generate the images, starting the QtUI
	/// for the Qt viewer.
	void generateImages(int evPerImage, int nImages);

	/// readSection() reads 1 section from VISUAL_DEF_FILENAME
	void readSection(G4String section);

	/// BeginOfRunAction() from BLManager::RunAction.
	/// Executes the "beginRun:" commands of the selected viewer.
	void BeginOfRunAction(const G4Run *run);

	/// EndOfRunAction() from BLManager::RunAction.
	/// Executes the "endRun:" commands of the selected viewer.
	void EndOfRunAction(const G4Run *run);
};

#endif // BLVISMANAGER_H

#endif //G4BL_VISUAL
