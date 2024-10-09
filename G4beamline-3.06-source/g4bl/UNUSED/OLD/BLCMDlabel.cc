//	BLCMDlabel.cc
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

#include <vector>

#include "G4Text.hh"
#include "G4Polymarker.hh"

#include "BLCommand.hh"
#include "BLCallback.hh"
#include "BLCoordinates.hh"
#include "BLVisManager.hh"

/**	class BLCMDlabel implements the label command to display labels
 *	or markers.
 *
 *	callback() draws the labels, which works for visualization with
 *	no beam. But that call-back happens too early when beam is run,
 *	so EndOfRunAction() calls it.
 **/
class BLCMDlabel : public BLCommand, public BLManager::RunAction, BLCallback {
	G4String text;
	G4String color;
	G4double size;
	G4String coordinates;
	G4String file;
	static std::vector<G4Polymarker*> markers;
	static std::vector<G4Text*> texts;
	static bool registered;
public:
	BLCMDlabel();

	G4String commandName() { return "label"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs() {
		argString(text,"text","The text of the label; empty => "
							"circular marker.");
		argString(color,"color","The color (1,1,1).");
		argDouble(size,"size","The size in pixels.");
		argString(coordinates,"coordinates","Coordinates (global).");
		argString(file,"file","File to read for points.");
	}

	// BeginOfRunAction() from BLManager::RunAction.
	virtual void BeginOfRunAction(const G4Run *run) { }

	// EndOfRunAction() from BLManager::RunAction.
	virtual void EndOfRunAction(const G4Run *run) {
		callback(4);
	}

	void callback(int type) {
		BLVisManager *vm = BLVisManager::getObject();
		if(!vm) return;
		for(unsigned i=0; i<markers.size(); ++i)
			vm->Draw(*markers[i]);
		for(unsigned i=0; i<texts.size(); ++i)
			vm->Draw(*texts[i]);
	}
};
std::vector<G4Polymarker*> BLCMDlabel::markers;
std::vector<G4Text*> BLCMDlabel::texts;
bool BLCMDlabel::registered = false;
BLCMDlabel defineLabel;

BLCMDlabel::BLCMDlabel()
{
	registerCommand(BLCMDTYPE_OTHER);
	setSynopsis("Display labels or markers visually in 3-D space.");
	setDescription("Each positional argument is X,Y,Z of the position; "
	"as many positions as desired can be used. The same text is displayed "
	"at each position. Empty text means use a circular marker.\n\n"
	"The marker or label always faces the camera.\n\n"
	"The default height of text is 12 pixels; the default diameter of "
	"markers is 5 pixels; some viewers have rather narrow limits on "
	"sizes.\n\n"
	"This command is not placed into the geometry.\n\n"
	"Coordinates can be either centerline or global (the default).\n\n"
	"The file argument is most useful when text is empty; as only the "
	"first three columns are read, a BLTrackFile can be used.");
	text = "";
	color = "1,1,1";
	size = -1;
	coordinates = "global";
	file = "";
}

int BLCMDlabel::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	// restore defaults
	text = "";
	color = "1,1,1";
	size = -1;
	coordinates = "global";
	file = "";

	int retval = handleNamedArgs(namedArgs);

	if(argv.size() == 0 && file.size() == 0) return retval;

	G4Polymarker *m = 0;
	G4Text *t = 0;
	if(text == "") {
		if(size <= 0) size = 5;
		m = new G4Polymarker;
		m->SetMarkerType(G4Polymarker::dots);
		m->SetScreenSize(size);
		m->SetVisAttributes(getVisAttrib(color));
	} else {
		if(size <= 0) size = 12;
		t = new G4Text(text);
		t->SetLayout(G4Text::left);
		t->SetScreenSize(size);
		t->SetVisAttributes(getVisAttrib(color));
	}

	BLCoordinateType coord = BLCoordinates::getCoordinateType(coordinates);

	for(unsigned i=0; i<argv.size(); ++i) {
		std::vector<G4double> p = getList(argv[i],",");
		if(p.size() != 3) {
			G4Exception("label","Invalid label/marker position",
							JustWarning, argv[i]);
			continue;
		}
		G4Point3D pos(p[0],p[1],p[2]);
		if(coord == BLCOORD_CENTERLINE) {
			G4ThreeVector tmp(pos[0],pos[1],pos[2]);
			G4ThreeVector pos3v;
			BLCoordinates::getGlobalAnywhere(tmp,pos3v);
			pos[0] = pos3v[0];
			pos[1] = pos3v[1];
			pos[2] = pos3v[2];
		}
		if(m) {
			m->push_back(pos);
		} else if(t) {
			t->SetPosition(pos);
			texts.push_back(new G4Text(*t));
		}
	}

	if(file != "") {
		FILE *in = fopen(file.c_str(),"r");
		if(!in) {
			G4Exception("label","Cannot read file",
							JustWarning, file);
		} else {
			char buf[4096];
			while(fgets(buf,sizeof(buf),in) != 0) {
				if(buf[0] == '#') continue;
				G4ThreeVector pos;
				if(sscanf(buf,"%lf%lf%lf",&pos[0],&pos[1],
					&pos[2]) != 3) continue;
				if(m) {
					m->push_back(pos);
				} else if(t) {
					t->SetPosition(pos);
					texts.push_back(new G4Text(*t));
				}
			}
			fclose(in);
		}
	}

	if(m) markers.push_back(m);
	if(t) delete t;

	if(retval == 0 && !registered) {
		//BLManager::getObject()->registerRunAction(this,false);
		BLManager::getObject()->registerCallback(this,4);
		registered = true;
	}

	print("");

	return retval;
}
