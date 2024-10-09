//	BLMarkers.hh
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

#include "G4Polymarker.hh"
#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Text.hh"
#include "G4Polyline.hh"

#include "BLCommand.hh"

#define MARKER_SIZE 5		/* pixels */
#define MARKER_COLOR "1,0,1"	/* magenta */

/**	class BLMarkers draws a set of markers.
 *	It does nothing if not in visualization mode (except collect the 
 *	points of the markers).
 *
 *	This class is mostly used for debugging -- sometimes a visual set
 *	of marked points is very useful. It is also used by the "label"
 *	command.
 *
 *	NOTE: this class MUST be created via new (so it still exists
 *	when its EndOfRunAction() is called). The constructor registers it
 *	to be drawn appropriately. You can create an instance and call 
 *	clear() and addMarker() during any stage of G4beamline; drawing
 *	occurs whenever the viewer is refreshed (at startup, and at end of
 *	run).
 *
 *	Multiple instances of BLMarkers can be used; they should probably
 *	be given different colors.
 **/
class BLMarkers {
	G4Polymarker markers;
	G4String text;
	bool line;
	static std::vector<BLMarkers*> allMarkers;
	static std::vector<G4Text*> allTexts;
	static std::vector<G4Polyline*> allLines;
public:
	/// Constructor. color is the usual color string "1,1,1" for white,
	/// "1,0,0" for red, etc.; size is in pixels. The defaults are
	/// magenta and a modest size.
	BLMarkers(G4String color=MARKER_COLOR, G4int size=MARKER_SIZE) :
					markers(), text(), line(false) { 
		markers.SetVisAttributes(BLCommand::getVisAttrib(color));
		markers.SetMarkerType(G4Polymarker::circles);
		markers.SetScreenDiameter(size);
		markers.SetFillStyle(G4VMarker::filled);
		allMarkers.push_back(this);
	}

	/// Destructor.
	virtual ~BLMarkers() { 
		for(std::vector<BLMarkers*>::iterator it=allMarkers.begin();
					it != allMarkers.end(); ++it) {
			if(*it == this) {
				allMarkers.erase(it);
				break;
			}
		}
	}

	/// clear() will remove all markers.
	void clear() { markers.clear(); }

	/// size() returns the number of markers present.
	int size() const { return markers.size(); }

	/// setText() will cause it to display text at each point.
	void setText(const G4String _text) { text = _text; }

	/// setLine() will cause it to draw a polyline through all points, in 
	/// the order they were added.
	void setLine() { line = true; }

	/// addMarker() adds a marker, global coordinates.
	void addMarker(const G4ThreeVector point) {
		if(text == "") {
			markers.push_back(point);
		} else {
			G4Text *t = new G4Text(text);
			t->SetLayout(G4Text::left);
			t->SetScreenSize(markers.GetScreenDiameter());
			t->SetVisAttributes(markers.GetVisAttributes());
			t->SetPosition(point);
			allTexts.push_back(t);
		}
	}

	// displayAll() will display all instances of BLMarkers.
	static void displayAll() {
		G4VVisManager* vm = G4VVisManager::GetConcreteInstance();
		if(!vm) return;
		// handle lines, converting from BLMarker to G4Polyline
		for(std::vector<BLMarkers*>::iterator it=allMarkers.begin();
				it != allMarkers.end(); /* no incr */) {
			if(!(*it)->line) { ++it; continue;}
			G4Polyline *pl = new G4Polyline();
			pl->SetVisAttributes((*it)->markers.GetVisAttributes());
			for(unsigned j=0; j<(*it)->markers.size(); ++j)
				pl->push_back((*it)->markers[j]);
			allLines.push_back(pl);
			it = allMarkers.erase(it);
		}
		for(unsigned i=0; i<allMarkers.size(); ++i) {
			vm->Draw(allMarkers[i]->markers);
		}
		for(unsigned i=0; i<allTexts.size(); ++i) {
			vm->Draw(*allTexts[i]);
		}
		for(unsigned i=0; i<allLines.size(); ++i) {
			vm->Draw(*allLines[i]);
		}
	}
};
