//	BLKillTrack.hh
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

#ifndef BLKILLTRACK_HH
#define BLKILLTRACK_HH

#include "BLParam.hh"

/**	class BLKillTrack - class to kill tracks.
 *
 *	Any BLElement that wants to kill tracks that enter a given physical 
 *	volume (e.g. a box command with kill=1) should register an
 *	instance of this class with the BLManager, associated with the
 *	physical volume.
 **/
class BLKillTrack : public BLManager::SteppingAction {
	bool verbose;
	G4String name;
public:
	BLKillTrack(G4String &_name) : BLManager::SteppingAction(), name(_name)
	{ 
		verbose = Param.getInt("steppingVerbose") != 0;
	}
	void UserSteppingAction(const G4Step *step) {
		G4Track *track = step->GetTrack();
		track->SetTrackStatus(fStopAndKill);
		if(verbose) printf("Track killed by '%s' with kill=1\n",
			name.c_str());
	}
};

#endif // BLKILLTRACK_HH
