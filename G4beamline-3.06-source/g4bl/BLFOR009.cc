//	BLFOR009.cc
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

#include <map>
#include <vector>

#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

#include "BLFOR009.hh"
#include "BLWriteAsciiFile.hh"

struct Track {
	G4ThreeVector pos;
	G4double time;
	G4ThreeVector momentum;
	G4ThreeVector Bfield;
	G4ThreeVector Efield;
	int PDGid;
	int eventId;
	int trackId;
	int parentId;
	G4double weight;
	bool isEqual(Track &t) {
		return eventId == t.eventId &&  trackId == t.trackId &&
		   fabs(time-t.time) < 1E-6 && (pos-t.pos).mag2() < 1E-6 &&
		   (momentum-t.momentum).mag2() < 1E-12;
	}
};

struct Region {
	std::vector<Track> track;
	double Z;
	Region(double v) : track() { Z=v; }
	void addTrack(Track &t) { 
		// ensure a single center track (for MPI)
		if(t.eventId == 0 && track.size() > 0) {
			if(!t.isEqual(track[0])) {
			    G4Exception("BLFOR009 Output",
			    	"Different Center Tracks in for009.dat",
				JustWarning,"");
			}
			return;
		}
		track.push_back(t); 
	}
};


BLFOR009::BLFOR009(G4String _filename, G4String title, G4String _mode)
{
	filename = _filename;
	status = BLTF_ERROR;
	file = 0;
	mode = _mode;
	map = new RegionMap();

	if(mode != "w") {
		fprintf(stderr,"BLFOR009: only write mode is implemented\n");
		return;
	}

	file = BLWriteAsciiFile::fopen(filename);
	if(!file) {
		fprintf(stderr,"BLFOR009: Cannot open file '%s' for writing\n",
			filename.c_str());
		return;
	}
	status = BLTF_OK;
	fprintf(file,"#%s\n",title.c_str());
	fprintf(file,"#Units are ns, meters, GeV/c, Tesla, and V/m\n");
	fprintf(file,"#IEVT IPNUM IPTYP IPFLG JSRG T X Y Z Px Py Pz Bx By Bz Weight Ex Ey Ez SARC POLx POLy POLz\n");

}

BLFOR009::~BLFOR009()
{
	close();
	// no attempt to delete map (too complicated, and unnecessary)
}

void BLFOR009::close()
{
	if(status != BLTF_OK || !file) return;

	RegionMap::iterator i;
	int jsrg;
	for(i=map->begin(),jsrg=1; i!=map->end(); ++i,++jsrg) {
		int ignored = 0;
		int n = i->second->track.size();
		bool hasReference = false;
		for(int k=0; k<n; ++k) {
			Track t(i->second->track[k]);
			int iptyp=0;
			switch(t.PDGid) {
			case -11:	iptyp =  1;	break;
			case  11:	iptyp = -1;	break;
			case -13:	iptyp =  2;	break;
			case  13:	iptyp = -2;	break;
			case -211:	iptyp = -3;	break;
			case  211:	iptyp =  3;	break;
			case -321:	iptyp = -4;	break;
			case  321:	iptyp =  4;	break;
			case -2212:	iptyp = -5;	break;
			case  2212:	iptyp =  5;	break;
			default:
				++ignored;
				continue;
			}
			if(t.eventId == 0) hasReference = true;
			fprintf(file,"%d %d %d %d %d"
					" %.6g"
					" %.6g %.6g %.6g"
					" %.6g %.6g %.6g"
					" %.6g %.6g %.6g"
					" %.4f"
					" %.6g %.6g %.6g"
					" 0   0 0 0\n",
				t.eventId,t.trackId,iptyp,0,jsrg,
				t.time/second,
				t.pos[0]/meter,t.pos[1]/meter,t.pos[2]/meter,
				t.momentum[0]/GeV,t.momentum[1]/GeV,t.momentum[2]/GeV,
				t.Bfield[0]/tesla,t.Bfield[1]/tesla,t.Bfield[2]/tesla,
				t.weight,
				t.Efield[0]/(volt/meter),t.Efield[1]/(volt/meter),t.Efield[2]/(volt/meter)
			);
		}
		if(!hasReference) 
			G4Exception("FOR009.DAT output","No reference",
							JustWarning, "");
		printf("FOR009.DAT Region %d z=%.1f  %d tracks,  %d "
				"unsupported particles ignored\n",
						jsrg,i->second->Z,n,ignored);
	}
	
	BLWriteAsciiFile::fclose(file);
	file = 0;
}

BLTrackFileStatus BLFOR009::write(G4ThreeVector& pos, G4double time,
			G4ThreeVector& momentum, G4ThreeVector& Bfield,
			G4ThreeVector& Efield, int PDGid, int eventId,
			int trackId, int parentId, G4double weight)
{
	if(status != BLTF_OK)
		return status;

	double z = pos[2];
	Region *m = (*map)[z];
	if(!m) {
		m = new Region(z);
		(*map)[z] = m;
	}

	// eventId=-1 is the center particle, and eventId=0 is a real event.
	// in FOR009.DAT center must be event 0.
	if(eventId == 0) {
		// yes, one line per region is printed.
		printf("BLFOR009: beam event 0 omitted (0 must be reference)\n");
		return BLTF_OK;
	} else if(eventId == -1) {
		eventId = 0;
	}

	Track t;
	t.pos = pos;
	t.time = time;
	t.momentum = momentum;
	t.Bfield = Bfield;
	t.Efield = Efield;
	t.PDGid = PDGid;
	t.eventId = eventId;
	t.trackId = trackId;
	t.parentId = parentId;
	t.weight = weight;

	m->addTrack(t);

	return BLTF_OK;
}

BLTrackFileStatus BLFOR009::read(G4ThreeVector& pos, G4double& time,
			G4ThreeVector& momentum, G4ThreeVector& Bfield,
			G4ThreeVector& Efield, int& PDGid, int& eventId,
			int& trackId, int& parentId)
{
	return BLTF_ERROR;
}

BLTrackFileStatus BLFOR009::read(G4ThreeVector& pos, G4double& time,
			G4ThreeVector& momentum, G4ThreeVector& Bfield,
                        G4ThreeVector& Efield, int& PDGid, int& eventId,
			int& trackId, int& parentId, G4double &weight)
{
	return BLTF_ERROR;
}
