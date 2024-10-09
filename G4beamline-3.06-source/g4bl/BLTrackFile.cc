//	BLTrackFile.cc
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

#include <string.h>

#include "globals.hh" // for G4Exception()
#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

#include "BLTrackFile.hh"
#include "BLWriteAsciiFile.hh"

BLTrackFile::BLTrackFile()
{
	status = BLTF_DUMMY;
	mode = "";
	file = 0;
	version = 1;
	unit = mm;
}

BLTrackFile::BLTrackFile(G4String filename, G4String comment, G4String _mode,
								int _version)
{
	status = BLTF_ERROR;
	mode = _mode;
	version = _version;
	unit = mm;
	if(mode == "w")
		file = BLWriteAsciiFile::fopen(filename);
	else
		file = fopen(filename.c_str(),"r");
	if(!file || (mode != "r" && mode != "w")) {
		fprintf(stderr,"BLTrackFile: cannot open '%s' in mode '%s'",
					filename.c_str(),mode.c_str());
		return;
	}
	status = BLTF_OK;
	if(mode == "w") {
	    if(version == 1) {
		fprintf(file,"#BLTrackFile %s\n",comment.c_str());
		fprintf(file,"#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight\n");
		fprintf(file,"#mm mm mm MeV/c MeV/c MeV/c ns - - - - -\n");
	    } else if(version == 2) {
		fprintf(file,"#BLTrackFile2 %s\n",comment.c_str());
		fprintf(file,"#x y z Px Py Pz t PDGid EventID TrackID ParentID Weight Bx By Bz Ex Ey Ez ProperTime PathLength PolX PolY PolZ InitX InitY InitZ InitT InitKE\n");
		fprintf(file,"#mm mm mm MeV/c MeV/c MeV/c ns - - - - - T T T MV/m MV/m MV/m ns mm - - - mm mm mm ns MeV\n");
	    } else {
		G4Exception("BLTrackFile output","Invalid version.",
				FatalException, "");
	    }
	} else if(mode == "r") {
		char line[512];
		if(!fgets(line,sizeof(line),file)) {
			status = BLTF_ENDOFFILE;
			return;
		}
		G4String s(line);
		if(s.substr(0,12) != "#BLTrackFile") {
			G4Exception("BLTrackFile input","Possible wrong format",
				JustWarning, "Assumed OK");
			rewind(file);
		}
	}
}

BLTrackFile::~BLTrackFile()
{
	if(file && mode == "w") BLWriteAsciiFile::fclose(file);
	if(file && mode == "r") fclose(file);
	file = 0;
	status = BLTF_ERROR;
}

BLTrackFileStatus BLTrackFile::write(const G4ThreeVector& pos, G4double time,
			const G4ThreeVector& momentum, int PDGid, int eventId,
			int trackId, int parentId, G4double weight)
{
	if(status == BLTF_DUMMY) return BLTF_OK;
	if(status != BLTF_OK || mode != "w") return BLTF_ERROR;

	// left justified to avoid pesky initial spaces in column 1
	fprintf(file,"%.6g %.6g %.6g ",pos[0]/unit,pos[1]/unit,pos[2]/unit);
	fprintf(file,"%.6g %.6g %.6g ",momentum[0]/MeV,momentum[1]/MeV,
						momentum[2]/MeV);
	fprintf(file,"%.6g %d %d %d %d %.6g\n",time/ns,PDGid,eventId,
						trackId,parentId,weight);
	return status;
}

BLTrackFileStatus BLTrackFile::write2(const G4ThreeVector& pos, G4double time,
			const G4ThreeVector& momentum, int PDGid, int eventId,
			int trackId, int parentId, G4double weight,
			const G4ThreeVector& bfield,
			const G4ThreeVector& efield,
			G4double properTime, G4double pathLength,
			const G4ThreeVector& polarization,
			const G4ThreeVector& initialPos,
			G4double initialT, G4double initialKE)
{
	if(status == BLTF_DUMMY) return BLTF_OK;
	if(status != BLTF_OK || mode != "w") return BLTF_ERROR;

	// left justified to avoid pesky initial spaces in column 1
	fprintf(file,"%.6g %.6g %.6g ",pos[0]/unit,pos[1]/unit,pos[2]/unit);
	fprintf(file,"%.6g %.6g %.6g ",momentum[0]/MeV,momentum[1]/MeV,
						momentum[2]/MeV);
	fprintf(file,"%.6g %d %d %d %d %.6g\n",time/ns,PDGid,eventId,
						trackId,parentId,weight);
	fprintf(file,"%.6g %.6g %.6g ",bfield[0]/tesla,bfield[1]/tesla,
						bfield[2]/tesla);
	fprintf(file,"%.6g %.6g %.6g ",efield[0]/(megavolt/meter),
			efield[1]/(megavolt/meter),efield[2]/(megavolt/meter));
	fprintf(file,"%.6g ",properTime/ns);
	fprintf(file,"%.6g ",pathLength/unit);
	fprintf(file,"%.6g %.6g %.6g ",polarization[0],polarization[1],
						polarization[2]);
	fprintf(file,"%.6g %.6g %.6g ",initialPos[0]/unit,initialPos[1]/unit,
						initialPos[2]/unit);
	fprintf(file,"%.6g\n",initialT/ns);
	fprintf(file,"%.6g\n",initialKE/MeV);
	return status;
}

BLTrackFileStatus BLTrackFile::read(G4ThreeVector& pos, G4double& time,
			G4ThreeVector& momentum, int& PDGid, int& eventId,
			int& trackId, int& parentId, G4double &weight)
{
	if(status == BLTF_DUMMY) return BLTF_ENDOFFILE;
	if(status != BLTF_OK || mode != "r") return BLTF_ERROR;

	char line[1024];
	do {
		if(!fgets(line,sizeof(line),file)) {
			status = BLTF_ENDOFFILE;
			return status;
		}
		if(line[0] == '#' && strstr(line,"cm cm cm"))
			unit = cm;
	} while(line[0] == '#' || line[0] == '\n' || line[0] == '\r' || 
							line[0] == '\0');

	// Tom Roberts 2010-04-05
	// Read into doubles, and then convert. This permits event numbers 1.0e7
	// and 10000001 to be read as accurately as possible.
	double data[12];

	//Ajit Kurup 2008-06-02
	//Read in PDGid, eventId,trackId and parentId as integers
	//since the float->int conversion doesn't work properly for
	//large numbers

	data[11] = 1.0;	// in case weight is missing from file
	if(sscanf(line,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
				data+0,data+1,data+2,data+3,data+4,data+5,
				data+6,data+7,data+8,data+9,data+10,
				data+11) < 11) {
		printf("BLTrackFile invalid input line (file abandoned): '%s'\n",
				line);
		fflush(stdout);
		status = BLTF_ERROR;
		return status;
	}

	pos[0] = data[0]*unit;
	pos[1] = data[1]*unit;
	pos[2] = data[2]*unit;
	momentum[0] = data[3]*MeV;
	momentum[1] = data[4]*MeV;
	momentum[2] = data[5]*MeV;
	time = data[6]*ns;
	PDGid = (int)data[7];
	eventId = (int)data[8];
	trackId = (int)data[9];
	parentId = (int)data[10];
	weight = data[11];

	return status;
}

BLTrackFileStatus BLTrackFile::read2(G4ThreeVector& pos, G4double& time,
			G4ThreeVector& momentum, int& PDGid, int& eventId,
			int& trackId, int& parentId, G4double &weight,
			G4ThreeVector& bfield, G4ThreeVector& efield,
			G4double& properTime, G4double& pathLength,
			G4ThreeVector& polarization, G4ThreeVector& initialPos,
			G4double& initialT, G4double& initialKE)
{
	if(status == BLTF_DUMMY) return BLTF_ENDOFFILE;
	if(status != BLTF_OK || mode != "r") return BLTF_ERROR;

	char line[1024];
	do {
		if(!fgets(line,sizeof(line),file)) {
			status = BLTF_ENDOFFILE;
			return status;
		}
		if(line[0] == '#' && strstr(line,"cm cm cm"))
			unit = cm;
	} while(line[0] == '#' || line[0] == '\n' || line[0] == '\r' || 
							line[0] == '\0');

	double data[28];
	int n = sscanf(line,"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
				data+0,data+1,data+2,data+3,data+4,data+5,
				data+6,data+7,data+8,data+9,data+10,data+11,
				data+12,data+13,data+14,data+15,data+16,
				data+17,data+18,data+19,data+20,data+21,
				data+22,data+23,data+24,data+25,data+26,
				data+27);
	if(n != 12 && n != 28) {
		G4Exception("BLTrackFile input","Invalid data format",
				JustWarning, "File Abandoned.");
		status = BLTF_ERROR;
		return status;
	}

	pos[0] = data[0]*unit;
	pos[1] = data[1]*unit;
	pos[2] = data[2]*unit;
	momentum[0] = data[3]*MeV;
	momentum[1] = data[4]*MeV;
	momentum[2] = data[5]*MeV;
	time = data[6]*ns;
	PDGid = (int)data[7];
	eventId = (int)data[8];
	trackId = (int)data[9];
	parentId = (int)data[10];
	weight = data[11];
	if(n == 28) {
		bfield[0] = data[12]*tesla;
		bfield[1] = data[13]*tesla;
		bfield[2] = data[14]*tesla;
		efield[0] = data[15]*(megavolt/meter);
		efield[1] = data[16]*(megavolt/meter);
		efield[2] = data[17]*(megavolt/meter);
		properTime = data[18]*ns;
		pathLength = data[19]*unit;
		polarization[0] = data[20];
		polarization[1] = data[21];
		polarization[2] = data[22];
		initialPos[0] = data[23]*unit;
		initialPos[1] = data[24]*unit;
		initialPos[2] = data[25]*unit;
		initialT = data[26]*ns;
		initialKE = data[27]*MeV;
	} else {
		bfield[0] = 0.0;
		bfield[1] = 0.0;
		bfield[2] = 0.0;
		efield[0] = 0.0;
		efield[1] = 0.0;
		efield[2] = 0.0;
		properTime = 0.0;
		pathLength = 0.0;
		polarization[0] = 0.0;
		polarization[1] = 0.0;
		polarization[2] = 0.0;
		initialPos[0] = 0.0;
		initialPos[1] = 0.0;
		initialPos[2] = 0.0;
		initialT = 0.0;
		initialKE = 0.0;
	}

	return status;
}

