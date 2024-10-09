//	BLTrackNTuple.cc
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

#include <ctype.h>

#include "BLAssert.hh"
#include "BLTrackNTuple.hh"
#include "BLNTuple.hh"

static const int STANDARD=12, TRACE=18, EXTENDED=28;
static const char *fieldNames[] = {
	"x", "y", "z", "Px", "Py", "Pz", "t", "PDGid", "EventID", "TrackID",
	"ParentID", "Weight",
	"Bx", "By", "Bz", "Ex", "Ey", "Ez",
	"ProperTime", "PathLength", "PolX", "PolY", "PolZ", 
	"InitX", "InitY", "InitZ", "InitT", "InitKE"
};

BLTrackNTuple::BLTrackNTuple()
{
	ntuple = 0;
	require = "";
	eval = 0;
	noSingles = 0;
	nFields = 0;
	coordinateType = BLCOORD_CENTERLINE;
}

BLTrackNTuple *BLTrackNTuple::create(G4String format, G4String category,
	G4String name, G4String filename, BLCoordinateType _coordinateType,
	G4String _require, int _noSingles)
{
	BLTrackNTuple *t = new BLTrackNTuple();
	t->require = _require;
	t->noSingles = _noSingles;
	t->coordinateType = _coordinateType;
	t->nFields = STANDARD;

	if(t->require != "") t->eval = new BLEvaluator();

	// get nFields and fields based on format
	for(G4String::size_type i=0; i<format.size(); ++i)
		format[i] = tolower(format[i]);
	G4String::size_type j=format.find("trace");
	if(j != format.npos) {
		t->nFields = TRACE;
		format.erase(j,format.npos);
	}
	j=format.find("extended");
	if(j != format.npos) {
		t->nFields = EXTENDED;
		format.erase(j,format.npos);
	}
	j=format.find("for009");
	if(j != format.npos) {
		t->nFields = TRACE;
	}
	G4String fields=fieldNames[0];
	for(int i=1; i<t->nFields; ++i) {
		fields += ":";
		fields += fieldNames[i];
	}

	t->ntuple = BLNTuple::create(format,category,name,fields,filename);
	if(!t->ntuple) {
		delete t;
		t = 0;
	}
	return t;
}

bool BLTrackNTuple::appendTrack(const G4Track *track)
{
	return appendTrack(track,track->GetGlobalTime(),track->GetPosition(),
			track->GetMomentum(),track->GetProperTime(),
			track->GetTrackLength());
}

bool BLTrackNTuple::appendTrack(const G4Track *track, double time, 
				G4ThreeVector position, G4ThreeVector momentum,
				double properTime, double trackLength)
{
	G4RunManager* runmgr = G4RunManager::GetRunManager();
	const G4Event* event = runmgr->GetCurrentEvent();
	int evId = event->GetEventID();
	G4ThreeVector globalPosition = position;

	// transform to desired coordinates, if available
	BLCoordinates *coord = (BLCoordinates *)track->GetUserInformation();
	if(coord && coord->isValid()) {
		coord->getCoords(coordinateType,position);
		momentum = coord->getRotation() * momentum;
	} else {
		printf("BLTrackNTuple::SteppingAction: track has no "
			"BLCoordinates object\n");
	}

	double data[EXTENDED];
	BLAssert((unsigned int)nFields <= sizeof(data)/sizeof(data[0]));
	data[0] = position[0]/mm;		// x (mm)
	data[1] = position[1]/mm;		// y (mm)
	data[2] = position[2]/mm;		// z (mm)
	data[3] = momentum[0]/MeV;		// Px (MeV/c)
	data[4] = momentum[1]/MeV;		// Py (MeV/c)
	data[5] = momentum[2]/MeV;		// Pz (MeV/c)
	data[6] = time/ns;			// t (ns)
	data[7] = track->GetDefinition()->GetPDGEncoding();
	data[8] = evId;				// Event ID
	data[9] = BLManager::getObject()->getExternalTrackID(track);
	data[10] = BLManager::getObject()->getExternalParentID(track);
	data[11] = track->GetWeight();		// Weight
	if(nFields > STANDARD) {
		G4double field[6], point[4];
		point[0] = globalPosition[0];
		point[1] = globalPosition[1];
		point[2] = globalPosition[2];
		point[3] = time;
		BLGlobalField::getObject()->GetFieldValue(point,field);
		G4ThreeVector B(field[0],field[1],field[2]);
		G4ThreeVector E(field[3],field[4],field[5]);
		if(coord && coord->isValid()) {
			B = coord->getRotation() * B;
			E = coord->getRotation() * E;
		}
		data[12] = B[0]/tesla;
		data[13] = B[1]/tesla;
		data[14] = B[2]/tesla;
		data[15] = E[0]/(megavolt/meter);
		data[16] = E[1]/(megavolt/meter);
		data[17] = E[2]/(megavolt/meter);
	}
	if(nFields > TRACE) {
		G4ThreeVector pol, vertexPosition;
		pol = track->GetPolarization();
		if(coord && coord->isValid()) {
			pol = coord->getRotation() * pol;
		}
		vertexPosition = track->GetVertexPosition();
		data[18] = properTime/ns;
		data[19] = trackLength/mm;
		data[20] = pol[0];
		data[21] = pol[1];
		data[22] = pol[2];
		data[23] = vertexPosition[0]/mm;
		data[24] = vertexPosition[1]/mm;
		data[25] = vertexPosition[2]/mm;
		data[26] = (time - track->GetLocalTime())/ns;
		data[27] = track->GetVertexKineticEnergy()/MeV;
	}

	// implement require
	if(eval) {
		for(int i=0; i<nFields; ++i)
			eval->setVariable(fieldNames[i],data[i]);
		double v=eval->evaluate(require.c_str());
		if(!eval->isOK())
			G4Exception("BLTrackNTuple",
				"Invalid require expression",
				FatalException,require);
		if(fabs(v) < 1.0e-16) return false;
	}

	if(noSingles) 
		ntuple->doCallbacks(data,nFields);
	else
		ntuple->appendRow(data,nFields); // calls doCallbacks()

	return true;
}

G4String BLTrackNTuple::getFormatList()
{
	G4String s = BLNTuple::getFormatList();
	s += " Extended asciiExtended";
	if(s.find("root") != s.npos) s += " rootExtended";
	return s;
}
