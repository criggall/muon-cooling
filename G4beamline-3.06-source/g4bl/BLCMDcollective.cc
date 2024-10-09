//	BLCMDcollective.cc
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

#include <stdio.h>
#include <vector>

#include "BLCommand.hh"
#include "BLCollectiveComputation.hh"
#include "BLRunManager.hh"
#include "BLNTuple.hh"
#include "BLParam.hh"
#include "BLGlobalField.hh"
#include "BLTime.hh"
#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

/**	class BLCMDcollective - example collective computation.
 *
 *	This class simply computes the mean and sigma of the time stepping
 *	of BLRunManager, and enters them into an NTuple. It optinally generates
 *	an NTuple of field values at specified points in global coordinates.
 **/
class BLCMDcollective : public BLCommand, public BLCollectiveComputation {
	G4double deltaT;
	G4String format;
	G4String filenamePrefix;
	BLNTuple *timeNtuple;
	std::vector<G4ThreeVector> point;
	std::vector<BLNTuple*> fieldNtuple;
	BLRunManager *runManager;
	BLGlobalField *globalField;
	int nstep;
	long prevTime;
	long startTime;
public:
	BLCMDcollective();

	G4String commandName() { return "collective"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	void setupFieldPoint(G4String v);

	virtual void beginCollectiveTracking(std::vector<BLTrackData>& v);

	virtual void collectiveStep(std::vector<BLTrackData>& v);

	virtual void endCollectiveTracking(std::vector<BLTrackData>& v);
};

BLCMDcollective defaultCollectiveCommand;

BLCMDcollective::BLCMDcollective() : point(), fieldNtuple()
{
	registerCommand(BLCMDTYPE_OTHER);
	setSynopsis("Monitor collective computation");
	setDescription("This command computes the means and sigmas related "
		"to the time stepping in BLRunManager (global coordinates), "
		"generating a TimeStep NTuple. If the simulation has "
		"multiple bunches, this NTuple combines them all (and is thus "
		"almost useless).\n\n"
		"This command can also generate field NTuple-s at specified "
		"points in x,y,z -- unnamed parameters should be 'x,y,z' "
		"values for monitoring E and B fields (global coordinates).\n\n"
		"NOTE: This command must come AFTER other commands that "
		"compute collective fields; otherwise stale field values "
		"will be used from the previous time step. If deltaT is set "
		"> 0.0, this command will put the RunManager into collective "
		"mode and set its deltaT; otherwise the previous commands "
		"should do that, and this command won't modify deltaT.\n\n"
		"This command is not placed into the geometry.");

	// initialize class variables here
	deltaT = -1.0*ns;
	format = "";
	filenamePrefix = "";
	timeNtuple = 0;
	runManager = 0; // not available, yet -- see command(...)
	globalField = 0; // not available, yet -- see command(...)
	nstep = 0;
	prevTime = BLTime::timems();
	startTime = BLTime::time();
}

int BLCMDcollective::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	runManager = BLRunManager::getObject();
	globalField = BLGlobalField::getObject();

	handleNamedArgs(namedArgs);

	runManager->registerCollectiveComputation(this);
	if(deltaT > 0.0) {
		runManager->setCollectiveMode(true);
		runManager->setDeltaT(deltaT);
	}

	// setup the NTuple-s
	timeNtuple = BLNTuple::create(format, "NTuples", "TimeStep",
	    "t:Ntracks:meanDT:sigmaDT:meanX:sigmaX:meanY:sigmaY:meanZ:sigmaZ:cpuTime",
	    filenamePrefix+"TimeStep");
	assert(timeNtuple != 0);

	print("");

	for(unsigned i=0; i<argv.size(); ++i)
		setupFieldPoint(argv[i]);

	return 0;
}

void BLCMDcollective::defineNamedArgs()
{
	argDouble(deltaT,"deltaT","Time step (-1 ns).",ns);
	argString(format,"format","Format of NTuples.");
	argString(filenamePrefix,"filenamePrefix","Prefix of NTuple filenames.");
}	

void BLCMDcollective::setupFieldPoint(G4String v)
{
	std::vector<double> vv = getList(v,",");
	if(vv.size() != 3) {
		printError("collective: invalid field point '%s'",
								v.c_str());
		return;
	}
	point.push_back(G4ThreeVector(vv[0],vv[1],vv[2]));
	fieldNtuple.push_back(BLNTuple::create(format, "Field", v,
			"x:y:z:t:Bx:By:Bz:Ex:Ey:Ez",filenamePrefix+v));
	printf("                           Sample Fields: %s\n",v.c_str());
}

void BLCMDcollective::beginCollectiveTracking(std::vector<BLTrackData>& v)
{
	printf("beginCollectiveTracking %d tracks\n",(int)v.size());
	fflush(stdout);

	nstep = 0;
	prevTime = BLTime::timems();
	startTime = BLTime::time();
}

void BLCMDcollective::collectiveStep(std::vector<BLTrackData>& v)
{
#define fmax(A,B) ((A)>(B) ? (A) : (B))
	++nstep;
	long time = BLTime::timems();
	float cpuTime = (float)(time-prevTime)/1000.0f;
	prevTime = time;

	// compute means and sigmas of DT, X, Y, and Z.
	G4double stepTime = runManager->getStepTime();
	G4double sumDT=0.0, sumDT2=0.0, sumX=0.0, sumX2=0.0, 
		 sumY=0.0, sumY2=0.0, sumZ=0.0, sumZ2=0.0;
	int active=0;
	G4double meanZ=std::strtod("nan()",0);
	for(unsigned i=0; i<v.size(); ++i) {
		// get the track pointer, and verify track is alive
		G4Track *track = v[i].track;
		G4TrackStatus trackStatus = track->GetTrackStatus();
		if(trackStatus != fAlive && trackStatus != fStopButAlive)
			continue;
		++active;
		G4ThreeVector position = track->GetPosition();
		G4double time = track->GetGlobalTime();
		G4double dt = time - stepTime;
		sumDT += dt;
		sumDT2 += dt*dt; 
		sumX += position[0];
		sumX2 += position[0]*position[0];
		sumY += position[1];
		sumY2 += position[1]*position[1];
		sumZ += position[2];
		sumZ2 += position[2]*position[2];
	}
	if(active >= 1) {
	    G4double meanDT = sumDT/active;
	    G4double sigmaDT = sqrt(fmax(sumDT2/active - meanDT*meanDT,0.0));
	    G4double meanX = sumX/active;
	    G4double sigmaX = sqrt(fmax(sumX2/active - meanX*meanX,0.0));
	    G4double meanY = sumY/active;
	    G4double sigmaY = sqrt(fmax(sumY2/active - meanY*meanY,0.0));
	    meanZ = sumZ/active;
	    G4double sigmaZ = sqrt(fmax(sumZ2/active - meanZ*meanZ,0.0));
	    double data[11];
	    data[0] = stepTime;
	    data[1] = active;
	    data[2] = meanDT;
	    data[3] = sigmaDT;
	    data[4] = meanX;
	    data[5] = sigmaX;
	    data[6] = meanY;
	    data[7] = sigmaY;
	    data[8] = meanZ;
	    data[9] = sigmaZ;
	    data[10] = cpuTime;
	    assert(timeNtuple != 0);
	    timeNtuple->appendRow(data,sizeof(data)/sizeof(data[0]));

	    for(unsigned i=0; i<point.size(); ++i) {
		G4ThreeVector p=point[i];
		double pp[4],ff[6];
		pp[0] = p[0];
		pp[1] = p[1];
		pp[2] = p[2];
		pp[3] = stepTime;
		globalField->GetFieldValue(pp,ff);
		double data[10];
		data[0] = p[0]/mm;
		data[1] = p[1]/mm;
		data[2] = p[2]/mm;
		data[3] = stepTime/ns;
		data[4] = ff[0]/tesla;
		data[5] = ff[1]/tesla;
		data[6] = ff[2]/tesla;
		data[7] = ff[3]/(megavolt/meter);
		data[8] = ff[4]/(megavolt/meter);
		data[9] = ff[5]/(megavolt/meter);
		assert(i < fieldNtuple.size());
		fieldNtuple[i]->appendRow(data,sizeof(data)/sizeof(data[0]));
	    }
	}
	printf("CollectiveStep %d t=%.3f  %d tracks (%d active)  meanZ=%.3f\n",
		nstep,runManager->getStepTime(),(int)v.size(),active,meanZ);
	fflush(stdout);
}

void BLCMDcollective::endCollectiveTracking(std::vector<BLTrackData>& v)
{
	long time = BLTime::time() - startTime;
	printf("endCollectiveTracking %d tracks  %ld seconds\n",(int)v.size(),
								time);
	fflush(stdout);
}
