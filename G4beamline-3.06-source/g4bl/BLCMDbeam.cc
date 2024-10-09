//	BLCMDbeam.cc
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
#define _USE_MATH_DEFINES
#include <math.h>

#include "G4RunManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "CLHEP/Units/SystemOfUnits.h"
using namespace CLHEP;

#include "BLAssert.hh"
#include "BLBeam.hh"
#include "BLParam.hh"
#include "BLGroup.hh"
#include "BLCoordinates.hh"
#include "BLNTuple.hh"
#include "BLTrackFile.hh"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif

const G4double UNDEFINED = -3.7e21;
const G4int ALL_EVENTS = 0x7FFFFFFF; // used to indicate an unset argument
// TrackFields is for Root input
const char TrackFields[] =
    "x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight:PolX:PolY:PolZ";
const unsigned NTrackFields = 12;
const unsigned NTrackFields2 = 15;

/**	class BLCMDbeam implements the beam command.
 **/
class BLCMDbeam : public BLBeam, public BLCommand {
	enum BeamType { NONE, GAUSSIAN, RECTANGULAR, ASCII, ELLIPSE, ROOT };
	BeamType type;
	G4String particle;
	G4int eventsGenerated;
	G4int nEvents;
	G4int eventID;
	G4int firstEvent;
	G4int lastEvent;
	G4double beamX;
	G4double beamY;
	G4double beamZ;
	G4double maxR;
	G4String rotation;
	G4RotationMatrix *rotationMatrix;
	G4ThreeVector position;
	G4int renumber;
	G4double weight;
	G4int secondaryTrackID;
	// Gaussian beam arguments:
	G4double meanMomentum;
	G4double sigmaX;
	G4double sigmaY;
	G4double sigmaZ;
	G4double sigmaXp;
	G4double sigmaYp;
	G4double sigmaP;
	G4double sigmaT;
	G4double sigmaE;
	G4double meanXp;
	G4double meanYp;
	G4double meanT;
	G4String polarization;
	// rectangular beam arguments:
	G4double beamHeight;
	G4double beamWidth;
	// Histo beam arguments:
	G4String filename;
	G4String directory;
	G4int uid;
	G4String name;
	BLNTuple *ntuple;
	double *data;
	// ASCII beam arguments:
	G4String format;
	// internal variables
	G4int index;
	G4ParticleGun *particleGun;
	G4ParticleDefinition *particleDefinition;
	BLTrackFile *trackFile;
	G4int nvar;
	G4int prevEventID;
	G4ThreeVector polarization3V;
	// Random number: Gaussian if sigma>0; flat if sigma<0 (-sigma is
	// the halfwidth).
	G4double myrand(G4double mean, G4double sigma) {
		if(sigma >= 0.0) return sigma*CLHEP::RandGauss::shoot() + mean;
		return mean+sigma-2.0*sigma*G4UniformRand();
	}
public:
	/// Constructor.
	BLCMDbeam();

	// (accept the default copy constructor)

	/// commandName() returns "beam".
	virtual G4String commandName() { return "beam"; }
	
	/// command() implements the beam command.
	virtual int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// printBeam() will print a description.
	virtual void printBeam();

	/// defineNamedArgs() defines the named arguments for this command.
	virtual void defineNamedArgs();

	/// getNEvents() returns the # events to process.
	virtual int getNEvents() const { return nEvents; }

	/// init() will initialize internal variables.
	virtual void init();

	/// generateReferenceParticle() generates the reference particle.
	virtual bool generateReferenceParticle(G4Event *event)
		{ return false; }

	/// nextBeamEvent() generates the next beam event.
	virtual bool nextBeamEvent(G4Event *event);

	/// summary() will print a summary, if necessary.
	virtual void summary() { }
};

BLCMDbeam defineBeam;

BLCMDbeam::BLCMDbeam() : rotation() , position()
{
	registerCommand(BLCMDTYPE_BEAM);
	setSynopsis("Define the Beam.");
	setDescription("The beam command is: beam type arg1=v1 ...\n\n"
		"Types are: gaussian, rectangular, ellipse, ascii, and root.\n\n"
		"Gaussian beams are randomly generated to emanate from "
		"beamX,beamY,beamZ with the given "
		"sigmas; negative sigma means flat with |sigma| as "
		"halfwidth.\n\n"
		"Rectangular beams are randomly generated to emanate "
		"from the rectangle beamHeight by beamWidth centered at "
		"beamX,beamY,beamZ. \n\n"
		"Ellipse beams are randomly generated on "
		"the ellipses in (X,Xp), (Y,Yp), (T,E), with meanE determined "
		"from meanP and the sigmas used as half-widths; tracks\n"
		"are generated on the ellipse with uniform density when "
		"plotted with scales such that the ellipse is a circle.\n\n"
		"ASCII beams are read from a file using the format "
		"specified; the formats supported are: BLTrackFile2. "
		"The value of polarization in the file overwrites the value "
		"from the command line. "
		"Note the original BLTrackFile format is supported, but the "
		"additional fields are not present. Note also that properTime, "
		"pathLength, initialPosition, initialT, and initialKE are for "
		"THIS simulation, and their values in the file are ignored.\n\n"
		"Root beams are read from a .root file using the "
		"TNtuple named directory/name in the file. It must have the "
		"same fields as used in BLTrackFile2 format. Note that "
		"EventIDs greater than 16,777,216 will be rounded.\n\n"
		"When reading a file (ascii or root), "
		"beamX and beamY are added to input tracks; if beamZ is set "
		"it will overwrite the z of the track, but if it is not set "
		"the z of the track in the file is kept.\n\n"
		"All coordinates are centerline coordinates.\n\n"
		"Multiple beam commands can be given, and they will generate\n"
		"events in the order they appear in input.file.\n\n"
		"Events are generated starting at firstEvent, until either\n"
		"nEvents have been generated or lastEvent would be "
		"exceeded.\n\n"
		"This command places itself into the geometry.\n\n"
		"Beam tracks with momentum=0 will take one step as a stoppd "
		"particle (as desired).\n\n"
		"For gaussian, rectangular, and ellipse beams, the beam "
		"particle can be given as either a particle name or "
		"its integer PDGid. Some common beam particle names are: "
		"proton, anti_proton, pi+, pi-, mu+, mu-, e+, e-, kaon+, "
		"kaon-, kaon0, nu_e, anti_nu_e.\n"
		"See the User's Guide for a complete list of particle "
		"names; ions can also be specified by integer PDGid.\n\n"
		"Ions are supported as beam particles, but only fully-ionized "
		"ions, and not for Root input (a Root NTuple field cannot hold "
		"the PDGid). "
		"Ions are specified by an integer PDGid=100ZZZAAA0, where ZZZ "
		"is the 3-digit charge (# protons), and AAA is the 3-digit "
		"atomic number. Note you must use physics processes "
		"appropriate for ions.");
	// initial default values:
	type = NONE;
	particle = "mu+";
	eventsGenerated = 0;
	nEvents = ALL_EVENTS;
	eventID = ALL_EVENTS;
	firstEvent = -1;
	lastEvent = ALL_EVENTS;
	meanMomentum = 200.0*MeV;
	beamX = 0.0;
	beamY = 0.0;
	beamZ = UNDEFINED;
	maxR = 1.0*kilometer;
	renumber = 0;
	weight = 1.0;
	secondaryTrackID = 1001;
	rotationMatrix = 0;
	sigmaX = 0.0;
	sigmaY = 0.0;
	sigmaZ = 0.0;
	sigmaXp = 0.0;
	sigmaYp = 0.0;
	sigmaP = 0.0;
	sigmaT = 0.0;
	sigmaE = UNDEFINED;
	meanXp = 0.0;
	meanYp = 0.0;
	meanT = 0.0;
	polarization = "";
	beamHeight = 0.0;
	beamWidth = 0.0;
	filename = "";
	directory = "";
	uid = 0;
	name = "";
	ntuple = 0;
	data = 0;
	format = "BLTrackFile";
	index = 0;
	particleDefinition = 0;
	particleGun = 0;
	trackFile = 0;
	nvar = 0;
	prevEventID = -9999;
	polarization3V = G4ThreeVector(0.0,0.0,0.0);
}

int BLCMDbeam::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	BLCMDbeam *b = new BLCMDbeam(*this);

	if(argv.size() < 1) { // default
		b->type = GAUSSIAN;
		b->nEvents = 1;
	} else if(argv[0] == "gaussian" || argv[0] == "gauss") {
		b->type = GAUSSIAN;
		b->nEvents = 1;
	} else if(argv[0] == "rectangular" || argv[0] == "rect") {
		b->type = RECTANGULAR;
		b->nEvents = 1;
	} else if(argv[0] == "ascii" || argv[0] == "ASCII") {
		b->type = ASCII;
		b->nEvents = ALL_EVENTS;
	} else if(argv[0] == "ellipse" || argv[0] == "ELLIPSE") {
		b->type = ELLIPSE;
		b->nEvents = 1;
	} else if(argv[0] == "root") {
		b->type = ROOT;
		b->nEvents = ALL_EVENTS;
	} else {
		printError("beam: invalid type '%s'",argv[0].c_str());
	}

	b->rotation = "";
	int retval = b->handleNamedArgs(namedArgs);

	if(b->firstEvent != -1 && b->lastEvent != ALL_EVENTS)
		b->nEvents = b->lastEvent - b->firstEvent + 1;

	if(b->rotation != "") {
		b->rotationMatrix = stringToRotationMatrix(b->rotation);
	} else {
		b->rotationMatrix = new G4RotationMatrix();
	}
	// as usual, the order seems backward because (C R C^-1) C = C R.
	*b->rotationMatrix = *BLCoordinates::getCurrentRotation() * 
							*b->rotationMatrix;
	G4ThreeVector local(b->beamX,b->beamY,
					(b->beamZ==UNDEFINED ? 0.0 : b->beamZ));
	BLCoordinates::getCurrentGlobal(local,b->position);

	// ensure beam is within the world
	BLGroup::getWorld()->setMinWidth(fabs(b->position[0])*2.0);
	BLGroup::getWorld()->setMinHeight(fabs(b->position[1])*2.0);
	BLGroup::getWorld()->setMinLength(fabs(b->position[2])*2.0);

	BLManager::getObject()->registerBeam(b);

	// handle polarization
	if(b->polarization != "") {
		std::vector<G4double> polar = getList(b->polarization,",");
		if(polar.size() == 3)
			b->polarization3V = G4ThreeVector(polar[0],polar[1],
								polar[2]);
		else
			G4Exception("beam command","Invalid Polarization",
				FatalException, "Not 3-vector");
	}

	b->printBeam();

	return retval;
}

void BLCMDbeam::defineNamedArgs()
{
	argString(particle,"particle","Beam particle name (default=mu+).");
	argInt(nEvents,"nEvents","Number of events to process (default=1 "
		"for generating events, default=ALL for reading files) "
		"set to lastEvent-firstEvent+1 if both are set.");
	argInt(firstEvent,"firstEvent","First event # to process (default "
		"is the next sequential eventID, 1 if none)");
	argInt(lastEvent,"lastEvent","Last  (highest) event # to process");
	argDouble(beamX,"beamX","Beam location in X (mm)");
	argDouble(beamY,"beamY","Beam location in Y (mm)");
	argDouble(beamZ,"beamZ","Beam location in Z (mm)");
	argDouble(maxR,"maxR","Beam maximum radius (mm)");
	argString(rotation,"rotation","Rotation of the beam");
	argInt(renumber,"renumber","If nonzero, renumber events sequentially."); 
	argDouble(weight,"weight","Weight for events, overwritten by value from input file (1.0)."); 
	argInt(secondaryTrackID,"secondaryTrackID","The next TrackID for secondaries (1001)."); 
	argDouble(meanMomentum,"meanMomentum",
				"Gaussian Beam mean momentum (MeV/c)");
	argDouble(meanMomentum,"meanP","Synonym for meanMomentum.");
	argDouble(meanMomentum,"P","Synonym for meanMomentum.");
	argDouble(sigmaX,"sigmaX","Gaussian Beam sigma in X (mm)");
	argDouble(sigmaY,"sigmaY","Gaussian Beam sigma in Y (mm)");
	argDouble(sigmaZ,"sigmaZ","Gaussian Beam sigma in Z (mm)");
	argDouble(sigmaXp,"sigmaXp","Gaussian Beam sigma in dxdz (slope)");
	argDouble(sigmaYp,"sigmaYp","Gaussian Beam sigma in dydz (slope)");
	argDouble(sigmaP,"sigmaP","Gaussian Beam sigma in P (MeV/c)");
	argDouble(sigmaT,"sigmaT","Gaussian Beam sigma in T (ns)");
	argDouble(sigmaE,"sigmaE","Elliptical Beam sigma in E (MeV)");
	argDouble(meanXp,"meanXp","Gaussian Beam mean in Xp (slope)");
	argDouble(meanYp,"meanYp","Gaussian Beam mean in Yp (slope)");
	argDouble(meanT,"meanT","Gaussian Beam mean in T (ns)");
	argString(polarization,"polarization","Polarization 3-vector (0,0,0)");
	argDouble(beamHeight,"beamHeight","Rectangular Beam height (mm)");
	argDouble(beamWidth,"beamWidth","Rectangular Beam width (mm)");
	argString(filename,"filename","input file name");
	argString(filename,"file","synonym for filename.");
	argString(directory,"directory","Root-file directory of NTuple");
	argString(directory,"category","Deprecated synonym for directory.");
	argInt(uid,"uid","HistoScope uid of NTuple");
	argString(name,"name","Root name of NTuple.");
	argString(format,"format","ASCII file format (Default=BLTrackFile)");
	argDouble(meanXp,"beamXp","Synonym for meanXp.");
	argDouble(meanYp,"beamYp","Synonym for meanYp.");
	argDouble(beamX,"x","Synonym for beamX.");
	argDouble(beamY,"y","Synonym for beamY.");
	argDouble(beamZ,"z","Synonym for beamZ.");
}

void BLCMDbeam::printBeam()
{
	switch(type) {
	case NONE:
		printf("beam    NONE    Invalid beam command\n");
		break;
	case GAUSSIAN:
		printf("beam    GAUSSIAN particle=%s ",particle.c_str());
		printf("nEvents=%d ",nEvents);
		printf("firstEvent=%d ",firstEvent);
		printf("lastEvent=%d ",lastEvent);
		printf("beamX=%.1f ",beamX);
		printf("beamY=%.1f ",beamY);
		printf("beamZ=%.1f ",beamZ);
		printf("maxR=%.1f ",maxR);
		printf("\n\t\t");
		printf("meanMomentum=%.1f ",meanMomentum);
		printf("weight=%1f ",weight);
		printf("\n\t\t");
		printf("sigmaX=%.1f ",sigmaX);
		printf("sigmaY=%.1f ",sigmaY);
		printf("sigmaZ=%.1f ",sigmaZ);
		printf("sigmaXp=%.5f ",sigmaXp);
		printf("sigmaYp=%.5f ",sigmaYp);
		printf("\n\t\t");
		printf("sigmaP=%.1f ",sigmaP);
		printf("sigmaT=%.3f ",sigmaT);
		printf("meanXp=%.5f ",meanXp);
		printf("meanYp=%.5f ",meanYp);
		printf("meanT=%.3f ",meanT);
		printf("\n");
		break;
	case RECTANGULAR:
		printf("beam    RECTANGULAR particle=%s ",particle.c_str());
		printf("nEvents=%d ",nEvents);
		printf("firstEvent=%d ",firstEvent);
		printf("lastEvent=%d ",lastEvent);
		printf("beamX=%.1f ",beamX);
		printf("beamY=%.1f ",beamY);
		printf("beamZ=%.1f ",beamZ);
		printf("maxR=%.1f ",maxR);
		printf("\n\t\t");
		printf("meanMomentum=%.1f ",meanMomentum);
		printf("weight=%1f ",weight);
		printf("\n\t\t");
		printf("beamHeight=%.1f ",beamHeight);
		printf("beamWidth=%.1f ",beamWidth);
		printf("sigmaXp=%.5f ",sigmaXp);
		printf("sigmaYp=%.5f ",sigmaYp);
		printf("\n\t\t");
		printf("sigmaP=%.1f ",sigmaP);
		printf("sigmaT=%.3f ",sigmaT);
		printf("meanXp=%.5f ",meanXp);
		printf("meanYp=%.5f ",meanYp);
		printf("meanT=%.3f ",meanT);
		printf("\n");
		break;
	case ROOT:
		printf("beam    ROOT   ");
		printf("nEvents=%d ",nEvents);
		printf("firstEvent=%d ",firstEvent);
		printf("lastEvent=%d ",lastEvent);
		printf("beamZ=%.1f ",beamZ);
		printf("maxR=%.1f ",maxR);
		printf("renumber=%d ",renumber);
		printf("weight=%1f ",weight);
		printf("\n\t\t");
		printf("filename=%s ",filename.c_str());
		printf("directory=%s ",directory.c_str());
		printf("name=%s ",name.c_str());
		printf("\n\t\t");
		printf("\n");
		break;
	case ASCII:
		printf("beam    ASCII   ");
		printf("nEvents=%d ",nEvents);
		printf("firstEvent=%d ",firstEvent);
		printf("lastEvent=%d ",lastEvent);
		printf("beamZ=%.1f ",beamZ);
		printf("maxR=%.1f ",maxR);
		printf("renumber=%d ",renumber);
		printf("weight=%1f ",weight);
		printf("\n\t\t");
		printf("filename=%s ",filename.c_str());
		printf("format=%s ",format.c_str());
		printf("\n\t\t");
		printf("\n");
		break;
	case ELLIPSE:
		printf("beam    ELLIPSE particle=%s ",particle.c_str());
		printf("nEvents=%d ",nEvents);
		printf("firstEvent=%d ",firstEvent);
		printf("lastEvent=%d ",lastEvent);
		printf("beamX=%.1f ",beamX);
		printf("beamY=%.1f ",beamY);
		printf("beamZ=%.1f ",beamZ);
		printf("maxR=%.1f ",maxR);
		printf("\n\t\t");
		printf("meanMomentum=%.1f ",meanMomentum);
		printf("weight=%1f ",weight);
		printf("\n\t\t");
		printf("sigmaX=%.1f ",sigmaX);
		printf("sigmaY=%.1f ",sigmaY);
		printf("sigmaXp=%.5f ",sigmaXp);
		printf("sigmaYp=%.5f ",sigmaYp);
		printf("\n\t\t");
		printf("sigmaE=%.1f ",sigmaE);
		printf("sigmaT=%.3f ",sigmaT);
		printf("meanXp=%.5f ",0.0);
		printf("meanYp=%.5f ",0.0);
		printf("meanT=%.3f ",0.0);
		printf("\n");
		break;
	}
}

void BLCMDbeam::init()
{
	eventsGenerated = 0;
	index = 0;
	if(particleDefinition != 0) return;

	if(isdigit(particle(0)) || particle(0) == '-') {
		G4int pdgid = atoi(particle.c_str());
		particleDefinition = G4ParticleTable::GetParticleTable()->
							FindParticle(pdgid);
		if(!particleDefinition && pdgid > 1000000000) { // Ions
		    G4int Z = (pdgid/10000) % 1000;
		    G4int A = (pdgid/10) % 1000;
		    G4double excitationEnergy = 0.0*keV;
		    particleDefinition = G4ParticleTable::GetParticleTable()->
				GetIonTable()->GetIon(Z,A,excitationEnergy);
		} else if(!particleDefinition && pdgid < -1000000000) {
		    G4Exception("beam command","UnknownParticle",FatalException,
						"Anti-ions not supported");
		}
	} else {
		particleDefinition = G4ParticleTable::GetParticleTable()->
					FindParticle(particle);
	}
	if(!particleDefinition)
		G4Exception("beam command","UnknownParticle",FatalException,
						"Unknown particle type");
	particleGun = new G4ParticleGun(1);
	particleGun->SetParticleDefinition(particleDefinition);

	switch(type) {
	case NONE:
		break;
	case GAUSSIAN:
	case RECTANGULAR:
		break;
	case ELLIPSE:
		if(sigmaE == UNDEFINED) sigmaE = sigmaP;
		break;
	case ROOT:
		{ ntuple = BLNTuple::read("Root",directory,name,TrackFields,
				  				filename);
		  if(ntuple == 0) break;
		  nvar = ntuple->getNData();
		  if(nvar == NTrackFields) {
		  	printf("beam: Root TNtuple %s:%s/%s\n",
			    filename.c_str(),directory.c_str(),name.c_str());
		  } else if(nvar == NTrackFields2) {
		  	printf("beam: Root TNtuple %s:%s/%s with polarization\n"
			    , filename.c_str(),directory.c_str(),name.c_str());
		  } else {
		  	delete ntuple;
			G4Exception("beam command", "Invalid Root input",
				JustWarning, "NTuple abandoned");
			nvar = 0;
			break;
		  }
		  data = new double[nvar];
		}
		break;
	case ASCII:
		if(format == "BLTrackFile") {
			trackFile = new BLTrackFile(filename,"","r",2);
		} else {
			printError("beam: Unknown format '%s', known formats: "
					"BLTrackFile.",format.c_str());
		}
		break;
	}
}

bool BLCMDbeam::nextBeamEvent(G4Event *event)
{
	G4double mass = particleDefinition->GetPDGMass();
	G4ThreeVector pos;
	G4ThreeVector direction;
	G4double time = 0.0;
	G4double momentum = 0.0;
	G4double ke = 0.0;
	G4int PDGid, trackID, parentID;

	BLManager *manager = BLManager::getObject();

	// default eventID -- changed when reading a file unless renumber!=0
	eventID = manager->getEventID();
	if(eventsGenerated == 0 && firstEvent != -1)
		eventID = firstEvent;
	if(eventID <= 0) {
		eventID = 1;
		printf("EventID <= 0 skipped\n");
	}

	if(++eventsGenerated > nEvents) return false;

	// loop - must not generate any event with r > maxR
	// beam pos = (0,0,0), beam nominal direction = (0,0,1)
	// rotation and offset come later
	trackID = -1;
	parentID = 0;
	bool again=false;
	do {
	    again = false;
	    switch(type) {
	    case NONE:
		return false;
	    case GAUSSIAN:
		setRandomSeedToGenerate(eventID);
		// meanX=beamX, meanY=beamY (implicitly)
		pos[0] = myrand(0.0,sigmaX);
		pos[1] = myrand(0.0,sigmaY);
		pos[2] = myrand(0.0,sigmaZ);
		direction[0] = myrand(meanXp,sigmaXp);
		direction[1] = myrand(meanYp,sigmaYp);
		direction[2] = 1.0/sqrt(1.0 + direction[0]*direction[0] +
						direction[1]*direction[1]);
		direction[0] *= direction[2];
		direction[1] *= direction[2];
		momentum = myrand(meanMomentum,sigmaP);
		time = myrand(meanT,sigmaT);
		trackID = 1;
		break;
	    case RECTANGULAR:
		setRandomSeedToGenerate(eventID);
		pos[0] = beamWidth*G4UniformRand()-beamWidth/2.0;
		pos[1] = beamHeight*G4UniformRand()-beamHeight/2.0;
		pos[2] = myrand(0.0,sigmaZ);
		direction[0] = myrand(meanXp,sigmaXp);
		direction[1] = myrand(meanYp,sigmaYp);
		direction[2] = 1.0/sqrt(1.0 + direction[0]*direction[0] +
						direction[1]*direction[1]);
		direction[0] *= direction[2];
		direction[1] *= direction[2];
		momentum = myrand(meanMomentum,sigmaP);
		time = myrand(meanT,sigmaT);
		trackID = 1;
		break;
	    case ROOT:
		{ do {
			if(ntuple == 0 || data == 0 || nvar == 0) return false;
			if(!ntuple->readRow(data,nvar)) return false;
		  	if(renumber == 0) eventID = (G4int)data[8];
		  } while(renumber == 0 && eventID < firstEvent);
		  particleDefinition = G4ParticleTable::GetParticleTable()->
			  			FindParticle((G4int)data[7]);
		  // Ions not supported here (float cannot hold their PDGid)
		  if(!particleDefinition) {
			char tmp[64];
			sprintf(tmp,"Invalid PDGid=%d -- Track abandoned",
								(int)data[7]);
			G4Exception("beam command","Invalid PDGid",
							JustWarning, tmp);
			again = true;
			continue;
		  }
		  mass = particleDefinition->GetPDGMass();
		  particleGun->SetParticleDefinition(particleDefinition);
		  pos[0] = data[0]*mm;
		  pos[1] = data[1]*mm;
		  pos[2] = data[2]*mm;
		  data[3] *= MeV;
		  data[4] *= MeV;
		  data[5] *= MeV;
		  momentum = sqrt(data[3]*data[3]+data[4]*data[4]+data[5]*data[5]);
		  direction[0] = data[3]/momentum;
		  direction[1] = data[4]/momentum;
		  direction[2] = data[5]/momentum;
		  time = data[6]*ns;
		  trackID = (G4int)data[9];
		  parentID = (G4int)data[10];
		  if(nvar >= NTrackFields) weight = data[11];
		  if(nvar >= NTrackFields2)
		  	polarization3V = G4ThreeVector(data[12],data[13],
								data[14]);
		}
		if(beamZ != UNDEFINED) pos[2] = 0.0;
		break;
	    case ASCII:
		if(format == "BLTrackFile") {
			BLAssert(trackFile != 0);
			do {
				int tmpID=0;
				G4ThreeVector bfield,efield,pol,initialPos;
				G4double properTime,pathLength,initialT,
								initialKE;
				BLTrackFileStatus stat=trackFile->read2(pos,
					time,direction,PDGid,tmpID,trackID,
					parentID,weight,bfield,efield,
					properTime,pathLength,pol,initialPos,
					initialT,initialKE);
				if(pol.mag2()>0.0) polarization3V = pol.unit();
				if(stat == BLTF_ERROR)
					G4Exception("beam command",
						"Invalid BLTrackFile input",
						JustWarning, "File abandoned");
				if(stat != BLTF_OK) {
					delete trackFile;
					trackFile = 0;
					return false;
				}
		  		if(renumber == 0) eventID = tmpID;
		  	} while(renumber == 0 && eventID < firstEvent);
		  	particleDefinition = G4ParticleTable::GetParticleTable()
							->FindParticle(PDGid);
			if(!particleDefinition && PDGid > 1000000000) { // Ions
			    G4int Z = (PDGid/10000) % 1000;
			    G4int A = (PDGid/10) % 1000;
		    	    G4double excitationEnergy = 0.0*keV;
			    particleDefinition = G4ParticleTable::
			       GetParticleTable()->GetIonTable()->GetIon(Z,A,excitationEnergy);
			} else if(!particleDefinition && PDGid < -1000000000) {
			    G4Exception("beam command","UnknownParticle",
			    	FatalException, "Anti-ions not supported");
			}
		  	if(!particleDefinition) {
			    char tmp[64];
			    sprintf(tmp,"Invalid PDGid=%d -- Track abandoned",
									PDGid);
			    G4Exception("beam command","Invalid PDGid",
							JustWarning, tmp);
			    again = true;
			    continue;
		  	}
		  	mass = particleDefinition->GetPDGMass();
		  	particleGun->SetParticleDefinition(particleDefinition);
			momentum = sqrt(direction[0]*direction[0]+
					direction[1]*direction[1]+
					direction[2]*direction[2]);
			direction[0] /= momentum;
			direction[1] /= momentum;
			direction[2] /= momentum;
		} else {
			G4Exception("beam command","Invalid Format",
					FatalException, "");
		}
		if(beamZ != UNDEFINED) pos[2] = 0.0;
		break;
	    case ELLIPSE:
		setRandomSeedToGenerate(eventID);
		// meanX=beamX, meanY=beamY (implicitly)
		double a = 2.0*M_PI*G4UniformRand();
		double b = 2.0*M_PI*G4UniformRand();
		double c = 2.0*M_PI*G4UniformRand();
		pos[0] = sigmaX*cos(a);
		pos[1] = sigmaY*cos(b);
		pos[2] = myrand(0.0,sigmaZ);
		direction[0] = sigmaXp*sin(a);
		direction[1] = sigmaYp*sin(b);
		direction[2] = 1.0/sqrt(1.0 + direction[0]*direction[0] +
						direction[1]*direction[1]);
		direction[0] *= direction[2];
		direction[1] *= direction[2];
		double E = sqrt(meanMomentum*meanMomentum + mass*mass) +
								sigmaE*sin(c);
		momentum = sqrt(E*E - mass*mass);
		time = sigmaT*cos(c);
		trackID = 1;
		break;
	    }
	    if(manager->skipEvent(eventID)) {
		if(manager->getSteppingVerbose())
			printf("eventcuts: EventID %d skipped\n",eventID);
	    	++eventID;
		again = true;
		continue;
	    }
	} while(again || sqrt(pos[0]*pos[0]+pos[1]*pos[1]) > maxR);

	BLAssert(trackID >= 0);

	if(eventID > lastEvent) return false;

	// apply rotation and offset
	if(rotationMatrix) {
		direction = *rotationMatrix * direction;
		pos = *rotationMatrix * pos;
		polarization3V = *rotationMatrix * polarization3V;
	}
	pos += position;

	ke = sqrt(momentum*momentum + mass*mass) - mass;
	particleGun->SetParticleTime(time);
	particleGun->SetParticlePosition(pos);
	particleGun->SetParticleEnergy(ke);
	particleGun->SetParticleMomentumDirection(direction);
	if(particleDefinition->GetPDGSpin() != 0.0)
		particleGun->SetParticlePolarization(polarization3V);
	else
		particleGun->SetParticlePolarization(G4ThreeVector(0.,0.,0.));
	particleGun->GeneratePrimaryVertex(event);
	event->SetEventID(eventID);
	event->GetPrimaryVertex()->SetWeight(weight);
	if(eventID != prevEventID) {
		setRandomSeedToTrack(eventID);
		manager->clearTrackIDMap();
		manager->setNextSecondaryTrackID(secondaryTrackID);
		prevEventID = eventID;
	}

	manager->setPrimaryTrackID(trackID,parentID);
	if(trackID >= secondaryTrackID) {
		G4Exception("beam command","Large Primary TrackID",JustWarning,
			"Confusion with secondary tracks is likely");
	}
	if(polarization3V.mag2() > 0.0 && 
	   !manager->getPhysics()->isSpinTrackingEnabled()) {
		G4Exception("beam command","Nonzero Polarization",JustWarning,
			"spinTracking=0 in physics command.");
	}
	return true;
}
