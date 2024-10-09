//	BLCMDcosmicraybeam.cc
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

#include "BLAssert.hh"
#include "BLBeam.hh"
#include "BLGroup.hh"
#include "BLParam.hh"

#include <CLHEP/Random/RandGeneral.h>

const G4double UNDEFINED = -3.7e21;

/**	class BLCMDcosmicraybeam -- implement a cosmic-ray muon "beam"
 *
 *	Generates a muon beam that fills the beam box with muons distributed
 *	like cosmic rays. After the simulation is complete, an estimate of
 *	the exposure time is printed (for apparatus at sea level). The
 *	primary vertices are on the "celestial sphere", which should be outside
 *	the apparatus; all muons generated will intersect the beam box (unless
 *	scattered or absorbed or decayed beforehand), so CPU time is not
 *	wasted on tracks missing the apparatus. Define the beam box 
 *	appropriately.
 *
 *	Note that the z axis is vertical, increasing DOWNWARD. The "nominal"
 *	beam muons propagate in the +z direction. All primary vertices will
 *	be in the z<0 half of the celestial sphere, and all muons will have
 *	a direction with positive z component; no beam muon is more than 70
 *	degrees from heading vertically downward.
 *
 *	For now, this is just mu+.
 *
 *	Method:
 *	First, an initial vertex is selected uniformly in the beam box
 *	(beamWidth,beamHeight,beamLength), and a muon is generated with
 *	appropriate momentum and angular distribution. This track is
 *	extended to the center of the beam box, and hits counts the number
 *	of hits inside the beam box at its midplane. Then if radius>0 the
 *	track is extended backwards to where it intersects the
 *	"celestial sphere" -- the sphere with R=radius centered at the
 *	center of the beam box. If radius<=0 then the primary vertices will
 *	be on the midplane of the beam box.
 **/
class BLCMDcosmicraybeam : public BLBeam, public BLCommand {
	G4double meanMomentum;
	G4double beamZ;
	G4double beamWidth;
	G4double beamHeight;
	G4double beamLength;
	G4ParticleGun *muonGun;
	G4ParticleDefinition *muon;
	G4double radius;
	G4long hits;
	G4double sterradians;
	G4double hitsPerM2PerSecPerSterrad;
	G4String particle;
	G4int evNumber;
public:
	/// Constructor.
	BLCMDcosmicraybeam();

	/// commandName() returns "cosmicraybeam".
	G4String commandName() { return "cosmicraybeam"; }
	
	/// command() implements the cosmicraybeam command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// printBeam() will print a description.
	void printBeam();

	/// defineNamedArgs() defines the named arguments for this command.
	void defineNamedArgs();

	/// getNEvents() returns the # events to process.
	int getNEvents() const { return nEvents; }

	/// init() will initialize internal variables.
	void init();

	/// generateReferenceParticle() generates the reference particle.
	bool generateReferenceParticle(G4Event *event);

	/// nextBeamEvent() generates the next beam event.
	bool nextBeamEvent(G4Event *event);

	/// summary() will print a summary.
	void summary();

	/// cosmicRayMuonMomentum() returns a random momentum value distributed
	/// like the muons from cosmic rays. Average momentum is 3.0*GeV/c, 
	/// cutoff is 120*GeV/c. Return value is in geant4 units (MeV).
	/// Data from Kremer et al, Phys. Rev. Lett., Vol 83 no 21, p4241 (1999)
	double cosmicRayMuonMomentum();

	/// cosmicRayMuonAngle() returns a random angle distributed like the 
	/// polar angle of cosmic ray muons.
	/// The return value is the polar angle from vertical (radians). It is
	/// cut off at 70 degrees.
	/// Note this is the distribution for muons of ~3 GeV/c, which is the
	/// average momentum. In fact, lower energy muons have a steeper
	/// distribution and higher ones have a flatter distribution.
	/// But this is a reasonable approximation.
	/// Particle Data Group, Review of Particle Properties, 2002.
	/// Section 23.3.1.
	double cosmicRayMuonAngle();
};

BLCMDcosmicraybeam defineCosmicRayBeam;

BLCMDcosmicraybeam::BLCMDcosmicraybeam() : BLBeam()
{
	registerCommand(BLCMDTYPE_BEAM);
	setSynopsis("Define a Cosmic-Ray muon 'beam'.");
	setDescription("The muon beam is nominally headed in the +Z direction,\n"
		"implying that +Z is physically DOWN.\n"
		"The beam intersects a box defined by beamWidth, beamHeight,\n"
		"and beamLength, centered at X=Y=0 and beamZ. For each\n"
		"event a point is selected randomly within this box, angles\n"
		"theta and phi and the muon momentum are generated according\n"
		"to a fit to their sea-level distributions, the track is\n"
		"extended backwards to the 'celestial sphere', and that is\n"
		"the initial beam position for the event. The muon flux\n"
		"through the rectangle at Z=beamZ is used to display an\n"
		"estimate of the sea-level exposure time for the run.\n\n"
		"This command places itself into the geometry.");

	meanMomentum = 3.0*GeV;
	beamZ = 0.0;
	beamWidth = 0.0;
	beamHeight = 0.0;
	beamLength = 0.0;
	muonGun = 0;
	muon = 0;
	radius = 0.0;
	hits = 0;
	sterradians = 0.0;
	hitsPerM2PerSecPerSterrad = 0.0;
	particle = "mu+";
	evNumber = 0;
}

int BLCMDcosmicraybeam::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 0) {
		printError("Invalid cosmicraybeam command");
		return -1;
	}

	handleNamedArgs(namedArgs);

	// ensure our muons are within the world.
	BLGroup::getWorld()->setMinWidth(radius*2.0);
	BLGroup::getWorld()->setMinHeight(radius*2.0);
	BLGroup::getWorld()->setMinLength(radius*2.0);

	BLManager::getObject()->registerBeam(this);

	printBeam();

	printf("********************************************************\n"
	       "*                                                      *\n"
	       "*  NOTE: This command needs to be completely re-done.  *\n"
	       "*                                                      *\n"
	       "********************************************************\n");

	return 0;
}

void BLCMDcosmicraybeam::defineNamedArgs()
{
	argInt(nEvents,"nEvents","Number of events to process");
	argDouble(beamZ,"beamZ","Beam location in Z (mm)",mm);
	argDouble(radius,"radius","Radius of celestial sphere (mm)",mm);
	argDouble(beamHeight,"beamHeight","Rectangular Beam height (mm)",mm);
	argDouble(beamWidth,"beamWidth","Rectangular Beam width (mm)",mm);
	argDouble(beamLength,"beamLength","Rectangular Beam length (mm)",mm);
}

void BLCMDcosmicraybeam::printBeam()
{
	printf("cosmicraybeam mu+ particle=%s ",particle.c_str());
	printf("nEvents=%d ",nEvents);
	printf("beamZ=%.1f ",beamZ);
	printf("radius=%.3f ",radius);
	printf("\n\t\t");
	printf("beamHeight=%.1f ",beamHeight);
	printf("beamWidth=%.1f ",beamWidth);
	printf("beamLength=%.3f ",beamLength);
	printf("\n");
}

void BLCMDcosmicraybeam::init()
{
	if(muon != 0) return;

	muon = G4ParticleTable::GetParticleTable()->FindParticle(particle);
	muonGun = new G4ParticleGun(1);
	muonGun->SetParticleDefinition(muon);
	evNumber = 0;
}

bool BLCMDcosmicraybeam::generateReferenceParticle(G4Event *event)
{
	return false;
}

bool BLCMDcosmicraybeam::nextBeamEvent(G4Event *event)
{
	for(;;) {
		if(++evNumber > nEvents) return false;
		if(!BLManager::getObject()->skipEvent(evNumber)) break;
	}

	setRandomSeedToGenerate(evNumber);

	G4double mass = muon->GetPDGMass();
	G4ThreeVector position;
	G4ThreeVector direction;
	G4double time = 0.0;
	G4double momentum = 0.0;

	// work in local coordinates (relative to the beam box)
	G4double x = beamWidth*G4UniformRand() - beamWidth/2.0;
	G4double y = beamHeight*G4UniformRand() - beamHeight/2.0;
	G4double z = beamLength*G4UniformRand() - beamLength/2.0;
	momentum = cosmicRayMuonMomentum();
	double theta = cosmicRayMuonAngle();
	double phi = 2.0*pi*G4UniformRand();
	direction[0] = momentum*sin(theta)*cos(phi);
	direction[1] = momentum*sin(theta)*sin(phi);
	direction[2] = momentum*cos(theta);

	// propagate to z=0 (center of the beam box)
	x -= direction[0]/direction[2]*z;
	y -= direction[1]/direction[2]*z;
	z = 0.0;

	// count midplane hits in the beam box (for normalization in summary())
	if(fabs(x) < beamWidth/2.0 && fabs(y) < beamHeight/2.0)
		++hits;

	// project backwards (negative z) onto the celestial sphere
	if(radius > 0.0) {
		G4double xp = direction[0]/direction[2];
		G4double yp = direction[1]/direction[2];
		G4double a = 1.0 + xp*xp + yp*yp;
		G4double b = 2.0 * (x*xp + y*yp);
		G4double c = x*x + y*y - radius*radius;
		z = (-b - sqrt(b*b-4.0*a*c)) / (2.0 * a);
		position[0] = x + xp*z;
		position[1] = y + yp*z;
		position[2] = z + beamZ;
	} else {
		position[0] = x;
		position[1] = y;
		position[2] = beamZ;
	}

	G4double ke = sqrt(momentum*momentum + mass*mass) - mass;
	muonGun->SetParticleTime(time);
	muonGun->SetParticlePosition(position);
	muonGun->SetParticleEnergy(ke);
	muonGun->SetParticleMomentumDirection(direction);
	muonGun->GeneratePrimaryVertex(event);

	setRandomSeedToTrack(evNumber);
	return true;
}

void BLCMDcosmicraybeam::summary()
{
	G4double area = beamWidth*beamHeight/meter2;
	G4double hitsPerSec = hitsPerM2PerSecPerSterrad * area * sterradians;
	printf("CosmicRayBeam: %ld hits in midplane, area %.3f meter^2  "
		" %.3f sterradians\n",
		hits,area,sterradians);
	if(hitsPerSec > 0.0)
		printf("               %.1f hits/sec, estimated exposure"
			" %.1f sec\n",hitsPerSec,hits/hitsPerSec);
}

double BLCMDcosmicraybeam::cosmicRayMuonMomentum()
{
	// Data from Kremer et al, Phys. Rev. Lett., Vol 83 no 21, p4241 (1999).
	// values are lower bin edge, bin average, mu+ rate, mu- rate
	// (laid out this silly way so verification with the paper is easy)
	// NOTE: units are GeV/c, and counts/(GeV/c m^2 sr s)
	static double vals[] = {
		0.0,	0.0,	0.0,	0.0,
		0.2,	0.25,	14.0,	11.0,
		0.3,	0.35,	16.8,	13.6,
		0.4,	0.47,	17.2,	14.4,
		0.55,	0.62,	16.6,	13.5,
		0.70,	0.78,	15.6,	13.3,
		0.85,	0.92,	14.8,	12.1,
		1.0,	1.1,	13.0,	11.0,
		1.2,	1.3,	12.0,	10.1,
		1.4,	1.5,	10.2,	8.7,
		1.6,	1.84,	9.1,	7.3,
		2.1,	2.49,	6.6,	5.2,
		2.94,	3.49,	4.12,	3.38,
		4.12,	4.78,	2.53,	1.98,
		5.5,	6.21,	1.61,	1.25,
		7.0,	8.37,	0.90,	0.69,
		10.0,	12.42,	0.389,	0.309,
		15.5,	18.85,	0.138,	0.108,
		23.0,	26.68,	0.063,	0.046,
		31.1,	36.69,	0.028,	0.019,
		43.6,	51.47,	0.0099,	0.0071,
		61.1,	72.08,	0.0036,	0.0030,
		85.6,	100.96,	0.0014,	0.0012,
		120.0,	120.0,	0.0,	0.0}; // cutoff at 120 GeV/c
	const int nvals = sizeof(vals)/sizeof(vals[0]);
	const int nbins = nvals/4 - 1;
	const int npdf=256;
	static double pdf[npdf];
	static double pmax = vals[4*nbins];
	static bool init=true;

	if(init) {
		// RandGeneral needs equal-sized bins for pdf[]
		// it returns a value in the range [0,1)
		hitsPerM2PerSecPerSterrad = 0.0;
		for(int i=0,ibin=0; i<npdf; ++i) {
			double p = (i+0.5)*pmax/npdf;
			while(p >= vals[4*ibin+5]) ++ibin;
			BLAssert(ibin <= nbins);
			double f = (p - vals[4*ibin+1]) /
						(vals[4*ibin+5]-vals[4*ibin+1]);
			BLAssert(0.0 <= f && f <= 1.0);
			pdf[i] = (1.0-f)*(vals[4*ibin+2]+vals[4*ibin+3]) +
					f*(vals[4*ibin+6]+vals[4*ibin+7]);
			hitsPerM2PerSecPerSterrad += pdf[i] * pmax/npdf;
		}
		init = false;
	}

	CLHEP::RandGeneral generator(pdf,npdf); // BUG in RandGeneral - cannot use new
	return generator.shoot() * pmax * GeV;
}

double BLCMDcosmicraybeam::cosmicRayMuonAngle()
{
	const int npdf=128;
	static double pdf[npdf];
	const double thetamax = 70.0*deg;
	static bool init = true;

	if(init) {
		// RandGeneral needs equal-sized bins for pdf[]
		// it returns a value in the range [0,1)
		sterradians = 0.0;
		G4double dtheta = thetamax / npdf;
		for(int i=0; i<npdf; ++i) {
			// Particle Data Group, Review of Particle Properties,
			// 2002. Section 23.3.1.
			double c = cos(dtheta*i);
			pdf[i] = c*c;
			sterradians += 2.0*pi*c*c*sin(dtheta*i)*dtheta;
		}
		init = false;
	}

	CLHEP::RandGeneral generator(pdf,npdf); // BUG in RandGeneral - cannot use new
	return generator.shoot() * thetamax;
}

