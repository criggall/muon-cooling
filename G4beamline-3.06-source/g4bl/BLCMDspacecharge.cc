//	BLCMDspacecharge.cc
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

#ifdef G4BL_FFTW

#include <stdio.h>

#include "G4UImanager.hh"

#include "BLCommand.hh"
#include "BLAssert.hh"
#include "BLCollectiveComputation.hh"
#include "BLRunManager.hh"
#include "BLManager.hh"
#include "BLGlobalField.hh"
#include "BLNTuple.hh"
#include "BLParam.hh"
#include "BLElementField.hh"
#include "BLTime.hh"

#include "PoissonConvolve3D.hh"

const G4double UNSET=-1.0e99;

// parameters to control grid sizing (fractions of grid size):
const G4double MINGRID=0.5;		// percentile is kept > this
const G4double MAXGRID=0.667;		// percentile is kept < this
const G4double SUPPRESSGRID=0.85;	// charge > this is suppressed

/**	class Bunch represents a bunch.
 *	It maintains the trajectory for a reference particle and computes the
 *	field for all particles in the bunch related to that reference.
 *
 *	A Bunch is created in PreUserTrackingAction() when tracking a
 *	reference particle. The tracks in the bunch are added in
 *	beginCollectiveTracking(). The fields from the Bunch are computed 
 *	in computeField(), which is called by collectiveStep(), and they
 *	are then used for tracking via addFieldValue().
 *
 *	The grid for PoissonConvolve3D is always centered at the reference
 *	particle, and moves with it between time steps. No solution is
 *	possible if any charge is at or very close to the boundary, so the
 *	specified percentile of the charge distribution is used to determine
 *	the size of the grid, keeping the percentile of the distribution
 *	between MINGRID and MAXGRID of the grid size. Any charge located 
 *	>SUPPRESSGRID of the grid size is suppressed (in x, y, or z).
 *
 *	The challenge is to provide suitably accurate boundary conditions
 *	without taking too much CPU. For a perfectly Gaussian bunch that can
 *	be done by applying an approximation grid 4 or 8 times smaller than 
 *	the Poisson3D grid. For a cylindrical bunch that does not work.
 *	So for now it is left up to the user to specify the approximation
 *	grid size. Unlike the Poisson3D grid, the approximation grid can be
 *	any size. An attempt is made to merge cells with tiny charge (<1% of 
 *	the total) with its neighbor, but this does not always reduce the
 *	number of approximation charges very much.
 *
 *	This version requires the reference particle's momentum to be parallel
 *	to the z axis.
 **/
class Bunch {
	static const double extrapolateFrac; // fraction of deltaT to extrap.
	static const double minDeltaT; // min dt between reference steps
	struct Point {	// one point in the reference trajectory
		G4ThreeVector position;
		G4ThreeVector momentum;
		double time;
		Point(double t, G4ThreeVector p, G4ThreeVector m)
			{ time=t; position=p; momentum=m; }
	};
	/// class Approx does a 3-d histogram of the charge to generate
	/// an approximation for PoissonConvolve3D. Within each bin (voxel) it
	/// use the mean position and the total charge.
	class Approx {
		int nx, ny, nz;
		double totalC;
		double dx, dy, dz;
		double xMin, yMin, zMin;
		double *xx, *yy, *zz, *cc;
		int index(int ix, int iy, int iz) 
			{ return ix+nx*(iy+ny*iz); }
		std::vector<Charge> v;
	public:
		Approx(int _nx, int _ny, int _nz) : v() 
		{	nx=_nx; ny=_ny; nz=_nz;
			xx = new double[nx*ny*nz];
			yy = new double[nx*ny*nz];
			zz = new double[nx*ny*nz];
			cc = new double[nx*ny*nz];
		}
		void reset(Bunch *b) {
			totalC = 0.0;
			dx = b->bound_x*2.0/nx; xMin = -b->bound_x;
			dy = b->bound_y*2.0/ny; yMin = -b->bound_y;
			dz = b->bound_z*2.0/nz; zMin = -b->bound_z;
			int nBins3 = nx*ny*nz;
			for(int i=0; i<nBins3; ++i)
				xx[i] = yy[i] = zz[i] = cc[i] = 0.0;
			v.clear();
		}
		void putCharge(double x, double y, double z, double c);
		void computeApprox();
		void reduce();
		std::vector<Charge> getVector() const { return v; }
	};
	friend class Approx;
	std::vector<Point> point;
	unsigned index;
	double prevTime;
	G4ParticleDefinition *definition;
	double charge;				// charge of each macro-particle
	double maxBeta;				// max beta in bunch frame
	double bound_x, bound_y, bound_z;	// max distance particle to ref
						// (lab frame, global coords)
	int nParticles;				// # particles in the bunch
	double percentile;			// %tile used for grid sizing
	int minActive;				// Min # particles in bunch
	PoissonConvolve3D *poisson;
	bool goodGuess;
	bool fixedGrid;
	Approx approx;
	std::vector<G4Track*> tracks;		// tracks in the bunch
	//
	BLRunManager *runManager;
	static int verbose;
	static int useApproximationOnly;
public:
	Bunch(G4ParticleDefinition *def, int nx, int ny, int nz,
		double dx, double dy, double dz, double _charge, 
		double _maxBeta, bool _fixedGrid, double _percentile,
		int nxApprox, int nyApprox, int nzApprox, int _minActive) :
		point(), approx(nxApprox,nyApprox,nzApprox), tracks()
		{ index=0; prevTime=-DBL_MAX; definition=def; 
		  charge=_charge; maxBeta=_maxBeta; fixedGrid=_fixedGrid;
		  bound_x=dx; bound_y=dy; bound_z=dz;
		  nParticles=0; percentile=_percentile; minActive=_minActive;
		  poisson = new PoissonConvolve3D(nx,ny,nz,-dx,dx,-dy,dy,-dz,dz);
		  goodGuess=false; runManager=BLRunManager::getObject();
		  approx.reset(this);
		}

	/// add a point to the reference trajectory.
	void addReferencePoint(const G4Track *track) {
		double time=track->GetGlobalTime();
		if(time >= prevTime+minDeltaT) { // valid for prevTime=-DBL_MAX
			point.push_back(Point(time, track->GetPosition(),
				track->GetMomentum()));
			prevTime = time;
		}
	}

	/// add a point to the reference trajectory.
	void addReferencePoint(double t, G4ThreeVector pos, G4ThreeVector mom)
		{ point.push_back(Point(t,pos,mom)); }

	/// nTrajectoryPoints() returns the # points in the reference trajectory
	int nTrajectoryPoints() { return point.size(); }

	/// getNParticles() returns # particles in bunch
	int getNParticles() const { return nParticles; }

	/// getParticleDefinition()
	G4ParticleDefinition *getDefinition() { return definition; }

	/// Add a track to the Bunch.
	/// returns false if this track is not part of this bunch.
	bool addTrack(G4Track *t) {
		if(!isInBunch(t)) return false;
		tracks.push_back(t);
		return true;
	}

	int size() { return tracks.size(); }
	int getNapprox() { return approx.getVector().size(); }
	int getMinActive() { return minActive; }
	G4ThreeVector getBound()
		{ return G4ThreeVector(bound_x,bound_y,bound_z); }

	/// Interpolate in the point[] vector to get position and 4-momentum of
	/// the reference for this bunch.
	/// Returns true if OK, false if time is after the valid range
	/// for this reference. If time is before the reference was created,
	/// its first two trajectory points are linearly extrapolated backwards.
	/// Because it must search for the enclosing points of the trajectory,
	/// this is most efficient when time does not jump very much (it starts
	/// the search at the place in the trajectory of the previous call).
	/// (Note: the const is a lie, but is required by BLElementField; only
	///  index is modified, which is essentially benign.)
	bool getReferencePos4Mom(double time, G4ThreeVector &position,
					G4LorentzVector &fourMomentum) const;

	/// init() initializes the computation, sizing the grid for the
	/// initial bunch, and computing the field for the initial time step.
	/// Called after all tracks are added to the bunch, before stepping.
	/// Returns false if failure.
	bool init();

	/// isInBunch() returns true if this track is in the bunch.
	/// Tracks are selected by these criteria:
	///     same particle as reference
	///     within bound_x,bound_y,bound_z of the bunch's reference traj.
	///     v/c in the bunch's beam frame less than maxBeta
	bool isInBunch(G4Track *track);

	/// All tracks in the bunch are added into the PoissonConvolve3D
	/// and the field is computed. returns false if current step time
	/// is beyond the range of the bunch's reference.
	bool computeField();

	/// resizes the grid, based on current tracks.
	/// returns true if grid was resized.
	bool resizeGrid(bool force=false);

	/// addFieldValue adds the field for this bunch to field[].
	void addFieldValue(const G4double point[4], G4double field[6]) const;

	/// testEfield() will test the units and value of the E field in the
	/// beam frame for a point charge.
	static void testEfield();

	/// testBfield() will test the B field in the lab frame for a bunch
	/// that is a line of moving charges, 1000 Amps.
	static void testBfield();

	/// setVerbose() sets the debug print level
	static void setVerbose(int v) { verbose=v; }

	/// setUseApproximationOnly() sets it
	static void setUseApproximationOnly(int v) { useApproximationOnly=v; }
};
const double Bunch::minDeltaT=1.0e-6*ns;
const double Bunch::extrapolateFrac=0.01;
int Bunch::verbose=0;
int Bunch::useApproximationOnly=0;

/**	class BLCMDspacecharge - space-charge computation
 *	This computation works within the Geant4 tracking framework, and does
 *	the computation in the various callbacks. The algorithm tracks
 *	individual particles, and computes their field using macro-particles
 *	that are co-located with the tracked particles. So it computes E and B
 *	fields from the macro-particles and adds this in to the BLGlobalField,
 *	which Geant4 uses for tracking. 
 *
 *	This computation uses a Poisson solver in the beam frame. The beam
 *	frame is determined from the reference particle's position and
 *	4-momentum.
 *
 *	It is best to terminate the simulation with "trackcuts maxTime=...",
 *	or minActive=.... If it is terminated because particles leave the
 *	world or hit an object with kill=1, when the # particles in the bunch
 *	gets very low PoissonConvolve3D may not converge, aborting the run.
 *	Fortunately the NTuples wre written, so this is OK, it just looks ugly.
 *
 *	Note that fewer than minActive particles in the bunch is treated as
 *	a NORMAL way to end the simulation, not an error. But messages
 *	beginning with "***" are printed to alert the user that this happened.
 *
 *	This version requires the reference particle's momentum to be parallel
 *	to the z axis.
 **/
class BLCMDspacecharge : public BLCommand, public BLCollectiveComputation,
		public BLManager::SteppingAction, 
		public BLManager::TrackingAction,
		public BLElementField {
	G4double deltaT;		// time step (ns)
	G4double charge;		// charge of Macroparticle-s
	G4int nx;
	G4int ny;
	G4int nz;
	G4double dx;
	G4double dy;
	G4double dz;
	G4double maxBeta;
	G4int verbose;			// verbosity conrol for debugging
	G4int ignoreFieldWhenTracking;
	G4int useApproximationOnly;
	G4int fixedGrid;
	double percentile;
	int nxApprox;
	int nyApprox;
	int nzApprox;
	int minActive;
	// misc. data:
	BLManager *manager;
	BLRunManager *runManager;
	int nstep;
	std::vector<Bunch> bunch;
	Bunch *currentBunch;
	bool tracking;
public:
	BLCMDspacecharge();

	G4String commandName() { return "spacecharge"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	virtual void beginCollectiveTracking(std::vector<BLTrackData>& v);

	virtual void collectiveStep(std::vector<BLTrackData>& v);

	virtual void endCollectiveTracking(std::vector<BLTrackData>& v);

	/// TrackingAction and SteppingAction are used only to get the
	/// points of each reference trajectory into bunch[].
	/// PreUserTrackingAction() creates each Bunch object, 
	/// seting currentBunch which is used by the others.
	virtual void PreUserTrackingAction(const G4Track *track);
        virtual void PostUserTrackingAction(const G4Track *track);
	virtual void UserSteppingAction(const G4Step *step);

	virtual void addFieldValue(const G4double point[4], G4double field[6])
		const;
};
BLCMDspacecharge defaultBLCMDspacecharge;

BLCMDspacecharge::BLCMDspacecharge() : BLCollectiveComputation(),
		BLManager::SteppingAction(), BLManager::TrackingAction(),
		BLElementField(), bunch()
{
	registerCommand(BLCMDTYPE_PHYSICS);
	setSynopsis("Beam-frame Green's function space charge computation");
	setDescription("This is a space charge computation for bunched beams. "
		"It uses a grid in the beam frame to solve Poisson's equation "
		"via a Green's function with infinite boundary conditions; "
		"the E field is boosted back to the lab frame E and B for "
		"tracking.\n\n"
		"Macro-particles are used to enable the simulation of larger "
		"bunches than can be feasibly simulated as individual "
		"particles; the macro-particles have zero radius, "
		"but are pro-rated into the nearest eight grid points when "
		"placed into the grid. This computation can handle up to about "
		"a million macro-particles, but 100,000 is more sensible for "
		"all but the simplest physical situations.\n\n"
		"The bunch is created from the beam tracks before tracking "
		"begins. There is one bunch for each reference particle. "
		"Particles in the bunch must be the same particle as the "
		"reference, must initially be within {dx,dy,dz} of the "
		"reference particle, and when boosted to the reference "
		"particle's rest frame must initially have beta < maxBeta.\n\n"
		"After boosting the particles to the beam frame, they are "
		"placed into the grid, pro-rating to the eight nearest grid "
		"points. The grid is dynamically re-sized to keep the "
		"99th percentile of the particles between 0.5 and 0.67 of "
		"the grid size. This maintains a reasonable balance between "
		"resolution of grid points within the bunch and covering all "
		"of the particles. Particles located at >85% of the grid size "
		"do not contribute to the field computation, but are tracked "
		"using the field of the rest of the bunch (and other bunches)."
		"\n\n"
		"The grid has {nx,ny,nz} points; there is a small "
		"computational advantage to using powers of 2, but any values "
		">1 can be used. For efficiency, the convolution of the "
		"Green's function with the charge grid is performed using "
		"FFTs; the grid is doubled in each dimension with the proper "
		"symmetry applied to the Green's function, so the cyclical "
		"convolution of the FFTs gives the proper potential with "
		"infinite boundary conditions.\n\n"
		"Outside the grid an approximation is used. An approximation "
		"grid is constructed, with the same size of the Poisson grid, "
		"but using {nxApprox,nyApprox,nzApprox} points. Particles are "
		"placed into this approximation grid, and the mean position "
		"is kept as well as the charge. Approximation grid points "
		"with less than 1% of the total charge are "
		"consolidated with their inner neighbors. The non-zero "
		"approximation grid points are treated as point charges when "
		"computing the potential outside the grid. For reasonably "
		"Gaussian bunches, {7,7,7} are reasonable values for the "
		"approximation grid sizes.\n\n"
		"The E field in the beam frame is computed via the derivatives "
		"of the linear interpolating function using the eight nearest "
		"grid points, and is boosted back to the lab "
		"frame E and B for tracking by the usual Geant4 routines.\n\n"
		"Bunch particles that get destroyed cease contributing to the "
		"bunch. As the bunch particles are selected during start-up, "
		"no additional particles are ever added to a bunch. This "
		"algorithm handles multiple bunches of any particle types.\n\n"
		"NOTE: For now, the reference MUST be parallel to the Z "
		"axis.");

	// initialize class variables here
	deltaT = -1.0;
	charge = -1.0;
	nx = 65;
	ny = 65;
	nz = 65;
	dx = -1.0;
	dy = -1.0; 
	dz = -1.0;
	maxBeta = 0.1;
	verbose = -1;
	ignoreFieldWhenTracking = 0;
	useApproximationOnly = 0;
	fixedGrid = 0;
	percentile = 99;
	nxApprox = 7;
	nyApprox = 7;
	nzApprox = 7;
	minActive = -95;
	manager = 0;	// not available yet -- set in command(...)
	runManager = 0;	// not available yet -- set in command(...)
	nstep = 0;
	currentBunch = 0;
	tracking = false;
}

int BLCMDspacecharge::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(manager != 0) {
		printError("spacecharge - can be used at most once.");
		return -1;
	}

	handleNamedArgs(namedArgs);

	manager = BLManager::getObject();
	runManager = BLRunManager::getObject();

	if(verbose < 0) verbose = Param.getInt("steppingVerbose");
	if(verbose < 0) verbose = 0;
	Bunch::setVerbose(verbose);
	Bunch::setUseApproximationOnly(useApproximationOnly);

	if(nxApprox < 0) nxApprox = nx;
	if(nyApprox < 0) nyApprox = ny;
	if(nzApprox < 0) nzApprox = nz;

	if(charge <= 0.0 || deltaT <= 0.0)
		printError("spacecharge: charge and deltaT must both be > 0");
	if(dx <= 0.0 || dy <= 0.0 || dz <= 0.0)
		printError("spacecharge: dx, dy, and dz must all be > 0");

	runManager->registerCollectiveComputation(this);
	runManager->setCollectiveMode(true);
	runManager->setDeltaT(deltaT);

	manager->registerSteppingAction(this);
	manager->registerTrackingAction(this);

	BLGlobalField::getObject()->addElementField(this); // infinite bounding
							   // box 
	print("");

	Bunch::testEfield();
	Bunch::testBfield();

	return 0;
}

void BLCMDspacecharge::defineNamedArgs()
{
	argDouble(deltaT,"deltaT","Time step (ns).",ns);
	argDouble(charge,"charge","Charge of macro particles (times particle charge).");
	argInt(nx,"nx","Number of grid points in x (65).");
	argInt(ny,"ny","Number of grid points in y (65).");
	argInt(nz,"nz","Number of grid points in z (65).");
	argDouble(dx,"dx","Max distance of particle to reference in x (mm)",mm);
	argDouble(dy,"dy","Max distance of particle to reference in y (mm)",mm);
	argDouble(dz,"dz","Max distance of particle to reference in z (mm)",mm);
	argInt(nxApprox,"nxApprox","# bins in x in approximation (7).");
	argInt(nyApprox,"nyApprox","# bins in y in approximation (7).");
	argInt(nzApprox,"nzApprox","# bins in z in approximation (7).");
	argDouble(maxBeta,"maxBeta","Max beta (v/c) of particle in beam frame (0.1).");
	argInt(verbose,"verbose","Non-zero for verbose prints (0).");
	argInt(ignoreFieldWhenTracking,"ignoreFieldWhenTracking","For testing only (0).");
	argInt(useApproximationOnly,"useApproximationOnly","For testing only (0).");
	argInt(fixedGrid,"fixedGrid","Nonzero prevents re-sizing the grid (0).");
	argDouble(percentile,"percentile","Percentile of charge distribution used for grid sizing (99).");
	argInt(minActive,"minActive","Minimum # active tracks in bunch; if "
	     "< 0 is % of initial bunch size (-95).");
}	

void BLCMDspacecharge::beginCollectiveTracking(std::vector<BLTrackData>& v)
{
	if(verbose >= 1) printf("BLCMDspacecharge::beginCollectiveTracking v.size=%d\n",(int)v.size());

	runManager->setDeltaT(deltaT);
	nstep = 0;

	// add tracks to Bunch-es
	std::map<G4ParticleDefinition*,int> orphans;
	int nactive=0, norphan=0;
	for(unsigned j=0; j<v.size(); ++j) {
		G4Track *track = v[j].track;
		G4TrackStatus trackStatus = track->GetTrackStatus();
		if(trackStatus != fAlive && trackStatus != fStopButAlive)
			continue;
		++nactive;
		bool match=false;
		for(unsigned i=0; i<bunch.size(); ++i) {
			if(bunch[i].addTrack(track))
				goto found;
			if(track->GetDefinition() == bunch[i].getDefinition())
				match = true;
		}
		++orphans[track->GetDefinition()];
		if(match) ++norphan;
found:		continue;
	}


	for(unsigned i=0; i<bunch.size(); ++i) {
		BLAssert(bunch[i].init());
		printf("Bunch %d:\n",i);
		printf("        Bunch has %d tracks, minActive=%d\n",
				    bunch[i].size(),bunch[i].getMinActive());
		G4ThreeVector b=bunch[i].getBound();
		printf("        Grid has nx,ny,nx=%d,%d,%d\n",nx,ny,nz);
		printf("        Initial Grid size: %.3f,%.3f,%.3f\n",
				    b.x(),b.y(),b.z());
		printf("        Initial Approximation has nx,ny,nx=%d,%d,%d "
				    " %d charges\n",nxApprox,nyApprox,nzApprox,
				    bunch[i].getNapprox());
	}

	if(orphans.size() > 0) {
		printf("spacecharge: orphan tracks:");
		std::map<G4ParticleDefinition*,int>::iterator i;
		char sep=' ';
		for(i=orphans.begin(); i!=orphans.end(); ++i) {
			printf("%c%d*%s",sep,i->second,
					i->first->GetParticleName().c_str());
			sep = ',';
		}
		printf("\n");
	}
	if(norphan > nactive/1000)
		G4Exception("spacecharge","Too many orphan tracks",JustWarning,
			"More than 0.1% of matching tracks not in any Bunch");
}

void BLCMDspacecharge::collectiveStep(std::vector<BLTrackData>& v)
{
	for(unsigned i=0; i<bunch.size(); ++i) {
		if(!bunch[i].computeField()) {
			printf("***Ending computation\n");
			runManager->AbortRun();
			break;
		}
	}

	++nstep;
}

void BLCMDspacecharge::endCollectiveTracking(std::vector<BLTrackData>& v)
{
	printf("BLCMDspacecharge::endCollectiveTracking after %d time steps, v.size=%d\n",
		nstep,(int)v.size());
}

void BLCMDspacecharge::PreUserTrackingAction(const G4Track *track)
{
	tracking = true;
	currentBunch = 0;
	if(manager->getState() == REFERENCE && 
	   track->GetDefinition()->GetPDGCharge() != 0.0) {
	      bunch.push_back(Bunch(track->GetDefinition(),nx,ny,nz,dx,dy,dz,
  				charge,maxBeta,fixedGrid,percentile,
				nxApprox,nyApprox,nzApprox,minActive));
	      currentBunch = &bunch[bunch.size()-1];
	      if(verbose >= 1) printf("spacecharge: starting Reference\n");
	}
}

void BLCMDspacecharge::PostUserTrackingAction(const G4Track *track)
{
	tracking = false;
	if(currentBunch == 0) return;

	if(currentBunch->nTrajectoryPoints() < 2)
		G4Exception("BLCMDspacecharge","Too few reference steps",
			FatalException, "Reference particle took < 2 steps");
	currentBunch = 0;

	if(verbose >= 1) printf("spacecharge: ending Reference\n");
}

void BLCMDspacecharge::UserSteppingAction(const G4Step *step)
{
	if(currentBunch == 0) return;

	G4Track *track = step->GetTrack();

	/// require the reference be parallel to the z axis
	G4ThreeVector direction = track->GetMomentumDirection();
	if(direction.perp() > 1.0e-9)
		G4Exception("BLCMDspacecharge","Invalid Bunch momentum",
			FatalException, "Reference not parallel to z axis.");

	currentBunch->addReferencePoint(track);
}

void BLCMDspacecharge::addFieldValue(const G4double point[4], G4double field[6]) const
{
	if(manager->getState() != BEAM) return;
	if(ignoreFieldWhenTracking != 0 && tracking) return;

	for(unsigned i=0; i<bunch.size(); ++i) {
		bunch[i].addFieldValue(point,field);
	}
}

bool Bunch::getReferencePos4Mom(double time,
		G4ThreeVector &position, G4LorentzVector &fourMomentum) const
{
	// find the current entry in point[] -- return false if beyond the end.
	while(time < point[index].time && index > 0)
		--((Bunch*)this)->index; // discard const
	while(index+1 < point.size()-1 && time > point[index+1].time) {
		++((Bunch*)this)->index; // discard const
	}
	if(index+1 >= point.size()) return false;
	if(!(index+1 < point.size()))
		printf("index=%d point.size=%ld time=%.3f last_time=%.3f\n",
			index,(long)point.size(),time,point[point.size()-1].time);
	BLAssert(index+1 < point.size());

	// Because minDeltaT was imposed while filling point[], we know the
	// points are distinct.
	double dt = point[index+1].time - point[index].time;
	BLAssert(dt >= minDeltaT);

	// interpolate linearly. Works correctly even if time < point[0].time.
	double f=(time-point[index].time)/dt;
	position = (1.0-f)*point[index].position + f*point[index+1].position;
	G4ThreeVector momentum = (1.0-f)*point[index].momentum +
						f*point[index+1].momentum;
	fourMomentum.setVectM(momentum,definition->GetPDGMass());

	if(verbose >= 3) printf("Bunch::getReferencePos4Mom(%.4f) "
			"pos=%.3f,%.3f,%.3f  4-mom=%.3f,%.3f,%.3f,%.3f\n",
			time,position[0],position[1],position[2],
			fourMomentum[0],fourMomentum[1],fourMomentum[2],
			fourMomentum[3]);

	return true;
}

bool Bunch::init()
{
	if(minActive < 0)
		minActive = (tracks.size() * (-minActive))/100;

	resizeGrid(true);

	// get position and gamma of reference
	G4double time=runManager->getStepTime();
	G4ThreeVector refPos;
	G4LorentzVector ref4mom;
	if(!getReferencePos4Mom(time,refPos,ref4mom)) {
	    if(verbose >= 3) printf("Bunch::init returning false\n");
	    return false;
	}
	G4double gamma=ref4mom.gamma();

	poisson->updatePosition(-bound_x,bound_x,-bound_y,
					bound_y,-gamma*bound_z,gamma*bound_z);

	bool retval=computeField();

	return retval;
}

bool Bunch::isInBunch(G4Track *track)
{
	if(verbose >= 2) {
		G4ThreeVector pos=track->GetPosition();
		G4ThreeVector mom=track->GetMomentum();
		printf("isInBunch? TrackID=%d pos=%.3f,%.3f,%.3f 3-momentum="
			"%.3f,%.3f,%.3f\n", track->GetTrackID(),
			pos[0],pos[1],pos[2],mom[0],mom[1],mom[2]);
	}

	// track must be alive
	G4TrackStatus trackStatus = track->GetTrackStatus();
	if(trackStatus != fAlive && trackStatus != fStopButAlive) {
		if(verbose >= 3) printf("Not in bunch, not alive\n");
		return false;
	}

	///     same particle as reference
	if(track->GetDefinition() != definition) {
		if(verbose >= 3) printf("Not in bunch, wrong particle\n");
		return false;
	}

	///     within bound_x,bound_y,bound_z of the bunch's reference traj.
	G4ThreeVector refPos;
	G4LorentzVector ref4mom;
	if(!getReferencePos4Mom(track->GetGlobalTime(),refPos,ref4mom)) {
		if(verbose >= 3) printf("Not in bunch, time too large\n");
		return false;
	}
	G4ThreeVector p = track->GetPosition() - refPos;
	if(fabs(p.x()) > bound_x || fabs(p.y()) > bound_y || 
							fabs(p.z()) > bound_z) {
		if(verbose >= 3) printf("Not in bunch, too far away\n");
		return false;
	}

	///     v/c in the bunch's beam frame less than maxBeta
	G4ThreeVector boost=ref4mom.findBoostToCM();
	G4LorentzVector v;
	v.setVectM(track->GetMomentum(),track->GetDefinition()->GetPDGMass());
	if(v.boost(boost).beta() > maxBeta) {
		if(verbose >= 3)
			printf("Not in bunch, beta=%.4f > %.4f=maxBeta\n",
					v.boost(boost).beta(),maxBeta);
		return false;
	}

	if(verbose >= 3) printf("Is in bunch\n");
	return true;
}

bool Bunch::computeField()
{
	if(verbose >= 3) printf("Bunch::computeField entered\n");

	// get position and gamma of reference
	G4double time=runManager->getStepTime();
	G4ThreeVector refPos;
	G4LorentzVector ref4mom;
	if(!getReferencePos4Mom(time,refPos,ref4mom)) {
	    printf("***Bunch: Invalid Reference\n");
	    return false;
	}
	G4double gamma=ref4mom.gamma();
	if(verbose >= 3) 
	    printf("Bunch::computeField t=%.3f refPos=%.3f,%.3f,%.3f gamma=%.4f\n",
	    	time,refPos[0],refPos[1],refPos[2],gamma);


	// resize the grid, if necessary
	if(fixedGrid == 0) {
		if(resizeGrid()) {
			poisson->updatePosition(-bound_x,bound_x,-bound_y,
					bound_y,-gamma*bound_z,gamma*bound_z);
			goodGuess = false;
		}
	}

	poisson->zeroRhs();
	nParticles = 0;
	approx.reset(this);

	int nSuppress=0;
	for(unsigned i=0; i<tracks.size(); ++i) {
		G4Track *track = tracks[i];
		if(!track) continue;
		G4TrackStatus trackStatus = track->GetTrackStatus();
		if(trackStatus != fAlive && trackStatus != fStopButAlive) {
			tracks[i] = 0; // erase is expensive for large # tracks
			continue;
		}

		G4ThreeVector p=track->GetPosition();
		double dt=track->GetGlobalTime() - runManager->getStepTime();
		G4double dtExtrapolate=runManager->getDeltaT()*extrapolateFrac;
		if(fabs(dt) > dtExtrapolate) {
			G4ThreeVector tmp(dt * track->GetVelocity());
			p -= tmp;
		}
		p -= refPos;
		if(fabs(p.x()) > bound_x*SUPPRESSGRID ||
		   fabs(p.y()) > bound_y*SUPPRESSGRID ||
		   fabs(p.z()) > bound_z*SUPPRESSGRID) {
		    	if(verbose >= 3) printf("spacecharge: suppressed "
					"p=%.3f,%.3f,%.3f\n",p[0],p[1],p[2]);
			++nSuppress;
			continue;
		}
		double c=charge*definition->GetPDGCharge();
		if(!poisson->putCharge(p.x(),p.y(),gamma*p.z(),c/epsilon0)) {
			printf("***Not in Poisson grid; p=%.3f,%.3f,%.3f\n",
		    					p.x(),p.y(),p.z());
			continue;
		}
		approx.putCharge(p.x(),p.y(),gamma*p.z(),c/epsilon0);
		++nParticles;
	}

	if(nParticles < minActive) {
		printf("***Bunch: Too few tracks; have %d, require %d\n",
					nParticles,minActive);
		return false;
	}

	if(nSuppress > 0) {
		if(verbose >= 2) 
			printf("Bunch::computeField suppressed %d tracks\n",
								nSuppress);
	}

	// Solve Poisson's equation in the beam frame.
	approx.computeApprox();
	poisson->setApproximation(approx.getVector());
	if(useApproximationOnly == 0) {
		if(verbose >= 3) printf("Bunch::computeField calling solve\n");
		poisson->solve();
	}

	goodGuess = true;

	if(verbose >= 3) printf("Bunch::computeField returning true\n");
	return true;
}

bool Bunch::resizeGrid(bool force)
{
	if(fixedGrid != 0) return false;

	// get position of reference
	G4double time=runManager->getStepTime();
	G4ThreeVector refPos;
	G4LorentzVector ref4mom;
	getReferencePos4Mom(time,refPos,ref4mom);

	// histogram x,y,z for tracks in Bunch
	// (histograms extend to twice the bounds.)
	const int NHIST=200;
	int histX[NHIST], histY[NHIST], histZ[NHIST];
	for(int i=0; i<NHIST; ++i)
		histX[i] = histY[i] = histZ[i] = 0;
	double dx=2.0*bound_x/NHIST;
	double dy=2.0*bound_y/NHIST;
	double dz=2.0*bound_z/NHIST;
	int n=0;
	for(unsigned i=0; i<tracks.size(); ++i) {
		G4Track *track = tracks[i];
		if(track == 0) continue;
		G4TrackStatus trackStatus = track->GetTrackStatus();
		if(trackStatus != fAlive && trackStatus != fStopButAlive) {
			tracks[i] = 0; // erase is expensive for large # tracks
			continue;
		}
		G4ThreeVector p=track->GetPosition();
		double dt=track->GetGlobalTime() - runManager->getStepTime();
		G4double dtExtrapolate=runManager->getDeltaT()*extrapolateFrac;
		if(fabs(dt) > dtExtrapolate) {
			G4ThreeVector tmp(dt * track->GetVelocity());
			p -= tmp;
		}
		p -= refPos;
		++n;
		int ix=(int)(fabs(p.x())/dx);
		if(ix > NHIST-1) ix=NHIST-1;
		++histX[ix];
		int iy=(int)(fabs(p.y())/dy);
		if(iy > NHIST-1) iy=NHIST-1;
		++histY[iy];
		int iz=(int)(fabs(p.z())/dz);
		if(iz > NHIST-1) iz=NHIST-1;
		++histZ[iz];
	}
	if(n == 0) return false;

	int need=(int)(percentile*n/100.0);
	double avg = (MINGRID+MAXGRID)/2.0;

	// find percentile in x, and set bound_x so it's in [MINGRID,MAXGRID]
	int t=0, ix=0, iy=0, iz=0;
	for(ix=0; ix<NHIST; ++ix) {
		t += histX[ix];
		if(t >= need) break;
	};
	double new_x = ix*dx;
	if(new_x < MINGRID*bound_x || new_x > MAXGRID*bound_x) {
		bound_x = new_x/avg;
		force = true;
	}

	// find percentile in y, and set bound_y so it's in [MINGRID,MAXGRID]
	t = 0;
	for(iy=0; iy<NHIST; ++iy) {
		t += histY[iy];
		if(t >= need) break;
	};
	double new_y = iy*dy;
	if(new_y < MINGRID*bound_y || new_y > MAXGRID*bound_y) {
		bound_y = new_y/avg;
		force = true;
	}

	// find percentile in z, and set bound_z so it's in [MINGRID,MAXGRID]
	t = 0;
	for(iz=0; iz<NHIST; ++iz) {
		t += histZ[iz];
		if(t >= need) break;
	};
	double new_z = iz*dz;
	if(new_z < MINGRID*bound_z || new_z > MAXGRID*bound_z) {
		bound_z = new_z/avg;
		force = true;
	}

	if(!force) return false; 

	if(verbose >= 2) printf("Bunch new grid size:%.3f,%.3f,%.3f\n",
		    				bound_x,bound_y,bound_z);
	return true;
}

void Bunch::addFieldValue(const G4double point[4], G4double field[6]) const
{
	// return if no particles in bunch
	if(nParticles == 0) return;

	// check extrapolation of reference position is OK
	G4double stepTime = runManager->getStepTime();
	G4double deltaT = runManager->getDeltaT();

	// get position of reference
	G4ThreeVector pos;
	G4LorentzVector ref4mom;
	if(!getReferencePos4Mom(point[3],pos,ref4mom))
		return;
	G4ThreeVector boost=ref4mom.findBoostToCM();
	double gamma=ref4mom.gamma();
	if(verbose >= 3) printf("Bunch::addFieldValue: t=%.3f refPos=%.3f,%.3f,%.3f point=%.3f,%.3f,%.3f\n",point[3],pos[0],pos[1],pos[2],point[0],point[1],point[2]);

	G4ThreeVector Ebeam;
	if(useApproximationOnly)
		Ebeam = poisson->approxE(point[0]-pos[0],
				point[1]-pos[1],gamma*(point[2]-pos[2]));
	else
		Ebeam = poisson->getE(point[0]-pos[0],
				point[1]-pos[1],gamma*(point[2]-pos[2]));

	// boost of EM fields from Jefimenko, AJP _64_(5), p618 (1996), 
	// eq12-17 with xyz->zxy and B'=0; boost[] is beta[]=v[]/c_light;
	// There is an additional minus sign in boost[] from findBoostToCM().
	// Reference is parallel to the z axis.
	double f[6];
	f[0] = gamma*boost[2]*Ebeam[1]/c_light;		// Bx
	f[1] = -gamma*boost[2]*Ebeam[0]/c_light;	// By
	f[2] = 0.0;					// Bz
	f[3] = gamma*Ebeam[0];		// Ex
	f[4] = gamma*Ebeam[1];		// Ey
	f[5] = Ebeam[2];		// Ez

	if(verbose >= 3) 
		printf("Bunch::addFieldValue(%.3f,%.3f,%.3f,%.3f)\n"
			"    pos=%.3f,%.3f,%.3f boost=%.3f,%.3f,%.3f "
			"gamma=%.3f  Ebeam=%.3f,%.3f,%.3f\n"
			"    B=%.3f,%.3f,%.3f  E=%.3f,%.3f,%.3f\n",
			point[0]/mm,point[1]/mm,point[2]/mm,point[3]/ns,
			pos[0]/mm,pos[1]/mm,pos[2]/mm,
			boost[0],boost[1],boost[2],gamma,
			Ebeam[0]/(megavolt/meter),
			Ebeam[1]/(megavolt/meter),
			Ebeam[2]/(megavolt/meter),
			f[0]/tesla,f[1]/tesla,f[2]/tesla,
			f[3]/(megavolt/meter),f[4]/(megavolt/meter),
			f[5]/(megavolt/meter));

	field[0] += f[0];
	field[1] += f[1];
	field[2] += f[2];
	field[3] += f[3];
	field[4] += f[4];
	field[5] += f[5];
}

void Bunch::testEfield()
{
	/* UNITS TEST
	   E = k*Q/r/r,  k=8.987551787E9 (Newton*meter*meter/Coulomb/Coulomb)
	   And remember that 1 volt/meter == 1 Newton/Coulomb.
	   So a 1 Coulomb charge has E=8.99E3 MV/m at a distance of 1 meter.
	   1 coulomb = 6.24150E+18 e+.
	*/
	const double HW=100.0*mm; // half-width of grid
	PoissonConvolve3D *p = new PoissonConvolve3D(33,33,33,-HW,HW,-HW,HW,-HW,HW);
	p->zeroRhs();
	p->putCharge(0.0,0.0,0.0,6.24150E+18/epsilon0);
	p->defaultApproximation();
	p->solve();
	G4ThreeVector E = p->getE(1.0*meter,0.0,0.0);
	printf("spacecharge: testEfield: E field should be %.2f MV/m: %.2f\n",
			8987.55, E.x()/(megavolt/meter));
	BLAssert(fabs(E.x()/(megavolt/meter)-8987.55) < 0.1 && \
		fabs(E.y()/(megavolt/meter)) < 0.1 && \
		fabs(E.z()/(megavolt/meter)) < 0.1);
/***
	// check for continuity at the edge of the grid.
	printf("#x Phi Ex\n");
	for(double x=0.9*HW; x<1.1*HW; x+=0.01*HW) {
		E = p->getE(x,0.0,0.0);
		printf("%.2f %.5g %.1f\n",x,p->phi(x,0.0,0.0),
				E.x()/(megavolt/meter));
	}
***/
	delete p;
}

void Bunch::testBfield()
{
	/* Magnetic field test
	   A uniform current of 1000 Amps, at radius 0.1 m, has B=0.002 T.
	   1000 Amps = 0.001C charges 1 cm apart moving at 10000 m/s.
	*/
	double c=0.001*coulomb;
	double beta=10000*meter/second/c_light;
	double gamma=1.0/sqrt(1.0-beta*beta);
	double d=1.0*cm;
	PoissonConvolve3D *p = new PoissonConvolve3D(129,129,129,-1000.0,1000.0,
				-1000.0,1000.0,-gamma*1000.0,gamma*1000.0);
	p->zeroRhs();
	std::vector<Charge> v;
	for(double z=-50.0*d; z<50.01*d; z+=d) {
		v.push_back(Charge(0.0,0.0,gamma*z,c/epsilon0));
		p->putCharge(0.0,0.0,gamma*z,c/epsilon0);
	}
	p->setApproximation(v);
	p->solve();
	G4ThreeVector E = p->getE(0.1*meter,0.0,0.0);
	double By = gamma*beta*E[0]/c_light;	// By
	printf("spacecharge: testBfield: By should be 0.002 tesla: %.6f\n",
								By/tesla);
	BLAssert(By-0.002*tesla < 0.00005*tesla);
}

void Bunch::Approx::putCharge(double x, double y, double z, double c)
{
	// (Must add all charge somewhere, as it was put into the grid.)
	int ix=(int)floor((x-xMin)/dx);
	if(ix < 0 || ix >= nx) ix = nx/2;
	int iy=(int)floor((y-yMin)/dy);
	if(iy < 0 || iy >= ny) iy = ny/2;
	int iz=(int)floor((z-zMin)/dz);
	if(iz < 0 || iz >= nz) iz = nz/2;
	int i=index(ix,iy,iz);
	xx[i] += c*x;
	yy[i] += c*y;
	zz[i] += c*z;
	cc[i] += c;
	totalC += fabs(c);
}

void Bunch::Approx::computeApprox()
{
	reduce();
	v.clear();
	for(int ix=0; ix<nx; ++ix) {
		for(int iy=0; iy<ny; ++iy) {
			for(int iz=0; iz<nz; ++iz) {
				int i=index(ix,iy,iz);
				if(cc[i] == 0.0) continue;
				v.push_back(Charge(xx[i]/cc[i], yy[i]/cc[i],
							zz[i]/cc[i],cc[i]));
			}
		}
	}
}

void Bunch::Approx::reduce()
{
	// This routine finds bins in the 3-d histogram that contain less than
	// minFrac of the total charge, and moves their charge 1 bin towards
	// the center. This reduces the number of charges in the approximation.
	// For tests with gaussian beams and nx=ny=nz=9, the reduction was from
	// 124 to 26.
	const double minFrac=0.02;

	// ix<nx/2 region
	for(int ix=0; ix<nx/2; ++ix) {
		for(int iy=0; iy<ny; ++iy) {
			for(int iz=0; iz<nz; ++iz) {
				int i=index(ix,iy,iz);
				if(cc[i] == 0.0) continue;
				if(fabs(cc[i]) < totalC*minFrac) {
					int j=index(ix+1,iy,iz);
					xx[j] += xx[i];
					yy[j] += yy[i];
					zz[j] += zz[i];
					cc[j] += cc[i];
					xx[i] = yy[i] = zz[i] = cc[i] = 0.0;
				}
			}
		}
	}

	// ix>nx/2 region
	for(int ix=nx/2+1; ix<nx; ++ix) {
		for(int iy=0; iy<ny; ++iy) {
			for(int iz=0; iz<nz; ++iz) {
				int i=index(ix,iy,iz);
				if(cc[i] == 0.0) continue;
				if(fabs(cc[i]) < totalC*minFrac) {
					int j=index(ix-1,iy,iz);
					xx[j] += xx[i];
					yy[j] += yy[i];
					zz[j] += zz[i];
					cc[j] += cc[i];
					xx[i] = yy[i] = zz[i] = cc[i] = 0.0;
				}
			}
		}
	}

	// iy<ny/2 region
	for(int iy=0; iy<ny/2; ++iy) {
		for(int ix=0; ix<nx; ++ix) {
			for(int iz=0; iz<nz; ++iz) {
				int i=index(ix,iy,iz);
				if(cc[i] == 0.0) continue;
				if(fabs(cc[i]) < totalC*minFrac) {
					int j=index(ix,iy+1,iz);
					xx[j] += xx[i];
					yy[j] += yy[i];
					zz[j] += zz[i];
					cc[j] += cc[i];
					xx[i] = yy[i] = zz[i] = cc[i] = 0.0;
				}
			}
		}
	}

	// iy>ny/2 region
	for(int iy=ny/2+1; iy<ny; ++iy) {
		for(int ix=0; ix<nx; ++ix) {
			for(int iz=0; iz<nz; ++iz) {
				int i=index(ix,iy,iz);
				if(cc[i] == 0.0) continue;
				if(fabs(cc[i]) < totalC*minFrac) {
					int j=index(ix,iy-1,iz);
					xx[j] += xx[i];
					yy[j] += yy[i];
					zz[j] += zz[i];
					cc[j] += cc[i];
					xx[i] = yy[i] = zz[i] = cc[i] = 0.0;
				}
			}
		}
	}

	// iz<nz/2 region
	for(int iz=0; iz<nz/2; ++iz) {
		for(int ix=0; ix<nx; ++ix) {
			for(int iy=0; iy<ny; ++iy) {
				int i=index(ix,iy,iz);
				if(cc[i] == 0.0) continue;
				if(fabs(cc[i]) < totalC*minFrac) {
					int j=index(ix,iy,iz+1);
					xx[j] += xx[i];
					yy[j] += yy[i];
					zz[j] += zz[i];
					cc[j] += cc[i];
					xx[i] = yy[i] = zz[i] = cc[i] = 0.0;
				}
			}
		}
	}

	// iz>nz/2 region
	for(int iz=nz/2+1; iz<nz; ++iz) {
		for(int ix=0; ix<nx; ++ix) {
			for(int iy=0; iy<ny; ++iy) {
				int i=index(ix,iy,iz);
				if(cc[i] == 0.0) continue;
				if(fabs(cc[i]) < totalC*minFrac) {
					int j=index(ix,iy,iz-1);
					xx[j] += xx[i];
					yy[j] += yy[i];
					zz[j] += zz[i];
					cc[j] += cc[i];
					xx[i] = yy[i] = zz[i] = cc[i] = 0.0;
				}
			}
		}
	}
}

#else // G4BL_FFTW
int dummy_spacecharge=0;
#endif // G4BL_FFTW
