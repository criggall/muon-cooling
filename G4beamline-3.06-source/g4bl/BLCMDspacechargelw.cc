//	BLCMDspacechargelw.cc
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

const G4double UNSET=-1.0e99;
const char *FIELDS = "x:y:z:t:nStep:deltaT:Btot:Etot:deltaPos:deltaB:deltaE";

/**	struct TrajectoryPoint represents one point on the trajectory of a 
 *	MacroParticle. Global coordinates.
 *	pos={x,y,z,t}, v={Vx,Vy,Vz}  [Vx=dx/dt, ...].
 *
 *	A consistency check is available because successive points on the
 *	trajectory can be differenced to get the velocity, which can be
 *	compared to both points' values of v[].
 **/
struct TrajectoryPoint {
	G4double pos[4];
	G4double v[3];
	TrajectoryPoint(const G4double _pos[4], const G4double _v[3])
		{ pos[0]=_pos[0]; pos[1]=_pos[1]; pos[2]=_pos[2];
		  pos[3]=_pos[3]; v[0]=_v[0]; v[1]=_v[1]; v[2]=_v[2]; }
};

/**	class MacroParticle represents a spherical distribution of charge
 *	with uniform density and specified radius, following a specified 
 *	trajectory.
 *
 *	APPROXIMATION: the nonzero radius is used only to prevent infinities
 *	for field computations inside the macro-particle. The speed-of-light
 *	delay within the macro-particle is neglected -- IOW the entire macro
 *	particle is treated as if it were located at its center, but the
 *	charge of the macro-particle is scaled by r/radius when r<radius.
 *
 *	tmin and tmax default to -infinity and +infinity; they are used to
 *	limit the field of the macroparticle when it is created or destroyed
 *	(i.e. secondaries should set tmin, and particles that interact or
 *	decay should set tmax when they are killed).
 **/
class MacroParticle {
	G4double charge;	// total charge 
	G4double radius;	// radius of macro-particle
	G4int K;		// exponent of charge density
	G4double tmin;		// time of creation (or -infinity)
	G4double tmax;		// time of destruction (or +infinity)
	G4int minIndex;		// for pruning trajectory[]
	G4int nCalls;		//  "    "           "
	std::vector<TrajectoryPoint> trajectory;
public:
	static int verbose;
public:
	/// Constructor -- no trajectory point.
	MacroParticle(G4double _charge, G4double _radius, G4int _K=1) :
								trajectory()
		{ charge=_charge; radius=_radius; K=_K;  tmin=-DBL_MAX;
		  tmax=DBL_MAX; minIndex=0; nCalls=0; }

	/// constructor from track, includes trajectory point.
	MacroParticle(const G4Track *track, G4double _charge, G4double _radius,
						G4int _K=1) : trajectory() {
		charge = track->GetDefinition()->GetPDGCharge() * _charge;
		radius=_radius;
		K=_K;
		tmin=-DBL_MAX;
		tmax=DBL_MAX;
		minIndex=0;
		nCalls=0;
		addPoint(track);
	}

	/// setTmin() sets tmin (time of creation).
	void setTmin(G4double v) { tmin = v; }

	/// setTmax() sets tmax (time of destruction).
	void setTmax(G4double v) { tmax = v; }

	/// addPoint() adds a trajectory point to the MacroParticle.
	/// Points must be added in order of increasing time (e.g. in
	/// UserSteppingAction()).
	void addPoint(const G4Track *track) {
		G4double pos[4], vel[3];
		G4ThreeVector p=track->GetPosition();
		pos[0] = p[0];
		pos[1] = p[1];
		pos[2] = p[2];
		pos[3] = track->GetGlobalTime();
		p = track->GetMomentumDirection();
		G4double speed = track->GetVelocity();
		vel[0] = p[0] * speed;
		vel[1] = p[1] * speed;
		vel[2] = p[2] * speed;
		addPoint(pos,vel);
	}

	/// addPoint() adds a trajectory point to the MacroParticle.
	/// Points must be added in order of increasing time (e.g. in
	/// UserSteppingAction()).
	/// pos={x,y,z,t}  v={Vx,Vy,Vz}
	void addPoint(const G4double pos[4], const G4double v[3])
		{ trajectory.push_back(TrajectoryPoint(pos,v)); }

	/// getMostRecentTime() returns the time of the most recent
	/// point in trajectory[].
	G4double getMostRecentTime() {
		BLAssert(trajectory.size() >= 1);
		return trajectory[trajectory.size()-1].pos[3];
	}

	/// computeField() computes the E and B fields due to this macro
	/// particle at the specified point. point[] = {x,y,z,t},
	/// field[]={Bx,By,Bz,Ex,Ey,Ez}; global coordinates.
	/// The MacroParticle must have at least one TrajectoryPoint.
	void computeField(const G4double point[4], G4double field[6]) {
		field[0]=field[1]=field[2]=field[3]=field[4]=field[5]=0.0;
		addField(point,field);
	}

	/// addField() is just like computeField(), except the value
	/// for this MacroParticle is added to the current value of field[].
	void addField(const G4double point[4], G4double field[6]);

	/// pruneTrajectory() will remove unneeded entries in trajectory[]
	/// It requires 100 calls to addfield() between prunes, leaves 10
	/// extra entries (just in case of error in minIndex), and always
	/// leaves at least 20 entries in trajectory[].
	void pruneTrajectory();

	/// testUnits will test the E and B field units.
	/// User validates by looking at stdout (correct and computed values
	/// are printed).
	static void testUnits();
};
int MacroParticle::verbose=0;

/**	class BLCMDspacechargelw - Lienard-Wiechert space-charge computation
 *	This computation works within the Geant4 tracking framework, and does
 *	the computation in the various callbacks. The algorithm tracks
 *	individual particles, and computes their field using macro-particles
 *	that are co-located with the tracked particles. So it computes E and B
 *	fields from the macro-particles and adds this in to the BLGlobalField,
 *	which Geant4 uses for tracking. The field from each macro-particle is
 *	computed via Lienard-Wiechert potentials. When necessary, the
 *	trajectories of other macro-particles are linearly extrapolated from
 *	their current positions (because of the retarded time, this is needed
 *	only needed for macro-particles that are close to the particle currently
 *	being tracked). To handle particle creation and destruction, each
 *	macro-particle has time limits, and does not contribute any field if
 *	its trajectory is evaluated outside those limits. The expensive N^2
 *	computation is done in the beginTrack() callback, not in the much more
 *	frequently called addElementField(). 
 *
 *	For a given particle being tracked, this version first evaluates the
 *	other MacroParticles' field at its position at the start of the
 *	time step; it then linearly extrapolates the particle's trajectory
 *	to the next time step, and evaluates the other MacroParticles' field
 *	there; addElementField() then interpolates between those values
 *	linearly in global time. As the interpolation is linear, it gives a
 *	reasonable value even a short time outside the time step.
 *
 **/
class BLCMDspacechargelw : public BLCommand, public BLCollectiveComputation,
		public BLManager::SteppingAction, 
		public BLManager::TrackingAction,
		public BLElementField {
	G4double deltaT;		// time step (ns)
	G4double radius;		// radius of MacroParticle-s
	G4double charge;		// charge of Macroparticle-s
	G4int trackTwice;		// non-zero to track particle in 
					// beginTrack; 0 to extrapolate linearly
	G4int verbose;			// verbosity conrol for debugging
	G4int K;			// exponent of MP density
	G4int ignoreFieldWhenTracking;
	bool tracking;
	std::vector<MacroParticle> mp;	// parallels vector<BLTrackData>
	// data for EM field computation:
	int index;
	G4double pos0[4], pos1[4];	// for time interpolation
	G4double field0[6], field1[6];	// for time interpolation
	// monitor data:
	struct Monitor {
		double pos[4];
		double field[6];
		double deltaT;
		Monitor() { pos[0]=pos[1]=pos[2]=pos[3]=UNSET;
		    field[0]=field[1]=field[2]=field[3]=field[4]=field[5]=UNSET;
		    deltaT=UNSET; }
		Monitor(G4double p[4], G4double f[6], G4double dT) {
			pos[0]=p[0]; pos[1]=p[1]; pos[2]=p[2]; pos[3]=p[3];
			field[0]=f[0]; field[1]=f[1]; field[2]=f[2];
			field[3]=f[3]; field[4]=f[4]; field[5]=f[5];
			deltaT=dT;
		}
	};
	std::vector<Monitor> mon;
	BLNTuple *ntuple;
	// misc. data:
	BLManager *manager;
	BLRunManager *runManager;
	BLPhysics *physics;
	int nstep;
	std::vector<int> restoreList;
	G4double rangeCut;
public:
	BLCMDspacechargelw();

	G4String commandName() { return "spacechargelw"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	void setStochasticsOff(G4Track *track);

	void restoreStochastics(G4Track *track);

	virtual void beginCollectiveTracking(std::vector<BLTrackData>& v);

	virtual void beginTrack(std::vector<BLTrackData>& v, int _index);

	virtual void collectiveStep(std::vector<BLTrackData>& v);

	virtual void endCollectiveTracking(std::vector<BLTrackData>& v);

	virtual void UserSteppingAction(const G4Step *step);

	virtual void PreUserTrackingAction(const G4Track *track);
        virtual void PostUserTrackingAction(const G4Track *track);

	virtual void addFieldValue(const G4double point[4], G4double field[6])
		const;
};
BLCMDspacechargelw defaultBLCMDspacechargelw;

BLCMDspacechargelw::BLCMDspacechargelw() : BLCollectiveComputation(), BLManager::SteppingAction(),
		BLManager::TrackingAction(), BLElementField(), mp(), mon(),
		restoreList()
{
	registerCommand(BLCMDTYPE_PHYSICS);
	setSynopsis("Lienard-Wiechert space charge computation");
	setDescription("This is a space charge computation that uses "
		"macro-particles to simulate more particles than is feasible "
		"to track individually. Each macro-particle is tracked as a "
		"single particle, but its charge is multiplied by the "
		"macro-particle charge when computing the field. The radius "
		"of the macro-particle is used to avoid the singularity from "
		"a point charge; outside the radius the macro-particle is "
		"treated as a point charge; inside the radius the point-charge "
		"field is multiplied by (r/radius)^K (radius and K are "
		"parameters).\n\n"
		"The trajectory of every particle is kept, and when computing "
		"the field at a point, the intersection of the point's past "
		"lightcone with the trajectory is used to determine the field "
		"from the macro-particle; there is a loop over all particles "
		"except the one currently being tracked. This computation "
		"scales as N^2, where N is the number of macro-particles; that "
		"makes it computationally infeasible for more than a few "
		"hundred macro-particles. But for the particles used, it is "
		"correct to within the following approximations: a) using "
		"macro-particles, b) linearly interpolating between steps, "
		"c) omitting the radiation term in the L-W potential.\n\n"
		"The fields are computed using eq. 63.8-9 (p 162) of Landau "
		"and Lifshitz, _Classical_Theory_of_Fields_, ignoring the "
		"radiation term. The computed fields are used in the usual "
		"Geant4 tracking. Particle creation and destruction are "
		"handled properly.\n\n"
		"This algorithm is primarily intended to test other space "
		"charge algorithms.");

	// initialize class variables here
	deltaT = -1.0;
	radius = -1.0;
	charge = -1.0;
	trackTwice = 0;
	verbose = -1;
	K = 1;
	ignoreFieldWhenTracking = 0;
	tracking = false;
	index = -1;
	ntuple = 0;
	manager = 0;	// not available yet -- set in command(...)
	runManager = 0;	// not available yet -- set in command(...)
	physics = 0;	// not available yet -- set in beginCollectiveTracking
	nstep = 0;
	rangeCut = 0.0;
}

int BLCMDspacechargelw::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	manager = BLManager::getObject();
	runManager = BLRunManager::getObject();

	handleNamedArgs(namedArgs);

	if(verbose < 0) verbose=Param.getInt("steppingVerbose");
	if(verbose < 0) verbose = 0;
	MacroParticle::verbose = verbose;

	if(radius <= 0.0 || charge <= 0.0 || deltaT <= 0.0)
		printError("spacechargelw: radius, charge, deltaT must all be > 0");

	if(K < 0 || K > 4)
		printError("spacechargelw: K must be one of {0,1,2,3,4}");

	runManager->registerCollectiveComputation(this);
	runManager->setCollectiveMode(true);
	runManager->setDeltaT(deltaT);

	manager->registerSteppingAction(this);
	manager->registerTrackingAction(this);

	BLGlobalField::getObject()->addElementField(this); // infinite bounding
							   // box 
	print("");

	MacroParticle::testUnits();

	return 0;
}

void BLCMDspacechargelw::defineNamedArgs()
{
	argDouble(deltaT,"deltaT","Time step (ns).",ns);
	argDouble(radius,"radius","Radius of macro-particles (mm).",mm);
	argDouble(charge,"charge","Charge of macro-particles (times particle charge).");
	argInt(trackTwice,"trackTwice","0=linear extrapolation, 1=track (0)");
	argInt(verbose,"verbose","Non-zero for verbose prints (0).");
	argInt(K,"K","Exponent for macro-particle density (1).");
	argInt(ignoreFieldWhenTracking,"ignoreFieldWhenTracking","For testing only (0).");
}	

void BLCMDspacechargelw::setStochasticsOff(G4Track *track)
{
	restoreList.clear();

	if(!physics->getStochasticsEnabled()) return;

	rangeCut = physics->getRangeCut();

	G4UImanager* UI = G4UImanager::GetUIpointer();
	// turn off delta rays (and other particle production) by
	// setting production cut very high in energy by setting
	// the required range very large.
	UI->ApplyCommand("/run/setCut 1 parsec");
	// Turn off enegry-loss straggling
	UI->ApplyCommand("/process/eLoss/fluct false");

	// run through list of processes, disabling stochastic ones
	G4ProcessManager *pm = track->GetDefinition()->GetProcessManager();
	G4ProcessVector *pv=pm->GetProcessList();
	for(int i=0; i<pv->size(); ++i) {
		if(physics->isStochasticProcess((*pv)[i])) {
			if(pm->GetProcessActivation(i))
				restoreList.push_back(i);
			pm->SetProcessActivation(i,false);
		}
	}
}

void BLCMDspacechargelw::restoreStochastics(G4Track *track)
{
	if(!physics->getStochasticsEnabled()) return;

	G4UImanager* UI = G4UImanager::GetUIpointer();
	static char setCmd[64] = { 0 };
	if(setCmd[0] == '\0') sprintf(setCmd,"/run/setCut %g mm",rangeCut/mm);
	UI->ApplyCommand(setCmd);
	UI->ApplyCommand("/process/eLoss/fluct true");

	// run through list of disabled processes, restoring them
	G4ProcessManager *pm = track->GetDefinition()->GetProcessManager();
	for(unsigned int j=0; j<restoreList.size(); ++j) {
		pm->SetProcessActivation(restoreList[j],true);
	}
}

void BLCMDspacechargelw::beginCollectiveTracking(std::vector<BLTrackData>& v)
{
	physics = manager->getPhysics();

	// initialize mp[] and mon[]
	mp.clear();
	mon.clear();
	for(unsigned i=0; i<v.size(); ++i) {
		mp.push_back(MacroParticle(v[i].track,charge,radius,K));
		mon.push_back(Monitor());
	}

	runManager->setDeltaT(deltaT);

	nstep = 0;
	MacroParticle::verbose = verbose;

	ntuple = BLNTuple::create("","NTuples","CollectiveMonitor",FIELDS,"");
}

void BLCMDspacechargelw::beginTrack(std::vector<BLTrackData>& v, int _index)
{
	if(manager->getState() != BEAM) return;

	index = _index;
	BLAssert(index >= 0 && (unsigned)index<v.size());

	G4Track *track = v[index].track;
	if(verbose==1 || verbose == 2) 
		printf("spacechargelw::beginTrack(%d)\n",index);

	// if this is a new (secondary) particle, add it to mp[] and mon[]
	if((unsigned)index >= mp.size()) {
		BLAssert((unsigned)index == mp.size());
		mp.push_back(MacroParticle(track,charge,radius,K));
		mp[index].setTmin(track->GetGlobalTime());
		mon.push_back(Monitor());
		if(verbose >= 3)
			printf("spacechargelw::beginTrack(%d) added to mp[]\n",index);
	}

	// compute field of all other MacroParticle-s at this track's position
	G4ThreeVector p=track->GetPosition();
	pos0[0] = p[0];
	pos0[1] = p[1];
	pos0[2] = p[2];
	pos0[3] = track->GetGlobalTime();
	field0[0]=field0[1]=field0[2]=field0[3]=field0[4]=field0[5]=0.0;
	for(unsigned i=0; i<mp.size(); ++i) {
		if(i == (unsigned)index) continue;
		mp[i].addField(pos0,field0);
	}
	if(verbose >= 3)
		printf("spacechargelw::beginTrack(%4d) pos0:%.3f,%.3f,%.3f,%.3f field0:%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
			index,pos0[0],pos0[1],pos0[2],pos0[3],
			field0[0]/tesla,field0[1]/tesla,field0[2]/tesla,
			field0[3]/(megavolt/meter),field0[4]/(megavolt/meter),
			field0[5]/(megavolt/meter));

	if(trackTwice) {
		// track with stochastics off, full external field, and
		// constant space charge field (at track's current position).
		int isave = index;
		pos1[0] = pos0[0]; pos1[1] = pos0[1]; pos1[2] = pos0[2];
		pos1[3] = pos0[3] + deltaT; // values for a constant field
		G4Track *temp = new G4Track(*track);
		manager->setState(SPECIAL);
		setStochasticsOff(temp);
		if(manager->getSteppingVerbose() != 0)
			printf("spacechargelw: pre-tracking begun\n");
		runManager->processOneTrack(temp); // will be suspended
		if(manager->getSteppingVerbose() != 0)
			printf("spacechargelw: pre-tracking ended\n");
		restoreStochastics(temp);
		manager->setState(BEAM);
		p = temp->GetPosition();
		pos1[0] = p[0];
		pos1[1] = p[1];
		pos1[2] = p[2];
		pos1[3] = temp->GetGlobalTime();
		delete temp;
		index = isave;
		// protect against immediate track death
		if(pos1[3]-pos0[3] < 1.0E-4*ns) pos1[3] += 1.0E-4*ns;
	} else {
		// Linearly extrapolate this track to the next time step.
		p = track->GetMomentumDirection();
		G4double speed = track->GetVelocity();
		pos1[0] = pos0[0] + p[0]*speed*deltaT;
		pos1[1] = pos0[1] + p[1]*speed*deltaT;
		pos1[2] = pos0[2] + p[2]*speed*deltaT;
		pos1[3] = pos0[3] + deltaT;
	}

	// compute field of all other MacroParticle-s at extrapolated position
	// (this can be at a global time beyond the known trajectory of some of 
	//  the other macro-particles; they are linearly extrapolated inside
	//  their addfield()).
	field1[0]=field1[1]=field1[2]=field1[3]=field1[4]=field1[5]=0.0;
	for(unsigned i=0; i<mp.size(); ++i) {
		if(i == (unsigned)index) continue;
		mp[i].addField(pos1,field1);
	}
	if(verbose >= 3)
		printf("                        pos1:%.3f,%.3f,%.3f,%.3f field1:%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
			pos1[0],pos1[1],pos1[2],pos1[3],
			field1[0]/tesla,field1[1]/tesla,field1[2]/tesla,
			field1[3]/(megavolt/meter),field1[4]/(megavolt/meter),
			field1[5]/(megavolt/meter));

	// monitor the differences
	// (mon[index] contains the previous extrapolated values)
	if(mon[index].pos[0] != UNSET && fabs(pos0[3]-mon[index].pos[3]) < 0.002*ns) {
		G4ThreeVector deltaPos, deltaB, deltaE;
		for(int i=0; i<3; ++i) {
			deltaPos[i] = pos0[i] - mon[index].pos[i];
			deltaB[i] = field0[i] - mon[index].field[i];
			deltaE[i] = field0[i+3] - mon[index].field[i+3];
		}
		if(verbose >= 3)
		    printf("deltaPos=%.6f,%.6f,%.6f deltaB=%.6f,%.6f,%.6f "
		    	"deltaE=%.5f,%.6f,%.6f\n",
		    	deltaPos[0]/mm,deltaPos[1]/mm,deltaPos[2]/mm,
			deltaB[0]/tesla,deltaB[1]/tesla,deltaB[2]/tesla,
			deltaE[0]/(megavolt/meter),deltaE[1]/(megavolt/meter),
			deltaE[2]/(megavolt/meter));
		double data[11];
		data[0] = pos0[0];	// x
		data[1] = pos0[1];	// y
		data[2] = pos0[2];	// z
		data[3] = pos0[3];	// t
		data[4] = nstep;	// nStep
		data[5] = mon[index].deltaT;	// deltaT
		data[6] = sqrt(field0[0]*field0[0]+field0[1]*field0[1]+
					field0[2]*field0[2]); // Btot
		data[7] = sqrt(field0[3]*field0[3]+field0[4]*field0[4]+
					field0[5]*field0[5]); // Etot
		data[8] = deltaPos.mag(); // deltaPos
		data[9] = deltaB.mag();	// deltaB
		data[10] = deltaE.mag();	// deltaE
		ntuple->appendRow(data,11);
	}
	// put the new extrapolated values into mon[]
	mon[index] = Monitor(pos1,field1,deltaT);

	// Now that we have done the expensive N^2 computation, BLCMDspacechargelw::addField
	// will interpolate between pos0[] and pos1[] to approximate the field
	// at the trajectory of this particle with the corresponding value on 
	// the linearly-extrapolated trajectory.
}

void BLCMDspacechargelw::collectiveStep(std::vector<BLTrackData>& v)
{
	if(manager->getState() != BEAM) return;

	// TODO: update radius here
	// TODO: update deltaT here

	// prune MacroParticle trajectories
	for(unsigned i=0; i<mp.size(); ++i)
		mp[i].pruneTrajectory();
}

void BLCMDspacechargelw::endCollectiveTracking(std::vector<BLTrackData>& v)
{
}

void BLCMDspacechargelw::UserSteppingAction(const G4Step *step)
{
	if(manager->getState() != BEAM) return;	// includes SPECIAL

	BLAssert(index>=0 && (unsigned)index<mp.size());

	// add trajectory point unless tiny increment in time
	if(step->GetTrack()->GetGlobalTime()-mp[index].getMostRecentTime() >=
								deltaT/10.0) {
		mp[index].addPoint(step->GetTrack());
		if(verbose >= 3)
			printf("spacechargelw::UserSteppingAction index=%d point added\n",index);
	} else {
		if(verbose >= 3)
			printf("spacechargelw::UserSteppingAction index=%d NO POINT - dt too small\n",index);
	}
}

void BLCMDspacechargelw::PreUserTrackingAction(const G4Track *track)
{
	if(manager->getState() != BEAM) return;

	if(verbose >= 3)
		printf("spacechargelw::PreUserTrackingAction index=%d\n",index);

	BLAssert(index>=0 && (unsigned)index<mp.size());

	tracking = true;
}

void BLCMDspacechargelw::PostUserTrackingAction(const G4Track *track)
{
	if(manager->getState() != BEAM) return;

	BLAssert(index>=0 && (unsigned)index<mp.size());

	tracking = false;

	// add this point to the trajectory (no call to UserSteppingAction
	// when track is suspended).
	if(track->GetGlobalTime()-mp[index].getMostRecentTime() > deltaT/100.0)
		mp[index].addPoint(track);

	// set Tmax for tracks that are killed
	G4TrackStatus trackStatus = track->GetTrackStatus();
	if(trackStatus != fAlive && trackStatus != fStopButAlive && trackStatus != fSuspend) {
		mp[index].setTmax(track->GetGlobalTime());
		if(verbose >= 3)
			printf("spacechargelw::PostUserTrackingAction index=%d track killed\n",index);
	} else {
		if(verbose >= 3)
			printf("spacechargelw::PostUserTrackingAction index=%d\n",index);
	}

	index = -1;
}

void BLCMDspacechargelw::addFieldValue(const G4double point[4], G4double field[6]) const
{
	if(manager->getState() != BEAM) return;
	if(ignoreFieldWhenTracking != 0 && tracking) return;

	double fthis[6];
	fthis[0]=fthis[1]=fthis[2]=fthis[3]=fthis[4]=fthis[5]=0.0;

	if(index >= 0) {
		// interpolate field between pos0 and pos1 for this track
		BLAssert(pos1[3]-pos0[3] > 1.0e-8*ns);
		G4double f = (point[3]-pos0[3])/(pos1[3]-pos0[3]);
		fthis[0] = (1.0-f)*field0[0] + f*field1[0];
		fthis[1] = (1.0-f)*field0[1] + f*field1[1];
		fthis[2] = (1.0-f)*field0[2] + f*field1[2];
		fthis[3] = (1.0-f)*field0[3] + f*field1[3];
		fthis[4] = (1.0-f)*field0[4] + f*field1[4];
		fthis[5] = (1.0-f)*field0[5] + f*field1[5];
		if(verbose >= 3)
			printf("spacechargelw::addFieldValue index=%d f=%.4f point=%.3f,%.3f,%.3f,%.3f\n",
				index,f,point[0],point[1],point[2],point[3]);
	} else {
		// not during tracking, sum all macro-particles
		for(unsigned i=0; i<mp.size(); ++i) {
		    ((BLCMDspacechargelw*)this)->mp[i].addField(point,fthis);
		}
	}
	if(verbose >= 3)
		printf("spacechargelw::addFieldValue pos:%.3f,%.3f,%.3f,%.3f field:%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
			point[0],point[1],point[2],point[3],
			fthis[0]/tesla,fthis[1]/tesla,fthis[2]/tesla,
			fthis[3]/(megavolt/meter),fthis[4]/(megavolt/meter),
			fthis[5]/(megavolt/meter));

	field[0] += fthis[0];
	field[1] += fthis[1];
	field[2] += fthis[2];
	field[3] += fthis[3];
	field[4] += fthis[4];
	field[5] += fthis[5];
}


void MacroParticle::addField(const G4double p[4], G4double field[6])
{
	// p[] is the point at which the field is to be evaluated.
	// q[] is a point on the trajectory, v[] is its velocity at q.
	// pq[] = p - q.
	double q[4], v[3], pq[4];
	pq[0] = pq[1] = pq[2] = pq[3] = 0.0;

	if(charge == 0.0) return;
	BLAssert(trajectory.size() > 0);

	/* find the TrajectoryPoint closest before the lightcone from p
	   Use first point (i==0) if before trajectory began. */
	int i;
	for(i=(int)trajectory.size()-1; i>=0; --i) {
		const double *qq = trajectory[i].pos;
		pq[0] = p[0] - qq[0];
		pq[1] = p[1] - qq[1];
		pq[2] = p[2] - qq[2];
		pq[3] = p[3] - qq[3];
		double d = sqrt(pq[0]*pq[0]+pq[1]*pq[1]+pq[2]*pq[2]);
		if(i == 0 || d/c_light < pq[3])
			break;
	}
	q[0] = trajectory[i].pos[0];
	q[1] = trajectory[i].pos[1];
	q[2] = trajectory[i].pos[2];
	q[3] = trajectory[i].pos[3];
	v[0] = trajectory[i].v[0];
	v[1] = trajectory[i].v[1];
	v[2] = trajectory[i].v[2];

	// update trajectory pruning info
	if(i < minIndex) minIndex = i;
	++nCalls;

	/* now solve for the time of the intersection of the lightcone from p 
	   with the trajectory, assuming v is constant; t is the offset 
	   from q[3]. The other root is the unphysical advanced field. */
	double v2 = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];
	double c2 = c_light*c_light;
	double a = c2 - v2;
	double b = 2.0*(pq[0]*v[0]+pq[1]*v[1]+pq[2]*v[2]-c2*pq[3]);
	double c = c2*pq[3]*pq[3]-(pq[0]*pq[0]+pq[1]*pq[1]+pq[2]*pq[2]);
	double discriminant=b*b-4.0*a*c;
	// if discriminant ~ 0, p is deep inside this MacroParticle -> no field
	if(discriminant < 0.001) {
		if(verbose >= 3) printf("addField: discriminant < 0.001\n");
		return;
	}
	// note that a>0 so -sqrt(discriminant) gives the retarded root
	double t = (-b - sqrt(discriminant))/(2.0*a);

	/* update q[] and pq[] to the actual intersection */
	q[0] += t*v[0];
	q[1] += t*v[1];
	q[2] += t*v[2];
	q[3] += t;
	pq[0] = p[0] - q[0];
	pq[1] = p[1] - q[1];
	pq[2] = p[2] - q[2];
	pq[3] = p[3] - q[3];

	// apply creation/destruction limits.
	if(q[3] < tmin || q[3] > tmax) {
		if(verbose >= 3) printf("addField: t=%.4f tmin=%.4f tmax=%.4f\n",
					q[3],tmin,tmax);
		return;
	}

	// r = distance from q to p
	double r = sqrt(pq[0]*pq[0]+pq[1]*pq[1]+pq[2]*pq[2]);
	// require accuracy in the lightcone: 1 ppm or 0.01 micron
	double tolerance=(r>10.0*mm ? r/1000000.0 : 0.00001*mm);
	double light_cone_error=fabs(r-c_light*pq[3]);
	BLAssert(light_cone_error < tolerance);

	/* effective charge */
	double f = 1.0;
	if(r < radius) {
		double a=r/radius;
		switch(K) {
		case 0:	f = 1.0;	break;
		case 1:	f = a;		break;
		case 2:	f = a*a;	break;
		case 3:	f = a*a*a;	break;
		case 4:	f = a*a*a*a;	break;
		}
	}
	if(f < 1.0e-4) return;
	double ec = f * charge;

	/* get E and B fields  - Landau and Lifshitz, _Classical_Theory_of_
	   _Fields_, eq. 63.8-9 p 162. Ignore radiation term.  
	   Units done by hand: the usual Coulomb constant is divided by 1000,
	   and B is divided by c_light.
	*/
	double k = r - (pq[0]*v[0]+pq[1]*v[1]+pq[2]*v[2])/c_light;
	double E[3];
	k = ec * e_SI * 8.987551787E6 * (1.0-v2/c_light/c_light)/k/k/k;
	E[0] = k * (pq[0]-v[0]*r/c_light);
	E[1] = k * (pq[1]-v[1]*r/c_light);
	E[2] = k * (pq[2]-v[2]*r/c_light);
	field[0] += (pq[1]*E[2] - pq[2]*E[1])/r/c_light;
	field[1] += (pq[2]*E[0] - pq[0]*E[2])/r/c_light;
	field[2] += (pq[0]*E[1] - pq[1]*E[0])/r/c_light;
	field[3] += E[0];
	field[4] += E[1];
	field[5] += E[2];
}

void MacroParticle::pruneTrajectory()
{
	const int safety=10;	// extra entries left in trajectory[]
	const int minKeep=20;	// min # entries to keep
	const int minCalls=100;	// min # calls to addfield() between prunes

	if(nCalls < minCalls || minIndex <= safety+1) return;

	minIndex -= safety;
	if(trajectory.size()-(unsigned)minIndex < (unsigned)minKeep) return;

	for(int i=0; i<minIndex; ++i)
		trajectory.erase(trajectory.begin());

	if(verbose >= 3) printf("MacroParticle:prune removed %d newSize=%ld tMin=%.6f\n",
			minIndex,(long)trajectory.size(),trajectory[0].pos[3]);

	minIndex = trajectory.size();
	nCalls = 0;
}

void MacroParticle::testUnits()
{
	printf("MacroParticle::testUnits()\n");
	G4double pos[4]={0,0,0,0};
	G4double v[3]={0,0,0};
	G4double field[6]; 

	/* UNITS TEST
	   E = k*Q/r/r,  k=8.987551787E9 (Newton*meter*meter/Coulomb/Coulomb)
	   And remember that 1 volt/meter == 1 Newton/Coulomb.
	   So a 1 Coulomb charge has E=8.99E3 MV/m at a distance of 1 meter.
	*/
	MacroParticle q(coulomb,0.1*mm);
	pos[0]=pos[1]=pos[2]=pos[3]=0.0;
	v[0]=v[1]=v[2]=0.0;
	q.addPoint(pos,v);
	pos[0] = 1.0*meter;
	q.computeField(pos,field);
	printf("E field test: Ex should be 8988 MV/m   pos=%.0f,%.0f,%.0f  E=%.4g,%.4g,%.4g MV/m\n",
		pos[0]/mm,pos[1]/mm,pos[2]/mm,
		field[3]/(megavolt/meter),field[4]/(megavolt/meter),
		field[5]/(megavolt/meter));

	pos[0] = 0.0;
	pos[1] = 2.0*meter;
	q.computeField(pos,field);
	printf("E field test: Ey should be 2247 MV/m   pos=%.0f,%.0f,%.0f  E=%.4g,%.4g,%.4g MV/m\n",
		pos[0]/mm,pos[1]/mm,pos[2]/mm,
		field[3]/(megavolt/meter),field[4]/(megavolt/meter),
		field[5]/(megavolt/meter));


	pos[1] = 0.0;
	pos[2] = -1.0*meter;
	q.computeField(pos,field);
	printf("E field test: Ez should be -8988 MV/m   pos=%.0f,%.0f,%.0f  E=%.4g,%.4g,%.4g MV/m\n",
		pos[0]/mm,pos[1]/mm,pos[2]/mm,
		field[3]/(megavolt/meter),field[4]/(megavolt/meter),
		field[5]/(megavolt/meter));


	/* Magnetic field test
	   A uniform current of 1000 Amps, at radius 0.1 m, has B=0.002 T.
	   1000 Amps = 0.001C charges 1 cm apart moving at 10000 m/s.
	*/
	int N=8001;		// # charges
	double dz=10*mm;	// their interval along z
	double charge = 0.001*coulomb;	// their charge
	double speed=10000*meter/second; // their speed along z
	std::vector<MacroParticle> array;
	for(int i=0; i<N; ++i) {
		array.push_back(MacroParticle(charge,0.1*mm));
		pos[0]=pos[1]=pos[2]=pos[3]=0.0;
		pos[2] = (i-N/2) * dz;
		v[0]=v[1]=0.0;
		v[2] = speed;
		array[i].addPoint(pos,v);
	}
	field[0]=field[1]=field[2]=field[3]=field[4]=field[5]=0.0;
	pos[0]=pos[1]=pos[2]=pos[3]=0.0;
	pos[0] = 0.1*meter;
	for(int i=0; i<N; ++i) {
		double f[6];
		array[i].computeField(pos,f);
		for(int j=0; j<6; ++j) field[j] += f[j];
	}
	printf("B field test: By should be 0.002 T  pos=%.0f,%.0f,%.0f  B=%.6g,%.6g,%.6g T\n",
		pos[0]/mm,pos[1]/mm,pos[2]/mm,
		field[0]/tesla,field[1]/tesla,field[2]/tesla);
	field[0]=field[1]=field[2]=field[3]=field[4]=field[5]=0.0;
	pos[0]=pos[1]=pos[2]=pos[3]=0.0;
	pos[1] = 0.2*meter; // further away
	for(int i=0; i<N; ++i) {
		double f[6];
		array[i].computeField(pos,f);
		for(int j=0; j<6; ++j) field[j] += f[j];
	}
	printf("B field test: Bx should be -0.001 T  pos=%.0f,%.0f,%.0f  B=%.6g,%.6g,%.6g T\n",
		pos[0]/mm,pos[1]/mm,pos[2]/mm,
		field[0]/tesla,field[1]/tesla,field[2]/tesla);
}
