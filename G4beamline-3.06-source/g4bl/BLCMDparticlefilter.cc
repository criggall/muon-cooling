//	BLCMDparticlefilter.cc
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

#include <set>

#include "G4UserLimits.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Decay.hh"

#include "BLParam.hh"
#include "BLElement.hh"
#include "BLManager.hh"
#include "BLEvaluator.hh"

/**	class BLCMDparticlefilter - implements a ParticleFilter
 *	Each placement of this element generates a physical volume that will
 *	force the decay of any particle named in the "decay" argument, and
 *	will kill any particle named in its "kill" argument, or not named in
 *	its "keep" argument if that is non-empty.
 *
 *	Note that the standard Decay process is disabled for all particles 
 *	named in the "decay" agument, so they will only decay when they
 *	enter a ParticleFilter volume.
 *
 *	Note also the G4VProcess of the ParticleFilterPlacement is only
 *	applied to particles listed in the decay argument.
 *
 *	There is a bug: the decay happens after the first step in the volume;
 *	so maxStep is arranged to ensure two steps occur in the volume. This
 *	appears to be in the way Geant4 handles physical volumes while
 *	tracking....
 **/
class BLCMDparticlefilter : public BLElement {
	G4double radius;
	G4double innerRadius;
	G4double height;
	G4double width;
	G4double length;
	G4String material;
	G4String color;
	G4VSolid *solid;
	G4double maxStep;
	G4String decay;
	G4String kill;
	G4String keep;
	G4int nWait;
	G4int referenceWait;
	G4String require;
	G4int steppingVerbose;
	std::set<G4String> decaySet;
	std::set<G4String> killSet;
	std::set<G4String> keepSet;
	void setup();
	friend class ParticleFilterPlacement;
	static int instance;
public:
	/// Default constructor.
	BLCMDparticlefilter();

	/// Destructor.
	virtual ~BLCMDparticlefilter();

	/// clone()
	BLElement *clone() { return new BLCMDparticlefilter(*this); }

	/// Copy constructor.
	BLCMDparticlefilter(const BLCMDparticlefilter& r);

	/// commandName() returns "particlefilter".
	G4String commandName() { return "particlefilter"; }

	/// command() implements the particlefilter command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() will construct the particlefilter.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns this element's Length along the Z axis.
	virtual G4double getLength() { return length; }

	/// getWidth() returns this element's Width along the X axis.
	virtual G4double getWidth() { return width; }

	/// getHeight() returns this element's height along the Y axis.
	virtual G4double getHeight() { return height; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	virtual G4bool isOK() { return true; }

	/// generatePoints() from BLElement
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isOutside() from BLElement
	G4bool isOutside(G4ThreeVector &local, G4double tolerance);
};
int BLCMDparticlefilter::instance = 0;
BLCMDparticlefilter defaultParticleFilter;


/**	class ParticleFilterPlacement is one placement of a ParticleFilter
 *	element. It inherits from G4Decay which does the heavy lifting,
 *	all we do here is decide whether or not to decay. We also kill
 *	tracks in BLManager::SteppingAction, because you cannot do that inside
 *	PostStepGetPhysicalInteractionLength().
 *
 *	Because geometry construction time is too early to manipulate
 *	physics processes, it registers a callback with BLManager to
 *	do it before the reference particle.
 **/
class ParticleFilterPlacement : public G4Decay, public BLManager::SteppingAction,
				public BLManager::TrackingAction, public BLCallback {
	BLCMDparticlefilter *filter;
	G4VPhysicalVolume *vol;
	unsigned long nDecayed;
	unsigned long nKilled;
	int nWait;
	int referenceWait;
public:
	ParticleFilterPlacement(BLCMDparticlefilter *_pf,
				G4VPhysicalVolume *_pv, G4String _name);

	void callback(int type);

	G4bool IsApplicable(const G4ParticleDefinition&);
	G4double AtRestGetPhysicalInteractionLength(const G4Track& track,
					G4ForceCondition* condition);
	G4double PostStepGetPhysicalInteractionLength( const G4Track& track,
		G4double previousStepSize, G4ForceCondition* condition);
	void UserSteppingAction(const G4Step *step);
	void PreUserTrackingAction(const G4Track *track);
	void PostUserTrackingAction(const G4Track *track) { }
};


BLCMDparticlefilter::BLCMDparticlefilter() : BLElement(), decaySet(), killSet(),
								keepSet()
{
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("Will kill particles from a list, or force particles to decay.");
	setDescription("A particlefilter will force the decay of certain particles\n"
		"when they enter the physical volume of the element.\n"
		"The list of affected particles is in the 'decay' argument,\n"
		"and the normal Decay process is disabled for them.\n"
		"In addition, the 'kill' argument is a list of particles to\n"
		"kill when they enter the element, and the 'keep' argument\n"
		"will kill all particles not named, if it is not empty.\n\n"
		"require is an expression invloving track parameters that will "
		"kill the track if it evaluates to zero (use a comparison "
		"operator). The variables available are:\n"
		"    x,y,z,t,Px,Py,Pz,Ptot,PDGid,EventID,TrackID,ParentID\n"
		"Units are mm, ns, MeV/c.\n\n"
		"If nWait > 1, particles will not be killed until they hit "
		"this element nwait times; this can be used to limit the "
		"number of revolutions around a ring. Decays are unaffected "
		"by nWait. referenceWait does the same for the reference "
		"particle\n"
		"The element can be either a cylinder or a box: set length and "
		"radius, or set length and width and height.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.");
	
	// initial field values
	radius = -1.0;
	innerRadius = 0.0;
	height = 0.0;
	width = 0.0;
	length = 1.0*mm;
	material = "";
	color = "1,1,1";
	solid = 0;
	maxStep = -99.0;
	decay = "";
	kill = "";
	keep = "";
	nWait = 1;
	referenceWait = 1;
	require = "";
	steppingVerbose = 0;
}

BLCMDparticlefilter::~BLCMDparticlefilter()
{
}

BLCMDparticlefilter::BLCMDparticlefilter(const BLCMDparticlefilter& r) : BLElement(r),
		decaySet(r.decaySet), killSet(r.killSet), keepSet(r.keepSet)
{
	radius = r.radius;
	innerRadius = r.innerRadius;
	height = r.height;
	width = r.width;
	length = r.length;
	material = r.material;
	color = r.color;
	solid = 0;
	maxStep = r.maxStep;
	decay = r.decay;
	kill = r.kill;
	keep = r.keep;
	nWait = r.nWait;
	referenceWait = r.referenceWait;
	require = r.require;
	steppingVerbose = r.steppingVerbose;
}

int BLCMDparticlefilter::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("particlefilter: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultParticleFilter.handleNamedArgs(namedArgs);
	}

	BLCMDparticlefilter *t = new BLCMDparticlefilter(defaultParticleFilter);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	t->setup();

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	return retval;
}

void BLCMDparticlefilter::setup()
{
	// get set of particle names into decaySet
	std::vector<G4String> v=splitString(decay,',',true);
	for(unsigned i=0; i<v.size(); ++i) {
		if(v[i].size() == 0) continue;
		decaySet.insert(v[i]);
	}
/***
	unsigned int place, next;
	for(place=next=0; next<decay.size(); place=next+1) {
		next = decay.find(",",place);
		G4String p;
		if(next < decay.size())
			p = decay.substr(place,next-place);
		else
			p = decay.substr(place);
		if(p.size() == 0) break;
		decaySet.insert(p);
	}
***/
	// get set of particle names into killSet
	v=splitString(kill,',',true);
	for(unsigned i=0; i<v.size(); ++i) {
		if(v[i].size() == 0) continue;
		killSet.insert(v[i]);
	}
/***
	for(place=next=0; next<kill.size(); place=next+1) {
		next = kill.find(",",place);
		G4String p;
		if(next < kill.size())
			p = kill.substr(place,next-place);
		else
			p = kill.substr(place);
		if(p.size() == 0) break;
		killSet.insert(p);
	}
***/
	// get set of particle names into keepSet
	v=splitString(keep,',',true);
	for(unsigned i=0; i<v.size(); ++i) {
		if(v[i].size() == 0) continue;
		keepSet.insert(v[i]);
	}
/***
	for(place=next=0; next<keep.size(); place=next+1) {
		next = keep.find(",",place);
		G4String p;
		if(next < keep.size())
			p = keep.substr(place,next-place);
		else
			p = keep.substr(place);
		if(p.size() == 0) break;
		keepSet.insert(p);
	}
***/
	// if keep is nonempty, add decaySet to it (they will decay)
	if(keepSet.size() > 0) {
		std::set<G4String>::iterator i;
		for(i=decaySet.begin(); i!=decaySet.end(); ++i)
			keepSet.insert(*i);
	}
}

void BLCMDparticlefilter::defineNamedArgs()
{
	argDouble(radius,"radius","The radius of the cylindrical particlefilter (mm).");
	argDouble(innerRadius,"innerRadius","The inner radius of the cylindrical particlefilter (0 mm, solid).");
	argDouble(height,"height","The height of the rectangular particlefilter (mm).");
	argDouble(width,"width","The width of the rectangular particlefilter (mm).");
	argDouble(length,"length","The length of the particlefilter (mm).");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm).");
	argString(material,"material","The material (default=parent's).");
	argString(color,"color","The color of the particlefilter (white).");
	argString(decay,"decay","A comma-separated list of particle names to decay.",false);
	argString(kill,"kill","A comma-separated list of particle names to kill.",false);
	argString(keep,"keep","A comma-separated list of particle names to keep.",false);
	argInt(nWait,"nWait","Intersection # to do the kill (default = 1)");
	argInt(referenceWait,"referenceWait","Intersection # for reference (default = 1)");
	argString(require,"require","Expression which will kill the track if zero.",false);
	argInt(steppingVerbose,"steppingVerbose","Nonzero to display track kills.");
	argString(decay,"decays","Synonym for decay.",false);

	if(radius > 0.0)
		width = height = 2.0 * radius;
}

void BLCMDparticlefilter::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");
	steppingVerbose = Param.getInt("steppingVerbose");

	// ensure at least two steps will be taken inside the physical volume
	G4double minStep = (width<height ? width : height);
	minStep = (minStep<length ? minStep : length);
	minStep = minStep/2.0;
	maxStep = (minStep<maxStep ? minStep : maxStep);

	G4Material *mat;
	if(material != "")
		mat = getMaterial(material);
	else
		mat = parent->GetMaterial();

	G4String thisname = parentName+getName();

	if(!solid) {
		if(radius > 0.0) {
			solid = new G4Tubs(thisname+"Tubs", innerRadius, radius,
					length/2.0, 0.0, 2.0*pi);
		} else if(height > 0.0 && width > 0.0) {
			solid = new G4Box(thisname+"Box",width/2.0,
					height/2.0,length/2.0);
		} else {
			printError("particlefilter::construct %s INVALID - no "
				"radius or height&width",thisname.c_str());
			return;
		}
	}
	G4LogicalVolume *lv = new G4LogicalVolume(solid,mat, thisname+"LogVol");
	lv->SetVisAttributes(getVisAttrib(color));
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");
	lv->SetUserLimits(new G4UserLimits(maxStep));

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4VPhysicalVolume *vol = new G4PVPlacement(g4rot, relativePosition,
				lv,thisname,parent,false,0,surfaceCheck);

	char processname[64];
	sprintf(processname,"BLCMDparticlefilter%d",++instance);
	ParticleFilterPlacement *place = 
			new ParticleFilterPlacement(this,vol,processname);

	if(killSet.size() > 0 || keepSet.size() > 0 || require != "") {
	    BLManager::getObject()->registerTrackingAction(place);
	    BLManager::getObject()->registerSteppingAction(vol,place);
	}

	printf("BLCMDparticlefilter::Construct %s parent=%s\n",
		thisname.c_str(),parentName.c_str());
}

void BLCMDparticlefilter::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	if(radius > 0.0)
		generateTubs(npoints, innerRadius, radius, 0.0, 360.0*deg, length, v);
	else
		generateBox(npoints,width,height,length,v);
}

G4bool BLCMDparticlefilter::isOutside(G4ThreeVector &local, G4double tolerance)
{
	if(radius > 0.0) {
		G4double r = sqrt(local[0]*local[0]+local[1]*local[1]);
		return r > radius-tolerance || r < innerRadius+tolerance ||
			fabs(local[2]) > length/2.0-tolerance;
	} else {
		return fabs(local[0]) > width/2.0-tolerance ||
			fabs(local[1]) > height/2.0-tolerance ||
			fabs(local[2]) > length/2.0-tolerance;
	}
}

ParticleFilterPlacement::ParticleFilterPlacement(BLCMDparticlefilter *_pf,
		G4VPhysicalVolume *_pv, G4String _name) : 
					G4Decay(_name), BLCallback()
{
	filter=_pf;
	vol=_pv; 
	nDecayed = 0;
	nKilled = 0;

	// now is too early to manage the physics processes
	BLManager::getObject()->registerCallback(this,0);

	// also register a callback after tracking, to display statistics
	BLManager::getObject()->registerCallback(this,2);
}

void ParticleFilterPlacement::callback(int type)
{
	if(type == 0) {
	    // add ourselves to the process list for each applicable particle
	    // and remove the Decay process if necessary.
	    G4ParticleTable::G4PTblDicIterator *myParticleIterator =
			G4ParticleTable::GetParticleTable()->GetIterator();
	    myParticleIterator->reset();
	    while((*myParticleIterator)()) {
		G4ParticleDefinition *pd = myParticleIterator->value();
		if(IsApplicable(*pd)) {
			G4ProcessManager *pmgr = pd->GetProcessManager();
again:			G4ProcessVector *pv = pmgr->GetProcessList();
			int n = pv->size();
			for(int i=0; i<n; ++i) {
				G4VProcess *process = (*pv)[i];
				if(process->GetProcessName().find("Decay") !=
							G4String::npos) {
					pmgr->RemoveProcess(process);
					goto again;
				}
			}
			pmgr->AddProcess(this);
			pmgr->SetProcessOrdering(this,idxPostStep);
			pmgr->SetProcessOrdering(this,idxAtRest);
//for(int i=0; i<pv->size(); ++i) printf("ProcessName='%s'\n",(*pv)[i]->GetProcessName().c_str());
		}
	    }
	} else if(type == 2) {
		printf("particlefilter %s: %lu decays  %lu killed\n",
			filter->getName().c_str(),nDecayed,nKilled);
	}
}

G4bool ParticleFilterPlacement::IsApplicable(const G4ParticleDefinition& pd)
{
	if(filter->decaySet.count(pd.GetParticleName()) > 0)
		return G4Decay::IsApplicable(pd);
	return false;
}

G4double ParticleFilterPlacement::AtRestGetPhysicalInteractionLength(
			const G4Track& track, G4ForceCondition* condition)
{
	// we are only applied to particlies listed in decaySet.
	if(track.GetVolume() == vol || track.GetNextVolume() == vol) {
		++nDecayed;
		*condition = Forced;
		if(filter->steppingVerbose) 
			printf("particlefilter '%s' forced Decay\n",
					filter->getName().c_str());
		return DBL_MIN;
	}

	*condition = NotForced;
	return DBL_MAX;
}

G4double ParticleFilterPlacement::PostStepGetPhysicalInteractionLength(
			const G4Track& track, G4double previousStepSize,
			G4ForceCondition* condition)
{
	// we are only applied to particlies listed in decaySet.
	if(track.GetVolume() == vol || track.GetNextVolume() == vol) {
		++nDecayed;
		*condition = Forced;
		if(filter->steppingVerbose) 
			printf("particlefilter '%s' forced Decay\n",
					filter->getName().c_str());
		return DBL_MIN;
	}

	*condition = NotForced;
	return DBL_MAX;
}
void ParticleFilterPlacement::UserSteppingAction(const G4Step *step) 
{
	// handle all states

	G4Track *track = step->GetTrack();
	G4VPhysicalVolume *preVol = 
				step->GetPreStepPoint()->GetPhysicalVolume();
	G4VPhysicalVolume *postVol = 
				step->GetPostStepPoint()->GetPhysicalVolume();

	// only handle the step that enters the volume
	if(preVol == postVol || postVol != vol) return;

	// wait until track has hit us nWait times.
	if(BLManager::getObject()->getState() == BEAM) {
		if(--nWait > 0) return;
	} else if(BLManager::getObject()->getState() == VISUAL) {
		if(--nWait > 0) return;
	} else {
		if(--referenceWait > 0) return;
	}

	// we get called for every particle transiting our volume, so we
	// need to check for presence in killSet and keepSet, then require.
	G4ParticleDefinition *def = track->GetDefinition();
	if(filter->killSet.count(def->GetParticleName()) > 0 ||
			(filter->keepSet.size() > 0 &&
		         filter->keepSet.count(def->GetParticleName()) == 0)) {
		track->SetTrackStatus(fStopAndKill);
		++nKilled;
		if(filter->steppingVerbose) 
			printf("particlefilter '%s' killed track\n",filter->getName().c_str());
	} else if(filter->require != "") {
		static BLEvaluator *e=0;
		if(!e)
			e = new BLEvaluator();
		G4ThreeVector position = track->GetPosition();
		e->setVariable("x",position[0]);
		e->setVariable("y",position[1]);
		e->setVariable("z",position[2]);
		e->setVariable("t",track->GetGlobalTime());
		G4ThreeVector mom = track->GetMomentum();
		e->setVariable("Px",mom[0]);
		e->setVariable("Py",mom[1]);
		e->setVariable("Pz",mom[2]);
		e->setVariable("Ptot",sqrt(mom[0]*mom[0]+mom[1]*mom[1]+
							mom[2]*mom[2]));
		e->setVariable("PDGid",def->GetPDGEncoding());
		e->setVariable("EventID",BLManager::getObject()->
				getEventID());
		e->setVariable("TrackID",BLManager::getObject()->
				getExternalTrackID(track));
		e->setVariable("ParentID",BLManager::getObject()->
				getExternalParentID(track));
		if(e->evaluate(filter->require) == 0.0) {
			track->SetTrackStatus(fStopAndKill);
			++nKilled;
			if(filter->steppingVerbose) 
				printf("particlefilter '%s' killed track\n",filter->getName().c_str());
		}
	}
}

void ParticleFilterPlacement::PreUserTrackingAction(const G4Track *track)
{
	nWait = filter->nWait;
	referenceWait = filter->referenceWait;
}
