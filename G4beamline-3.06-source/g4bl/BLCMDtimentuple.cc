//	BLCMDtimentuple.cc

#include <map>
#include <vector>

#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4StepPoint.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"

#include "BLParam.hh"
#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLTrackNTuple.hh"
#include "BLCoordinates.hh"

const char TrackFields[] =
	"x:y:z:Px:Py:Pz:t:PDGid:EventID:TrackID:ParentID:Weight";
const unsigned NTrackFields = 12;


/**	class BLCMDtimentuple - implements a time NTuple
 *	Each command of this class generates an NTuple of the beam at
 *	a specified global time. This is a "perfect" detector in that
 *	it does not perturb the beam at all, and
 *	intrinsically has the resolution of a float (it measures position
 *	and momentum and time).
 *
 *	The userSteppingAction() for the timentuple is called for every
 *	volume, and checks if the step brackets the time specified; if
 *	so it linearly interpolates to the specified time and enters
 *	the result into the NTuple.
 *
 *	Note that if a BLCoordinates instance is linked into the track,
 *	its centerline coordinates are used; otherwise global coordinates
 *	are used.
 **/
class BLCMDtimentuple : public BLCommand {
	G4double time;
	G4String name;
	G4String format;
	G4String filename;
	G4String require;
	G4String coordinates;
	G4int referenceParticle;
	BLCoordinateType coordinateType;
public:
	/// Default constructor.
	BLCMDtimentuple();

	/// Destructor.
	virtual ~BLCMDtimentuple();

	/// Copy constructor.
	BLCMDtimentuple(const BLCMDtimentuple& r);

	/// commandName() returns "timentuple".
	G4String commandName() { return "timentuple"; }

	/// command() implements the timentuple command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// help() prints help text.
	void help(bool detailed) {
		if(description[description.size()-2] == ':')
			description += BLTrackNTuple::getFormatList(); 
		BLCommand::help(detailed);
	}

	void setName(G4String _name) { name = _name; }
};
BLCMDtimentuple defaultTimeNTuple;

/**	class TimeNTuple implements an NTuple for a BLCMDtimentuple.
 **/
class TimeNTuple : public BLManager::SteppingAction {
	G4double time;
	BLTrackNTuple *ntuple;
	G4String require;
	BLCoordinateType coordinateType;
public:
	/// constructor.
	TimeNTuple(G4String _name, G4double _time, G4String format,
			G4String filename, G4String _require,
			BLCoordinateType _coordinateType) : 
						BLManager::SteppingAction()
	{ time = _time; 
	  require = _require;
	  coordinateType = _coordinateType;
	  ntuple = BLTrackNTuple::create(format,"TimeNTuple",_name,
	  		filename,coordinateType,require);
	}

	/// UserSteppingAction() from BLManager::SteppingAction.
	void UserSteppingAction(const G4Step *step);
};

BLCMDtimentuple::BLCMDtimentuple() : BLCommand()
{
	registerCommand(BLCMDTYPE_DATA);
	setSynopsis("Construct an NTuple of tracks at a specified time.");
	setDescription("A time NTuple generates an NTuple of every track at a\n"
		"specified global time.\n"
		"It uses a linear interpolation in the step that straddles "
		"the required time, so accuracy will suffer for large steps. "
		"The NTuple uses centerline coordinates, if available.\n\n"
		"This command is not placed into the geometry.\n\n"
		"The standard NTuple fields are:\n"
		"    x,y,z (mm)\n"
		"    Px,Py,Pz (MeV/c)\n"
		"    t (ns)\n"
		"    PDGid (11=e-, 13=mu-, 22=gamma, 211=pi+, 2212=proton, ...)\n"
		"    EventID (may be inexact above 16,777,215)\n"
		"    TrackID\n"
		"    ParentID (0 => primary particle)\n"
		"    Weight (defaults to 1.0)\n\n"
		"The following additional fields are appended for "
		"format=Extended, format=asciiExtended, and "
		"format=rootExtended:\n"
		"    Bx, By, Bz (Tesla)\n"
		"    Ex, Ey, Ez (Megavolts/meter)\n"
		"    ProperTime (ns)\n"
		"    PathLength (mm)\n"
		"    PolX, PolY, PolZ (polarization)\n"
		"    InitX, initY, InitZ (initial position, mm)\n"
		"    InitT (initial time, ns)\n"
		"    InitKE (MeV when track was created)\n\n"
		"Valid Formats (ignore case): ");
	
	time = 0.0;
	name = "";
	format = "";
	filename = "";
	require = "";
	coordinates = "Centerline";
	referenceParticle = 0;
	coordinateType = BLCOORD_CENTERLINE;
}

BLCMDtimentuple::BLCMDtimentuple(const BLCMDtimentuple& r) : BLCommand(r)
{
	time = r.time;
	name = r.name;
	format = r.format;
	filename = r.filename;
	G4String require;
	coordinates = r.coordinates;
	referenceParticle = r.referenceParticle;
	coordinateType = r.coordinateType;
}

BLCMDtimentuple::~BLCMDtimentuple()
{
}

int BLCMDtimentuple::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("timentuple: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return handleNamedArgs(namedArgs);
	}

	BLCMDtimentuple *t = new BLCMDtimentuple(defaultTimeNTuple);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	t->coordinateType = BLCoordinates::getCoordinateType(t->coordinates);

	// ascii->bltrackfile format, for accuracy and consistency of output
	for(unsigned i=0; i<t->format.size(); ++i)
		t->format[i] = tolower(t->format[i]);
	if(t->format == "ascii")
		t->format = "bltrackfile";

	t->print(argv[0]);

	TimeNTuple *tn = new TimeNTuple(t->name,t->time,t->format,t->filename,
						t->require,t->coordinateType);
	BLManager::getObject()->registerBeamStep(0,tn);
	if(t->referenceParticle != 0)
		BLManager::getObject()->registerReferenceParticleStep(0,tn);

	return retval;
}

void BLCMDtimentuple::defineNamedArgs()
{
	argDouble(time,"time","The global time of the sampling (ns).");
	argString(format,"format","The NTuple format (see above for list).");
	argString(filename,"filename","The filename of the NTuple.");
	argString(filename,"file","Synonym for filename.");
	argString(require,"require","Expression which must be nonzero to include the track (default=1)",false);
	argString(coordinates,"coordinates","Coordinates: global, centerline, or reference (default=c).");
	argInt(referenceParticle,"referenceParticle","Set to 1 to include the Reference Particle.");
}


void TimeNTuple::UserSteppingAction(const G4Step *step)
{
	if(BLManager::getObject()->getState() == SPECIAL) return;

	// only use reference coordinates when they are valid
	BLManagerState state = BLManager::getObject()->getState();
	if(coordinateType == BLCOORD_REFERENCE && state != BEAM) return;

	G4Track *track = step->GetTrack();

	// get basic physical-volume info
	G4StepPoint *prePoint = step->GetPreStepPoint();
	if(!prePoint) return;
	G4StepPoint *postPoint = step->GetPostStepPoint();
	if(!postPoint) return;

	// make a row only if the pre- and post-step times bracket desired time
	// and protect against ntuple==null
	G4double preTime = prePoint->GetGlobalTime();
	G4double postTime = postPoint->GetGlobalTime();
	if(preTime > time || postTime < time || !ntuple) return;
	
	G4ThreeVector position = track->GetPosition();
	G4double time = track->GetGlobalTime();
	G4ThreeVector momentum = track->GetMomentum();

	// Linear interpolation
	G4double dt = step->GetDeltaTime();
	G4double err = this->time - time;
	G4double properTime = track->GetProperTime();
	G4double trackLength = track->GetTrackLength();
	if(dt > 0.001*ns  && fabs(err) <= dt*1.5) {
		time += err;
		position += (err/dt)*step->GetDeltaPosition();
		momentum += (err/dt)*(postPoint->GetMomentum()-
						prePoint->GetMomentum());
		double delta_t = track->GetGlobalTime() - time;
		properTime -= delta_t*track->GetDefinition()->GetPDGMass() /
						track->GetTotalEnergy();
		trackLength -= delta_t*track->GetVelocity();
	}

	ntuple->appendTrack(track,time,position,momentum,properTime,
								trackLength);
}
