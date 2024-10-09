//	BLZStep.hh
//#define ZSTEPWORLD

#ifndef BLZSTEP_HH
#define BLZSTEP_HH

#include <vector>
#include <map>

#include "globals.hh" 
#include "G4VUserParallelWorld.hh" 
#include "G4Box.hh" 
#include "G4LogicalVolume.hh" 
#include "G4PVPlacement.hh" 
#include "G4TransportationManager.hh" 
#include "G4Navigator.hh" 


#ifdef ZSTEPWORLD

struct ZStep {
	G4double z;
	void *action; // always BLManager::ZSteppingAction 
	int when;
	G4VPhysicalVolume *pv;
	ZStep(G4double _z, void *a, int w)
		{ z=_z; action=a; when=w; pv=0; }
};

class BLZStep : public G4VUserParallelWorld {
	G4Navigator *navigator; // set in setupProcesses() -- delay needed
	std::vector<ZStep> &zStepVector;
	std::map<G4VPhysicalVolume*,int> &zStepMap;
	G4LogicalVolume *lvworld;
public: 
	// Constructor. Gets vector and map from BLManager.
	BLZStep(std::vector<ZStep> &v, std::map<G4VPhysicalVolume*,int> &m) : 
					G4VUserParallelWorld("ZStepWorld"),
					zStepVector(v), zStepMap(m) { 
		navigator = 0;
		lvworld = 0;
		printf("ZStepWorld initialized\n");
	}

	virtual ~BLZStep() { }

	/// getVolume() returns the current physical volume in ZStepWorld.
	G4VPhysicalVolume *getVolume() {
		G4TouchableHistory *th = navigator->CreateTouchableHistory();
		G4VPhysicalVolume *v = th->GetVolume();
		delete th;
		return v;
	}

	/// Construct() from G4VUserParallelWorld. Creates and places volumes
	/// for every Z step that BLManager has.
	virtual void Construct();

	/// setupProcesses() sets up the necessary physics processes
	void setupProcesses();

	/// placeVolume() places one volume into ZStepWorld
	G4VPhysicalVolume *placeVolume(G4double z);
};

#endif // ZSTEPWORLD
#endif // BLZSTEP_HH
