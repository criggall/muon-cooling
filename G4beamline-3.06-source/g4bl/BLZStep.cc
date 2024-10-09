//	BLZStep.cc

#include "G4TransportationManager.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParallelWorldScoringProcess.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"

#include "BLAssert.hh"
#include "BLZStep.hh"
#include "BLCoordinates.hh"
#include "BLCommand.hh"
#include "BLManager.hh"

#ifdef ZSTEPWORLD

const double THICK=0.001*mm;	// thickness of solids (along z)

void BLZStep::Construct() {
	printf("BLZStep::Construct() entered\n");

	// get the ZStepWorld physical volume
	G4VPhysicalVolume *world = GetWorld(); 
}

void BLZStep::setupProcesses()
{
	// get the ZStepWorld physical volume
	G4VPhysicalVolume *world = GetWorld(); 

	// get the navigator in ZStepWorld
	navigator = G4TransportationManager::GetTransportationManager()->
						GetNavigator("ZStepWorld");
	BLAssert(navigator != 0);

	// omit everything if no Z steps (bookends are 2 entries)
	if(zStepVector.size() <= 2) return;

	// construct the volumes in ZStepWorld
	// Note that zStepVector and zStepMap were initialized from the
	// objects with the same names in BLManager.
	lvworld = world->GetLogicalVolume(); 
	double prevZ=-DBL_MAX;
	G4VPhysicalVolume *prevV=0;
	for(unsigned i=0; i<zStepVector.size(); ++i) {
		if(zStepVector[i].action == 0) continue;
		if(fabs(zStepVector[i].z-prevZ) <= THICK) {
			zStepVector[i].pv = prevV;
			continue;
		}
		prevZ = zStepVector[i].z;
		prevV = placeVolume(prevZ);
		zStepVector[i].pv = prevV;
		zStepMap[prevV] = i;
	}

	// Add parallel world scoring process 
	G4ParticleTable::G4PTblDicIterator *myParticleIterator =
			G4ParticleTable::GetParticleTable()->GetIterator();
	G4ParallelWorldScoringProcess* proc =
			new G4ParallelWorldScoringProcess("ZStepLimiter"); 
	proc->SetParallelWorld("ZStepWorld"); 
	myParticleIterator->reset(); 
	while( (*myParticleIterator)() ) { 
	    G4ParticleDefinition* particle = myParticleIterator->value(); 
	    if (!particle->IsShortLived()) { 
		G4ProcessManager* pmanager = particle->GetProcessManager(); 
		pmanager->AddProcess(proc); 
		pmanager->SetProcessOrderingToLast(proc, idxAtRest); 
		pmanager->SetProcessOrdering(proc, idxAlongStep, 1); 
		pmanager->SetProcessOrderingToLast(proc, idxPostStep); 
	    } 
	} 
	printf("G4ParallelWorldScoringProcess added to all particles\n");
}


G4VPhysicalVolume *BLZStep::placeVolume(G4double z)
{
	BLAssert(lvworld != 0);

	static double worldSize=0.0;
	if(worldSize <= 0.0) { 
		G4Box *worldBox = dynamic_cast<G4Box*>(lvworld->GetSolid());
		BLAssert(worldBox != 0);
		worldSize = worldBox->GetXHalfLength();
		if(worldBox->GetYHalfLength() > worldSize)
			worldSize = worldBox->GetYHalfLength();
		if(worldBox->GetZHalfLength() > worldSize)
			worldSize = worldBox->GetZHalfLength();
		worldSize *= 4; // guaranteed larger in any orientation
	}

	char name[32];
	sprintf(name,"Z%.0f",z);
	G4VSolid *solid=0;
	G4ThreeVector cl(0,0,z+THICK/2.0), global;
	G4double rc=BLCoordinates::getRadiusCut(cl);
	G4RotationMatrix *rot=BLCoordinates::getGlobalRotation(cl);
	BLCoordinates::getGlobalAnywhere(cl,global);

	// if no radiusCut, use a G4Box that is bigger than the world
	// (extending outside the world is OK for this)
	if(rc <= 0.0)
		solid = new G4Box(name, worldSize, worldSize, THICK/2.0); 
	else
		solid = new G4Tubs(name,0.0,rc,THICK,0.0,360*deg);

	G4LogicalVolume *lv = new G4LogicalVolume(solid,0,name,0,0,0); 

/*** also place a visible object into the mass world
lv->SetMaterial(BLCommand::getMaterial("Vacuum"));
lv->SetVisAttributes(BLCommand::getVisAttrib("0,1,1"));
G4LogicalVolume *lvmass = BLManager::getObject()->getWorldPhysicalVolume()->GetLogicalVolume();
BLAssert(lvmass!=0);
new G4PVPlacement(rot, global, lv, name, lvmass, 0, 0); 
// cannot share LogicalVOlume between worlds
lv = new G4LogicalVolume(solid,0,name,0,0,0); 
***/

//printf("BLZStep '%s'  at z=%.3f mm\n",name,z);
	return new G4PVPlacement(rot, global, lv, name, lvworld, 0, 0,surfaceCheck); 
}

#else // ZSTEPWORLD
int zstepworld_dummy=0;
#endif // ZSTEPWORLD
