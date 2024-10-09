//	BLGlobalField.cc
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

#include <time.h>
#include <cmath>

#include "G4EqMagElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4PropagatorInField.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4SimpleRunge.hh"
#include "G4MagIntegratorDriver.hh"

#include "BLGlobalField.hh"
#include "BLGroup.hh"
#include "BLCommand.hh"
#include "BLParam.hh"

// define parameters:
BLSetParam minStep("minStep","0.01","Minimum step size (mm)");
BLSetParam deltaChord("deltaChord","3.0","Geant4 tracking parameter");
BLSetParam deltaOneStep("deltaOneStep","0.01","Geant4 tracking parameter");
BLSetParam deltaIntersection("deltaIntersection","0.1","Geant4 tracking parameter");
BLSetParam epsMin("epsMin","2.5e-7","Geant4 tracking parameter");
BLSetParam epsMax("epsMax","0.05","Geant4 tracking parameter");
BLSetParam fieldVoxels("fieldVoxels","200,200,200","Size of voxels for field computation (mm)");

bool BLGlobalField::spinTracking = false;
G4EqEMFieldWithSpin *BLGlobalField::spinEqRhs = 0;
BLGlobalField *BLGlobalField::object = 0;

/**	class FieldVoxels is a voxelized field computation.
 **/
class FieldVoxels {
	G4double dx, dy, dz;
	unsigned int nx, ny, nz;
	G4double minX, minY, minZ;
	// voxels is a single-index C array of std::vector-s.
	// use index() to index it as a 3-d array.
	std::vector<const BLElementField*> *voxels;
public:
	FieldVoxels(const std::vector<const BLElementField*>& fields,
				G4double _dx, G4double _dy, G4double _dz);

	// index returns -1 if outside of global bounding box
	int index(int ix, int iy, int iz) const {
			if(ix < 0 || ix >= nx) return -1;
			if(iy < 0 || iy >= ny) return -1;
			if(iz < 0 || iz >= nz) return -1;
			return ix+nx*(iy+ny*iz);
		}
	int index(const G4double point[3]) const {
			int ix = (int)floor((point[0]-minX)/dx);
			int iy = (int)floor((point[1]-minY)/dy);
			int iz = (int)floor((point[2]-minZ)/dz);
			return index(ix,iy,iz);
		}

	int limit(int i, int minv, int maxv) const {
			if(i < minv) return minv;
			if(i > maxv) return maxv;
			return i;
		}

	void addFieldValue(const G4double point[4], G4double field[6]) const;
};

BLGlobalField::BLGlobalField() : G4ElectroMagneticField(), fields()
{
	first = true;
	forceEfieldZero = false;
	forceBfieldZero = false;
	fieldVoxels = 0;

	clear();

	// Get transportation, field, and propagator  managers
	G4TransportationManager  *pTransportMgr= 
               G4TransportationManager::GetTransportationManager();
	G4FieldManager* fieldMgr = pTransportMgr->GetFieldManager();
	G4PropagatorInField *pFieldPropagator =
					pTransportMgr->GetPropagatorInField(); 

	// Need to SetFieldChangesEnergy to account for electric fields
	fieldMgr->SetFieldChangesEnergy(true);

	fieldMgr->SetDetectorField(this);

	// Construct equation of motion of particles through e.m. fields
	G4MagIntegratorStepper *pStepper = 0;
	if(spinTracking) {
		spinEqRhs = new G4EqEMFieldWithSpin(this);
		G4EquationOfMotion *equation = spinEqRhs;
		pStepper = new G4ClassicalRK4(equation,12);
	} else {
		spinEqRhs = 0;
		G4EquationOfMotion *equation = new G4EqMagElectricField(this);
		pStepper = new G4ClassicalRK4(equation,8);
	}
	G4double min_step = Param.getDouble("minStep");
	G4MagInt_Driver *fIntgrDriver = new G4MagInt_Driver(min_step, 
				pStepper, pStepper->GetNumberOfVariables() ); 
	G4ChordFinder *pChordFinder = new G4ChordFinder(fIntgrDriver); 
	fieldMgr->SetChordFinder(pChordFinder); 

	// Set accuracy parameters
	G4double delta_chord = Param.getDouble("deltaChord");
	pChordFinder->SetDeltaChord( delta_chord );
	G4double delta_onestep = Param.getDouble("deltaOneStep");
	fieldMgr->SetAccuraciesWithDeltaOneStep(delta_onestep);
	G4double delta_intersection = Param.getDouble("deltaIntersection");
	fieldMgr->SetDeltaIntersection(delta_intersection); 
	G4double eps_min = Param.getDouble("epsMin");
	G4double eps_max = Param.getDouble("epsMax");
	pFieldPropagator->SetMinimumEpsilonStep(eps_min);
	pFieldPropagator->SetMaximumEpsilonStep(eps_max);
//	G4cout << "Accuracy Parameters:" <<
//    		" MinStep=" << min_step <<
//   		" DeltaChord=" << delta_chord <<
//  		" DeltaOneStep=" << delta_onestep << G4endl;
//	G4cout << "                    " <<
//		" DeltaIntersection=" << delta_intersection <<
//    		" EpsMin=" << eps_min <<
//   		" EpsMax=" << eps_max <<  G4endl;
	fieldMgr->SetChordFinder(pChordFinder);

	// set object
	object = this;
}

BLGlobalField::~BLGlobalField()
{
	clear();
}

BLGlobalField *BLGlobalField::getObject()
{
	if(!object) new BLGlobalField();
	return object;
}

void BLGlobalField::GetFieldValue(const G4double *point, G4double *field) const
{
	// NOTE: this routine often dominates the CPU time for tracking.
	// The test case is the MICE beamline + detector, which has 38 
	// entries in fields[].
	// Using bounding boxes sped it up by a factor of 8 (13 ev/sec => 
	// 100 ev/sec).
	// Using the simple array fp[] instead of fields[] directly sped
	// it up by an additional factor of two (200 ev/sec).
	// Making this a friend of BLElementField and putting the code for
	// isInBoundingBox() inline here didn't help at all.
	// 9/20/2012:
	// Using voxels for a system with 400 magnets sped it up by a factor
	// of 8 (!). So that is now the algorithm, unless there are less than 3
	// fields or the Parameter fieldVoxels=0,0,0.

	field[0] = field[1] = field[2] = field[3] = field[4] = field[5] = 0.0;

	// protect against Geant4 bug that calls us with point[] NaN.
	if(std::isnan(point[0])) return;

	if(first) {
		// clip all bounding boxes to the world
		//@BLGroup *world = BLGroup::getWorld();
		//@double wx = world->getWidth()/2.0;
		//@double wy = world->getHeight()/2.0;
		//@double wz = world->getLength()/2.0;
		G4VSolid *s = BLManager::getObject()->getWorldPhysicalVolume()->
				GetLogicalVolume()->GetSolid();
		G4Box *b = dynamic_cast<G4Box*>(s);
		BLAssert(b != 0);
		double wx = b->GetXHalfLength();
		double wy = b->GetYHalfLength();
		double wz = b->GetZHalfLength();
		bool change=false;
		for(size_t i=0; i<fields.size(); ++i) {
			BLElementField *p = (BLElementField *)fields[i];
			G4double xMin,xMax,yMin,yMax,zMin,zMax;
			p->getBoundingBox(xMin,xMax,yMin,yMax,zMin,zMax);
			if(xMin < -wx) { xMin = -wx; change=true; }
			if(xMax > wx)  { xMax = wx;  change=true; }
			if(yMin < -wy) { yMin = -wy; change=true; }
			if(yMax > wy)  { yMax = wy;  change=true; }
			if(zMin < -wz) { zMin = -wz; change=true; }
			if(zMax > wz)  { zMax = wz;  change=true; }
			if(change) p->setBoundingBox(xMin,xMax,yMin,yMax,zMin,zMax);
		}
		if(change) G4Exception("BLGlobalField",
			"EM Field Extends Outside World",
			JustWarning,
			"May give inaccurate tracking near world boundary.");
		if(fields.size() > 2) { // no voxels for 0,1,2 fields
			std::vector<G4double> list =
				BLCommand::getList(Param.getString("fieldVoxels"),",");
			if(list.size() != 3)
				G4Exception("BLGLobalField",
						"Invalid fieldVoxels",
						FatalException,"");
			if(list[0] != 0.0) // 0 => no voxels
				fieldVoxels = new FieldVoxels(fields,
						list[0],list[1],list[2]);
		}
		first = false;
	}

	if(fieldVoxels) {
		fieldVoxels->addFieldValue(point,field);
	} else {
		for(size_t i=0; i<fields.size(); ++i) {
			const BLElementField *p = fields[i];
			if(p->isInBoundingBox(point)) {
				p->addFieldValue(point,field);
			}
		}
	}

	if(forceBfieldZero) {
		field[0] = field[1] = field[2] = 0.0;
	}
	if(forceEfieldZero) {
		field[3] = field[4] = field[5] = 0.0;
	}
}

void BLGlobalField::clear()
{
	std::vector<const BLElementField*>::iterator i;
	for(i=fields.begin(); i!=fields.end(); ++i) {
		delete *i;
	}
	fields.clear();

	first = true;
}

void BLGlobalField::setSpinTracking(bool v)
{
	if(object) G4Exception("BLGlobalField","Invalid setSpinTracking() call",
				FatalException,"");
	spinTracking = v;
}

FieldVoxels::FieldVoxels(const std::vector<const BLElementField*>& fields,
				G4double _dx, G4double _dy, G4double _dz)
{
	dx = _dx;
	dy = _dy;
	dz = _dz;
	voxels = 0;
	BLAssert(dx>0.0 && dy>0.0 && dz>0.0);

	// get global bounding box around all bounding boxes
	minX = DBL_MAX, minY = DBL_MAX, minZ = DBL_MAX; // class variables
	double maxX = -DBL_MAX, maxY = -DBL_MAX, maxZ = -DBL_MAX; // local
	for(size_t m=0; m<fields.size(); ++m) {
		G4double xMin,xMax,yMin,yMax,zMin,zMax;
		fields[m]->getBoundingBox(xMin,xMax,yMin,yMax,zMin,zMax);
		if(xMin < minX) minX = xMin;
		if(xMax > maxX) maxX = xMax;
		if(yMin < minY) minY = yMin;
		if(yMax > maxY) maxY = yMax;
		if(zMin < minZ) minZ = zMin;
		if(zMax > maxZ) maxZ = zMax;
	}

	// clip to World
	G4VSolid *s = BLManager::getObject()->getWorldPhysicalVolume()->
				GetLogicalVolume()->GetSolid();
	G4Box *b = dynamic_cast<G4Box*>(s);
	BLAssert(b != 0);
	double wx = b->GetXHalfLength();
	double wy = b->GetYHalfLength();
	double wz = b->GetZHalfLength();
	if(minX < -wx) minX = -wx;
	if(maxX > wx) maxX = wx;
	if(minY < -wy) minY = -wy;
	if(maxY > wy) maxY = wy;
	if(minZ < -wz) minZ = -wz;
	if(maxZ > wz) maxZ = wz;

	// move voxel boundaries off even values
	minX -= 0.1*mm;
	maxX += 0.1*mm;
	minY -= 0.1*mm;
	maxY += 0.1*mm;
	minZ -= 0.1*mm;
	maxZ += 0.1*mm;

	// determine nx,ny,nz and update dx,dy,dz
	nx = (int)(maxX-minX)/dx + 1;
	dx = (maxX-minX)/nx;
	ny = (int)(maxY-minY)/dy + 1;
	dy = (maxY-minY)/ny;
	nz = (int)(maxZ-minZ)/dz + 1;
	dz = (maxZ-minZ)/nz;

	// allocate voxels[] (this initializes all of its std::vector-s)
	if(nx>50000 || ny>50000 || nz>50000 || nx*ny*nz > 100000000) {
		char tmp[256];
		sprintf(tmp,"nx=%u ny=%u nz=%u",nx,ny,nz);
		G4Exception("BLGLobalField", "Too many fieldVoxels",
						FatalException,tmp);
	}
	voxels = new std::vector<const BLElementField*>[nx*ny*nz];
	BLAssert(voxels);

	// fill each voxel with the BLElementField-s that intersect it
	for(size_t m=0; m<fields.size(); ++m) {
		G4double xMin,xMax,yMin,yMax,zMin,zMax;
		G4int ixMin,ixMax,iyMin,iyMax,izMin,izMax;
		fields[m]->getBoundingBox(xMin,xMax,yMin,yMax,zMin,zMax);
		ixMin = limit((int)floor((xMin-minX)/dx),0,nx-1);
		ixMax = limit((int)floor((xMax-minX)/dx),0,nx-1);
		iyMin = limit((int)floor((yMin-minY)/dy),0,ny-1);
		iyMax = limit((int)floor((yMax-minY)/dy),0,ny-1);
		izMin = limit((int)floor((zMin-minZ)/dz),0,nz-1);
		izMax = limit((int)floor((zMax-minZ)/dz),0,nz-1);
		for(int i=ixMin; i<=ixMax; ++i) {
		    for(int j=iyMin; j<=iyMax; ++j) {
			for(int k=izMin; k<=izMax; ++k) {
			    int ind = index(i,j,k);
			    BLAssert(ind >= 0);
			    voxels[ind].push_back(fields[m]);
			}
		    }
		}
	}

	// scan voxels and print statistics/summary
	int nVoxels=nx*ny*nz;
	int nMax=0, nG5=0, nG10=0, nG20=0;
	for(int i=0; i<nVoxels; ++i) {
		int n = voxels[i].size();
		if(n > nMax) nMax = n;
		if(n > 5) ++nG5;
		if(n > 10) ++nG10;
		if(n > 20) ++nG20;
	}
	printf("fieldVoxels: nx=%d ny=%d nz=%d, %d voxels, %ld fields\n",
				nx,ny,nz,nVoxels,fields.size());
	printf("fieldVoxels: max field count is %d fields,"
				" # voxels >5: %d, >10: %d, >20: %d\n",
				nMax,nG5,nG10,nG20);
	if(nVoxels > 1000000 || nMax > 3) {
		char tmp[256];
		sprintf(tmp,"nVoxels = %d, max fields/voxel = %d",
							nVoxels,nMax);
		G4Exception("BLGlobalField","Check number of field voxels",
				JustWarning, tmp);
	}
}

void FieldVoxels::addFieldValue(const G4double point[4],G4double field[6]) const
{
	int j = index(point);
	if(j < 0) return; // outside global bounding box

	std::vector<const BLElementField*> &p = voxels[j];
	size_t n = p.size();
	for(size_t i=0; i<n; ++i) {
		if(!p[i]->isInBoundingBox(point)) continue;
		p[i]->addFieldValue(point,field);
	}
}
