//	BLGlobalField.hh
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

#ifndef BLGLOBALFIELD_HH
#define BLGLOBALFIELD_HH

#include <vector>
#include "G4ElectroMagneticField.hh"
#include "G4EqEMFieldWithSpin.hh"

#include "BLElementField.hh"

class FieldVoxels;	// declared and implemented in BLGlobalField.cc

/**	BLGlobalField - handles the global ElectroMagnetic field
 *
 *	There is a single BLGlobalField object.
 *
 *	The field from each individual beamline element is given by a
 *	BLElementField object. Any number of overlapping BLElementField
 *	objects can be added to the global field. Any element command that
 *	represents an element with an EM field must add the appropriate
 *	BLElementField to the global BLGlobalField object.
 *
 *	For efficiency, as the search through the bounding boxes can be
 *	expensive with hundreds of individual BLFieldElement-s, 3-d
 *	voxels are used to speed up the search.
 **/
class BLGlobalField : public G4ElectroMagneticField {
protected:
	std::vector<const BLElementField *> fields;
	mutable bool first;
	bool forceEfieldZero;
	bool forceBfieldZero;
	mutable FieldVoxels *fieldVoxels;
	static bool spinTracking;
	static G4EqEMFieldWithSpin *spinEqRhs;
	static BLGlobalField *object;
	BLGlobalField();
	~BLGlobalField();
	BLGlobalField(const BLGlobalField&);	// undefined
	BLGlobalField& operator=(const BLGlobalField&);	// undefined
public:
	/// getObject() returns the single BLGlobalField object.
	/// It is constructed, if necessary.
	static BLGlobalField *getObject();

	/// GetFieldValue() returns the field value at a given point[].
	/// field is really field[6]: Bx,By,Bz,Ex,Ey,Ez.
	/// point[] is in global coordinates: x,y,z,t.
	void GetFieldValue(const G4double *point, G4double *field) const;

	void getFieldValue(G4ThreeVector &pos, G4double t, G4ThreeVector &B,
							G4ThreeVector &E) {
		G4double point[4], field[6];
		point[0]=pos[0]; point[1]=pos[1]; point[2]=pos[2],point[3]=t;
		GetFieldValue(point,field);
		B.setX(field[0]); B.setY(field[1]); B.setZ(field[2]);
		E.setX(field[3]); E.setY(field[4]); E.setZ(field[5]);
	}

	/// DoesFieldChangeEnergy() returns true.
	G4bool DoesFieldChangeEnergy() const { return true; }

	/// addElementField() adds the BLElementField object for a single
	/// element to the global field.
	void addElementField(const BLElementField* f)
		{ fields.push_back(f); }

	/// clear() removes all BLElementField-s from the global object,
	/// and destroys them. Used before the geometry is completely
	/// re-created.
	void clear();

	/// zeroEfield() sets the flag that forces the Efield to be zero
	/// everywhere. Used for ICOOL-style reference particle tracking.
	void zeroEfield(bool v) { forceEfieldZero = v; }

	/// zeroBfield() sets the flag that forces the Bfield to be zero
	/// everywhere. Used for ICOOL-style reference particle tracking.
	void zeroBfield(bool v) { forceBfieldZero = v; }

	/// setSpinTracking() sets a flag. MUST be called before the
	/// BLGlobalField object is created.
	static void setSpinTracking(bool v=true);

	/// setSpinAnomaly() sets the value of the spin anomaly.
	static void setSpinAnomaly(double v) 
		{ if(spinEqRhs != 0) spinEqRhs->SetAnomaly(v); }
};

#endif // BLGLOBALFIELD_HH
