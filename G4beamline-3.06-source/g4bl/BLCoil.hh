//	BLCoil.hh
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

#ifndef BLCOIL_HH
#define BLCOIL_HH

#include <map>

#include "globals.hh"
#include "BLCommand.hh"

/**	BLCoil creates and implements a coil (part of a solenoid)
 *
 *	A BLCoil consists only of the parameters describing the coil geometry
 *	and a field map for unit current density. It is not a BLElement, and
 *	cannot be placed geometrically -- see BLSolenoid for that.
 *
 *	Cylindrical symmetry around z is assumed. The field is approximated
 *	by computing the field from nSheets infinitesimally-thin current
 *	sheets interpolating from innerRadius to OuterRadius, each length
 *	long. The solution to Maxwell's Equations for each sheet is used.
 *
 *	The field has been verified to be correct within 0.2% for both a
 *	long solenoid and for a pair of Helmholtz coils.
 *
 *	Effect of nSheets on accuracy: consider the coil innerRadius=330.0*mm,
 *	outerRadius=505.0*mm, length=167.0*mm. Compute the maximum error
 *	for 1000 random points, excluding the region within 2 cm of the coils;
 *	use 99 sheets as the "correct" value. In the region 0<=r<2*outerRadius,
 *	0<z<10*outerRadius, for 0.1% accuracy one needs 11 sheets, and for 1%
 *	one needs 4 sheets. In the region 0<=r<innerRadius/2, same z, for 0.1%
 *	accuracy one needs 5 sheets, and for 1% one needs 2 sheets. Note that
 *	inside the coils this approximation of using thin sheets is not very
 *	good; but interpolating on a grid and coordinating sheet positions
 *	with the grid is probably about as good as one can do -- at least the
 *	values will be smooth (which is not true for the bare computation --
 *	there's a singularity at every sheet, and Bz jumps from +infinity
 *	to -infinity across the sheet). 
 *
 *	Note that Bz does not decrease to 0.1% of its central value until
 *	z>10*(innerRadius+outerRadius)/2; it decreases to 1% at ~5 times that
 *	average radius.
 *
 *
 *	BLCoil can also read a mapFile defining the field (rather than doing
 *	the above computation). For instance, it could be from opera-2d.
 *	The mapFile gives tables of Bz and Br on a grid in the z-r plane.
 *	mapFile format:
 *		lines beginning with # are comments, and are ignored.
 *		The first line contains the parameters in name=value form:
 *			z0	initial value of z (mm)
 *			dz	grid interval along z (mm)
 *			nz	number of grid points along z
 *			dr	grid interval along r (mm)
 *			nr	number of grid points along r
 *			reflectZ set nonzero to reflect in the z=0 plane
 *				(inverting Br)
 *			current	the current for this field (A/mm^2)
 *		All parameters are required. If reflectZ is nonzero then
 *		z0 must be 0.0 (only z>=0 is contained in the file).
 *		The point with local coordinates (z=0,r=0) is what gets placed.
 *
 *		Then follow the two maps in any order. Each map begins
 *		with a line giving its name (Bz or Br) in column 1, followed
 *		by nz lines each containing nr floats (Tesla). The first
 *		value in each line is for r=0. Values are separated by one
 *		or more spaces (the first value on each line may but need not
 *		be preceeded by spaces).
 *	The main reason to use a mapFile is to accurately include the effects
 *	of iron. That implies that the map may not be accurate for values
 *	of current other than that for which the map was generated.
 *
 *	See BLCMDfieldmap for a more general EM field read from a map-file.
 *
 *
 *	BLCoil was separated from BLSolenoid so multiple solenoids with the
 *	same geometry (coil) could share a single (large) field map for the
 *	coil -- the B field is scaled by current/currentOfMap, where 
 *	currentOfMap is 1.0 for the computed field, and is the current 
 *	parameter from a mapFile.
 *
 *	This class does not use BLFieldMap becuase it predates it, and
 *	because the map must include the parameters so construct() will
 *	know whether or not the map corresponds to the current values.
 *	You can, of course, build a solenoid out of a BLTubs and a
 *	BLFieldMap.
 **/
class BLCoil {
	static std::map<G4String,BLCoil*> mapCoil;
	G4String name;
	G4double innerRadius;
	G4double outerRadius;
	G4double length;
	G4int nSheets;
	G4String material;
	G4double tolerance;
	G4double maxR;
	G4double minZ;
	G4double maxZ;
	G4double dR;
	G4double dZ;
	G4int nR;
	G4int nZ;
	G4String filename;
	G4String mapFile;
	G4int exactComputation;
	float *mapBr;
	float *mapBz;
	G4double sheetR0;
	G4double sheetDR;
	G4double norm;
	G4bool goodMap;
	G4bool reflectZ;
	G4double currentOfMap;
	friend class BLCMDcoil;
	friend class BLCMDsolenoid;
	friend class SolenoidField;
public:
	/// Default constructor.
	BLCoil();

	/// Complex Constructor. All arguments that default to 0 are related
	/// to generating the field map; those that are 0 will be determined
	/// by tolerance. To construct the map the minimum info is tolerance
	/// and maxR; if either of them is 0, no map will be generated or used.
	/// NOTE: without a map, field values inside the coils will fluctuate
	/// wildly, because points on the sheets are singular; map construction
	/// is designed to minimize this. Without a map, nSheets must be given.
	/// tolerance is a fraction of the best value at r=0,z=0; values
	/// smaller than about 0.001 are unrealistic (and expensive!).
	BLCoil(G4String _name, G4double _innerRadius, G4double _outerRadius,
		G4double _length, G4int _nSheets=0, G4double _tolerance=0.0,
		G4double _maxR=0.0, G4double _maxZ=0.0, G4double _dR=0.0,
		G4double _dZ=0.0, G4int _nR=0, G4int _nZ=0, 
		G4String _filename="") {
			name=_name; innerRadius=_innerRadius;
			outerRadius=_outerRadius; length=_length;
			nSheets=_nSheets; tolerance=_tolerance; maxR=_maxR;
			maxZ=_maxZ; dR=_dR; dZ=_dZ; nR=_nR; nZ=_nZ; 
			filename = _filename;
			exactComputation = 0;
			mapBr=0; mapBz=0; sheetR0=0.0; sheetDR=0.0; norm=0.0; 
			goodMap=false; 
			if(nSheets > 0) { 
				sheetDR=(outerRadius-innerRadius)/nSheets;
				sheetR0=innerRadius + sheetDR/2.0;
			}
			if(filename == "") filename = name + ".dat";
		}

	/// Destructor.
	virtual ~BLCoil();

	/// Copy constructor.
	BLCoil(const BLCoil& r);

	/// printCoil() will print a description
	void printCoil();

	/// find() returns a pointer to the named BLCoil.
	static BLCoil *find(G4String name);

	/// addField() adds the field value for the given point[] into field[].
	/// Note that point[] is relative to the center of the coil.
	/// Units for current are Amp/mm^2.
	void addField(const G4double point[4], G4double field[6],
					G4double currentAmpPerMm2) const;

	/// generateFieldMap() will generate the field map (for execution
	/// speedup). If just maxR and tolerance were given, all other
	/// map parameters will be automatically determined so that the
	/// estimated error is less than the specified tolerance.
	/// The field map is cached in a disk file which is automatically
	/// recreated if anything changes (may take several minutes).
	void generateFieldMap();

	/// readMapFile() will read the mapFile.
	void readMapFile();

	/// displayMapValues() will display the 4 map values around a point.
	void displayMapValues(G4double r, G4double z);
private:
	/// getUnitField() will return Br and Bz for unit current.
	void getUnitField(G4double r, G4double z, G4double& Br, G4double& Bz);

	/// determineMaxR() will determine maxR based on tolerance.
	/// Uses nSheets=99 to determine center field, and walks out
	/// far enough so (Bz/Bz0)<0.5*tolerance.
	void determineMaxR();

	/// determineMaxZ() will determine maxZ based on tolerance.
	/// Uses nSheets=99 to determine on-axis field, and walks out
	/// far enough so (Bz/Bz0)<0.5*tolerance.
	void determineMaxZ();

	/// determineNsheets() will determine nSheets, given maxZ and 
	/// tolerance. Increases nSheets until maxError < tolerance/4
	/// (the other 3/4 of tolerance is for mapping errors). Checks
	/// just a handful of points selected to be representative.
	void determineNsheets();

	/// estimateMapError() will estimate the maximum error in mapping,
	/// given maxZ, nSheets, dR, dZ, nR, nZ, and an allocated but
	/// uninitialized fieldMap. It generates 1000 random points and
	/// returns the maximum error; quits early if tolerance is
	/// exceeded. Error is (currentMap - best)/best(r=0,z=0), where 
	/// best is computed with no map and nSheets=99. For each point
	/// it generates the 4 surrounding values of fieldMap (this is faster
	/// than generating the entire map first).
	/// Returns the maximum error of either Bz or Br, normalized to
	/// Bz(r=0,z=0,nSheets=99).
	G4double estimateMapError();

	/// computeError() computes the error for a single point.
	/// Error is (currentMap - best)/best(r=0,z=0), where best should be
	/// an identical BLCoil with no map and nSheets=99. If makeMap is true,
	/// it generates the 4 surrounding values of the fieldMap.
	G4double computeError(G4double r, G4double z, BLCoil &best,
							G4bool makeMap);

	/// getSheetField() returns the B field from an infinitesimally-thin
	/// sheet of current. Arguments:
	/// r and z are the cylindrical coordinates of the point in question
	/// relative to center of the coil (axis along z), mm.
	/// Br and Bz return the computed B field at that point (geant4 units)
	/// sheetR is the radius of the current sheet (mm)
	/// sheetLen is the length of the current sheet (mm)
 	/// currentAmpPerMm is the current of the sheet (Amp per linear mm)
	void getSheetField(G4double r, G4double z, G4double& Br, G4double& Bz,
		G4double sheetR, G4double sheetLen, G4double currentAmpPerMm)
									const;
};

#endif // BLCOIL_HH
