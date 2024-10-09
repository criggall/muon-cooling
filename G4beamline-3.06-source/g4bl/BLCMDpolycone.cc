//	BLCMDpolycone.cc
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

#include <vector>

#include "G4VisAttributes.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLAssert.hh"
#include "BLGroupElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLKillTrack.hh"

/**	BLCMDpolycone implements a polycone of material, axis along Z.
 *
 **/
class BLCMDpolycone : public BLGroupElement {
	G4String innerRadius;
	G4String outerRadius;
	G4String z;
	G4double initialPhi;
	G4double finalPhi;
	G4String material;
	G4String color;
	G4int kill;
	G4double maxStep;
	G4Polycone *polycone;
	G4double length;
	G4double maxRadius;
	std::vector<G4double> irs;
	std::vector<G4double> ors;
	std::vector<G4double> zs;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDpolycone();

	/// Destructor.
	virtual ~BLCMDpolycone() { }

	/// Copy constructor.
	BLCMDpolycone(const BLCMDpolycone& r);

	/// clone()
	BLElement *clone() { return new BLCMDpolycone(*this); }

	/// commandName() returns "polycone".
	G4String commandName() { return "polycone"; }

	/// command() implements the polycone command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// setup() will handle the vector arguments.
	void setup();

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	virtual G4VSolid *getSolid();

	/// construct() - construct the tube/cylinder.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the tube/cylinder.
	G4double getLength() { return length; }

	/// getWidth() returns the width of the tube/cylinder.
	G4double getWidth() { return 2.0*maxRadius; }

	/// getHeight() returns the height of the tube/cylinder.
	G4double getHeight() { return 2.0*maxRadius; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,zs[0]);
		if(index == 1) return G4ThreeVector(0.0,0.0,zs[zs.size()-1]);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// isOutside() from BLElement.
	bool isOutside(G4ThreeVector &local, G4double tolerance) 
		{ BLAssert(polycone!=0);  return polycone->Inside(local) != kInside; }

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v)
		{ v.clear();
		  for(int n=0; n<npoints*5; ++n)
		  	v.push_back(polycone->GetPointOnSurface());
		}

	/// isWithin() from BLGroupElement.
	bool isWithin(G4ThreeVector &local, G4double tolerance) 
	    { BLAssert(polycone!=0);  return polycone->Inside(local) != kOutside; }
};

BLCMDpolycone defaultPolycone;	// default object

// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDpolycone::BLCMDpolycone() : BLGroupElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a polycone with axis along z");
	setDescription("This is a direct interface to G4Polycone.\n"
		"For a solid polycone, omit innerRadius and it will be\n"
		"filled with zeroes. The number of entries in z, innerRadius,\n"
		"and outerRadius must be the same. Note that a polycone is\n"
		"placed at its z=0,r=0 point, which need not be its\n"
		"geometric center.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.");

	// provide initial values for fields
	innerRadius = "";
	outerRadius = "";
	z = "";
	initialPhi = 0.0;
	finalPhi = 360.0*deg;
	material = "";
	color = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	polycone = 0;
	length = 0.0;
	maxRadius = 0.0;
}

// Copy constructor - be sure to use the copy constructor BLGroupElement(r)
BLCMDpolycone::BLCMDpolycone(const BLCMDpolycone& r) : BLGroupElement(r),
					irs(r.irs), ors(r.ors), zs(r.zs)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	innerRadius = r.innerRadius;
	outerRadius = r.outerRadius;
	z = r.z;
	initialPhi = r.initialPhi;
	finalPhi = r.finalPhi;
	material = r.material;
	color = r.color;
	kill = r.kill;
	maxStep = r.maxStep;
	polycone = r.polycone;
	length = r.length;
	maxRadius = r.maxRadius;
}

int BLCMDpolycone::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("polycone: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultPolycone.handleNamedArgs(namedArgs);
	}

	BLCMDpolycone *t = new BLCMDpolycone(defaultPolycone);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	t->setup();

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	return retval;
}

void BLCMDpolycone::defineNamedArgs()
{
	argString(innerRadius,"innerRadius","Comma-separated list of inner radii (mm)",false);
	argString(outerRadius,"outerRadius","Comma-separated list of outer radii (mm)",false);
	argString(z,"z","Comma-separated list of z positions (mm)",false);
	argDouble(initialPhi,"initialPhi","The initial Phi value (deg; 0 for all)",deg);
	argDouble(finalPhi,"finalPhi","The final Phi value (deg; 360 for all)",deg);
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)",mm);
	argString(material,"material","The material of the polycone");
	argString(color,"color","The color of the polycone (''=invisible)");
	argInt(kill,"kill","Set nonzero to kill every track that enters.");
}

void BLCMDpolycone::setup()
{
	// parse list arguments into the vectors
	zs = getList(z,',');
	if(zs.size() == 0)
		printError("polycone: invalid z list");
	irs = getList(innerRadius,",");
	if(irs.size() == 0 && innerRadius.size() > 0)
		printError("polycone: invalid innerRadius list");
	ors = getList(outerRadius,",");
	if(ors.size() == 0)
		printError("polycone: invalid outerRadius list");

	// if innerRadius=="", fill irs[] with zeroes
	if(irs.size() == 0) {
		for(unsigned i=0; i<zs.size(); ++i)
			irs.push_back(0.0);
	}

	// check valididty of the lists
	if(zs.size()!=ors.size() || zs.size()!=irs.size())
		printError("polycone: inconsistent array sizes");

	// compute object size and check z order
	G4double prevZ=-DBL_MAX, minZ=+DBL_MAX, maxZ=-DBL_MAX, maxR=0.0;
	bool ok=true;
	for(unsigned i=0; i<zs.size(); ++i) {
		double z=zs[i];
		if(z < prevZ) ok = false;
		prevZ = z;
		if(minZ > z) minZ = z;
		if(maxZ < z) maxZ = z;
		double r=ors[i];
		if(maxR < r) maxR = r;
	}
	if(!ok) printError("polycone: invalid order of z values");

	maxRadius = maxR;
	length = 2.0*(fabs(minZ)>fabs(maxZ) ? fabs(minZ) : fabs(maxZ));
/***
	unsigned int place, next;
	for(place=next=0; next<z.size(); place=next+1) {
		next = z.find(",",place);
		G4String p;
		if(next < z.size())
			p = z.substr(place,next-place);
		else
			p = z.substr(place);
		if(p.size() == 0) break;
		G4double z = atof(p.c_str());
		zs.push_back(z);
		if(minZ > z) minZ = z;
		if(maxZ < z) maxZ = z;
	}
	// get innerRadius positions into irs[]
	for(place=next=0; next<innerRadius.size(); place=next+1) {
		next = innerRadius.find(",",place);
		G4String p;
		if(next < innerRadius.size())
			p = innerRadius.substr(place,next-place);
		else
			p = innerRadius.substr(place);
		if(p.size() == 0) break;
		irs.push_back(atof(p.c_str()));
	}
	// get outerRadius positions into ors[]
	for(place=next=0; next<outerRadius.size(); place=next+1) {
		next = outerRadius.find(",",place);
		G4String p;
		if(next < outerRadius.size())
			p = outerRadius.substr(place,next-place);
		else
			p = outerRadius.substr(place);
		if(p.size() == 0) break;
		G4double r = atof(p.c_str());
		ors.push_back(r);
		if(maxR < r) maxR = r;
	}
***/
}

G4VSolid *BLCMDpolycone::getSolid()
{
	if(!polycone) {
		polycone = new G4Polycone("Polycone", 
				initialPhi,finalPhi-initialPhi,zs.size(),
				&zs[0],&irs[0],&ors[0]);
	}
	return polycone;
}

void BLCMDpolycone::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	getSolid(); // sets polycone

	G4Material *mat;
	if(material != "")
		mat = getMaterial(material);
	else
		mat = parent->GetMaterial();

	G4LogicalVolume *lv = new G4LogicalVolume(polycone,mat, thisname+"LogVol");
	lv->SetVisAttributes(getVisAttrib(color));
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");
	lv->SetUserLimits(new G4UserLimits(maxStep));

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4VPhysicalVolume *pv = new G4PVPlacement(g4rot, relativePosition,lv,
					thisname, parent,false,0,surfaceCheck);

	if(kill)
		BLManager::getObject()->
			registerSteppingAction(pv,new BLKillTrack(thisname));

	// get globalRotation and globalPosition
	G4RotationMatrix *globalRotation = 0;
	if(relativeRotation && parentRotation) {
		globalRotation = 
		    new G4RotationMatrix(*parentRotation * *relativeRotation);
	} else if(relativeRotation) {
		globalRotation = relativeRotation;
	} else if(parentRotation) {
		globalRotation = parentRotation;
	}
	G4ThreeVector globalPosition(relativePosition + parentPosition);
	if(parentRotation)
		globalPosition = *parentRotation * relativePosition +
				parentPosition;

	constructChildren(lv,thisname,globalRotation,globalPosition);

	printf("BLCMDpolycone::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f kill=%d\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2],
		globalPosition[2]-getLength()/2.0,
		globalPosition[2]+getLength()/2.0,
		kill);
}

