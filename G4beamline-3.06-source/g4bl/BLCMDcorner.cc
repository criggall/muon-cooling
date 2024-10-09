//	BLCMDcorner.cc
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

#include "G4UserLimits.hh"
#include "G4StepPoint.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4Track.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"

#include "BLParam.hh"
#include "BLManager.hh"
#include "BLElement.hh"
#include "BLGroup.hh"
#include "BLCoordinates.hh"
#include "BLCoordinateTransform.hh"


/**	class BLCMDcorner - implements a corner in the centerline.
 *	
 *	In addition to adding a corner to BLCoordinates, this class can be
 *	a BLElement that adds the corner angle to every track that enters it.
 **/
class BLCMDcorner : public BLElement {
	G4double radius;
	G4double height;
	G4double width;
	G4double length;
	G4String material;
	G4String color;
	G4VSolid *solid;
	G4double maxStep;
	G4double z;
	G4String rotation;
	G4RotationMatrix rotationMatrix;
	G4RotationMatrix centerlineRotation;
	G4double radiusCut;
	friend class BLCornerTrackRotation;
public:
	/// Default constructor.
	BLCMDcorner();

	/// Destructor.
	~BLCMDcorner();

	/// Copy constructor.
	BLCMDcorner(const BLCMDcorner& r);

	/// clone() - invalid.
	BLElement *clone() { return 0; }

	/// commandName() returns "corner".
	G4String commandName() { return "corner"; }

	/// command() implements the corner command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() will construct the corner.
	/// Used for normal placements of a Corner object.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns this element's Length along the Z axis.
	G4double getLength() { return length; }

	/// getWidth() returns this element's Width along the X axis.
	G4double getWidth() { return width; }

	/// getHeight() returns this element's height along the Y axis.
	G4double getHeight() { return height; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// isOutside() from BLElement.
	bool isOutside(G4ThreeVector &local, G4double tolerance)
		{ return true; }

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v)
		{ v.clear(); }
};

BLCMDcorner defaultCorner;

/**	class BLCornerTrackRotation implements a track rotation for a BLCMDcorner.
 **/
class BLCornerTrackRotation : public BLManager::SteppingAction {
	G4String name;
	G4VPhysicalVolume *thisVol;
	G4RotationMatrix rotation;
public:
	/// constructor.
	BLCornerTrackRotation(G4String _name, G4VPhysicalVolume *pv,
						G4RotationMatrix& rot)
	{ name=_name; thisVol = pv; rotation = rot; }

	/// UserSteppingAction() from BLManager::SteppingAction.
	void UserSteppingAction(const G4Step *step);
};

BLCMDcorner::BLCMDcorner() : BLElement(), rotation()
{
	registerCommand(BLCMDTYPE_LAYOUT);
	setSynopsis("Implement a corner in the centerline.");
	setDescription("The centerline is bent by a rotation.\n"
		"Every track that enters the volume also gets rotated.\n"
		"The z value is for the corner and the front face of the "
		"volume (if any).\n"
		"If the corner is paired with a bending magnet or other "
		"mechanism to bend the beam, no volume should be used.\n\n"
		"NOTE: This command is self-placing, do not use the place\n"
		"command; it also affects all following placements, and\n"
		"it cannot be issued inside a group.\n"
		"If radius=height=width=0 then no volume is associated\n"
		"with the corner, and a bending magnet should be placed\n"
		"nearby to bend the particles around the corner. Normally "
		"the bending magnet is placed before the corner, and is "
		"rotated by half the bend angle. For a sector bend, it's "
		"usually best to use the cornerarc command rather than this "
		"one.\n"
		"\nNOTE: all placements before this command must have z\n"
		"values before the corner, and all placements after this\n"
		"command must have z values after the corner. Note also\n"
		"that the angle is limited to 90 degrees.\n\n"
		"Note that the radiusCut is important to reduce or eliminate "
		"ambiguities in the global to centerline coordinate transform. "
		"It can also be used to 'shield' the beamline to prevent "
		"particles from taking unusual paths around the outside "
		"of beamline elements.\n\n"
		"This command places itself into the geometry.");
	
	// initial field values
	radius = 0.0;
	height = 0.0;
	width = 0.0;
	length = 1.0*mm;
	material = "Vacuum";
	color = "";
	solid = 0;
	maxStep = -99.0;
	z = 0.0;
	rotation = "";
	radiusCut = 0.0;
}

BLCMDcorner::~BLCMDcorner()
{
	if(solid) delete solid;
}

BLCMDcorner::BLCMDcorner(const BLCMDcorner& r) : BLElement(r), rotation(r.rotation)
{
	radius = r.radius;
	height = r.height;
	width = r.width;
	length = r.length;
	material = r.material;
	color = r.color;
	solid = 0;
	maxStep = r.maxStep;
	z = r.z;
	rotation = r.rotation;
	radiusCut = r.radiusCut;
}
int BLCMDcorner::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() >= 1 && argv[0] == "default") {
		return defaultCorner.handleNamedArgs(namedArgs);
	}

	if(BLGroup::getWorld() != BLGroup::getCurrent()) {
		printError("corner: cannot be issued inside a group.");
		return -1;
	}

	BLCMDcorner *t = new BLCMDcorner(defaultCorner);
	if(argv.size() >= 1) t->setName(argv[0]);
	t->z = -DBL_MAX;
	int retval = t->handleNamedArgs(namedArgs);

	if(argv.size() != 1 && (t->radius > 0.0 || t->height > 0.0)) {
		printError("corner: Invalid command, must have name");
		return -1;
	}

	if(t->z == -DBL_MAX) {
		printError("corner: z must be specified");
		return -1;
	}
	if(t->maxStep <= 0.0) t->maxStep = Param.getDouble("maxStep");
	if(t->radiusCut <= 0.0) t->radiusCut = 
					BLCoordinates::getCurrentRadiusCut();
	if(t->radiusCut <= 0.0)
		G4Exception("corner command","No radius cut",JustWarning, "");

	// get the rotationMatrix and the current centerline rotation
	t->rotationMatrix = *stringToRotationMatrix(t->rotation);
	t->centerlineRotation = *BLCoordinates::getCurrentRotation();

	// check for acute angle
	G4ThreeVector tmp(0.0,0.0,1.0);
	tmp = t->rotationMatrix * tmp;
	if(tmp[2] < 0.0)
		printError("corner: angle must be <= 90 degrees"); 

	G4ThreeVector cl_coord(0.0,0.0,t->z+t->length/2.0);

	BLGroup::getWorld()->placeElement(t,cl_coord,t->getName());
	BLCoordinates::corner(t->z,&t->rotationMatrix,t->radiusCut);

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(t->getName());

	return retval;
}

void BLCMDcorner::defineNamedArgs()
{
	argDouble(z,"z","The centerline Z of the corner (mm).",mm);
	argString(rotation,"rotation","The rotation of the corner (see above).");
	argDouble(radiusCut,"radiusCut","The radius cut for this following segment (mm default=previous).",mm);
	argDouble(radius,"radius","The radius of the circular corner volume (mm).",mm);
	argDouble(height,"height","The height of the rectangular corner volume (mm).",mm);
	argDouble(width,"width","The width of the rectangular corner volume (mm).",mm);
	argDouble(length,"length","The length of the corner volume (mm).",mm);
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm).",mm);
	argString(material,"material","The material of the corner volume.");
	argString(color,"color","The color of the corner volume (''=invisible).");

	if(radius > 0.0) width = height = radius*2.0;
}

void BLCMDcorner::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	if(!solid) {
		if(radius > 0.0) {
			solid = new G4Tubs(thisname+"Tubs", 0.0, radius,
					length/2.0, 0.0, 2.0*pi);
			height = width = 2.0*radius;
		} else if(height > 0.0 && width > 0.0) {
			solid = new G4Box(thisname+"Box",width/2.0,
					height/2.0,length/2.0);
		} else {
			// corner for a bending magnet
			return;
		}
	}
	G4Material *mat = getMaterial(material);
	G4LogicalVolume *lv = new G4LogicalVolume(solid,mat,thisname+"LogVol");
	lv->SetVisAttributes(getVisAttrib(color));
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");
	lv->SetUserLimits(new G4UserLimits(maxStep));

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4VPhysicalVolume *pv = new G4PVPlacement(g4rot,relativePosition,lv,
					thisname,parent,false,0,surfaceCheck);

	// BLCornerTrackRotation needs the rotation wrt the global coordinate
	// axes, so rotate our "local" rotationMatrix to the global axes
	G4RotationMatrix rot = centerlineRotation *
				rotationMatrix *
				centerlineRotation.inverse();

	BLCornerTrackRotation *cr = new BLCornerTrackRotation(thisname,pv,rot);
	BLManager::getObject()->registerSteppingAction(pv,cr);

	printf("BLCMDcorner::Construct %s parent=%s relZ=%.1f globZ=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		parentPosition[2]);
}


void BLCornerTrackRotation::UserSteppingAction(const G4Step *step)
{
	// process all states

	// get basic physical-volume info
	G4StepPoint *prePoint = step->GetPreStepPoint();
	if(!prePoint) return;
	G4VPhysicalVolume *preVol = prePoint->GetPhysicalVolume();
	if(!preVol) return;
	G4StepPoint *postPoint = step->GetPostStepPoint();
	if(!postPoint) return;
	G4VPhysicalVolume *postVol = postPoint->GetPhysicalVolume();
	if(!postVol) return;
	
	// return if not entering thisVol
	if(preVol == postVol || postVol != thisVol) return;

	// get more info
	G4Track *track = step->GetTrack();

	// rotate the track
	G4ThreeVector dir(track->GetMomentumDirection());
	dir = rotation * dir;
	((G4Track *)track)->SetMomentumDirection(dir);
	postPoint->SetMomentumDirection(dir);
	((G4Step *)step)->UpdateTrack();
}
