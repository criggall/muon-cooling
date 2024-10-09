//	BLCMDgenericbend.cc
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

#include "G4VisAttributes.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"

#include "BLGroupElement.hh"
#include "BLElementField.hh"
#include "BLGlobalField.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLEngeFunction.hh"
#include "BLTune.hh"
#include "BLKillTrack.hh"

#include "Randomize.hh"
#define rand G4UniformRand

const G4double FRINGE_ACCURACY=1.0e-4;


/**	BLCMDgenericbend implements a generic bending magnet with a box field.
 *
 *	The magnetic field is uniform throught the box, along the local Y 
 *	axis, modified by a fringe field computation.
 **/
class BLCMDgenericbend : public BLGroupElement {
	G4double fieldWidth;
	G4double fieldHeight;
	G4double fieldLength;
	G4double ironWidth;
	G4double ironHeight;
	G4double ironLength;
	G4double By;
	G4String fieldMaterial;
	G4String fieldColor;
	G4String ironMaterial;
	G4String ironColor;
	G4int kill;
	G4double maxStep;
	G4String fringe;
	G4double fringeFactor;
	G4int openAperture;
	G4Box *fieldBox;
	G4Box *ironBox;
	G4VSolid *ironSolid;
	friend class GenericBendField;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDgenericbend();

	/// Destructor.
	virtual ~BLCMDgenericbend() { }

	/// Copy constructor.
	BLCMDgenericbend(const BLCMDgenericbend& r);

	/// clone()
	BLElement *clone() { return new BLCMDgenericbend(*this); }

	/// commandName() returns "genericbend".
	G4String commandName() { return "genericbend"; }

	/// command() implements the genericbend command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	/// construct() - construct the generic bending magnet
	void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the bend
	G4double getLength() 
		{ return ironLength>fieldLength ? ironLength : fieldLength; }

	/// getWidth() returns the width of the bend
	G4double getWidth() 
		{ return ironWidth; }

	/// getHeight() returns the height of the bend
	G4double getHeight() 
		{ return ironHeight; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,-getLength()/2.0);
		if(index == 1) return G4ThreeVector(0.0,0.0,getLength()/2.0);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return openAperture==0 || getNChildren()==0; }

	/// isOutside() from BLElement.
	bool isOutside(G4ThreeVector &local, G4double tolerance);

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isWithin() from BLGroupElement.
	bool isWithin(G4ThreeVector &local, G4double tolerance);
};

BLCMDgenericbend defaultGenericBend;	// default object

/**	GenericBendField represents one placement of a generic bending magnet.
 *
 **/
class GenericBendField : public BLElementField {
	G4double halfwidth;
	G4double halfheight;
	G4double halflength;
	G4double *By;
	BLCoordinateTransform global2local;
	G4RotationMatrix rotation;
	G4double fringeMaxZ;
	G4double fringeDepth;
	BLEngeFunction enge;
public:
	/// constructor. 
	GenericBendField(BLCoordinateTransform& _global2local, BLCMDgenericbend *bend);

	/// addFieldValue() adds the field for this solenoid into field[].
	/// point[] is in global coordinates.
	void addFieldValue(const G4double point[4], G4double field[6]) const;
};


// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDgenericbend::BLCMDgenericbend() : BLGroupElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a generic bending magnet.");
	setDescription("The field region is a box with By specified.\n"
		"A fringe field computation based on the method of COSY\n"
		"INFINITY is included by default, extending the field in a "
		"rectangle extending the straight aperture along the local z.\n"
		"This is first order only, and assumes the magnet is "
		"infinitely wide; the fringe field extends outside of "
		"the magnet aperture only in a region extending the "
		"aperture in x and y. As the fringe field is first order "
		"only, it is slightly non-Maxwellian. It is calculated using "
		"Enge functions.\n\n"
		"Note that there is no field inside the 'iron'; this can "
		"result in gross tracking errors for particles in the iron, "
		"and implies that kill=1 is desirable.\n\n"
		"By default, the aperture is filled with a box volume "
		"of the fieldMaterial; this prevents placing any object "
		"inside the aperture. With openAperture=1 no aperture volume "
		"is used, and objects can be placed into the parent volume "
		"that are inside the aperture.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.\n\n"
		"If ironLength<=0, ironWidth<=0, or ironHeight<=0, no iron "
		"is used.");

	// provide initial values for fields
	fieldWidth = 0.0;
	fieldHeight = 0.0;
	fieldLength = 0.0;
	ironWidth = 0.0;
	ironHeight = 0.0;
	ironLength = 0.0;
	By = 0.0;
	fieldMaterial = "Vacuum";
	fieldColor = "";
	ironMaterial = "Fe";
	ironColor = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	fringe = "";
	fringeFactor = 1.0;
	openAperture = 0;
	fieldBox = 0;
	ironBox = 0;
	ironSolid = 0;
}

// Copy constructor - be sure to use the copy constructor BLGroupElement(r)
BLCMDgenericbend::BLCMDgenericbend(const BLCMDgenericbend& r) : BLGroupElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	fieldWidth = r.fieldWidth;
	fieldHeight = r.fieldHeight;
	fieldLength = r.fieldLength;
	ironWidth = r.ironWidth;
	ironHeight = r.ironHeight;
	ironLength = r.ironLength;
	BLTune::copyTunableArg(&By,&r.By);
	fieldMaterial = r.fieldMaterial;
	fieldColor = r.fieldColor;
	ironMaterial = r.ironMaterial;
	ironColor = r.ironColor;
	kill = r.kill;
	maxStep = r.maxStep;
	fringe = r.fringe;
	fringeFactor = r.fringeFactor;
	openAperture = r.openAperture;
	fieldBox = r.fieldBox;
	ironBox = r.ironBox;
	ironSolid = r.ironSolid;
}

int BLCMDgenericbend::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("genericbend: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultGenericBend.handleNamedArgs(namedArgs);
	}

	BLCMDgenericbend *t = new BLCMDgenericbend(defaultGenericBend);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	getMaterial(t->fieldMaterial);
	getMaterial(t->ironMaterial);

	t->print(argv[0]);

	return retval;
}

void BLCMDgenericbend::defineNamedArgs()
{
	argDouble(fieldWidth,"fieldWidth","The width of the field region (mm)");
	argDouble(fieldHeight,"fieldHeight","The height of the field region (mm)");
	argDouble(fieldLength,"fieldLength","The length of the field region (mm)");
	argDouble(ironWidth,"ironWidth","The width of the iron region (mm)");
	argDouble(ironHeight,"ironHeight","The height of the iron region (mm)");
	argDouble(ironLength,"ironLength","The length of the iron region (mm)");
	argTunable(By,"By","The magnetic field (Tesla)",tesla);
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)");
	argString(fieldMaterial,"fieldMaterial","The material of the field region.");
	argString(fieldColor,"fieldColor","The color of the field region.");
	argString(ironMaterial,"ironMaterial","The material of the iron region.");
	argString(ironColor,"ironColor","The color of the iron region.");
	argInt(kill,"kill","Set nonzero to kill particles hitting the iron.");
	argString(fringe,"fringe","Fringe field computation, set to 0 to disable,"
		" or a comma-separated list of 6 Enge function parameters.");
	argDouble(fringeFactor,"fringeFactor","Fringe depth factor (1.0).");
	argInt(openAperture,"openAperture","Set nonzero to omit the aperture volume.",false);
}

void BLCMDgenericbend::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)

{
	G4String thisname = parentName+getName();

	if(!fieldBox) {
	    double len = (ironLength>fieldLength ? ironLength : fieldLength);
	    fieldBox = new G4Box(thisname+"GenericBendField",
					fieldWidth/2.0, fieldHeight/2.0,
					len/2.0+0.1);
	}
	G4LogicalVolume *lvfield=0;
	if(!openAperture) {
	    G4Material *mat = getMaterial(fieldMaterial);
	    lvfield = new G4LogicalVolume(fieldBox,mat,
					thisname+"LogVol");
	    lvfield->SetVisAttributes(getVisAttrib(fieldColor));
	    if(maxStep <= 0.0) maxStep = Param.getDouble("maxStep");
	    lvfield->SetUserLimits(new G4UserLimits(maxStep));
	}

	G4LogicalVolume *lviron=0;
	if(ironHeight > 0.0 && ironWidth > 0 && ironLength > 0.0) {
	    if(!ironBox)
		ironBox = new G4Box(thisname+"GenericBendIron",
					ironWidth/2.0, ironHeight/2.0,
					ironLength/2.0);
	    if(!ironSolid)
		ironSolid = new G4SubtractionSolid(thisname+"GenericBendSolid",
					ironBox, fieldBox);
	    G4Material *mat = getMaterial(ironMaterial);
	    lviron = new G4LogicalVolume(ironSolid,mat,
					thisname+"LogVol");
	    lviron->SetVisAttributes(getVisAttrib(ironColor));
	    lviron->SetUserLimits(new G4UserLimits(maxStep));
	}

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	if(lvfield)
		new G4PVPlacement(g4rot,relativePosition,lvfield,thisname,
					parent,false,0,surfaceCheck);
	G4PVPlacement *pv=0;
	if(lviron)
		pv = new G4PVPlacement(g4rot,relativePosition,lviron,
					thisname, parent,false,0,surfaceCheck);

	if(kill && pv)
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

	BLCoordinateTransform global2local(globalRotation,globalPosition);

	G4double zmin = globalPosition[2]-getLength()/2.0;
	G4double zmax = globalPosition[2]+getLength()/2.0;

	GenericBendField *p = new GenericBendField(global2local,this);
	BLGlobalField::getObject()->addElementField(p);

	if(lvfield != 0) {
		constructChildren(lvfield,getName(),globalRotation,globalPosition);
	} else if(getNChildren() != 0) {
		printError("genericbend: cannot be a parent with openAperture=1");
	}

	printf("BLCMDgenericbend::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2], zmin,zmax);
}

GenericBendField::GenericBendField(BLCoordinateTransform& _global2local,
				BLCMDgenericbend *bend) :
				BLElementField(), rotation(), enge(ENGE_BEND)
{
	halfwidth = bend->fieldWidth/2.0;
	halfheight = bend->fieldHeight/2.0;
	halflength = bend->fieldLength/2.0;
	By = &bend->By;
	global2local = _global2local;
	rotation = global2local.getRotation().inverse();
	fringeDepth = bend->fringeFactor * halfheight * 2.0;

	if(bend->fringe != "") {
		std::vector<G4double> v = BLCommand::getList(bend->fringe,',');
		if(v.size() == 1 && v[0] == 0.0)
			enge.set(0,0,0,0,0,0);
		else if(v.size() == 6)
			enge.set(v[0],v[1],v[2],v[3],v[4],v[5]);
		else
			BLCommand::printError("Invalid fringe value '%s'\n",
						bend->fringe.c_str());
	}

	if(enge.getType() == ENGE_BLOCK) {
		fringeMaxZ = halflength;
	} else {
		for(int i=0; i<1000; ++i) {
			fringeMaxZ = i*halfwidth/5.0 + halflength;
			if(enge((fringeMaxZ-halflength)/fringeDepth) <
							FRINGE_ACCURACY)
				break;
		}
	}

	// set global bounding box
	G4double local[4], global[4];
	local[3] = 0.0;
	for(int i=0; i<2; ++i) {
		local[0] = (i==0 ? -1.0 : 1.0) * halfwidth;
		for(int j=0; j<2; ++j) {
			local[1] = (j==0 ? -1.0 : 1.0) * halfheight;
			for(int k=0; k<2; ++k) {
				local[2] = (k==0 ? -1.0 : 1.0) * fringeMaxZ;
				global2local.getGlobal(local,global);
				setGlobalPoint(global);
			}
		}
	}
}

void GenericBendField::addFieldValue(const G4double point[4], G4double field[6])
								const
{
	G4ThreeVector global(point[0],point[1],point[2]);
	G4ThreeVector local;

	global2local.getLocal(local,global);
	if(fabs(local[0]) > halfwidth || fabs(local[1]) > halfheight ||
	   fabs(local[2]) > fringeMaxZ)
		return;

	/* Enge function as applied in Cosy Infinity */
	double fringeZ = (fabs(local[2])-halflength)/fringeDepth;
	G4double e = enge(fringeZ);
	G4double e1 = enge.prime(fringeZ)/fringeDepth;
	G4double y = local[1];
	G4double B0 = *By;
/* 1,2,3 makes almost no visible difference in tracking; 
   the fields do differ visibly. */
#define COMPUTE_TERMS 1
#if COMPUTE_TERMS==1
	G4ThreeVector B(0.0,B0*e,B0*e1*y);
#elif COMPUTE_TERMS==2
	G4double e2 = enge.second(fringeZ)/fringeDepth/fringeDepth;
	G4double e3 = enge.third(fringeZ)/fringeDepth/fringeDepth/fringeDepth;
	G4double y2 = y*y;
	G4double y3 = y*y2;
	G4ThreeVector B(0.0,B0*(e-0.5*e2*y2),B0*(e1*y-y3*e3/6.0));
#elif COMPUTE_TERMS==3
	G4double e2 = enge.second(fringeZ)/fringeDepth/fringeDepth;
	G4double e3 = enge.third(fringeZ)/fringeDepth/fringeDepth/fringeDepth;
	G4double e4 = enge.fourth(fringeZ)/fringeDepth/fringeDepth/fringeDepth/fringeDepth;
	G4double e5 = enge.fifth(fringeZ)/fringeDepth/fringeDepth/fringeDepth/fringeDepth/fringeDepth;
	G4double y2 = y*y;
	G4double y3 = y*y2;
	G4double y4 = y*y3;
	G4double y5 = y*y4;
	G4ThreeVector B(0.0,B0*(e-0.5*e2*y2+y4*e4/24.0),
					B0*(e1*y-y3*e3/6.0+y5*e5/120.0));
#else
#error invalid value of COMPUTE_TERMS
#endif
	if(local[2] < 0.0) B[2] = -B[2];

	if(global2local.isRotated())
		B = rotation * B;

	field[0] += B[0];
	field[1] += B[1];
	field[2] += B[2];
}

G4bool BLCMDgenericbend::isOutside(G4ThreeVector &local, G4double tolerance)
{
	// if any iron dimension <= 0, always returns true.
	if(openAperture)
		return fabs(local[0]) > ironWidth/2.0-tolerance ||
			fabs(local[1]) > ironHeight/2.0-tolerance ||
			fabs(local[2]) > ironLength/2.0-tolerance ||
			fabs(local[0]) < fieldWidth/2.0+tolerance ||
			fabs(local[1]) < fieldHeight/2.0+tolerance;
	return fabs(local[0]) > ironWidth/2.0-tolerance ||
		fabs(local[1]) > ironHeight/2.0-tolerance ||
		fabs(local[2]) > ironLength/2.0-tolerance;
}

void BLCMDgenericbend::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	v.clear();
	// if any iron dimension <= 0, generates no points (no solid to check)
	if(ironHeight <= 0.0 || ironWidth <= 0.0 || ironHeight <= 0.0)
		return;
	if(openAperture) {
		// last two points are 0,0,l and 0,0,-l.
		generateBox(0,ironWidth,ironHeight,ironLength,v);
		v.pop_back();
		v.pop_back();
		G4double w=fieldWidth/2.0, h=fieldHeight/2.0, l=fieldLength/2.0;
		// the 8 inside corners
		v.push_back(G4ThreeVector(w,h,l));
		v.push_back(G4ThreeVector(w,h,-l));
		v.push_back(G4ThreeVector(w,-h,l));
		v.push_back(G4ThreeVector(w,-h,-l));
		v.push_back(G4ThreeVector(-w,h,l));
		v.push_back(G4ThreeVector(-w,h,-l));
		v.push_back(G4ThreeVector(-w,-h,l));
		v.push_back(G4ThreeVector(-w,-h,-l));
		// the 12 inside edge midpoints
		v.push_back(G4ThreeVector(0,h,l));
		v.push_back(G4ThreeVector(0,h,-l));
		v.push_back(G4ThreeVector(0,-h,l));
		v.push_back(G4ThreeVector(0,-h,-l));
		v.push_back(G4ThreeVector(w,0,l));
		v.push_back(G4ThreeVector(w,0,-l));
		v.push_back(G4ThreeVector(-w,0,l));
		v.push_back(G4ThreeVector(-w,0,-l));
		v.push_back(G4ThreeVector(w,h,0));
		v.push_back(G4ThreeVector(w,-h,0));
		v.push_back(G4ThreeVector(-w,h,0));
		v.push_back(G4ThreeVector(-w,-h,0));
		// the 4 inside face centers
		v.push_back(G4ThreeVector(w,0,0));
		v.push_back(G4ThreeVector(-w,0,0));
		v.push_back(G4ThreeVector(0,h,0));
		v.push_back(G4ThreeVector(0,-h,0));
		// random points on the surfaces
		G4double wi=ironWidth/2.0, hi=ironHeight/2.0, li=ironLength/2.0;
		while(v.size() < (unsigned)npoints) {
again:			G4double x=0,y=0,z=0;
			switch(v.size() % 10) {
			case 0: x=wi;			// left outer face
				y=rand()*ironHeight-hi;
				z=rand()*ironLength-li;
				break;
			case 1: x=-wi;			// right outer face
				y=rand()*ironHeight-hi;
				z=rand()*ironLength-li;
				break;
			case 2: x=rand()*ironWidth-wi;	// upper outer face
				y=hi;
				z=rand()*ironLength-li;
				break;
			case 3: x=rand()*ironWidth-wi;	// lower outer face
				y=-hi;
				z=rand()*ironLength-li;
				break;
			case 4: x=rand()*ironWidth-wi;	// downstream face
				y=rand()*ironHeight-hi;
				z=li;
				goto check;
			case 5: x=rand()*ironWidth-wi;	// upstream face
				y=rand()*ironHeight-hi;
				z=-li;
				// omit points in the aperture
check:				if(w < wi-10.0*mm && h < hi-10.0*mm && 
				   (fabs(x) < w && fabs(y) < h)) goto again;
				break;
			case 6: x=w;			// left inner face
				y=rand()*fieldHeight-h;
				z=rand()*ironLength-li;
				break;
			case 7: x=-w;			// right inner face
				y=rand()*fieldHeight-h;
				z=rand()*ironLength-li;
				break;
			case 8: x=rand()*fieldWidth-w;	// upper inner face
				y=h;
				z=rand()*ironLength-li;
				break;
			case 9: x=rand()*fieldWidth-w;	// lower inner face
				y=-h;
				z=rand()*ironLength-li;
				break;
			}
			v.push_back(G4ThreeVector(x,y,z));
		}
	} else {
		generateBox(npoints,ironWidth,ironHeight,ironLength,v);
	}
}

G4bool BLCMDgenericbend::isWithin(G4ThreeVector &local, G4double tolerance)
{
	if(ironHeight <= 0.0 || ironWidth <= 0.0 || ironHeight <= 0.0)
		return false;
	if(openAperture)
		return fabs(local[0]) < ironWidth/2.0+tolerance &&
			fabs(local[1]) < ironHeight/2.0+tolerance &&
			fabs(local[2]) < ironLength/2.0+tolerance &&
			fabs(local[0]) > fieldWidth/2.0-tolerance &&
			fabs(local[1]) > fieldHeight/2.0-tolerance;
	return fabs(local[0]) < ironWidth/2.0+tolerance &&
		fabs(local[1]) < ironHeight/2.0+tolerance &&
		fabs(local[2]) < ironLength/2.0+tolerance;
}
