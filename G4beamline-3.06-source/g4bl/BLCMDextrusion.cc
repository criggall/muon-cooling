//	BLCMDextrusion.cc
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
#include "G4ExtrudedSolid.hh"
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

/**	BLCMDextrusion implements an extrusion of material, axis along Z.
 *
 **/
class BLCMDextrusion : public BLGroupElement {
	G4double length;
	G4String vertices;
	G4double scale1;
	G4double scale2;
	G4String material;
	G4String color;
	G4int kill;
	G4double maxStep;
	G4ExtrudedSolid *extrusion;
	std::vector<G4TwoVector> polygon;
	G4double height;
	G4double width;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDextrusion();

	/// Destructor.
	virtual ~BLCMDextrusion() { }

	/// Copy constructor.
	BLCMDextrusion(const BLCMDextrusion& r);

	/// clone()
	BLElement *clone() { return new BLCMDextrusion(*this); }

	/// commandName() returns "extrusion".
	G4String commandName() { return "extrusion"; }

	/// command() implements the extrusion command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// setup() will handle the vertices, and determine height and width.
	void setup();

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	G4VSolid *getSolid();

	/// construct() - construct the extrusion
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the extrusion
	G4double getLength() { return length; }

	/// getWidth() returns the width of the extrusion
	G4double getWidth() { return width; }

	/// getHeight() returns the height of the extrusion
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
		{ BLAssert(extrusion!=0);  return extrusion->Inside(local) != kInside; }

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v)
		{ v.clear();
		  for(int n=0; n<npoints*5; ++n) {
			G4ThreeVector p=extrusion->GetPointOnSurface();
			if(extrusion->Inside(p) == kSurface)
		  		v.push_back(p);
		  }
		}

	/// isWithin() from BLGroupElement.
	bool isWithin(G4ThreeVector &local, G4double tolerance) 
	    { BLAssert(extrusion!=0);  return extrusion->Inside(local) != kOutside; }
};

BLCMDextrusion defaultExtrusion;	// default object

// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDextrusion::BLCMDextrusion() : BLGroupElement(), polygon()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a solid extrusion with axis along z");
	setDescription("This is a basic interface to G4ExtrudedSolid.\n"
		"A simple polygon in the X-Y plane is extruded along z, "
		"with optional scales in XY at the two ends (which generate "
		"a linear scaling along z).\n\n"
		"The polygon must be simple (no two sides intersect, no two "
		"vertices are equal). The vertices are listed starting from "
		"any vertex and traversing the polygon in either direction "
		"without lifting the pencil from the paper (Geant4 requires "
		"the traversal to be clockwise but this element internally "
		"reverses it if required). For an N-sided polygon give N "
		"vertices -- a side will be added from last to "
		"first to close the polygon; N is determined by "
		"counting the entries in the vertices argument.\n\n"
		"Note that while you cannot make an extrusion with a hole, you "
		"can make such an object in two parts or by placing a "
		"daughter volume in this one.\n\n"
		"Note the position placed is x=0,y=0,z=0, which is centered "
		"along z, but need not be near the center of the polygon in "
		"XY.\n\n"
		"With scale1!=scale2 this is not really an extrusion; by "
		"making one of them 0.001 or so, you can construct a sharp "
		"apex.\n\n"
		"Any x or y value in vertices can be an expression using "
		"double constants and the usual C operators and functions.\n\n"
		"This element must be placed (via the place command), and "
		"children can be placed inside it.");

	// provide initial values for fields
	length = 0.0;
	vertices = "";
	scale1 = 1.0;
	scale2 = 1.0;
	material = "";
	color = "1,1,1";
	kill = 0;
	maxStep = -1;
	extrusion = 0;
	height = 0.0;
	width = 0.0;
}

// Copy constructor
BLCMDextrusion::BLCMDextrusion(const BLCMDextrusion& r) : BLGroupElement(r),
							polygon(r.polygon)
{
	length = r.length;
	vertices = r.vertices;
	scale1 = r.scale1;
	scale2 = r.scale2;
	material = r.material;
	color = r.color;
	kill = r.kill;
	maxStep = r.maxStep;
	extrusion = 0;
	height = r.height;
	width = r.width;
}

int BLCMDextrusion::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	if(argv.size() != 1) {
		printError("extrusion: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultExtrusion.handleNamedArgs(namedArgs);
	}

	BLCMDextrusion *t = new BLCMDextrusion(defaultExtrusion);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	t->setup();

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	t->print(argv[0]);

	return retval;
}

void BLCMDextrusion::defineNamedArgs()
{
	argDouble(length,"length","Length of the extrusion (mm).",1.0,"",false);
	argString(vertices,"vertices","List of vertices of the XY polygon "
		"(mm): 'x0,y0;x1,y1;...'; a line from last to first "
		"is added. A 2 mm square is: "
		"'-1,-1;-1,1;1,1;1,-1'",false);
	argDouble(scale1,"scale1","The XY scale at the upstream (-z) end (1.0).",1.0,"",false);
	argDouble(scale2,"scale2","The XY scale at the downstream (+z) end (1.0).",1.0,"",false);
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)",mm);
	argString(material,"material","The material of the extrusion");
	argString(color,"color","The color of the extrusion (''=invisible)");
	argInt(kill,"kill","Set nonzero to kill every track that enters.");
	argString(vertices,"vertexes","Synonym for vertices.",false);
}

void BLCMDextrusion::setup()
{
	width = 0.0;
	height = 0.0;

	// get vertices into polygon[]
	polygon.clear();
	std::vector<G4String> v=splitString(vertices,";",true);
	for(unsigned i=0; i<v.size(); ++i) {
		if(v[i].size() == 0) continue;
		std::vector<G4double> p=getList(v[i],",");
		if(p.size() != 2) {
			printError("Syntax error in vertices");
			polygon.clear();
			break;
		}
		G4TwoVector point(p[0],p[1]);
		polygon.push_back(point);
		if(width < fabs(p[0])) width = fabs(p[0]);
		if(height < fabs(p[1])) height = fabs(p[1]);
	}
	if(polygon.size() < 3) {
		printError("extrusion: polygon has Fewer than three vertices");
		return;
	}

	// Compute the signed area to determine if it is traversed clockwise or 
	// ccw; invert the order if ccw because G4ExtrudedSolid requires cw.
	// (This is known as the "surveyor's formula".)
	double area=0.0;
	polygon.push_back(polygon[0]); // extra point for simplicity
	for(unsigned i=0; i<polygon.size()-1; ++i) {
		area += polygon[i].x()*polygon[i+1].y() -
			polygon[i+1].x()*polygon[i].y(); // omit the 1/2
	}
	polygon.pop_back(); // remove extra point
	if(area > 0.0) {	// invert the order of the vertices
		std::vector<G4TwoVector> tmp;
		for(unsigned i=0; i<polygon.size(); ++i)
			tmp.push_back(polygon[polygon.size()-1-i]);
		for(unsigned i=0; i<polygon.size(); ++i)
			polygon[i] = tmp[i];
	}

	width *= 2 * (scale1>scale2 ? scale1 : scale2);
	height *= 2 * (scale1>scale2 ? scale1 : scale2);
	return;
}

G4VSolid *BLCMDextrusion::getSolid()
{
	if(!extrusion) {
		G4TwoVector offset(0.0,0.0);
		extrusion = new G4ExtrudedSolid("Extrusion", 
				polygon,length/2.0,offset,scale1,offset,scale2);
	}
	return extrusion;
}

void BLCMDextrusion::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	getSolid(); // sets extrusion

	G4Material *mat=0;
	if(material.size() == 0) 
		mat = parent->GetMaterial();
	else
		mat = getMaterial(material);
	G4LogicalVolume *lv = new G4LogicalVolume(extrusion,mat, thisname+"LogVol");
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

	printf("BLCMDextrusion::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f kill=%d\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2],
		globalPosition[2]-getLength()/2.0,
		globalPosition[2]+getLength()/2.0,
		kill);
}

