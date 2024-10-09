//	BLCMDtessellatedsolid.cc
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

#include <set>

#include "G4VisAttributes.hh"
#include "G4TessellatedSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Color.hh"
#include "G4UserLimits.hh"
#include "G4Material.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include "G4VVisManager.hh"
#include "G4Polymarker.hh"

#include "BLCommandAlias.hh"
#include "BLGroupElement.hh"
#include "BLParam.hh"
#include "BLManager.hh"
#include "BLKillTrack.hh"

#define MARKER_SIZE 5	/* pixels */

/**	BLCMDtessellatedsolid implements a tessellatedsolid of material.
 *
 **/
class BLCMDtessellatedsolid : public BLGroupElement, 
						public BLManager::RunAction  {
	G4String filename;
	G4String material;
	G4String color;
	G4int kill;
	G4double maxStep;
	G4int debug;
	G4double height;
	G4double width;
	G4double length;
	G4double Zmin, Zmax;
	G4TessellatedSolid *tessellatedsolid;
	G4VPhysicalVolume *physVol;
public:
	/// Default constructor. Defines the command, args, etc.
	BLCMDtessellatedsolid();

	/// Destructor.
	virtual ~BLCMDtessellatedsolid() { }

	/// Copy constructor.
	BLCMDtessellatedsolid(const BLCMDtessellatedsolid& r);

	/// clone()
	BLElement *clone() { return new BLCMDtessellatedsolid(*this); }

	/// commandName() returns "tessellatedsolid".
	G4String commandName() { return "tessellatedsolid"; }

	/// command() implements the tessellatedsolid command.
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	/// defineNamedArgs() defines the named arguments for the command.
	void defineNamedArgs();

	virtual G4VSolid *getSolid();

	/// construct() - construct the tessellatedsolid.
	virtual void construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition);

	/// getLength() returns the length of the tessellatedsolid.
	G4double getLength() { return length; }

	/// getWidth() returns the outer radius of the tessellatedsolid.
	G4double getWidth() { return width; }

	/// getHeight() returns the outer radius of the tessellatedsolid.
	G4double getHeight() { return height; }

	/// getSurveyPoint() returns points in LOCAL coordinates.
	G4ThreeVector getSurveyPoint(int index) {
		if(index == 0) return G4ThreeVector(0.0,0.0,Zmin);
		if(index == 1) return G4ThreeVector(0.0,0.0,Zmax);
		throw "UNIMPLEMENTED";
	}

	/// isOK() returns true.
	G4bool isOK() { return true; }

	/// isOutside() from BLElement.
	bool isOutside(G4ThreeVector &local, G4double tolerance);

	/// generatePoints() from BLElement.
	void generatePoints(int npoints, std::vector<G4ThreeVector> &v);

	/// isWithin() from BLGroupElement.
	bool isWithin(G4ThreeVector &local, G4double tolerance);

	/// BeginOfRunAction() from BLManager::RunAction.
	void BeginOfRunAction(const G4Run *run) { }

	/// EndOfRunAction() from BLManager::RunAction.
	void EndOfRunAction(const G4Run *run);
};

BLCMDtessellatedsolid defaultTessellatedSolid;	// default object
BLCommandAlias aliasTessellatedSolid("tess",defaultTessellatedSolid);

// Default constructor - be sure to use the default constructor BLGroupElement()
BLCMDtessellatedsolid::BLCMDtessellatedsolid() : BLGroupElement()
{
	// register the commandName(), and its synopsis and description.
	registerCommand(BLCMDTYPE_ELEMENT);
	setSynopsis("construct a tessellatedsolid.");
	setDescription("A tessellatedsolid is defined by facets on its "
	"surface. Each facet is either a triangle or a planar quadrilateral, "
	"and the set of surfaces must be closed.\n\n"
	"The surface definition can be read from a file or from the input.file."
	"Leading whitespace is removed, and '*' represents a series of 0 or "
	"more non-blank chars. "
	"The first word in the line describes its content:\n"
	" #*     comment (ignored).\n"
	" tess*  also ignored.\n"
	" v*     a vertex containing 3 doubles (x,y,z in local coords).\n"
	" V*     a vertex containing an integer and 3 doubles (index,x,y,z).\n"
	" f*     a facet, containing 3 or 4 indexes into the vertex array.\n"
	" end*   end of input (or EOF).\n\n"
	"The only requirement on order is that the vertices mentioned in a "
	"facet line must already be present in the vertex array. Each vertex "
	"is appended to the end of the vertex array; for VERTEX lines, large "
	"gaps in the sequence of indices will waste memory. The first vertex "
	"is entry 0.\n\n"
	"This element must be placed (via the place command), and "
	"children can be placed inside it.\n\n");

	// provide initial values for fields
	filename = "";
	material = "";
	color = "1,1,1";
	kill = 0;
	maxStep = -1.0;
	debug = 0;
	height = 0.0;
	width = 0.0;
	length = 0.0;
	Zmin = DBL_MAX;
	Zmax = -DBL_MAX;
	tessellatedsolid = 0;
	physVol = 0;
}

// Copy constructor - be sure to use the copy constructor BLGroupElement(r)
BLCMDtessellatedsolid::BLCMDtessellatedsolid(const BLCMDtessellatedsolid& r) : BLGroupElement(r)
{
	// copy fields one at a time (transfers default values from the
	// default object to this new object).
	filename = r.filename;
	material = r.material;
	color = r.color;
	kill = r.kill;
	maxStep = r.maxStep;
	debug = r.debug;
	height = r.height;
	width = r.width;
	length = r.length;
	Zmin = r.Zmin;
	Zmax = r.Zmax;
	tessellatedsolid = 0;
	physVol = 0;
}

int BLCMDtessellatedsolid::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	std::vector<G4ThreeVector> vtx;
	int nFacets = 0;

	if(argv.size() != 1) {
		printError("tessellatedsolid: Invalid command, must have name");
		return -1;
	}

	if(argv[0] == "default") {
		return defaultTessellatedSolid.handleNamedArgs(namedArgs);
	}

	BLCMDtessellatedsolid *t = new BLCMDtessellatedsolid(defaultTessellatedSolid);
	t->setName(argv[0]);
	int retval = t->handleNamedArgs(namedArgs);

	if(t->maxStep < 0.0) t->maxStep = Param.getDouble("maxStep");

	// check material exists
	if(t->material.size() > 0) getMaterial(t->material);

	// read definition
	FILE *in=0;
	if(t->filename != "") {
		in = fopen(t->filename.c_str(),"r");
		if(!in) {
			printError("tessellatedsolid: cannot open '%s'",
							t->filename.c_str());
			retval = -1;
			t->print(argv[0]);
			return -1;
		}
	}
	G4TessellatedSolid *ts = new G4TessellatedSolid(argv[0]+"Solid");
	for(;;) {
		G4String *line=0;
		if(in != 0) {
			char tmp[1024];
			if(fgets(tmp,sizeof(tmp),in)) {
				static G4String tmpline("");
				tmpline = tmp;
				line = &tmpline;
			}
		} else {
			line = BLCommand::getNextCommand();
		}
		if(line == 0) break;
		// remove initial whitespace
		G4String::size_type i = line->find_first_not_of(" \t\r\n\v");
		if(i == line->npos) i = line->size();
		line->erase(0,i);
		// skip empty line
		if(line->size() < 1) continue;
		// get start of data
		i = line->find_first_of(" \t");
		if(i == line->npos) i = 0;
		// get first char
		char c = line->c_str()[0];
		// switch on line type
		if(c == '#' || c == '\r' || c == '\n') {
			continue;
		} else if(line->find("end") == 0) {
			break;
		} else if(line->find("tess") == 0) {
			continue;
		} else if(c == 'v') {
			G4String l = *line;
			for(int j=0; (j=l.find(',',j))!=l.npos; )
				l.replace(j,1," ");
			double x,y,z;
			if(sscanf(l.c_str()+i,"%lf%lf%lf",&x,&y,&z) != 3)
				goto syntax;
			vtx.push_back(G4ThreeVector(x,y,z));
			if(2.0*fabs(x) > width) width = 2.0*fabs(x);
			if(2.0*fabs(y) > height) height = 2.0*fabs(y);
			if(2.0*fabs(z) > length) length = 2.0*fabs(z);
			if(z < t->Zmin) t->Zmin = z;
			if(z > t->Zmax) t->Zmax = z;
		} else if(c == 'V') {
			G4String l = *line;
			for(int j=0; (j=l.find(',',j))!=l.npos; )
				l.replace(j,1," ");
			int j;
			double x,y,z;
			if(sscanf(l.c_str()+i,"%d%lf%lf%lf",&j,&x,&y,&z)!=4)
				goto syntax;
			vtx[j] = G4ThreeVector(x,y,z);
			if(2.0*fabs(x) > width) width = 2.0*fabs(x);
			if(2.0*fabs(y) > height) height = 2.0*fabs(y);
			if(2.0*fabs(z) > length) length = 2.0*fabs(z);
			if(z < Zmin) Zmin = z;
			if(z > Zmax) Zmax = z;
		} else if(c == 'f') {
			G4String l = *line;
			for(int j=0; (j=l.find(',',j))!=l.npos; )
				l.replace(j,1," ");
			int v1,v2,v3,v4;
			int j = sscanf(l.c_str()+i,"%d%d%d%d",&v1,&v2,&v3,
									&v4);
			if(v1 < 0 || v1 >= vtx.size()) goto syntax;
			if(v2 < 0 || v2 >= vtx.size()) goto syntax;
			if(v3 < 0 || v3 >= vtx.size()) goto syntax;
			if(j == 3) {
				G4TriangularFacet *f = new G4TriangularFacet(
					vtx[v1],vtx[v2],vtx[v3],ABSOLUTE);
				ts->AddFacet(f);
			} else if(j == 4) {
				if(v4 < 0 || v4 >= vtx.size()) goto syntax;
				G4QuadrangularFacet *f = 
					new G4QuadrangularFacet(vtx[v1],vtx[v2],
						vtx[v3],vtx[v4],ABSOLUTE);
				ts->AddFacet(f);
			} else goto syntax;
			++nFacets;
		} else {
syntax:			printError("tessellatedsolid: syntax error '%s'",
							line->c_str());
			retval = -1;
			break;
		}
	}
	if(in) fclose(in);
	ts->SetSolidClosed(true);
	t->tessellatedsolid = ts;

	t->print(argv[0]);
	if(nFacets > 0)
		printf("                           nVertexes=%ld nFacets=%d\n",
							vtx.size(),nFacets);

	if(t->debug) BLManager::getObject()->registerRunAction(t,false);

	return retval;
}

void BLCMDtessellatedsolid::defineNamedArgs()
{
	argString(filename,"filename","File to read for definition.");
	argDouble(maxStep,"maxStep","The maximum stepsize in the element (mm)");
	argString(material,"material","The material of the tessellatedsolid");
	argString(color,"color","The color of the tessellatedsolid (''=invisible)");
	argInt(kill,"kill","Set nonzero to kill every track that enters.");
	argInt(debug,"debug","Set nonzero to display markers on surface.");
	argString(filename,"file","Synonym for filename.");
}

G4VSolid *BLCMDtessellatedsolid::getSolid()
{
	return tessellatedsolid;
}

void BLCMDtessellatedsolid::construct(G4RotationMatrix *relativeRotation,
			G4ThreeVector relativePosition, 
			G4LogicalVolume *parent, 
			G4String parentName,
			G4RotationMatrix *parentRotation,
			G4ThreeVector parentPosition)
{
	G4String thisname = parentName+getName();

	BLAssert(tessellatedsolid != 0);

	G4Material *mat;
	if(material != "")
		mat = getMaterial(material);
	else
		mat = parent->GetMaterial();

	G4LogicalVolume *lv = new G4LogicalVolume(tessellatedsolid,mat, thisname+"LogVol");
	lv->SetVisAttributes(getVisAttrib(color));
	if(maxStep < 0.0) maxStep = Param.getDouble("maxStep");
	lv->SetUserLimits(new G4UserLimits(maxStep));

	// geant4 rotation convention is backwards from g4beamline
	G4RotationMatrix *g4rot = 0;
	if(relativeRotation)
		g4rot = new G4RotationMatrix(relativeRotation->inverse());

	G4VPhysicalVolume *pv = new G4PVPlacement(g4rot,
					relativePosition,lv,thisname,
					parent,false,0,surfaceCheck);
	physVol = pv;

	// handle kill parameter
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

	printf("BLCMDtessellatedsolid::Construct %s parent=%s relZ=%.1f globZ=%.1f\n"
			"\tzmin=%.1f zmax=%.1f\n",
		thisname.c_str(),parentName.c_str(),relativePosition[2],
		globalPosition[2],
		globalPosition[2]-getLength()/2.0,
		globalPosition[2]+getLength()/2.0);
}

G4bool BLCMDtessellatedsolid::isOutside(G4ThreeVector &local, G4double tolerance)
{
	//@ Note tolerance is ignored
	return tessellatedsolid->Inside(local) == kOutside;
}

struct Compare {
	bool operator()(const G4ThreeVector& a, const G4ThreeVector& b) {
		if(a.x() < b.x()) return true;
		if(a.x() > b.x()) return false;
		if(a.y() < b.y()) return true;
		if(a.y() > b.y()) return false;
		if(a.z() < b.z()) return true;
		return false;
	}
};

void BLCMDtessellatedsolid::generatePoints(int npoints, std::vector<G4ThreeVector> &v)
{
	std::set<G4ThreeVector,Compare> vertexSet;

	int nf = tessellatedsolid->GetNumberOfFacets();
	for(int i=0; i<nf; ++i) {
		G4VFacet *f = tessellatedsolid->GetFacet(i);
		int nv = f->GetNumberOfVertices();
		for(int j=0; j<nv; ++j)
			vertexSet.insert(f->GetVertex(j));
	}

	v.clear();

	std::set<G4ThreeVector,Compare>::iterator it;
	for(it=vertexSet.begin(); it!=vertexSet.end(); ++it) {
		v.push_back(*it);
		--npoints;
	}
	while(npoints-- > 0) {
		G4ThreeVector p = tessellatedsolid->GetPointOnSurface();
		v.push_back(p);
	}
}

G4bool BLCMDtessellatedsolid::isWithin(G4ThreeVector &local, G4double tolerance)
{
	//@ Note tolerance is ignored
	return tessellatedsolid->Inside(local) == kInside;
}

void BLCMDtessellatedsolid::EndOfRunAction(const G4Run *run)
{
	static G4Polymarker markers; // must be persistent
	G4VVisManager* pvm = G4VVisManager::GetConcreteInstance();
	if(!pvm) return;
	
	const G4RotationMatrix *rot = physVol->GetObjectRotation();
	G4ThreeVector trans = physVol->GetObjectTranslation();

	std::vector<G4ThreeVector> v;
	generatePoints(1000,v);
	for(unsigned i=0; i<v.size(); ++i) {
		G4ThreeVector p = *rot * v[i] + trans;
		markers.push_back(p);
	}

	markers.SetMarkerType(G4Polymarker::circles);
	markers.SetScreenSize(MARKER_SIZE);
	markers.SetFillStyle(G4VMarker::filled);
	G4VisAttributes va(G4Colour(1.,0.,0.));  // red
	markers.SetVisAttributes(&va);
	pvm->Draw(markers);
}
