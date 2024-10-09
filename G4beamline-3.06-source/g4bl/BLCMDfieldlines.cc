//	BLCMDfieldlines.cc
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

#ifdef G4BL_VISUAL

#include <stdio.h>
#include <vector>

#include "BLAssert.hh"
#include "BLCommand.hh"
#include "BLManager.hh"
#include "BLGroup.hh"
#include "BLVisManager.hh"
#include "BLGlobalField.hh"
#include "BLMarkers.hh"

extern void g4bl_exit(int);

class BLCMDfieldlines : public BLCommand, public BLCallback {
	G4double t;
	G4String centerStr;
	G4double radius;
	G4int nLines;
	G4double dl;
	G4String color;
	G4double minField;
	G4int maxPoints;
	G4int subdivide;
	G4int N;
	G4int exit;
	G4int square;
	G4int Efield;
	G4int forever;
	G4ThreeVector center;
	std::vector<G4ThreeVector> points;
	std::vector<G4ThreeVector> argPoints;
	BLMarkers *marker;
	BLVisManager* visManager;
	BLGlobalField *field;
	int ic;		// grid index of center (=N/2)
	double dx;	// grid cell size (in x and in y)
	double scale;
	G4RotationMatrix *rot;
	double worldHalfHeight;
	double worldHalfWidth;
	double worldHalfLength;
public:
	BLCMDfieldlines();
	//% BLCMDfieldlines(const BLCMDfieldlines &r);

	G4String commandName() { return "fieldlines"; }
	
	int command(BLArgumentVector& argv, BLArgumentMap& namedArgs);

	void defineNamedArgs();

	bool isInWorld(G4ThreeVector &pos);

	void generatePoints();

	void generateOneFieldLine(G4ThreeVector &point);

	void callback(int type);

	double bmag(int i, int j);

	G4ThreeVector point(int i, int j);

	int excludeRadiusSq(double bmag);
};

BLCMDfieldlines defaultFieldlinesCommand;	// registers the command, and holds
					// default values of the arguments.

BLCMDfieldlines::BLCMDfieldlines() : BLCommand(), BLCallback(), 
							points(), argPoints()
{
	registerCommand(BLCMDTYPE_OTHER);
	setSynopsis("Display magnetic field lines.");
	setDescription("Field lines are drawn starting within a circle "
	"specified as center and radius (global coordinates); the plane of "
	"the circle is normal to the B field at its center. Field lines "
	"are distributed within the circle "
	"with a density inversely proportional to |B|. While it is attempted "
	"to keep their spacing as uniform as possible, there are both "
	"ambiguity and randomness involved in placing the lines within the "
	"circle. Lines are placed within the circle from its center outward, "
	"and if |B|<minField then no more lines are placed.\n\n"
	"nLines is only approximate, and the actual number of lines drawn will "
	"be within a factor of 2 of the value. Asking for fewer than 10 or "
	"more than 1000 lines is likely to be ineffective. Unnamed parameters "
	"can contain specific x,y,z values of points to start a field line.\n\n"
	"Lines are drawn in both directions starting from the plane of the "
	"circle, and each half-line stops when it again reaches that plane. "
	"Line drawing also stops "
	"whenever |B| is less than minField or when it leaves the world. If "
	"you are interested in field lines far outside the magnets, you may "
	"need to add some object with a large value of x, y, and/or z, in "
	"order to expand the world volume. "
	"With dl=1 and subdivide=10, the accuracy of their meeting is usually "
	"better than 0.1 mm after several tens of meters.\n\n"
	"This command does nothing if not in visualization mode. For best "
	"results, use the Open Inventor viewer, and give your magnets a "
	"transparency of about 0.3 (e.g. color=1,1,1,0.3); if necessary, use "
	"the right-button menu to set the DrawStyle/TransparencyType to "
	"'sorted object add'. With field lines in 3-d, you will want the "
	"ability to zoom, rotate, and move the image interactively.\n\n"
	"This command is not placed into the geometry.\n\n");

	t = 0.0*ns;
	centerStr = "0,0,0";
	radius = 100.0*mm;
	nLines = 100;
	dl = 10.0*mm;
	color="1,1,1";
	minField = 0.001*tesla;
	maxPoints = 10000;
	subdivide = 10;
	N = 128;
	exit = 0;
	square = 0;
	Efield = 0;
	forever = 0;
	marker = 0;
	visManager = 0; // set in callback(4)
	field = 0;      // set in callback(4)
	ic = N/2;
	dx = 0.0;
	scale = ic/4;
	rot = 0;	// set in callback(4)
	worldHalfHeight = worldHalfWidth = worldHalfLength = 0.0;
}

/* %
BLCMDfieldlines::BLCMDfieldlines(const BLCMDfieldlines &r)
{
	t = r.t;
	centerStr = r.centerStr;
	radius = r.radius;
	nLines = r.nLines;
	dl = r.dl;
	color=r.color;
	minField = r.minField;
	maxPoints = r.maxPoints;
	subdivide = r.subdivide;
	N = r.N;
	exit = r.exit;
	square = r.square;
	Efield = r.Efield;
	forever = r.forever;
	visManager = 0; // set in callback(4)
	field = 0;      // set in callback(4)
	ic = N/2;
	dx = 0.0;
	scale = ic/4;
	rot = 0;	// set in callback(4)
	worldHalfHeight = worldHalfWidth = worldHalfLength = 0.0;
}
% */

int BLCMDfieldlines::command(BLArgumentVector& argv, BLArgumentMap& namedArgs)
{
	handleNamedArgs(namedArgs);

	argPoints.clear();
	for(unsigned i=0; i<argv.size(); ++i) {
		std::vector<double> list = getList(argv[i],",");
		if(list.size() != 3) {
			printError("Invalid point: %s",argv[i].c_str());
			continue;
		}
		argPoints.push_back(G4ThreeVector(list[0],list[1],list[2]));
	}

	if(N < 10 || N > 16383)
		printError("Invalid N: < 10 or > 16383");
	ic = N/2;
	dx = radius/ic;
	scale = ic/4;

	std::vector<double> tmp = getList(centerStr,",");
	if(tmp.size() != 3)
		printError("Invald value for center '%s'",centerStr.c_str());
	else
		center = G4ThreeVector(tmp[0],tmp[1],tmp[2]);

	if(Efield)
		minField = minField/tesla*(megavolt/meter);

	print("");

	BLManager::getObject()->registerCallback(new BLCMDfieldlines(*this),4);

	return 0;
}

void BLCMDfieldlines::defineNamedArgs()
{
	argDouble(t,"t","Time at which field lines are plotted (0 ns).",ns);
	argString(centerStr,"center","Center of circle (x,y,z) to start lines (mm, global).");
	argDouble(radius,"radius", "Radius of circle to start lines (mm).",mm);
	argInt(nLines,"nLines","Approximate number of field lines to plot (100).");
	argDouble(dl,"dl","Interval between points plotted (10 mm).",mm);
	argString(color,"color","Color of field lines (white=\"1,1,1\").");
	argDouble(minField,"minField","Minimum B field (0.001 tesla)",tesla);
	argInt(maxPoints,"maxPoints","Max # points plotted in a line (10000).");
	argInt(subdivide,"subdivide","# field integration points between plotted points (10).");
	argInt(N,"N","# grid points in x and y (128).");
	argInt(exit,"exit","Deprecated and ignored.");
	argInt(square,"square","Set nonzero to start from square rather than circle (0).");
	argInt(Efield,"Efield","Set nonzero to draw E field (0); minField is in MegaVolts/meter.");
	argInt(forever,"forever","Set nonzero to draw lines until maxPoints is reached or |B|<minField, not stopping at the initial plane.");
}

bool BLCMDfieldlines::isInWorld(G4ThreeVector &pos)
{
	return fabs(pos[0]) <= worldHalfWidth && 
		fabs(pos[1]) <= worldHalfHeight &&
		fabs(pos[2]) <= worldHalfLength;
}

void BLCMDfieldlines::callback(int type)
{
	// initialize variables that must be delayed.
	visManager = BLVisManager::getObject();
	if(!visManager) return;
	field = BLGlobalField::getObject();
	BLGroup *world = BLGroup::getWorld();
	worldHalfHeight = world->getHeight()/2.0;
	worldHalfWidth = world->getWidth()/2.0;
	worldHalfLength = world->getLength()/2.0;

	if(radius > 0.0 && nLines > 0) {
		// get rotation from Z axis to B at center.
		G4ThreeVector B, E, axis(0.0,0.0,1.0);
		field->getFieldValue(center,t,B,E);
		if(Efield) B = E;
		if(B.mag() < minField) return;
		axis = axis.cross(B.unit());
		if(axis.mag() > 1.0e-6)
		    rot = new G4RotationMatrix(axis.unit(),acos(B.unit()[2]));
		else
		    rot = new G4RotationMatrix();
		G4ThreeVector test = *rot * G4ThreeVector(0.0,0.0,1.0);
		BLAssert((test-B.unit()).mag()<1.0e-6 || (test+B.unit()).mag()<1.0e-6);

		fflush(stdout);
		fprintf(stderr,"\n");

		// generate the points from which field lines are drawn
		// loop over scale to get within a factor of 2 of nLines.
		scale = (ic/4)*bmag(ic,ic);
		for(int i=0; i<15; ++i) {
		    fprintf(stderr,"fieldlines: trying for ~ %d field lines...",
									nLines);
		    points.clear();
		    generatePoints();
		    int n = points.size();
		    fprintf(stderr," got %d\n",n);
		    if(n <= 0) return;
		    if(n < nLines/2)
			scale = 0.75*scale;
		    else if(n > nLines*2)
			scale = 1.2*scale;
		    else
			break;
		}
	}

	// add argPoints
	for(unsigned i=0; i<argPoints.size(); ++i)
		points.push_back(argPoints[i]);

	// generate the field lines.
	fprintf(stderr,"fieldlines: generating %d field lines...",
						(int)points.size());
	for(unsigned i=0; i<points.size(); ++i) {
		generateOneFieldLine(points[i]);
	}
	fprintf(stderr," done.\n\n");
}

void BLCMDfieldlines::generateOneFieldLine(G4ThreeVector &point)
{
	BLAssert(subdivide > 0);

	double step = dl/subdivide;
	G4ThreeVector B, E, norm;
	G4ThreeVector pos=point;
	if(!isInWorld(pos)) return;
	field->getFieldValue(pos,t,B,E);
	if(Efield) B = E;
	if(B.mag() < minField) return;
	// norm is the normal to the plane perpendicular to B at this point.
	norm = B.unit();

	marker = new BLMarkers(color);
	marker->setLine();

	// half-line on positive side of norm
	pos = point;
	marker->addMarker(pos);
	for(int n=0; n<maxPoints; ++n) {
		for(int j=0; j<subdivide; ++j) {
			if(!isInWorld(pos)) goto end1;
			if(forever == 0 && (pos-point)*norm < 0.0) break;
			field->getFieldValue(pos,t,B,E);
			if(Efield) B = E;
			if(B.mag() < minField) goto end1;
			pos += step*B.unit();
		}
		marker->addMarker(pos);
		if(forever == 0 && (pos-point)*norm < 0.0) break;
	}
end1:	marker = new BLMarkers(color);
	marker->setLine();

	// half-line on negative side of norm
	pos = point;
	marker->addMarker(pos);
	for(int n=0; n<maxPoints; ++n) {
		for(int j=0; j<subdivide; ++j) {
			if(!isInWorld(pos)) goto end2;
			if(forever == 0 && (pos-point)*norm > 0.0) break;
			field->getFieldValue(pos,t,B,E);
			if(Efield) B = E;
			if(B.mag() < minField) goto end2;
			pos -= step*B.unit();
		}
		marker->addMarker(pos);
		if(forever == 0 && (pos-point)*norm > 0.0) break;
	}
end2:	marker = 0; // they are automatically drawn
}

double BLCMDfieldlines::bmag(int i, int j) 
{
	G4ThreeVector p(point(i,j)),B, E;
	field->getFieldValue(p,t,B,E);
	if(Efield) B = E;
	return B.mag();
}

G4ThreeVector BLCMDfieldlines::point(int i, int j)
{
	return center + *rot * G4ThreeVector((i-ic)*dx,(j-ic)*dx,0.0);
}

int BLCMDfieldlines::excludeRadiusSq(double bmag)
{
	if(bmag < minField) return N*N;

	return (int)((scale/bmag)*(scale/bmag));
}

void BLCMDfieldlines::generatePoints()
{
/*	This routine places points for field lines into the circle in the
	plane normal to B at center.

	The basic algorithm is to cover the circle with an NxN grid, to
	place points in the grid, and then to exclude grid points within
	a given radius of the point; the radius depends on the value of B
	at the point just placed. The first point is placed at the center,
	and successive points are placed at the non-excluded grid point
	closest to the center.
*/
#define index(i,j) (i*N+j)
#define distance2(i1,j1,i2,j2) ((i1-i2)*(i1-i2) + (j1-j2)*(j1-j2))
#define exclude(I,J,D2) 						\
	for(int i=0; i<N; ++i) {				\
		for(int j=0; j<N; ++j) {			\
			if(distance2(i,j,I,J) <= D2)		\
				grid[index(i,j)] = 0;		\
		}						\
	}

	BLAssert(N >= 10 && N <= 16383); // ensure distance2() fits in an int

	unsigned char *grid = new unsigned char[N*N];

	// initialize grid(), excluding points outside radius
	int excludeSq = ic*ic;
	for(int i=0; i<N; ++i) {
		for(int j=0; j<N; ++j) {
			if(!square && distance2(i,j,ic,ic) > excludeSq)
				grid[index(i,j)] = 0;
			else
				grid[index(i,j)] = 1;
		}
	}

	// place as many points as will fit, starting with (ic,ic)
	int imin=ic, jmin=ic;
	do {
		// place this point
		G4ThreeVector p = point(imin,jmin);
		points.push_back(p);
		exclude(imin,jmin,excludeRadiusSq(bmag(i,j)));

		// find un-excluded point closest to (ic,ic)
		int d2min=N*N;
		imin = -1;
		for(int i=0; i<N; ++i) {
			for(int j=0; j<N; ++j) {
				if(grid[index(i,j)] == 0) continue;
				int d2=distance2(i,j,ic,ic);
				if(d2 < d2min) {
					d2min = d2;
					imin = i;
					jmin = j;
				}
			}
		}
	} while(imin >= 0);
}

#else // G4BL_VISUAL
int dummyBLCMDfieldlines = 0;
#endif // G4BL_VISUAL
