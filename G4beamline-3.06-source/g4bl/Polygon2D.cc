//	Polygon2D.cc - compute if point is inside a polygon
//
//#define TEST // un-comment to compile a main program for testing

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Polygon2D.h"

// from http://geomalgorithms.com/a03-_inclusion.html, 20160122
// Copyright 2000 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.
 
// TJR mods: added to class Polygon2D;  Point2D-s are float-s, not int-s;
// A tolerance was tested, but found unnecessary.

// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line (extremely rare for float-s)
//            <0 for P2  right of the line
//    See: Algorithm 1 "Area of Triangles and Polygons"
inline float Polygon2D::isLeft( Point2D P0, Point2D P1, Point2D P2 ) const
{
    float v = (P1.x - P0.x) * (P2.y - P0.y)
            - (P2.x - P0.x) * (P1.y - P0.y);
    return v;
}
//===================================================================

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
int Polygon2D::cn_PnPoly( Point2D P, Point2D* V, int n ) const
{
    int    cn = 0;    // the  crossing number counter

    // loop through all edges of the polygon
    for(int i=0; i<n; i++) {    // edge from V[i]  to V[i+1]
       if(((V[i].y <= P.y) && (V[i+1].y > P.y))     // an upward crossing
        || ((V[i].y > P.y) && (V[i+1].y <=  P.y))) { // a downward crossing
            // compute  the actual edge-ray intersect x-coordinate
            float vt = (float)(P.y  - V[i].y) / (V[i+1].y - V[i].y);
            if(P.x <  V[i].x + vt * (V[i+1].x - V[i].x)) // P.x < intersect
                 ++cn;   // a valid crossing of y=P.y right of P.x
        }
    }
    return (cn&1);    // 0 if even (out), and 1 if  odd (in)

}
//===================================================================


// wn_PnPoly(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
int Polygon2D::wn_PnPoly( Point2D P, Point2D* V, int n ) const
{
    int    wn = 0;    // the  winding number counter

    // loop through all edges of the polygon
    for(int i=0; i<n; i++) {   // edge from V[i] to  V[i+1]
        if(V[i].y <= P.y) {          // start y <= P.y
            if(V[i+1].y  > P.y)      // an upward crossing
                 if(isLeft( V[i], V[i+1], P) > 0)  // P left of  edge
                     ++wn;            // have  a valid up intersect
        }
        else {                        // start y > P.y (no test needed)
            if(V[i+1].y  <= P.y)     // a downward crossing
                 if(isLeft( V[i], V[i+1], P) < 0)  // P right of  edge
                     --wn;            // have  a valid down intersect
        }
    }
    return wn;
}
//===================================================================

void Polygon2D::setVertices(std::vector<Point2D> vtx)
{
	if(vertices != 0) delete[] vertices;

	nVertices = vtx.size();
	vertices = new Point2D[nVertices+1]; // duplicate vtx[0] at the end
	for(unsigned i=0; i<vtx.size(); ++i)
		vertices[i] = vtx[i];
	vertices[nVertices] = vtx[0];        // duplicate vtx[0] at the end
	minX = maxX = vtx[0].x;
	minY = maxY = vtx[0].y;
	for(unsigned i=0; i<vtx.size(); ++i) {
		if(minX > vertices[i].x) minX = vertices[i].x;
		if(maxX < vertices[i].x) maxX = vertices[i].x;
		if(minY > vertices[i].y) minY = vertices[i].y;
		if(maxY < vertices[i].y) maxY = vertices[i].y;
	}
}

#ifdef TEST

int main(int argc, char *argv[])
{
	Polygon2D poly;
	std::vector<Point2D> vtx;
	vtx.push_back(Point2D(0.0,0.0));
	vtx.push_back(Point2D(1.0,0.0));
	vtx.push_back(Point2D(1.0,1.0));
	vtx.push_back(Point2D(0.0,1.0));
	poly.setVertices(vtx);

	char line[2048];
	while(fgets(line,sizeof(line),stdin)) {
		float x,y;
		int i = sscanf(line,"%f%f",&x,&y);
		if(i == 2)
		    printf("isInside(%.3f,%.3f)=%d\n",x,y,poly.isInside(x,y));
		else
		    printf("Invalid\n");
	}
	return 0;
}

#endif // TEST
