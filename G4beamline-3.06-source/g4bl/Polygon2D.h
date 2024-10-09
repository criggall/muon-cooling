//	Polygon2D.h 

#include <vector>

struct Point2D {
	float x, y; 
	Point2D() : x(0.0), y(0.0) { }
	Point2D(float _x, float _y) : x(_x), y(_y) { }
};

/**	Class Polygon2D represents a 2-d polygon.
 *	Intended for testing if a point is inside the polygon.
 *
 *	The polygon MUST be simple (sides do not intersect).
 *	If no vertices have been set, isInside() returns true.
 **/
class Polygon2D {
	Point2D *vertices; // has nVertices+1 valid entries
	int nVertices;
	float minX, maxX, minY, maxY;
	inline float isLeft( Point2D P0, Point2D P1, Point2D P2 ) const;
	int cn_PnPoly( Point2D P, Point2D* V, int n ) const;
	int wn_PnPoly( Point2D P, Point2D* V, int n ) const;
public:
	Polygon2D() : vertices(0), nVertices(0) { }
	~Polygon2D() { if(vertices != 0) delete[] vertices; }

	// duplicates vtx[0] at the end.
	void setVertices(std::vector<Point2D> vtx);

	int getNvertices() const { return nVertices; }

	// returns true if Point2D(x,y) is inside the polygon, or if no
	// vertices have been set.
	// Note that a Point2D right on an edge may be reported as either
	// inside or outside; generally those on a lower or left edge will
	// be inside, and those on an upper or right edge will be outside.
	// As the resolution is a float, being on the edge is quite rare.
	bool isInside(float x, float y) const {
		if(vertices == 0) return true;
		if(x < minX || x > maxX || y < minY || y > maxY) return false;
		return wn_PnPoly(Point2D(x,y),vertices,nVertices) != 0;
	}
};
