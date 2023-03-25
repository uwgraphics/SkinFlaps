// File: insidePolygon.h
// Author: Court Cutting
// Purpose: 2D inside polygon test modified from Dan Sunday http://geomalgorithms.com/a03-_inclusion.html

// Copyright 2000 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

#ifndef __INSIDE_2D_POLYGON__
#define __INSIDE_2D_POLYGON__

#include <vector>
#include "Vec2f.h"
#include "Vec2d.h"
#include "Vec3f.h"

class insidePolygon
{
public:
	bool insidePolygon2f(const Vec2f &P, const std::vector<Vec2f> &polygon)
	{  // winding number test for a point in a polygon
		// Assumes closed polygon without repeat of first point.
		//      Input:   P = a point, V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
		//      Return:  wn = the winding number (=0 only when P is outside)
		int n = (int)polygon.size();
		int lastIndex = n - 1, wn = 0;    // the  winding number counter
		// loop through all edges of the polygon
		for (int i = 0; i<n; i++) {   // edge from V[lastIndex] to  V[1]
			if (polygon[lastIndex].Y <= P.Y) {          // start y <= P.y
				if (polygon[i].Y  > P.Y)      // an upward crossing
					if (isLeft(polygon[lastIndex], polygon[i], P) > 0)  // P left of  edge
						++wn;            // have  a valid up intersect
			}
			else {                        // start y > P.y (no test needed)
				if (polygon[i].Y <= P.Y)     // a downward crossing
					if (isLeft(polygon[lastIndex], polygon[i], P) < 0.0f)  // P right of  edge
						--wn;            // have  a valid down intersect
			}
			lastIndex = i;
		}
		return wn != 0;
	}

	bool insidePolygon2d(const Vec2d &P, const std::vector<Vec2d> &polygon)
	{  // winding number test for a point in a polygon
		// Assumes closed polygon without repeat of first point.
		//      Input:   P = a point, V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
		//      Return:  wn = the winding number (=0 only when P is outside)
		int n = (int)polygon.size();
		int lastIndex = n - 1, wn = 0;    // the  winding number counter
		// loop through all edges of the polygon
		for (int i = 0; i<n; i++) {   // edge from V[lastIndex] to  V[1]
			if (polygon[lastIndex].Y <= P.Y) {          // start y <= P.y
				if (polygon[i].Y  > P.Y)      // an upward crossing
					if (isLeft(polygon[lastIndex], polygon[i], P) > 0)  // P left of  edge
						++wn;            // have  a valid up intersect
			}
			else {                        // start y > P.y (no test needed)
				if (polygon[i].Y <= P.Y)     // a downward crossing
					if (isLeft(polygon[lastIndex], polygon[i], P) < 0.0f)  // P right of  edge
						--wn;            // have  a valid down intersect
			}
			lastIndex = i;
		}
		return wn != 0;
	}

	bool insidePolygon3f(const Vec3f &P, const std::vector<Vec3f> &polygon, int ignoreDimension)
	{  // winding number test for a point in a polygon
		// Assumes closed polygon without repeat of first point.
		// Input ignoreDimension is Cartesian dimension not considered in this 2d polygon inside test
		//      Input:   P = a point, V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
		//      Return:  wn = the winding number (=0 only when P is outside)
		dim0 = (ignoreDimension + 1) % 3;
		dim1 = (ignoreDimension + 2) % 3;
		int n = (int)polygon.size();
		int lastIndex = n - 1, wn = 0;    // the  winding number counter
		// loop through all edges of the polygon
		for (int i = 0; i<n; i++) {   // edge from V[lastIndex] to  V[1]
			if (polygon[lastIndex].xyz[dim1] <= P.xyz[dim1]) {          // start y <= P.y
				if (polygon[i].xyz[dim1]  > P.xyz[dim1])      // an upward crossing
					if (isLeftIgnore(polygon[lastIndex], polygon[i], P) > 0)  // P left of  edge
						++wn;            // have  a valid up intersect
			}
			else {                        // start y > P.y (no test needed)
				if (polygon[i].xyz[dim1] <= P.xyz[dim1])     // a downward crossing
					if (isLeftIgnore(polygon[lastIndex], polygon[i], P) < 0.0f)  // P right of  edge
						--wn;            // have  a valid down intersect
			}
			lastIndex = i;
		}
		return wn != 0;
	}

	insidePolygon(){}
	~insidePolygon(){}

private:
	int dim0, dim1;

	inline float isLeft(const Vec2f &P0, const Vec2f &P1, const Vec2f &P2)
	{  // tests if a point is Left|On|Right of an infinite line.
		//    Input:  three points P0, P1, and P2
		//    Return: >0 for P2 left of the line through P0 and P1
		//            =0 for P2  on the line
		//            <0 for P2  right of the line
		return ((P1.X - P0.X) * (P2.Y - P0.Y)
			- (P2.X - P0.X) * (P1.Y - P0.Y));
	}

	inline double isLeft(const Vec2d &P0, const Vec2d &P1, const Vec2d &P2)
	{  // tests if a point is Left|On|Right of an infinite line.
		//    Input:  three points P0, P1, and P2
		//    Return: >0 for P2 left of the line through P0 and P1
		//            =0 for P2  on the line
		//            <0 for P2  right of the line
		return ((P1.X - P0.X) * (P2.Y - P0.Y)
			- (P2.X - P0.X) * (P1.Y - P0.Y));
	}

	inline float isLeftIgnore(const Vec3f &P0, const Vec3f &P1, const Vec3f &P2)
	{  // tests if a point is Left|On|Right of an infinite line.
		//    Input:  three points P0, P1, and P2
		//    Return: >0 for P2 left of the line through P0 and P1
		//            =0 for P2  on the line
		//            <0 for P2  right of the line
		return ((P1.xyz[dim0] - P0.xyz[dim0]) * (P2.xyz[dim1] - P0.xyz[dim1])
			- (P2.xyz[dim0] - P0.xyz[dim0]) * (P1.xyz[dim1] - P0.xyz[dim1]));
	}

};
#endif // __INSIDE_2D_POLYGON__

