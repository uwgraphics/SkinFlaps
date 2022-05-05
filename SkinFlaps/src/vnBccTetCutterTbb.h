// File: vnBccTetCutterTbb.h
// Author: Court Cutting MD
// Date: 1/21/202
// Purpose: Defines a class that takes a nonselfintersecting (although may be self coincident) closed manifold triangulated surface solid
//    and produces a virtual noded bcc tetrahedral solid as described in: Molino N,  Bao Z, and Fedkiw R: http://physbam.stanford.edu/~fedkiw/papers/stanford2004-01.pdf
//    In contrast to the original serial version, this one uses Intel's Threading Building Blocks for speed.

#ifndef _VN_BCC_TET_CUTTER_TBB_
#define _VN_BCC_TET_CUTTER_TBB_

#include <Vec3f.h>
#include <Vec3d.h>
#include <Vec2d.h>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
// add threading building blocks
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"

class surgicalActions;  // for debug

class vnBccTetCutterTbb
{
public:
	bool makeFirstVnTets(materialTriangles* mt, vnBccTetrahedra* vbt, int maximumCubeGridDimension);  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	bool remakeVnTets(materialTriangles* mt);  // uses old vbt to get material coords of mt vertices and grid data, from which it makes new vbt.
	void setSurgicalActions(surgicalActions* sa) { _sa = sa; }  // for debug.  May nuke later
	vnBccTetCutterTbb(void);
	~vnBccTetCutterTbb(void);

private:
	surgicalActions* _sa;  // for debugging

	materialTriangles* _mt;
	vnBccTetrahedra* _vbt;
	std::vector<Vec3d> _vMatCoords;
public:
	struct arrayShort3Hasher {
		std::size_t operator()(const std::array<short, 3>& k) const
		{  // hash function
			long long lli = k[0];
			lli <<= 16;
			lli += k[1];
			lli <<= 16;
			lli += k[2];
			std::hash<long long> hash_funct;
			return hash_funct(lli);
		}
	};
	std::unordered_map<std::array<short, 3>, int, arrayShort3Hasher> _interiorNodes;
	typedef tbb::concurrent_vector< std::array<short, 3> > LocusVec;
	LocusVec _nodeLoci;
	struct planeLineCrossing {
		unsigned int solidRight : 1;
		unsigned int triangle : 31;
	};
	typedef std::multimap<double, planeLineCrossing> PLANE_LINE;
	struct triSegment2 {  // plane polygon segment.  These should be strung together to form a closed ring.
		Vec2d uv;	// 2D coord of triangle's entering edge intersection with plane
		int triangle;  // triangle immediately following this intersection in a polygon.
	};
	typedef std::multimap<std::pair<short, short>, std::vector<triSegment2> >::iterator PFIT;
	struct vnTetFace {
		std::vector<int> edgeTriangles;
		std::vector<int> interiorTriangles;
		std::array<unsigned int, 3> tetNodes[2];
		unsigned short interiorNodes;  // bit 0 is vertex 0 is interior through bit 2 as vertex 2 is interior. 0 is no interior vertices
		unsigned short set;  // planeSet from 0-5
	};
	typedef std::multimap<std::pair<short, short>, vnTetFace> TFMAP;
	struct bccPlane {
		short D;
		std::list<std::list<triSegment2> > polygons2;
		TFMAP vnTetFaces;	// surface intersected, possibly virtual noded, tet faces
	};
	tbb::concurrent_vector<bccPlane> _planeSetsTbb[6];

private:
	std::vector<bccTetCentroid> _vertexTetLoci;
	struct tetType {
		std::array<int, 4> tetNodes;
		bccTetCentroid centroid;
	};
	tbb::concurrent_vector<tetType> _tets;
	typedef tbb::concurrent_unordered_map<long long, std::list<vnTetFace*>, std::hash<long long> > UMSTF;
	tbb::concurrent_vector<UMSTF> _stfVec;
	std::atomic<long long> _surfaceTetFaceNumber;

	bool setupBccIntersectionStructures(int maximumGridDimension);
	void collectSurfaceTetCentersFaces();
	void createVirtualNodedSurfaceTets();
	void createSurfaceTetNodes();
	void fillNonVnTetCenter();
	void tetConnectedSurface(bccTetCentroid tc, std::set<int>& triangles, std::vector<int>& vertices);
};
#endif	// #ifndef _VN_BCC_TET_CUTTER_TBB_
