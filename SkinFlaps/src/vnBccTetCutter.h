#ifndef _VN_BCC_TET_CUTTER_
#define _VN_BCC_TET_CUTTER_

#include <Vec3f.h>
#include <Vec3d.h>
#include <Vec2d.h>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"

class vnBccTetCutter
{
public:
	bool makeFirstVnTets(materialTriangles *mt, vnBccTetrahedra *vbt, int maximumCubeGridDimension);  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	bool remakeVnTets(materialTriangles *mt);  // uses old vbt to get material coords of mt vertices and grid data, from which it makes new vbt.
	vnBccTetCutter(void);
	~vnBccTetCutter(void);

private:
	materialTriangles *_mt;
	vnBccTetrahedra *_vbt;
	std::vector<Vec3d> _vMatCoords;
	struct arrayShort3Hasher{
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
	std::unordered_map<std::array<short, 3>, long, arrayShort3Hasher> _interiorNodes;
	struct planeLineCrossing{
		unsigned long solidRight : 1;
		unsigned long triangle : 31;
	};
	typedef std::multimap<double, planeLineCrossing> PLANE_LINE;
	struct triSegment2{  // plane polygon segment.  These should be strung together to form a closed ring.
		Vec2d uv;	// 2D coord of triangle's entering edge intersection with plane
		long triangle;  // triangle immediately following this intersection in a polygon.
	};
	typedef std::multimap<std::pair<short, short>, std::list<triSegment2> >::iterator PFIT;
	struct vnTetFace{
		std::vector<int> edgeTriangles;
		std::vector<int> interiorTriangles;
		std::array<unsigned long, 3> tetNodes[2];
		unsigned short interiorNodes;  // bit 0 is vertex 0 is interior through bit 2 as vertex 2 is interior. 0 is no interior vertices
		unsigned short set;  // planeSet from 0-5
	};
	typedef std::multimap<std::pair<short, short>, vnTetFace> TFMAP;
	std::vector<bccTetCentroid> _vertexTetLoci;
	struct bccPlane{
		short D;
		std::list<std::list<triSegment2> > polygons2;
		TFMAP vnTetFaces;	// surface intersected, possibly virtual noded, tet faces
	};
	std::vector<std::unique_ptr<bccPlane> > _planeSets2[6];
	std::atomic<long long> _surfaceTetFaceNumber;
	std::unordered_map<long long, std::list<vnTetFace*> > _surfaceTetFaces;
	std::atomic<long> _surfaceTetNumber;
	bool setupBccIntersectionStructures(int maximumGridDimension);
	bool getPlanePolygons(const int planeSet, const int plane);
	void processIntersectionPlane(const int planeSet, const int plane);
	void makeConnectedComponentTetFaces(int planeSet, int planeNumber, std::pair<PFIT, PFIT> &face,
		std::vector<PLANE_LINE> &planeHorizLines, std::vector<PLANE_LINE>(&planeDiagonals)[2]);
	void createSurfaceTetFaces(std::pair<short, short> &face, std::list<std::list<triSegment2> > &facePolygons, bccPlane *planeData);
	void collectSurfaceTetCentersFaces();
	void createSurfaceTetNodes();
	void fillNonVnTetCenter();
	void tetConnectedSurface(bccTetCentroid tc, std::set<long> &triangles, std::vector<long> &vertices);
	std::pair<short, short> uvToPlaneFace(const Vec2d &v);
	void getFaceEdgePath(const Vec2d *startUv, std::pair<short, short> startFace, const Vec2d *endUv, std::pair<short, short> endFace, int triangle, int secondAxis,
		std::list<std::pair<short, short>> &triFaces, std::vector<PLANE_LINE> &horizLines, std::vector<PLANE_LINE>(&diagonals)[2]);
	void createVirtualNodedSurfaceTets();

	// COURT - may use these next lines for github post
//	bool getTriangleAdjacencies(int nTris, const int (*triangles)[3]);
/*	struct longPair	{  // COURT - nuke. Use pair of longs and sort on input.
		long lMin;
		long lMax;
		longPair()	{ lMin = 0; lMax = 0; }
		longPair(int int0, int int1) { lMin = std::min(int0, int1); lMax = std::max(int0, int1); }
		longPair(long l0, long l1) { lMin = std::min(l0, l1); lMax = std::max(l0, l1); }
		bool operator==(const longPair& l) const	{ return lMin == l.lMin&&lMax == l.lMax; }
		bool operator!=(const longPair& l) const	{ return lMin != l.lMin || lMax != l.lMax; }
		longPair& operator=(const longPair& lp)	{ this->lMin = lp.lMin; this->lMax = lp.lMax; return *this; }
	};
	struct longPairTest : public std::binary_function<longPair, longPair, bool>	{
		bool operator()(const longPair &e1, const longPair &e2) const
		{
			if (e1.lMin < e2.lMin)
				return true;
			else if (e1.lMin > e2.lMin)
				return false;
			else
				return (e1.lMax < e2.lMax);
		}
	};
	int nVerts,nTris;
	const float (*vtx)[3];
	const int (*tri)[3];
	std::vector<unsigned int> _adjTris; */
};

#endif	// #ifndef _VN_BCC_TET_CUTTER_
