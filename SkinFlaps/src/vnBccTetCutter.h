#ifndef _VN_BCC_TET_CUTTER_
#define _VN_BCC_TET_CUTTER_

#include <Vec3f.h>
#include <Vec2d.h>
#include "boundingBox.h"
#include <list>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
// #include <atomic>  // no MP yet
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"

class vnBccTetCutter
{
public:
	bool makeFirstVnTets(materialTriangles* mt, vnBccTetrahedra* vbt, int maximumCubeGridDimension);  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	bool subcutMacroTets(std::vector<int>& macroTets);
	vnBccTetCutter(void) { _macroTetCentroidSubset.clear(); }
	~vnBccTetCutter(void){}

private:
	materialTriangles* _mt;
	vnBccTetrahedra* _vbt;
	std::vector<Vec3f> _vMatCoords;
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
	struct bccTetCentroidHasher {
		std::size_t operator()(const std::array<unsigned short, 3>& k) const
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
	struct tetTriangles {
		int tetindx;
		std::vector<int> tris;
	};
	std::unordered_map<std::array<short, 3>, int, arrayShort3Hasher> _interiorNodes;
	std::unordered_map<bccTetCentroid, std::list<tetTriangles>, bccTetCentroidHasher> _centroidTriangles;
	typedef std::unordered_map<bccTetCentroid, std::list<tetTriangles>, bccTetCentroidHasher>::iterator ctIterator;
	typedef std::list<tetTriangles>::iterator ttIterator;
	std::vector<bccTetCentroid> _vertexTetCentroids;
	std::vector<bccTetCentroid> _macroTetCentroidSubset;  // centroids of all macrotets in subset for top down subcut
	int _gridSize[3];
//	std::atomic<int> _nSurfaceTets;  //future MP - NUKE?
	std::vector<std::vector<std::multimap<double, bool> > > evenXy, oddXy;  // Lines parallel to Z axis. First is Z intersect location along line, second is if solid interval start
	// the 0th entry in evenXY will always be empty in i and j.
	struct tetNodeIndex {
		tetTriangles *ttPtr;
		int nodeIndex;
	};
	std::unordered_map<std::array<short, 3>, std::list<tetNodeIndex>, arrayShort3Hasher> _exteriorNodeIndices;

	void createInteriorNodes();
	void assignExteriorTetNodes();
	int nearestRayPatchHit(const Vec3f& rayBegin, Vec3f rayEnd, const std::vector<int>& tris, float& distanceSq);  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
	bool nearestPatchPoint(const short(&gl)[4][3], const int tetIdx, const std::vector<int>& tris, Vec3f& closeP, float& distanceSq);
	bool closestPatchPoint(const Vec3f& P, const std::vector<int>& tris, Vec3f& closeP, float& distanceSq);
	bool pointInsidePatchVertex(const Vec3f& P, const int triangle, const int tIndex);
	void getConnectedComponents(const bccTetCentroid& tc, std::list<tetTriangles >& ct);  // for this centroid split its triangles into solid connected components
	bool setupBccIntersectionStructures(int maximumGridDimension);
	void inputTriangle(int tri);
	void fillNonVnTetCenter();

	friend class bccTetDecimator;
};
#endif	// #ifndef _VN_BCC_TET_CUTTER_
