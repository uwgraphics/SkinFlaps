#ifndef _VN_BCC_TET_CUTTER_OMP_
#define _VN_BCC_TET_CUTTER_OMP_

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

class vnBccTetCutter_omp
{
public:
	bool makeFirstVnTets(materialTriangles* mt, vnBccTetrahedra* vbt, int maximumCubeGridDimension);  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	vnBccTetCutter_omp(void) {}
	~vnBccTetCutter_omp(void){}

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
	struct zIntrsct {
		bool odd, solidBegin;
		int x;
		int y;
		double zInt;
	};
	struct tetTriangles {
		bccTetCentroid tc;
		int tetindx;
		std::vector<int> tris;
	};
	std::unordered_map<std::array<short, 3>, int, arrayShort3Hasher> _interiorNodes;

	std::unordered_map<bccTetCentroid, int, bccTetCentroidHasher> _centroidIndices;
	std::vector<tetTriangles> _tetTris;

	std::atomic<int> _nSurfaceTets;
	struct newTet{
		int tetIdx;
		bccTetCentroid tc;
		std::array<int, 4> tetNodes;
		std::vector<int> tris;
	};
	struct tetTris {
		int tetIdx;
		std::vector<int> tris;
	};
	std::unordered_map<bccTetCentroid, std::list<tetTris>, bccTetCentroidHasher> _centroidTriangles;  // replace this with next 2

	typedef std::unordered_map<bccTetCentroid, std::list<tetTriangles>, bccTetCentroidHasher>::iterator ctIterator;
	typedef std::list<tetTriangles>::iterator ttIterator;
	std::vector<bccTetCentroid> _vertexTetCentroids;
	int _gridSize[3];
	std::vector<std::vector<std::multimap<double, bool> > > evenXy, oddXy;  // Lines parallel to Z axis. First is Z intersect location along line, second is if solid interval start
	// the 0th entry in evenXY will always be empty in i and j.
	struct nodeTetSegment {
		std::vector<int> tetNodeTris;
		int tetIdx;
		int tetNodeIndex;
	};
//	std::unordered_map<std::array<short, 3>, std::list<tetNodeIndex>, arrayShort3Hasher> _exteriorNodeIndices;
	std::unordered_map<std::array<short, 3>, int, arrayShort3Hasher> nts_global;
	std::vector< std::list<nodeTetSegment> > nts_vec;
	std::vector< std::array<short, 3> > nts_locs;
	struct extNode {
		std::array<short, 3> loc;
		std::vector<std::pair<int, int> > tiPairs;  // first is tetrahedron number, second is its node index 0-3
	};

	void createInteriorNodes();
	void assignExteriorTetNodes(int exteriorNodeNum, std::vector<extNode>& eNodes);
	int nearestRayPatchHit(const Vec3f& rayBegin, Vec3f rayEnd, const std::vector<int>& tris, float& distanceSq);  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
	bool nearestPatchPoint(const short(&gl)[4][3], const int tetIdx, const std::vector<int>& tris, Vec3f& closeP, float& distanceSq);
	bool closestPatchPoint(const Vec3f& P, const std::vector<int>& tris, Vec3f& closeP, float& distanceSq);
	bool pointInsidePatchVertex(const Vec3f& P, const int triangle, const int tIndex);
	void getConnectedComponents(tetTriangles& tt, std::vector<newTet>& nt_vec, std::unordered_map<std::array<short, 3>, std::list<nodeTetSegment>, arrayShort3Hasher>& local_nts);  // for this centroid split its triangles into solid connected components
	bool setupBccIntersectionStructures(int maximumGridDimension);
	void inputTriangle(int tri, std::unordered_map<bccTetCentroid, std::vector<int>, bccTetCentroidHasher>& tc_loc, std::vector<zIntrsct>& zi_loc);
	void fillNonVnTetCenter();

	friend class bccTetDecimator;
};
#endif	// #ifndef _VN_BCC_TET_CUTTER_OMP_
