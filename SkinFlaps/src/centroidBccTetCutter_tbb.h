#ifndef _CENTROID_BCC_TET_CUTTER_TBB_
#define _CENTROID_BCC_TET_CUTTER_TBB_

#include <Vec3f.h>
#include <Vec2d.h>
#include <Vec3d.h>
#include "boundingBox.h"
#include <list>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
// add threading building blocks
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include <atomic>
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"

class centroidBccTetCutter_tbb
{
public:
	bool makeFirstVnTets(materialTriangles* mt, vnBccTetrahedra* vbt, int maximumCubeGridDimension);  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	bool remakeVnTets(materialTriangles* mt);  // use above setup with new material coord incisions made in mt
	centroidBccTetCutter_tbb(void) {}
	~centroidBccTetCutter_tbb(void){}

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
	tbb::concurrent_vector<zIntrsct> _zIntersects;

	struct tetTriangles {
		bccTetCentroid tc;
		int tetindx;
		std::vector<int> tris;
	};
	std::unordered_map<std::array<short, 3>, int, arrayShort3Hasher> _interiorNodes;

	static std::unordered_map<bccTetCentroid, int, bccTetCentroidHasher> _centroidIndices;
	static std::vector<tetTriangles> _tetTris;

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
//	std::unordered_map<bccTetCentroid, std::list<tetTris>, bccTetCentroidHasher> _centroidTriangles;  // replace this with next 2

	tbb::concurrent_unordered_map<bccTetCentroid, std::vector<int>, bccTetCentroidHasher> _centroidTriangles;  // careful does not support concurrent erase

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
//	std::vector< std::array<short, 3> > nts_locs;

	tbb::concurrent_vector< std::array<short, 3> > nts_locs;


	struct extNode {
		std::array<short, 3> loc;
		std::vector<std::pair<int, int> > tiPairs;  // first is tetrahedron number, second is its node index 0-3
	};

	void createInteriorNodes();
	bool setupBccIntersectionStructures(int maximumGridDimension);
	void inputTriangle(int tri, std::unordered_map<bccTetCentroid, std::vector<int>, bccTetCentroidHasher>& tc_loc, std::vector<zIntrsct>& zi_loc);
	void fillNonVnTetCenter();
	bool tetCutCore();
	void assignExteriorTetNodes(int exteriorNodeNum, std::vector<extNode>& eNodes);

	void getConnectedComponents(tetTriangles& tt, std::vector<newTet>& nt_vec, std::unordered_map<std::array<short, 3>, std::list<nodeTetSegment>, arrayShort3Hasher>& local_nts);  // for this centroid split its triangles into solid connected components
	bool isInsidePatch(const Vec3d& P, const std::vector<int>& tris, Vec3d& closestP);
	int nearestRayPatchHit(const Vec3d& rayBegin, Vec3d rayEnd, const std::vector<int>& tris, Vec3d& hitP, double& distanceSq);  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.

	friend class bccTetDecimator;
};
#endif	// #ifndef _CENTROID_BCC_TET_CUTTER_TBB_
