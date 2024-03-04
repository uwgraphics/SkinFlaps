#ifndef _VN_BCC_TET_CUTTER_TBB_
#define _VN_BCC_TET_CUTTER_TBB_

#include <Vec3f.h>
#include <Vec2d.h>
#include <Vec3d.h>
#include "boundingBox.h"
#include <list>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>
#include <algorithm>

#include "oneapi/tbb/concurrent_hash_map.h"
#include "oneapi/tbb/concurrent_vector.h"
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/parallel_for.h"

// #include <atomic>  // no MP yet
#include "materialTriangles.h"
#include "multiresBccTets.h"

class vnBccTetCutter_tbb
{
public:
	bool makeFirstVnTets(materialTriangles* mt, multiresBccTets* vbt, int maximumCubeGridDimension);  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	bool remakeVnTets(materialTriangles* mt);  // use above setup with new material coord incisions made in mt
	void createFirstMacroTets(materialTriangles* mt, multiresBccTets* vbt, const int nLevels, const int maximumDimensionMacroSubdivs);  // creates initial macro tet environment
	vnBccTetCutter_tbb(void) {}
	~vnBccTetCutter_tbb(void){}

private:
	materialTriangles* _mt;
	multiresBccTets* _vbt;
	std::vector<Vec3f> _vMatCoords;
	int _tetSubdivisionLevels;  // if not doing multiresolution, this will be 1.  If using binary subdivision macrotets, this will be the highest level subdivision.

	struct unsigned3 {
		std::array<unsigned short, 3> tc;
		unsigned short pad;
	};
	struct signed3 {
		std::array<short, 3> ss3;
		short pad;
	};
	union btHash {
		uint64_t ll;
		unsigned3 us3;
		signed3 ss3;
	};
	struct bccTetCentroidHasher {
		std::size_t operator()(const std::array<unsigned short, 3>& k) const
		{  // hash function
			btHash bh;
			bh.us3.tc = k;
			bh.us3.pad = 0;
			std::hash<long long> hash_funct;
			return hash_funct(bh.ll);
		}
	};
	struct arrayShort3Hasher {
		std::size_t operator()(const std::array<short, 3>& k) const
		{  // hash function
			btHash bh;
			bh.ss3.ss3 = k;
			bh.ss3.pad = 0;
			std::hash<long long> hash_funct;
			return hash_funct(bh.ll);
		}
	};
	std::unordered_map<std::array<short, 3>, int, arrayShort3Hasher> _interiorNodes;
	struct tetTriangles {
		bccTetCentroid tc;
		std::vector<int> tris;
	};
	std::unordered_set<bccTetCentroid, bccTetCentroidHasher> _surfaceCentroids;
//	static std::vector<tetTriangles> _tetTris;

	std::atomic<int> _nSurfaceTets;
	int _megatetSize, _meganodeSize;  // size of the largest macrotets at the beginning of the tet arrays, and size of largest macro nodes at beginning of node arrays
	int _firstNewExteriorNode;  // Index of first new exterior micronode to be created. This is _meganodeSize plus size of new interior micronodes.
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

//	typedef std::unordered_map<bccTetCentroid, std::list<tetTriangles>, bccTetCentroidHasher>::iterator ctIterator;
//	typedef std::list<tetTriangles>::iterator ttIterator;
	std::vector<bccTetCentroid> _vertexTetCentroids;
	struct boundingNodeTris {
		int node;
		std::vector<int> tris;
	};
	std::unordered_multimap<std::array<short, 3>, boundingNodeTris, arrayShort3Hasher> _boundingNodeData;


	struct zIntersectFlags{
		// will usually occupy 2 bytes:
		unsigned char macroNode : 1;
		unsigned char solidBegin : 1;
		unsigned char surfaceTri : 1;
		unsigned char odd : 1;
	};
	struct zIntrsct {
		zIntersectFlags flags;
		int x;
		int y;
		float zInt;
	};
	std::vector<std::vector<std::multimap<double, zIntersectFlags> > > evenXy, oddXy;  // Lines parallel to Z axis. First is Z intersect location along line, second 1 bit if solid interval start & 2 bit if border triangle intersect
	// the 0th entry in evenXY will always be empty in i and j.

	struct nodeTetSegment {
		std::vector<int> tetNodeTris;
		int tetIdx;
		int tetNodeIndex;
	};
	struct extNode {
		std::array<short, 3> loc;
		int node;
		std::vector<std::pair<int, int> > tiPairs;  // first is tetrahedron number, second is its node index 0-3
	};
	std::unordered_map<int, std::pair<int, int> > _decimatedNodes;  // first is id of node decimated, second are the two nodes of the edge first node splits.
	typedef std::unordered_map<int, std::pair<int, int> >::iterator DNIT;

	struct ss3HashCompare {
		static size_t hash(const std::array<short, 3>& ss3) {
			btHash bh;
			bh.ss3.ss3 = ss3;
			bh.ss3.pad = 0;
			std::hash<long long> hash_funct;
			return hash_funct(bh.ll);
		}
		//! True if strings are equal
		static bool equal(const std::array<short, 3>& x, const std::array<short, 3>& y) {
			return x == y;
		}
	};
	struct btcHashCompare {
		static size_t hash(const bccTetCentroid& tc) {
			btHash bh;
			bh.us3.tc = tc;
			bh.us3.pad = 0;
			std::hash<long long> hash_funct;
			return hash_funct(bh.ll);
		}
		//! True if strings are equal
		static bool equal(const bccTetCentroid& x, const bccTetCentroid& y) {
			return x == y;
		}
	};
	typedef oneapi::tbb::concurrent_hash_map<bccTetCentroid, std::vector<int>, btcHashCompare> CENTtris;
	CENTtris _centTris;
	oneapi::tbb::concurrent_vector<zIntrsct> _zIntr;
	typedef oneapi::tbb::concurrent_hash_map<std::array<short, 3>, std::list<nodeTetSegment>, ss3HashCompare> NTS_HASH;
	NTS_HASH _ntsHash;
	oneapi::tbb::concurrent_vector<newTet> _newTets;

	void createInteriorNodes();
	void createInteriorMicronodes();
	bool setupBccIntersectionStructures(int maximumGridDimension);
	void fillNonVnTetCenter();
	void fillInteriorMicroTets(std::vector<bccTetCentroid>& recutMacrotets);
	bool tetCutCore();
	void assignExteriorTetNodes(std::array<short, 3>& locus, std::list<nodeTetSegment>& tetNodeIds, oneapi::tbb::concurrent_vector<extNode>& eNodes);
	void getConnectedComponents(const tetTriangles& tt, oneapi::tbb::concurrent_vector<newTet>& nt_vec, NTS_HASH& local_nts);
	bool isInsidePatch(const Vec3d& P, const std::vector<int>& tris, Vec3d& closestP);
	int nearestRayPatchHit(const Vec3d& rayBegin, Vec3d rayEnd, const std::vector<int>& tris, Vec3d& hitP, double& distanceSq);  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
	void zIntersectTriangleTbb(Vec3d(&tri)[3], const bool surfaceTriangle, oneapi::tbb::concurrent_vector<zIntrsct>& zi_loc);
	void inputTriangleTetsTbb(const int& surfaceTriangle, CENTtris& centTris);
	int faceAdjacentMacrotet(const bccTetCentroid& macroTc, const int face, bccTetCentroid& tcAdj);
	void addCentroidMicronodesZ(const bccTetCentroid& tc);
	void decimateInteriorMicroTets(int firstInteriorMicroTet, std::vector<std::array<int, 3> > &boundingTris);
	void pack();

};
#endif	// #ifndef _VN_BCC_TET_CUTTER_TBB_
