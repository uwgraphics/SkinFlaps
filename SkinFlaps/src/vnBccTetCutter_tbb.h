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

// #include <atomic>
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"

class remapTetPhysics;  // forward declaration

class vnBccTetCutter_tbb
{
public:
	bool makeFirstVnTets(materialTriangles* mt, vnBccTetrahedra* vbt, int maximumCubeGridDimension);  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	void createFirstMacroTets(materialTriangles* mt, vnBccTetrahedra* vbt, const int nLevels, const int maximumDimensionMacroSubdivs);  // creates initial macro tet environment
	void addNewMultiresIncision();  // after have done createFirstMacroTets() and possibly made other incisions, this routine inputs new incision(s) and creates new tet structure.
	inline void setRemapTetPhysics(remapTetPhysics* rtp) { _rtp = rtp; }  // for use in surgical simulation project to reset spatial coords after a topo change.  Can be ignored elsewhere if desired.
	vnBccTetCutter_tbb(void) { _rtp = nullptr; }
	~vnBccTetCutter_tbb(void){}

private:
	materialTriangles* _mt;
	vnBccTetrahedra* _vbt;
	remapTetPhysics* _rtp;
	std::vector<Vec3f> _vMatCoords;  // material coordinates of surface vertices

	std::unordered_set<int> _vnTris;
	std::vector<bccTetCentroid> _vnCentroids;
	int _lastTriangleSize, _lastVertexSize;
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
	std::unordered_map<bccTetCentroid, std::vector<int>, bccTetCentroidHasher> _surfaceCentroids;
	std::atomic<int> _nSurfaceTets;
	int _meganodeSize;  // size of largest macro nodes at beginning of node arrays. _megatetSize size of the largest macrotets at the beginning of the tet arrays is now kept in _vbt.
	int _firstNewExteriorNode;  // Index of first new exterior micronode to be created. This is _meganodeSize plus size of new interior micronodes.
	struct tetTris {
		int tetIdx;
		std::vector<int> tris;
	};
	std::vector<tetTris> _surfaceTetTris;  // COURT could be expanded to handle megatet tris as well.
	// COURT with new data structure is next variable necessary since is filled with empty tris as well.  Could add above variable to pack().
	std::unordered_map<bccTetCentroid, tetTris, bccTetCentroidHasher> _megatetTetTris;  // megatets don't virtual node and may or may not have triangles passing through them.
	typedef std::unordered_map<bccTetCentroid, tetTris, bccTetCentroidHasher>::const_iterator MTTIT;

	std::vector<bccTetCentroid> _vertexTetCentroids;
	struct boundingNodeTris {
		int node;
		std::vector<tetTris*> megaTetTris;  // each is guaranteed to be a sorted vector
	};
	std::unordered_map<std::array<short, 3>, boundingNodeTris, arrayShort3Hasher> _boundingNodeData;
	struct megatetFace {
		std::array<int, 3> nodes;
		std::vector<int> *tris;  // must be sorted
	};
	std::vector<megatetFace> _megatetBounds;

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

	struct patch {
		std::vector<int> tris;
		int tetIndex;
		std::array<int, 4> tetNodes;
		bool noEdgeCut;
	};
	struct hole {  // a patch that bounds a hole in a solid
		patch* pptr;
		Vec3d outerPoint;
		int penetratedFace;
	};
	struct extNode {
		std::array<short, 3> loc;
		int node;
		std::vector<std::pair<int, int> > tiPairs;  // first is tetrahedron number, second is its node index 0-3
	};

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
	struct tetTriangles {
		bccTetCentroid tc;
		std::vector<int> tris;
	};
	typedef oneapi::tbb::concurrent_hash_map<bccTetCentroid, std::vector<int>, btcHashCompare> CENTtris;
	CENTtris _centTris;
	oneapi::tbb::concurrent_vector<zIntrsct> _zIntr;
	struct nodeTetSegment {
		std::vector<int> tetNodeTris;  // COURT convert to pointer
		int tetIdx;
		int tetNodeIndex;
	};
	typedef oneapi::tbb::concurrent_hash_map<std::array<short, 3>, std::list<nodeTetSegment>, ss3HashCompare> NTS_HASH;
	NTS_HASH _ntsHash;

	inline bool sortedVectorsIntersect(const std::vector<int>& v0, const std::vector<int>& v1) {
		auto i0 = v0.begin();
		auto i1 = v1.begin();
		while (i0 != v0.end()) {
			if (*i0 < *i1)
				++i0;
			else if (*i0 > *i1) {
				do {
					++i1;
				} while (i1 != v1.end() && *i0 > *i1);
				if (i1 == v1.end())
					return false;
				if (*i0 == *i1)
					return true;
			}
			else
				return true;
		}
		return false;
	}

	struct newTet {
		int tetIdx;
		bccTetCentroid tc;
		std::array<int, 4> tetNodes;
		std::vector<int> tris;
	};
	oneapi::tbb::concurrent_vector<newTet> _newTets;

	bool latticeTest();
	void macrotetRecutCore();
	void createInteriorNodes();
	void createInteriorMicronodes();
	bool setupBccIntersectionStructures(int maximumGridDimension);
	void fillNonVnTetCenter();
	void fillInteriorMicroTets(std::vector<bccTetCentroid>& recutMacrotets);
	void assignExteriorTetNodes(std::array<short, 3>& locus, std::list<nodeTetSegment>& tetNodeIds, oneapi::tbb::concurrent_vector<extNode>& eNodes);
	void getConnectedComponents(const tetTriangles& tt, oneapi::tbb::concurrent_vector<newTet>& nt_vec, NTS_HASH& local_nts);
	void processHoles(std::list<patch>& patches, std::list<hole>& holes, const Vec3d (&tetNodes)[4], int(&inodes)[4]);
	bool isHole(const patch* patch, const bccTetCentroid& tc, const Vec3d(&tetNodes)[4], Vec3d& patchPoint, int& tetFaceHit);
	int nearestRayPatchHit(const Vec3d& rayBegin, Vec3d rayEnd, const std::vector<int>& tris, Vec3d& hitP, double& distanceSq);  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
	void zIntersectTriangleTbb(Vec3d(&tri)[3], const bool surfaceTriangle, oneapi::tbb::concurrent_vector<zIntrsct>& zi_loc);
	void inputTriangleTetsTbb(const int& surfaceTriangle, CENTtris& centTris);
	void addCentroidMicronodesZ(const bccTetCentroid& tc);
	void linkMicrotetsToMegatets();
	void pack();

};
#endif	// #ifndef _VN_BCC_TET_CUTTER_TBB_
