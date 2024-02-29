#ifndef __REMAP_TET_PHYSICS__
#define __REMAP_TET_PHYSICS__

#include <array>
#include <vector>
#include <unordered_map>
#include "Vec3f.h"
#include "vnBccTetrahedra.h"

// forward declarations
class deepCut;

class remapTetPhysics
{
public:
	typedef std::array<unsigned short, 3> bccTetCentroid;
	void getOldPhysicsData(vnBccTetrahedra *oldVnbt);
	void remapNewPhysicsNodePositions(vnBccTetrahedra *newVnbt);  // done before new physics library made
	inline void clearVnTetVerts() {_newSurfaceTetLocs.clear();}
	inline void insertSurfaceTetVertex(int oldTet, int vertex) { vnTetVert vtv; vtv.vertex = vertex; _newSurfaceTetLocs.insert(std::make_pair(oldTet, vtv)); }
	remapTetPhysics();
	~remapTetPhysics();
	
private:
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
		inline std::size_t operator()(const std::array<unsigned short, 3>& k) const
		{  // hash function
			btHash bh;
			bh.us3.tc = k;
			bh.us3.pad = 0;
			std::hash<long long> hash_funct;
			return hash_funct(bh.ll);
		}
	};
	struct arrayShort3Hasher {
		inline std::size_t operator()(const std::array<short, 3>& k) const
		{  // hash function
			btHash bh;
			bh.ss3.ss3 = k;
			bh.ss3.pad = 0;
			std::hash<long long> hash_funct;
			return hash_funct(bh.ll);
		}
	};
	struct vnTetVert {
		int vertex;
		Vec3f loc;
	};
	std::unordered_multimap<int, vnTetVert> _oldSurfaceTetLocs, _newSurfaceTetLocs;
	std::vector<Vec3f> _oldNodePositions;
	std::unordered_multimap<bccTetCentroid, int, vnBccTetrahedra::bccTetCentroidHasher> _oldTetHash;
	std::vector<std::array<int, 4> > _oldTets;
	std::vector<std::array<short, 3> > _oldNodes;
};

#endif  // __REMAP_TET_PHYSICS__