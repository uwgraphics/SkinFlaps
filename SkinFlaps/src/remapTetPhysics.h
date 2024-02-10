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
	std::unordered_multimap<bccTetCentroid, int, bccTetCentroidHasher> _oldTetHash;  // bccTetCenter and index into _tetNodes
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
	std::unordered_map<std::array<short, 3>, int, arrayShort3Hasher> _oldNodeLocs;
	std::vector<Vec3f> _oldNodePositions;
	std::vector<std::array<int, 4> > _oldTets;
};

#endif  // __REMAP_TET_PHYSICS__