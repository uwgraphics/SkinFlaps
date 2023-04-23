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
	void restoreOldNodePositions(vnBccTetrahedra *newVnbt);  // done after new physics library made so new spatial position pointer set
	remapTetPhysics();
	~remapTetPhysics();
	
private:

	std::unordered_map<long long, Vec3f> _oldNPos;
	inline long long a3ToLl(const std::array<short, 3>& a3) {
		long long ll = a3[0] + 1;  // allow for -1
		ll <<= 16;
		ll += a3[1] + 1;
		ll <<= 16;
		ll += a3[2] + 1;
		return ll;
	}

//	std::vector<int> _oldFixedNodes;  // no longer using
	std::vector<Vec3f> _oldNodePositions;
	std::vector<std::array<int, 4> > _oldTets;
	std::vector<bccTetCentroid> _oldCentroids;

	struct bccTetCentroidHasher {
		std::size_t operator()(const bccTetCentroid& k) const
		{  // hash function
			long long ll = k[0];
			ll <<= 16;
			ll += k[1];
			ll <<= 16;
			ll += k[2];
			std::hash<long long> hash_funct;
			return hash_funct(ll);
		}
	};
	std::unordered_multimap<bccTetCentroid, int, bccTetCentroidHasher> _oldTetHash;  // bccTetCenter and index into _tetNodes

//	std::unordered_multimap<long long, int> _oldTetHash;
	std::vector<int> _oldVertexTets;
	std::vector<int> _newToOldNodes;
};

#endif  // __REMAP_TET_PHYSICS__