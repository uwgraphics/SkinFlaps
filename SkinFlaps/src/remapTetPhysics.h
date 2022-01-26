#ifndef __REMAP_TET_PHYSICS__
#define __REMAP_TET_PHYSICS__

#include <vector>
#include <unordered_map>
#include "Vec3f.h"
#include "vnBccTetrahedra.h"

// forward declarations
class deepCut;

class remapTetPhysics
{
public:
	void getOldPhysicsData(vnBccTetrahedra *oldVnbt);
	void remapNewPhysicsNodePositions(vnBccTetrahedra *newVnbt);  // done before new physics library made
	void restoreOldNodePositions(vnBccTetrahedra *newVnbt);  // done after new physics library made so new spatial position pointer set
	remapTetPhysics();
	~remapTetPhysics();
	
private:
	std::vector<long> _oldFixedNodes;
	std::vector<Vec3f> _oldNodePositions;
	std::vector<std::array<long, 4> > _oldTets;
	std::vector<bccTetCentroid> _oldCentroids;
	std::unordered_multimap<long long, long> _oldTetHash;
	std::vector<long> _oldVertexTets;
	std::vector<long> _newToOldNodes;
};

#endif  // __REMAP_TET_PHYSICS__