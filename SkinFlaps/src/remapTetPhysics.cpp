#include <assert.h>
#include <set>
#include "remapTetPhysics.h"

void remapTetPhysics::getOldPhysicsData(vnBccTetrahedra *oldVnbt)
{
	_oldNodePositions.clear();
	_oldNodePositions.reserve(oldVnbt->nodeNumber());
	for (int n = oldVnbt->nodeNumber(), i = 0; i < n; ++i)
		_oldNodePositions.push_back(oldVnbt->nodeSpatialCoordinate(i));
	//	std::unordered_map<std::array<short, 3>, int, arrayShort3Hasher> _oldNodeLocs;
	_oldNodeLocs.clear();
	_oldNodeLocs.reserve(oldVnbt->nodeNumber());
	for (int n = oldVnbt->nodeNumber(), i = 0; i < n; ++i)
		_oldNodeLocs.insert(std::make_pair(oldVnbt->_nodeGridLoci[i], i));

	_oldTets.clear();
	const auto otna = oldVnbt->getTetNodeArray();
	_oldTets.assign(otna.begin(), otna.end());  // just copy these as remakeNewTets() needs them.
	_oldTetHash.clear();
	_oldTetHash.reserve(oldVnbt->_tetHash.size());
	_oldTetHash.insert(oldVnbt->_tetHash.begin(), oldVnbt->_tetHash.end());
}

void remapTetPhysics::remapNewPhysicsNodePositions(vnBccTetrahedra *newVnbt)
{  // with new multires tet formulation all new physics nodes are no longer clones of old ones following decimation.
	// new low end nodes may be subnodes of larger tets.  Multithread?
	std::vector<bccTetCentroid> nodeCentroids;
	nodeCentroids.assign(newVnbt->_nodeGridLoci.size(), bccTetCentroid());
	for (int n = newVnbt->_tetNodes.size(), i = 0; i < n; ++i) {
		auto const &tn = newVnbt->_tetNodes[i];
		auto const& tc = newVnbt->_tetCentroids[i];
		for (int j = 0; j < 4; ++j)
			nodeCentroids[tn[j]] = tc;
	}
	auto newSpatialCoord = [&](int nodeIdx, int oldTet, Vec3f &bw) {
		auto& sc = newVnbt->_nodeSpatialCoords[nodeIdx];
		auto& tnOld = _oldTets[oldTet];
		sc = _oldNodePositions[tnOld[0]] * (1.0f - bw[0] - bw[1] - bw[2]);
		for (int i = 0; i < 3; i++)
			sc += _oldNodePositions[tnOld[i + 1]] * bw[i];
	};
	for (int n = newVnbt->_nodeGridLoci.size(), i = 0; i < n; ++i) {
		auto pr = _oldNodeLocs.equal_range(newVnbt->_nodeGridLoci[i]);
		if (pr.first == pr.second) {
			auto tc = nodeCentroids[i];
			int level = newVnbt->centroidLevel(tc);
			auto p2 = _oldTetHash.equal_range(tc);
			while(p2.first == p2.second && level < newVnbt->_tetSubdivisionLevels){
				tc = newVnbt->centroidUpOneLevel(tc);
				++level;
				p2 = _oldTetHash.equal_range(tc);
			}
			if (p2.first == p2.second) {  // look down
				tc = nodeCentroids[i];
				level = newVnbt->centroidLevel(tc);
				do {
					bccTetCentroid subC[8];
					newVnbt->subtetCentroids(tc, subC);
					--level;
					Vec3f bw, gl((short(&)[3]) * newVnbt->_nodeGridLoci[i].data());
					for (int j = 0; j < 8; ++j) {
						newVnbt->gridLocusToBarycentricWeight(gl, subC[j], bw);
						if (bw[0] < 0.0f || bw[0] > 1.0f || bw[1] < 0.0f || bw[1] > 1.0f || bw[2] < 0.0f || bw[2] > 1.0f || bw[0] + bw[1] + bw[2] > 1.0f)
							continue;
						p2 = _oldTetHash.equal_range(tc);
						if (p2.first == p2.second)
							continue;
						else if (std::distance(p2.first, p2.second) < 2) {
							newSpatialCoord(i, p2.first->second, bw);
							int junk = 0;
						}
						else
							int junk = 0;
					}
				} while (p2.first == p2.second && level > 1);
				int junk = 0;
			}
			else if (std::distance(p2.first, p2.second) < 2) {
				Vec3f bw, gl((short (&)[3]) *newVnbt->_nodeGridLoci[i].data());
				newVnbt->gridLocusToBarycentricWeight(gl, tc, bw);
				newSpatialCoord(i, p2.first->second, bw);
			}
			else {
				auto& tnUp = newVnbt->_tetNodes[pr.first->second];
				int junk = 0;
			}
		}
		else if (std::distance(pr.first, pr.second) < 2) {
			newVnbt->_nodeSpatialCoords[i] = _oldNodePositions[pr.first->second];
		}
		else {
			int junk = 0;
		}
	}
}

remapTetPhysics::remapTetPhysics()
{
}


remapTetPhysics::~remapTetPhysics()
{
}
