#include <assert.h>
#include <set>
#include <algorithm>
#include "materialTriangles.h"
#include "remapTetPhysics.h"

void remapTetPhysics::getOldPhysicsData(vnBccTetrahedra *oldVnbt)
{
	_oldNodePositions.clear();
	_oldNodePositions.reserve(oldVnbt->nodeNumber());
//	_oldNodeLocs.clear();
//	_oldNodeLocs.reserve(oldVnbt->nodeNumber());
	for (int n = oldVnbt->nodeNumber(), i = 0; i < n; ++i) {
//		const Vec3f& nsc = oldVnbt->nodeSpatialCoordinate(i);
		_oldNodePositions.push_back(oldVnbt->nodeSpatialCoordinate(i));
//		_oldNodeLocs.insert(std::make_pair(oldVnbt->_nodeGridLoci[i], nsc));
	}
	_oldTets.clear();
	_oldTets.assign(oldVnbt->_tetNodes.begin(), oldVnbt->_tetNodes.end());  // just copy these as remakeNewTets() needs them.
	_oldNodes.clear();
	_oldNodes.assign(oldVnbt->_nodeGridLoci.begin(), oldVnbt->_nodeGridLoci.end());
	_oldTetHash.clear();
	_oldTetHash = std::move(oldVnbt->_tetHash);
//	_oldTetHash.insert(oldVnbt->_tetHash.begin(), oldVnbt->_tetHash.end());
	_oldSurfaceTetLocs = std::move(_newSurfaceTetLocs);
	materialTriangles* mt = oldVnbt->getMaterialTriangles();
	for (auto& vtl : _oldSurfaceTetLocs)
		mt->getVertexCoordinate(vtl.second.vertex, vtl.second.loc.xyz);
}

void remapTetPhysics::remapNewPhysicsNodePositions(vnBccTetrahedra *newVnbt)
{  // with new multires tet formulation all new physics nodes are no longer clones of old ones following decimation.
	// new low end nodes may be subnodes of larger tets.  Multithread?
	materialTriangles *mt = newVnbt->getMaterialTriangles();
	for (auto& vtl : _newSurfaceTetLocs)
		mt->getVertexCoordinate(vtl.second.vertex, vtl.second.loc.xyz);
	std::vector<char> nodes(newVnbt->_nodeGridLoci.size(), 0x00);
	for (int n = newVnbt->_tetNodes.size(), i = 0; i < n; ++i) {
		auto tc = newVnbt->_tetCentroids[i];
		const auto& tn = newVnbt->_tetNodes[i];
		auto pr = _oldTetHash.equal_range(tc);
		int level = 1;
		bool sameSize = true;
		while (pr.first == pr.second && level <= newVnbt->_tetSubdivisionLevels) {
			sameSize = false;
			tc = newVnbt->centroidUpOneLevel(tc);
			++level;
			pr = _oldTetHash.equal_range(tc);
		}
		assert(level < 5);
		if (std::distance(pr.first, pr.second) > 1) {  // origin was virtual noded level 1 tet.
			assert(sameSize);  // new tet must also be level 1
			auto ntl = _newSurfaceTetLocs.find(i);
			assert(ntl != _newSurfaceTetLocs.end());
			int bestOldTet = -1;
			float dsq, minD = FLT_MAX;
			while (pr.first != pr.second) {
				auto otl = _oldSurfaceTetLocs.find(pr.first->second);
				assert(otl != _oldSurfaceTetLocs.end());
				dsq = (otl->second.loc - ntl->second.loc).length2();
				// unfortunately incisions can change triangles so they can't always be used for a match, but when they do dsq == 0
				if (dsq < 1e-5f) {
					bestOldTet = pr.first->second;
					break;
				}
				if (minD > dsq) {
					minD = dsq;
					bestOldTet = pr.first->second;
				}
				++pr.first;
			}
			auto& oN = _oldTets[bestOldTet];
			for (int j = 0; j < 4; ++j) {
				if (nodes[tn[j]])
					continue;
				nodes[tn[j]] = 1;
				newVnbt->_nodeSpatialCoords[tn[j]] = _oldNodePositions[oN[j]];
			}
		}
		else if (pr.first != pr.second) {
			auto& oN = _oldTets[pr.first->second];
			if (sameSize) {
				for (int j = 0; j < 4; ++j) {
					if (nodes[tn[j]])
						continue;
					nodes[tn[j]] = 1;
					newVnbt->_nodeSpatialCoords[tn[j]] = _oldNodePositions[oN[j]];
				}
			}
			else {
				for (int j = 0; j < 4; ++j) {
					if (nodes[tn[j]])
						continue;
					nodes[tn[j]] = 1;
					Vec3f nodeLocus((const short(&)[3]) * newVnbt->_nodeGridLoci[tn[j]].data());
					Vec3f bw, C;
					newVnbt->gridLocusToBarycentricWeight(nodeLocus, tc, bw);
					C = _oldNodePositions[oN[0]] * (1.0f - bw[0] - bw[1] - bw[2]);
					C += _oldNodePositions[oN[1]] * bw[0];
					C += _oldNodePositions[oN[2]] * bw[1];
					C += _oldNodePositions[oN[3]] * bw[2];
					newVnbt->_nodeSpatialCoords[tn[j]] = C;
				}
			}
		}
		else  // only doing further topo refinements during surgical simulation
			throw(std::logic_error("Logic error in remapTetPhysics.\n"));
	}
}

remapTetPhysics::remapTetPhysics()
{
}


remapTetPhysics::~remapTetPhysics()
{
}
