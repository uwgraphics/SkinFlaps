#include <assert.h>
#include <set>
#include <algorithm>
#include "materialTriangles.h"
#include "remapTetPhysics.h"

void remapTetPhysics::getOldPhysicsData(vnBccTetrahedra *oldVnbt)
{
	_oldNodePositions.clear();
	_oldNodePositions.reserve(oldVnbt->nodeNumber());
	for (int n = oldVnbt->nodeNumber(), i = 0; i < n; ++i)
		_oldNodePositions.push_back(oldVnbt->nodeSpatialCoordinate(i));
	_oldTetCentroids.clear();
	_oldTetCentroids.assign(oldVnbt->_tetCentroids.begin(), oldVnbt->_tetCentroids.end());
	_oldVertexTets.clear();
	_oldVertexTets.assign(oldVnbt->_vertexTets.begin(), oldVnbt->_vertexTets.end());
	_oldTets.clear();
	_oldTets.assign(oldVnbt->_tetNodes.begin(), oldVnbt->_tetNodes.end());  // just copy these as remakeNewTets() needs them.
	_oldNodes.clear();
	_oldNodes.assign(oldVnbt->_nodeGridLoci.begin(), oldVnbt->_nodeGridLoci.end());
	_oldTetHash.clear();
	_oldTetHash = std::move(oldVnbt->_tetHash);
	_oldVnTetTris = std::move(_newVnTetTris);
}

void remapTetPhysics::remapNewPhysicsNodePositions(vnBccTetrahedra *newVnbt)
{  // with new multires tet formulation all new physics nodes are no longer clones of old ones following decimation.
	// new low end nodes may be subnodes of larger tets.  Multithread?
	// Start known one to one tet correspondences and use any nodes to solve for vn multiplicities
	materialTriangles *mt = newVnbt->getMaterialTriangles();
	std::vector<char> nodes(newVnbt->_nodeGridLoci.size(), 0x00);
	auto processTet = [&](int tetIdx) ->bool {
		auto tc = newVnbt->_tetCentroids[tetIdx];
		const auto& tn = newVnbt->_tetNodes[tetIdx];
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
			int j;
			for (j = 0; j < 4; ++j) {
				if (nodes[tn[j]])
					break;
			}
			if (j < 4) { // prior one to one correspondence found
				float dsq, minD = FLT_MAX;
				int bestTet = -1;
				while (pr.first != pr.second) {
					auto& oN = _oldTets[pr.first->second];
					dsq = (newVnbt->nodeSpatialCoordinate(tn[j]) - _oldNodePositions[oN[j]]).length2();
					if (dsq < minD) {
						minD = dsq;
						bestTet = pr.first->second;
					}
					++pr.first;
				}
				auto& oN = _oldTets[bestTet];
				for (int k = 0; k < j; ++k)
					newVnbt->_nodeSpatialCoords[tn[k]] = _oldNodePositions[oN[k]];
				for (int k = j + 1; k < 4; ++k) {
					if (nodes[tn[k]])
						continue;
					nodes[tn[k]] = 1;
					newVnbt->_nodeSpatialCoords[tn[k]] = _oldNodePositions[oN[k]];
				}
				return true;
			}
			// no simple prior correspondence found
			auto ntl = _newVnTetTris.find(tetIdx);
			if (ntl == _newVnTetTris.end())
				throw(std::logic_error("Program error 3 in remapTetPhysics.\n"));
			std::set<int> newTetV;
			for (auto& tri : ntl->second) {
				int* tr = mt->triangleVertices(tri);
				for (int k = 0; k < 3; ++k)
					newTetV.insert(tr[k]);
			}
			int bestTet = -1;
			auto prStart = pr.first;
			struct oldTetV {
				int tetIdx;
				std::set<int> verts;
			};
			std::vector<oldTetV> oldTetVertices;
			while (prStart != pr.second) {
				auto otl = _oldVnTetTris.find(prStart->second);
				if (otl == _oldVnTetTris.end())
					throw(std::logic_error("Program error 4 in remapTetPhysics.\n"));
				oldTetVertices.push_back(oldTetV());
				oldTetVertices.back().tetIdx = prStart->second;
				auto& otv = oldTetVertices.back().verts;
				for (auto& triO : otl->second) {
					int* tr = mt->triangleVertices(triO);
					for (int k = 0; k < 3; ++k)
						otv.insert(tr[k]);
				}
				auto vit = newTetV.begin();
				for (; vit != newTetV.end(); ++vit) {
					if (otv.find(*vit) != otv.end()) {
						bestTet = prStart->second;
						break;
					}
				}
				if (bestTet > -1)
					break;
				++prStart;
			}
			if (bestTet < 0){  // least satisfying choice
				Vec3f newTetPos = { 0.0f, 0.0f, 0.0f };
				for (auto& nv : newTetV) {
					Vec3f tv;
					mt->getVertexCoordinate(nv, tv.xyz);
					newTetPos += tv;
				}
				newTetPos /= newTetV.size();
				float dsq, minD = FLT_MAX;
				for (auto& otv : oldTetVertices) {
					Vec3f oldTetPos = { 0.0f, 0.0f, 0.0f };
					for (auto& ov : otv.verts) {
						Vec3f tv;
						mt->getVertexCoordinate(ov, tv.xyz);
						oldTetPos += tv;
					}
					oldTetPos /= otv.verts.size();
					dsq = (newTetPos - oldTetPos).length2();
					if (minD > dsq) {
						minD = dsq;
						bestTet = otv.tetIdx;
					}
				}
			}
			auto& oN = _oldTets[bestTet];
			for (int j = 0; j < 4; ++j) {
				if (nodes[tn[j]])
					continue;
				if(prStart != pr.second)  // don't record the least satisfying correspondence. This means a vertex correspondence was found.
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
			throw(std::logic_error("Program error 5 in remapTetPhysics.\n"));
		return true;
	};
	// get node positions from vertexTet correspondence
	for (int n = _oldVertexTets.size(), i = 0; i < n; ++i) {
		if (newVnbt->_vertexTets[i] < 0)  // deleted vertex
			continue;
		auto tcNew = newVnbt->_tetCentroids[newVnbt->_vertexTets[i]];
		auto tcOld = _oldTetCentroids[_oldVertexTets[i]];
		const auto& nNew = newVnbt->_tetNodes[newVnbt->_vertexTets[i]];
		const auto& nOld = _oldTets[_oldVertexTets[i]];
		if (tcNew == tcOld) {
			for (int j = 0; j < 4; ++j){
				if(nodes[nNew[j]])
					continue;
				newVnbt->_nodeSpatialCoords[nNew[j]] = _oldNodePositions[nOld[j]];
				nodes[nNew[j]] = 1;
			}
			continue;
		}
		auto pr = _oldTetHash.equal_range(tcNew);
		int level = 1;
		while (pr.first == pr.second && level <= newVnbt->_tetSubdivisionLevels) {
			tcNew = newVnbt->centroidUpOneLevel(tcNew);
			++level;
			pr = _oldTetHash.equal_range(tcNew);
		}
		assert(level <= newVnbt->_tetSubdivisionLevels);
		if (std::distance(pr.first, pr.second) > 1) {  // this should be a macrotet which can't virtual node
			throw(std::logic_error("Program error 1 in remapTetPhysics.\n"));
		}
		else if (pr.first != pr.second) {
//			assert(_oldVertexTets[i] == pr.first->second);
			for (int j = 0; j < 4; ++j) {
				if (nodes[nNew[j]])
					continue;
				nodes[nNew[j]] = 1;
				Vec3f nodeLocus((const short(&)[3]) * newVnbt->_nodeGridLoci[nNew[j]].data());
				Vec3f bw, C;
				newVnbt->gridLocusToBarycentricWeight(nodeLocus, tcNew, bw);
				C = _oldNodePositions[nOld[0]] * (1.0f - bw[0] - bw[1] - bw[2]);
				C += _oldNodePositions[nOld[1]] * bw[0];
				C += _oldNodePositions[nOld[2]] * bw[1];
				C += _oldNodePositions[nOld[3]] * bw[2];
				newVnbt->_nodeSpatialCoords[nNew[j]] = C;
			}
		}
		else
			throw(std::logic_error("Program error 2 in remapTetPhysics.\n"));
	}
	// following 4 lines guaranteed tp be unique tets
	for(int i=0; i< newVnbt->_nMegatets; ++i)
		processTet(i);
	for (int n = newVnbt->_tetNodes.size(), i = newVnbt->_firstInteriorTet; i < n; ++i) 
		processTet(i);
	for (int i = newVnbt->_nMegatets; i < newVnbt->_firstInteriorTet; ++i)
		processTet(i);
	_oldNodePositions.clear();
	_oldTetCentroids.clear();
	_oldVertexTets.clear();
	_oldTets.clear();
	_oldNodes.clear();
	_oldTetHash.clear();
}

remapTetPhysics::remapTetPhysics()
{
}


remapTetPhysics::~remapTetPhysics()
{
}
