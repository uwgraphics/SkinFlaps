#include "bccTetDecimator.h"

#include <iostream>
#include <map>
#include <algorithm>
#include <unordered_set>
#include "vnBccTetCutter.h"
#include "materialTriangles.h"

void bccTetDecimator::createMacroTets(vnBccTetCutter* vbtc, int nLevels) {
	_vbtc = vbtc;
	--nLevels;
	int mult = (1 << nLevels);
	for (auto& n : _nodeGridLoci) {
		for (auto i = 0; i < 3; ++i)
			n[i] <<= nLevels;
	}
	_tetHash.clear();
	_tetHash.reserve(_tetCentroids.size() * 1.5);
	for (int n = _tetCentroids.size(), i = 0; i < n; ++i) {
		auto& tc = _tetCentroids[i];
		for (int j = 0; j < 3; ++j)
			tc[j] <<= nLevels;
		_tetHash.insert(std::make_pair(tc, i));
	}
	_unitSpacingInv *= mult;
	_unitSpacing = 1.0 / _unitSpacingInv;
	for(int i=0; i< 3; ++i)
		_gridSize[i] *= mult;  // Number of subdivisions in x,y, and z.
//	static Mat3x3f _barycentricInverses[6];
	for (auto& mc : _vbtc->_vMatCoords) {
		for (int i = 0; i < 3; ++i)
			mc[i] *= mult;
	}
	for (auto& vtc : vbtc->_vertexTetCentroids) {
		for (int i = 0; i < 3; ++i)
			vtc[i] *= mult;
	}
	for(auto &gs : vbtc->_gridSize)
		gs *= mult;

}

void bccTetDecimator::decimate2(int level) {
	if (level < 2)
		throw(std::logic_error("Tet decimation below level 2 is requested.\n"));
	_decimatedNodes.clear();
	decNodeConstraints.clear();
	int seed = 0, levelEnd = _tetNodes.size(), levelNow = 1;
	while (levelNow < level) {
		// to mark a tet for deletion its _tetNode[0] set -1 and its _tetCentroids[0] set USHRT_MAX.  Also is removed from _tetHash.
		while (seed < levelEnd) {
			if (_tetNodes[seed][0] > -1 && getLevel(_tetCentroids[seed]) == levelNow) {  // decimation candidate.  Do I need second test?
				bccTetCentroid subtets[8], macroTet = centroidUpOneLevel(_tetCentroids[seed]);
				subtetCentroids(macroTet, subtets);
				THIT stIt[8];
				int i;
				for (i = 0; i < 8; ++i) {
					if (subtets[i][0] == USHRT_MAX)
						break;
					auto pr = _tetHash.equal_range(subtets[i]);
					if (std::distance(pr.first, pr.second) != 1)  // must exist and not be virtual noded
						break;
					stIt[i] = pr.first;
				}
				if (i > 7) {  // macrotet found.  Decimate
					int hc, macroLevelUp;
					centroidHalfAxisSize(macroTet, hc, macroLevelUp);
					macroLevelUp <<= 1;
					//				bool up = (macroTet[hc] & macroLevelUp) == (macroTet[(hc + 2) % 3] & macroLevelUp) ? true : false;

					std::array<int, 4> nodes;
					nodes[0] = _tetNodes[stIt[0]->second][0];
					nodes[1] = _tetNodes[stIt[1]->second][1];
					nodes[2] = _tetNodes[stIt[2]->second][2];
					nodes[3] = _tetNodes[stIt[3]->second][3];

					short gl[4][3];  // debug
					centroidToNodeLoci(macroTet, gl);
					for (int j = 0; j < 4; ++j) {
						auto nl = _nodeGridLoci[nodes[j]];
						assert(nl[0] == gl[j][0] && nl[1] == gl[j][1] && nl[2] == gl[j][2]);
					}

					// now enter possibly decimated edge nodes
					_decimatedNodes.insert(std::make_pair(_tetNodes[stIt[0]->second][1], std::make_pair(nodes[0], nodes[1])));
					_decimatedNodes.insert(std::make_pair(_tetNodes[stIt[0]->second][2], std::make_pair(nodes[0], nodes[2])));
					_decimatedNodes.insert(std::make_pair(_tetNodes[stIt[0]->second][3], std::make_pair(nodes[0], nodes[3])));
					_decimatedNodes.insert(std::make_pair(_tetNodes[stIt[1]->second][2], std::make_pair(nodes[1], nodes[2])));
					_decimatedNodes.insert(std::make_pair(_tetNodes[stIt[1]->second][3], std::make_pair(nodes[1], nodes[3])));
					_decimatedNodes.insert(std::make_pair(_tetNodes[stIt[2]->second][3], std::make_pair(nodes[2], nodes[3])));

					// debug
					for (auto& dn : _decimatedNodes) {
						Vec3f n((short(&)[3]) * _nodeGridLoci[dn.first].data()), mean((short(&)[3]) * _nodeGridLoci[dn.second.first].data());
						mean += Vec3f((short(&)[3]) * _nodeGridLoci[dn.second.second].data());
						mean *= 0.5f;
						assert(n == mean);
					}

					// mark subtets deleted and add macroTet
					for (i = 0; i < 8; ++i) {
						_tetNodes[stIt[i]->second][0] = -1;
				//	Still need this for vertex tet assignment	_tetCentroids[stIt[i]->second][0] = USHRT_MAX;
						_tetHash.erase(stIt[i]->first);
					}
					_tetHash.insert(std::make_pair(macroTet, (int)_tetNodes.size()));
					_tetNodes.push_back(nodes);
					_tetCentroids.push_back(macroTet);
				}
			}
			++seed;
		}
		assert(seed == levelEnd);
		levelEnd = _tetNodes.size();
		++levelNow;
	}
	// reassign vertex tet and baryweight if original tet decimated
	for (int n = _vertexTets.size(), j, i = 0; i < n; ++i) {
		if (_tetNodes[_vertexTets[i]][0] < 0) {  // its tet was decimated, so could only have been a unique tet for its locus
			bccTetCentroid tc = _tetCentroids[_vertexTets[i]];
			Vec3f gridLocus;
			barycentricWeightToGridLocus(tc, _barycentricWeights[i], gridLocus);
			for (j = 0; j < 16; ++j) {  // will never get 16 levels high
				tc = centroidUpOneLevel(tc);
				auto tet = _tetHash.find(tc);
				if (tet != _tetHash.end()) {
					_vertexTets[i] = tet->second;
					break;
				}
			}
			if (j > 15)
				throw(std::logic_error("Program error in creating decimated tet heirarchy.\n"));
			gridLocusToBarycentricWeightLevel(gridLocus, tc, _barycentricWeights[i]);
		}
	}
	// pack decimated tets amd nodes
	std::vector<int> nodeMap, tetMap;
	nodeMap.assign(_nodeGridLoci.size(), -1);
	tetMap.assign(_tetNodes.size(), -1);
	int offset = 0;
	for (int n = _tetNodes.size(), i = 0; i < n; ++i) {
		auto& tn = _tetNodes[i];
		if (tn[0] > -1) {
			for (int j = 0; j < 4; ++j)
				nodeMap[tn[j]] = 1;
			_tetNodes[offset] = _tetNodes[i];
			_tetCentroids[offset] = _tetCentroids[i];
			tetMap[i] = offset;
			++offset;
		}
	}
	for (int n = _vertexTets.size(), i = 0; i < n; ++i)
		_vertexTets[i] = tetMap[_vertexTets[i]];
	tetMap.clear();
	_tetNodes.resize(offset);
	_tetCentroids.resize(offset);
	_tetHash.clear();
	_tetHash.reserve(offset);
	for (int i = 0; i < offset; ++i)
		_tetHash.insert(std::make_pair(_tetCentroids[i], i));
	offset = 0;
	for (int n = nodeMap.size(), i = 0; i < n; ++i) {
		if (nodeMap[i] > -1) {
			nodeMap[i] = offset;
			_nodeGridLoci[offset] = _nodeGridLoci[i];
			++offset;
		}
	}
	_nodeGridLoci.resize(offset);
	for (int n = _tetNodes.size(), i = 0; i < n; ++i) {
		auto& tn = _tetNodes[i];
		for (int j = 0; j < 4; ++j)
			tn[j] = nodeMap[tn[j]];
	}
	// binary tree structure for decimated nodes follows.  Only need leaves at the end.
	struct bNode {
		int node;
		DNIT nit;
		float depth;
	};
	std::list<bNode> parents, leaves;
	auto processFirstParent = [&]() {
		auto pn = parents.begin();
		bNode child;
		child.depth = pn->depth * 0.5f;
		child.node = pn->nit->second.first;
		child.nit = _decimatedNodes.find(child.node);
		if (child.nit == _decimatedNodes.end())
			leaves.push_back(child);
		else
			parents.push_back(child);
		child.node = pn->nit->second.second;
		child.nit = _decimatedNodes.find(child.node);
		if (child.nit == _decimatedNodes.end())
			leaves.push_back(child);
		else
			parents.push_back(child);
		parents.erase(pn);
	};
	for (auto dn = _decimatedNodes.begin(); dn != _decimatedNodes.end(); ++dn) {
		if (dn->first < 0)
			continue;
		int dNode = nodeMap[dn->first];
		if (dNode > -1) {  // used
			parents.clear();
			leaves.clear();
			// collect leaves of binary tree of decimated nodes
			bNode root;
			root.depth = 1.0f;
			root.nit = dn;
			root.node = dn->first;
			parents.push_back(root);
			while (!parents.empty())
				processFirstParent();
			auto dnc = decNodeConstraints.insert(std::make_pair(dNode, decimatedFaceNode()));
			auto& fn = dnc.first->second.faceNodes;
			auto& fp = dnc.first->second.faceParams;
			if (leaves.size() < 3) {
				fn.reserve(2);
				fn.push_back(nodeMap[leaves.front().node]);
				fn.push_back(nodeMap[leaves.back().node]);
				fp.assign(2, 0.5f);
			}
			else {
				std::map<int, float> dnVerts;
				for (auto& l : leaves) {
					auto pr = dnVerts.insert(std::make_pair(l.node, l.depth));
					if (!pr.second)
						pr.first->second += l.depth;
				}
				fn.reserve(dnVerts.size());
				fp.reserve(dnVerts.size());
				for (auto& dnv : dnVerts) {
					fn.push_back(nodeMap[dnv.first]);
					fp.push_back(dnv.second);
				}
			}
		}
	}
	_decimatedNodes.clear();
}

void bccTetDecimator::decimate(int level, int nSubtetsRequired, bool onlyInteriorTets) {
	if (level < 2)
		throw(std::logic_error("Tet decimation below level 2 is requested.\n"));
	decNodeConstraints.clear();
	_uniqueCornerNodes.clear();
	_decimatedNodes.clear();
	// next section generates all the super tet centroids generated at the requested level.
	// Later these could be loaded as part of the scene file or computed once at object load.
	// make level Cartesian spacing
	unsigned short levelBit = 2 << level;
	unsigned short halfLevelBit = levelBit >> 1, quarterLevelBit = levelBit >> 2;
	std::vector<bccTetCentroid> bigTets;
	bigTets.reserve(_tetNodes.size()/(levelBit << 1));
	int dims[3];
	for (int i = 0; i < 3; ++i)  // double since going directly to tet centroids.  COURT check later with higher levels.
		dims[i] = (_gridSize[i]) << 1;
	for (unsigned short i = 0; i < dims[0]; i += levelBit) {
		for (unsigned short j = 0; j < dims[1]; j += levelBit) {
			for (unsigned short k = 0; k < dims[2]; k += levelBit) {
				// make all possible tets at this level
				for (int c1 = 0; c1 < 3; ++c1) {
					int c2 = (c1 + 1) % 3, c0 = (c1 + 2) % 3;
					bccTetCentroid tc = { i, j, k };
					tc[c1] += halfLevelBit;
					for (int n = 0; n < 4; ++n) {
						bccTetCentroid c = tc;
						if(n&1)
							c[c0] += (n < 2 ? -quarterLevelBit : quarterLevelBit);
						else 
							c[c2] += (n < 2 ? -quarterLevelBit : quarterLevelBit);
						if (c[c0] < 65535 - halfLevelBit && c[c2] < 65535 - halfLevelBit) {
							bigTets.push_back(c);
						}
					}
				}
			}
		}
	}
	std::cout << "tet number was " << _tetNodes.size() << " and node number was " << _nodeGridLoci.size() << " before decimation and\n";
	for (int n = bigTets.size(), i = 0; i < n; ++i)
		recurseSubtets(bigTets[i], nSubtetsRequired, onlyInteriorTets ? _firstInteriorTet : 0);
	// reassign vertex tet and baryweight if original tet decimated
	for (int n = _vertexTets.size(), j, i = 0; i < n; ++i) {
		if (_tetNodes[_vertexTets[i]][0] < 0) {  // its tet was decimated, so could only have been a unique tet for its locus
			bccTetCentroid tc, oldCent = _tetCentroids[_vertexTets[i]];
			Vec3f gridLocus;
			barycentricWeightToGridLocus(oldCent, _barycentricWeights[i], gridLocus);
			for (j = 0; j < 16; ++j) {  // will never get 16 levels high
				centroidPromoteOneLevel(oldCent, tc);
				auto tet = _tetHash.find(tc);
				if (tet != _tetHash.end()) {
					_vertexTets[i] = tet->second;
					break;
				}
				oldCent = tc;
			}
			if (j > 15)
				throw(std::logic_error("Program error in creating decimated tet heirarchy.\n"));
			gridLocusToBarycentricWeightLevel(gridLocus, tc, _barycentricWeights[i]);
		}
	}
	// pack decimated tets amd nodes
	std::vector<int> nodeMap, tetMap;
	nodeMap.assign(_nodeGridLoci.size(), -1);
	tetMap.assign(_tetNodes.size(), -1);
	int offset = 0;
	for (int n = _tetNodes.size(), i = 0; i < n; ++i) {
		auto& tn = _tetNodes[i];
		if (tn[0] > -1) {
			for (int j = 0; j < 4; ++j)
				nodeMap[tn[j]] = 1;
			_tetNodes[offset] = _tetNodes[i];
			_tetCentroids[offset] = _tetCentroids[i];
			tetMap[i] = offset;
			++offset;
		}
	}
	for (int n = _vertexTets.size(), i = 0; i < n; ++i)
		_vertexTets[i] = tetMap[_vertexTets[i]];
	_tetNodes.resize(offset);
	_tetCentroids.resize(offset);
	_tetHash.clear();
	_tetHash.reserve(offset);
	for (int i = 0; i < offset; ++i)
		_tetHash.insert(std::make_pair(_tetCentroids[i], i));
	offset = 0;
	for (int n = nodeMap.size(), i = 0; i < n; ++i) {
		if (nodeMap[i] > -1) {
			nodeMap[i] = offset;
			_nodeGridLoci[offset] = _nodeGridLoci[i];
			++offset;
		}
	}
	_nodeGridLoci.resize(offset);
	for (int n = _tetNodes.size(), i = 0; i < n; ++i) {
		auto& tn = _tetNodes[i];
		for (int j = 0; j < 4; ++j)
			tn[j] = nodeMap[tn[j]];
	}
	std::cout << "after " << level << " level decimation with nSubtetsRequired at " << nSubtetsRequired << "\ntet number was " << _tetNodes.size() << " and node number was " << _nodeGridLoci.size() << "\n";
	_uniqueCornerNodes.clear();
	// binary tree structure for decimated nodes follows.  Only need leaves at the end.
	struct bNode {
		int node;
		DNIT nit;
		float depth;
	};
	std::list<bNode> parents, leaves;
	auto processFirstParent = [&]() {
		auto pn = parents.begin();
		bNode child;
		child.depth = pn->depth * 0.5f;
		child.node = pn->nit->second.first;
		child.nit = _decimatedNodes.find(child.node);
		if (child.nit == _decimatedNodes.end())
			leaves.push_back(child);
		else
			parents.push_back(child);
		child.node = pn->nit->second.second;
		child.nit = _decimatedNodes.find(child.node);
		if (child.nit == _decimatedNodes.end())
			leaves.push_back(child);
		else
			parents.push_back(child);
		parents.erase(pn);
	};
	for (auto dn = _decimatedNodes.begin(); dn != _decimatedNodes.end(); ++dn) {
		if (dn->first < 0)
			continue;
		int dNode = nodeMap[dn->first];
		if (dNode > -1) {  // used
			parents.clear();
			leaves.clear();
			// collect leaves of binary tree of decimated nodes
			bNode root;
			root.depth = 1.0f;
			root.nit = dn;
			root.node = dn->first;
			parents.push_back(root);
			while (!parents.empty())
				processFirstParent();
			auto dnc = decNodeConstraints.insert(std::make_pair(dNode, decimatedFaceNode()));
			auto& fn = dnc.first->second.faceNodes;
			auto& fp = dnc.first->second.faceParams;
			if (leaves.size() < 3) {
				fn.reserve(2);
				fn.push_back(nodeMap[leaves.front().node]);
				fn.push_back(nodeMap[leaves.back().node]);
				fp.assign(2, 0.5f);
			}
			else {
				std::map<int, float> dnVerts;
				for (auto& l : leaves) {
					auto pr = dnVerts.insert(std::make_pair(l.node, l.depth));
					if (!pr.second)
						pr.first->second += l.depth;
				}
				fn.reserve(dnVerts.size());
				fp.reserve(dnVerts.size());
				for (auto& dnv : dnVerts) {
					fn.push_back(nodeMap[dnv.first]);
					fp.push_back(dnv.second);
				}
			}
		}
	}
	_decimatedNodes.clear();
}

void bccTetDecimator::gridLocusToBarycentricWeightLevel(const Vec3f& gridLocus, const bccTetCentroid& tc, Vec3f& barycentricWeight)
{
	int level = getLevel(tc);
	assert(level > 1);
	int unitSpacing = 1;
	for (int i = 1; i < level; ++i)
		unitSpacing <<= 1;
	Vec3f B(gridLocus);
	// set barycentric coordinate within that tet
	int baryInv, hc = tc[0] & unitSpacing ? 0 : (tc[1] & unitSpacing ? 1 : 2);
	short xyz[3] = { tc[0] >> 1, tc[1] >> 1, tc[2] >> 1};
	if ((xyz[hc] & unitSpacing) == (xyz[(hc + 2) % 3] & unitSpacing)) {  // main axis below secondary
		baryInv = hc;
		// subtract grid locus of first point of tet
		B -= Vec3f((const short(&)[3]) xyz);
		B[hc] += (unitSpacing >> 1);
	}
	else {
		baryInv = hc + 3;
		B -= Vec3f((const short(&)[3])xyz);
		B[hc] -= (unitSpacing >> 1);
	}
	B[(hc + 1) % 3] += unitSpacing;
	B /= (float)unitSpacing;  // to use old one unit inverses
	barycentricWeight = _barycentricInverses[baryInv] * B;
	assert(barycentricWeight[0] > 0.0f && barycentricWeight[1] > 0.0f && barycentricWeight[2] > 0.0f && barycentricWeight[0] < 1.0f && barycentricWeight[1] < 1.0f && barycentricWeight[2] < 1.0f);
}


void bccTetDecimator::centroidPromoteOneLevel(const bccTetCentroid& tcLevelIn, bccTetCentroid& tc) {
	int levelBit = 1, levelX2, levelX4, level = getLevel(tcLevelIn);
	for (int i = 1; i < level; ++i)
		levelBit <<= 1;
	levelX2 = levelBit << 1;
	levelX4 = levelX2 << 1;
	int hc = -1;
	for (int i = 0; i < 3; ++i)
		if (tcLevelIn[i] & levelBit) {
			hc = i;
			break;
		}
	int c1 = (hc + 1) % 3, c2 = (hc + 2) % 3;
	tc = tcLevelIn;
	assert((tc[c1] & levelX2) != (tc[c2] & levelX2));
	// none of the 4 core subtets have the same hc as the supertet
	// if making tc[hc] a multiple of 4 with a one unit move creates a valid level up tet, is a center core subtet
	tc[hc] += (tc[hc] & levelX2) ? levelBit : -levelBit;
	if (tc[c1] & levelX2) {
		if ((tc[hc] & levelX4) != (tc[c2] & levelX4))  // valid level 2
			return;
	}
	if (tc[c2] & levelX2) {
		if ((tc[hc] & levelX4) != (tc[c1] & levelX4))  // valid level 2
			return;
	}
	// is corner subtet and not core so will have same hc axis in level 2
	tc[hc] = tcLevelIn[hc];
	tc[hc] += tc[hc] & levelX2 ? -levelBit : levelBit;
	if (tc[c1] & levelX2) {
		if (tc[c2] & levelX4) {
			tc[c1] += tc[c1] & levelX4 ? levelX2 : -levelX2;
		}
		else {
			tc[c1] += tc[c1] & levelX4 ? -levelX2 : levelX2;
		}
	}
	else {
		assert(tc[c2] & levelX2);
		if (tc[c1] & levelX4) {
			tc[c2] += tc[c2] & levelX4 ? levelX2 : -levelX2;
		}
		else {
			tc[c2] += tc[c2] & levelX4 ? -levelX2 : levelX2;
		}
	}
}

int bccTetDecimator::recurseSubtets(const bccTetCentroid &tc, const int &nNeeded, int firstInteriorTet) {

	if (tc[0] == 32 && tc[1] == 72 && tc[2] == 28)
		int junk = 0;

	int level, hc, c1, c2;
	for (level = 2; level < INT_MAX; level <<= 1) {
		if (tc[0] & level) {
			hc = 0;  c1 = 1;  c2 = 2;
			break;
		}
		else if (tc[1] & level) {
			hc = 1;  c1 = 2;  c2 = 0;
			break;
		}
		else if (tc[2] & level) {
			hc = 2;  c1 = 0;  c2 = 1;
			break;
		}
		else
			;
	}
	int halfLevel = level >> 1;
	// the cyclic geometry of the bcc tet is reflected in it's centroid algebra.
	bool up = (tc[c2] & (level << 1)) == (tc[hc] & (level<<1));
	int ntets = 0;
	bccTetCentroid c, ic;
	int subtets[8] = {-1,-1,-1,-1,-1,-1,-1,-1};  // none found. First 4 are corners 0-3. Last 4 are central cores
	auto findTet1 = [&]() ->int {  // processes level 1 tets


		if (c[0] == 32 && c[1] == 74 && c[2] == 31)
			int junk = 0;



		auto pr = _tetHash.equal_range(c);
		int d = (int)std::distance(pr.first, pr.second);
		if (d < 1)
			return -1;
		else if (d < 2) {
			if (pr.first->second < firstInteriorTet)
				return -1;
			else
				return pr.first->second;
		}
				else
			return -2;  // signals a virtual noded tet locus
	};
	auto getTet = [&](int subNum) ->bool {  // if true this is a level 1 virtual noded tet
		if(level > 2)
			subtets[subNum] = recurseSubtets(c, nNeeded, firstInteriorTet);
		else {
			subtets[subNum] = findTet1();
			if (subtets[subNum] < -1)
				return true;
		}
		if (subtets[subNum] > -1)
			++ntets;
		return false;
	};
	if (tc[hc] > 0){
		ic = tc;
		ic[hc] -= halfLevel;
		if (up) {
			if (tc[c1] > 0) {
				c = ic;
				c[c1] -= level;
				if (getTet(0))  // any level 1 tet that is virtual noded won't allow decimation into a level 2 so stop looking.  If above level 2 keep looking.
					return -2;
			}
			c = ic;
			c[c1] += level;
			if (getTet(1))
				return -2;
		}
		else {
			if (tc[c2] > 0) {
				c = ic;
				c[c2] -= level;
				if (getTet(2))
					return -2;
			}
			c = ic;
			c[c2] += level;
			if (getTet(3))
				return -2;
		}
	}
	ic = tc;
	ic[hc] += halfLevel;
	if (up) {
		c = ic;
		c[c2] += level;
		if (getTet(2))
			return -2;
		if (tc[c2] > 0) {
			c = ic;
			c[c2] -= level;
			if (getTet(3))
				return -2;
		}
	}
	else {
		if (tc[c1] > 0) {
			c = ic;
			c[c1] -= level;
			if (getTet(0))
				return -2;
		}
		c = ic;
		c[c1] += level;
		if (getTet(1))
			return -2;
	}
	// now have collected corner tets. Get the central 4 tet core around the hc axis.  Same for both up and down.
	if (tc[c1] > 0) {
		c = tc;
		c[c1] -= halfLevel;
		if (getTet(4))
			return -2;
	}
	c = tc;
	c[c1] += halfLevel;
	if (getTet(5))
		return -2;
	if (tc[c2] > 0) {
		c = tc;
		c[c2] -= halfLevel;
		if (getTet(6))
			return -2;
	}
	c = tc;
	c[c2] += halfLevel;
	if (getTet(7))
		return -2;
	// check for conditions blocking decimation
	if (level > 2) {
		for (int i = 0; i < 8; ++i)
			if (subtets[i] < -1)
				return -2;  // processing of lower levels done, but this level tet has virtual noding
	}
	if(ntets < nNeeded)
		return -1;
	// decimate this tet
	auto findMakeNewNode = [&](const std::array<short, 3>& locus) ->int {
		auto pr = _uniqueCornerNodes.insert(std::make_pair(locus, (int)_nodeGridLoci.size()));
		if (pr.second)
			_nodeGridLoci.push_back(locus);
		return pr.first->second;
	};
	std::array<int, 4>  dNodes = { -1, -1, -1, -1 };
	std::list<int> edgeSplitNodes[6];
	auto setESnode = [&](int esIdx, int node) {
		for (auto esx : edgeSplitNodes[esIdx])
			if(esx == node)
				return;
		edgeSplitNodes[esIdx].push_back(node);
	};
	if (subtets[0] > -1) {
		auto &tn = _tetNodes[subtets[0]];
		dNodes[0] = tn[0];
		setESnode(0, tn[1]);
		setESnode(1, tn[2]);
		setESnode(2, tn[3]);
	}
	else {  // since this tet empty can create a unique node for this level
		std::array<short, 3> newLoc;
		if(up)
			newLoc[hc] = (tc[hc] >> 1) - halfLevel;
		else
			newLoc[hc] = (tc[hc] >> 1) + halfLevel;
		newLoc[c1] = (tc[c1] >> 1) - level;
		newLoc[c2] = tc[c2] >> 1;
		dNodes[0] = findMakeNewNode(newLoc);
	}
	if (subtets[1] > -1) {
		auto& tn = _tetNodes[subtets[1]];
		dNodes[1] = tn[1];
		setESnode(0, tn[0]);
		setESnode(3, tn[2]);
		setESnode(4, tn[3]);
	}
	else {
		std::array<short, 3> newLoc;
		if (up)
			newLoc[hc] = (tc[hc] >> 1) - halfLevel;
		else
			newLoc[hc] = (tc[hc] >> 1) + halfLevel;
		newLoc[c1] = (tc[c1] >> 1) + level;
		newLoc[c2] = (tc[c2] >> 1);
		dNodes[1] = findMakeNewNode(newLoc);
	}
	if (subtets[2] > -1) {
		auto& tn = _tetNodes[subtets[2]];
		dNodes[2] = tn[2];
		setESnode(1, tn[0]);
		setESnode(3, tn[1]);
		setESnode(5, tn[3]);
	}
	else {
		std::array<short, 3> newLoc;
		if (up)
			newLoc[hc] = (tc[hc] >> 1) + halfLevel;
		else
			newLoc[hc] = (tc[hc] >> 1) - halfLevel;
		newLoc[c1] = (tc[c1] >> 1);
		if (up)
			newLoc[c2] = (tc[c2] >> 1) + level;
		else
			newLoc[c2] = (tc[c2] >> 1) - level;
		dNodes[2] = findMakeNewNode(newLoc);
	}
	if (subtets[3] > -1) {
		auto& tn = _tetNodes[subtets[3]];
		dNodes[3] = tn[3];
		setESnode(2, tn[0]);
		setESnode(4, tn[1]);
		setESnode(5, tn[2]);
	}
	else {
		std::array<short, 3> newLoc;
		if (up)
			newLoc[hc] = (tc[hc] >> 1) + halfLevel;
		else
			newLoc[hc] = (tc[hc] >> 1) - halfLevel;
		newLoc[c1] = (tc[c1] >> 1);
		if (up)
			newLoc[c2] = (tc[c2] >> 1) - level;
		else
			newLoc[c2] = (tc[c2] >> 1) + level;
		dNodes[3] = findMakeNewNode(newLoc);
	}
	assert(dNodes[0] > -1 && dNodes[1] > -1 && dNodes[2] > -1 && dNodes[3] > -1);
	for (int i = 0; i < 4; ++i)
		_uniqueCornerNodes.insert(std::make_pair(_nodeGridLoci[dNodes[i]], dNodes[i]));
	//  central core will also generate edgeSplitNodes[]
	if (up) {
				if (subtets[4] > -1) {
					auto& tn = _tetNodes[subtets[4]];
					setESnode(2, tn[0]);
					setESnode(1, tn[1]);
					setESnode(5, tn[2]);
					setESnode(0, tn[3]);
				}
				if (subtets[5] > -1) {
					auto &tn = _tetNodes[subtets[5]];
					setESnode(4, tn[0]);
					setESnode(3, tn[1]);
					setESnode(0, tn[2]);
					setESnode(5, tn[3]);
				}
				if (subtets[6] > -1) {
					auto &tn = _tetNodes[subtets[6]];
					setESnode(0, tn[0]);
					setESnode(5, tn[1]);
					setESnode(2, tn[2]);
					setESnode(4, tn[3]);
				}
				if (subtets[7] > -1) {
					auto& tn = _tetNodes[subtets[7]];
					setESnode(0, tn[0]);
					setESnode(5, tn[1]);
					setESnode(3, tn[2]);
					setESnode(1, tn[3]);
				}
			}
	else {
		if (subtets[4] > -1) {
			auto& tn = _tetNodes[subtets[4]];
			setESnode(1, tn[0]);
			setESnode(2, tn[1]);
			setESnode(0, tn[2]);
			setESnode(5, tn[3]);
		}
		if (subtets[5] > -1) {
			auto& tn = _tetNodes[subtets[5]];
			setESnode(3, tn[0]);
			setESnode(4, tn[1]);
			setESnode(5, tn[2]);
			setESnode(0, tn[3]);
		}
		if (subtets[6] > -1) {
			auto& tn = _tetNodes[subtets[6]];
			setESnode(5, tn[0]);
			setESnode(0, tn[1]);
			setESnode(1, tn[2]);
			setESnode(3, tn[3]);
		}
		if (subtets[7] > -1) {
			auto& tn = _tetNodes[subtets[7]];
			setESnode(5, tn[0]);
			setESnode(0, tn[1]);
			setESnode(4, tn[2]);
			setESnode(2, tn[3]);
		}
	}
	for (int i = 0; i < 8; ++i) {
		if (subtets[i] > -1) {
			_tetHash.erase(_tetCentroids[subtets[i]]);
			_tetNodes[subtets[i]][0] = -1;  // mark deleted. Remove corresponding _tetCentroid during pack.
		}
	}
	int superTet = (int)_tetNodes.size();
	_tetHash.insert(std::make_pair(tc, superTet));
	_tetNodes.push_back(dNodes);
	_tetCentroids.push_back(tc);
	for (int i = 0; i < 6; ++i) {
		for(auto esn : edgeSplitNodes[i]){
			std::pair<int, int> ip;
			if (i < 3) {
				ip.first = dNodes[0];
				ip.second = dNodes[i + 1];
			}
			else if (i < 5) {
				ip.first = dNodes[1];
				ip.second = dNodes[i - 1];
			}
			else {
				ip.first = dNodes[2];
				ip.second = dNodes[3];
			}
			_decimatedNodes.insert(std::make_pair(esn, ip));  // first is id of node decimated, second are the two nodes of the edge decimated node split.
		}
	}
	return superTet;
}

void bccTetDecimator::centroidToNodeLociLevels(const bccTetCentroid& centroid, short(&gridLoci)[4][3]) {
	int levelBit = 1, level = getLevel(centroid);
	for (int i = 1; i < level; ++i)
		levelBit <<= 1;
	int c1, c2, hc = centroid[0] & levelBit ? 0 : (centroid[1] & levelBit ? 1 : 2);
	c1 = hc < 2 ? hc + 1 : 0;
	c2 = hc > 0 ? hc - 1 : 2;
	auto tc = centroid;
	for (int i = 0; i < 3; ++i)
		tc[i] >>= 1;
	if (level > 1) // tc[hc] 0.5 not shifted out
		tc[hc] -= (levelBit >> 1);
	bool below01 = ((tc[hc] & levelBit) != (tc[c2] & levelBit));
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 4; ++j)
			gridLoci[j][i] = tc[i];
	}
	if (below01) {
		gridLoci[0][hc] += levelBit;
		gridLoci[1][hc] += levelBit;
		gridLoci[2][c2] -= levelBit;
		gridLoci[3][c2] += levelBit;
	}
	else {
		gridLoci[2][hc] += levelBit;
		gridLoci[3][hc] += levelBit;
		gridLoci[2][c2] += levelBit;
		gridLoci[3][c2] -= levelBit;
	}
	gridLoci[0][c1] -= levelBit;
	gridLoci[1][c1] += levelBit;
}

void bccTetDecimator::barycentricWeightToGridLocusLevels(const bccTetCentroid& tetCentroid, const Vec3f& barycentricWeight, Vec3f& gridLocus)
{
	short gridLoci[4][3];
	centroidToNodeLociLevels(tetCentroid, gridLoci);
	gridLocus = Vec3f((const short(&)[3])gridLoci[0]) * (1.0f - barycentricWeight.X - barycentricWeight.Y - barycentricWeight.Z);
	for (int i = 1; i < 4; ++i)
		gridLocus += Vec3f((const short(&)[3])gridLoci[i]) * barycentricWeight[i - 1];
}
