#include "multiresBccTets.h"

#include <iostream>
#include <map>
#include <algorithm>
#include <unordered_set>
#include "Mat2x2f.h"
//#include "vnBccTetCutter.h"
#include "materialTriangles.h"

/* void multiresBccTets::createMacroTets(materialTriangles* mt, const int nLevels, const int maximumDimensionMacroSubdivs) {  // takes current low res tets (level 1) and promotes them to level 4 for future subdivision 
	_mt = mt;
	_vbtc.makeFirstVnTets(_mt, this, maximumDimensionMacroSubdivs);
	int mult = (1 << (nLevels-1));
	for (auto& n : _nodeGridLoci) {
		for (auto i = 0; i < 3; ++i)
			n[i] *= mult;  // can be negative
	}
	_tetHash.clear();
	_tetHash.reserve(_tetCentroids.size() * 1.5);
	for (int n = _tetCentroids.size(), i = 0; i < n; ++i) {
		auto& tc = _tetCentroids[i];
		for (int j = 0; j < 3; ++j)
			tc[j] <<= (nLevels-1);
		_tetHash.insert(std::make_pair(tc, i));
	}
	_unitSpacingInv *= mult;
	_unitSpacing = 1.0 / _unitSpacingInv;

	for (int i = 0; i < 3; ++i)
		_gridSize[i] *= mult;  // Number of subdivisions in x,y, and z.
	//	static Mat3x3f _barycentricInverses[6];  // COURT fix
	//	for (auto& mc : _vbtc->_vMatCoords) {
	//		for (int i = 0; i < 3; ++i)
	//			mc[i] *= mult;
	//	}
	//	for (auto& vtc : vbtc->_vertexTetCentroids) {
	//		for (int i = 0; i < 3; ++i)
	//			vtc[i] *= mult;
	//	}
	//	for(auto &gs : vbtc->_gridSize)
	//		gs *= mult;

	// macrotets guaranteed not to virtual node.  Subcut any found at this stage.
	std::vector<char> vnTets;
	vnTets.assign(_tetCentroids.size(), 0x00);
	auto th = _tetHash.begin();
	auto hNext = th;
	do {
		++hNext;
		if (hNext == _tetHash.end() || th->first != hNext->first)
			th = hNext;
		else {
			vnTets[hNext->second] = 1;
			vnTets[th->second] = 1;
		}
	} while (th != _tetHash.end());
	std::vector<char> vnVerts;
	vnVerts.assign(_mt->numberOfVertices(), 0x00);
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i)
		if (vnTets[_vertexTets[i]])
			vnVerts[i] = 1;
	// get unique tet faces at the boundary of the vnTets to be removed in contact with tets that remain.
	std::set<std::array<int, 3> > boundingTris; // of vnTets
	std::set<int> boundingNodes, deadNodes;
	for (int n = vnTets.size(), i = 0; i < n; ++i) {
		if (!vnTets[i])
			continue;
		auto &tn = _tetNodes[i];
		for (int k, j = 0; j < 4; ++j) {
			deadNodes.insert(tn[j]);
			bccTetCentroid adjTc;
			int adjFace = faceAdjacentMacrotet(_tetCentroids[i], j, adjTc);
			if (adjFace < 0)
				continue;
			auto thit = _tetHash.find(adjTc);
			if (thit == _tetHash.end() || vnTets[thit->second])  // virtual noded
				continue;
			std::set<int> adjNodes;
			for (k = 0; k < 3; ++k)
				adjNodes.insert(_tetNodes[thit->second][(adjFace + k) & 3]);
			for (k = 0; k < 3; ++k) {
				if (adjNodes.find(tn[(i + k) & 3]) == adjNodes.end())
					break;
			}
			if (k < 3)
				continue;
			// this is a unique bounding tet face of the remaining solid. Order tri so pointing into the excised solid.
			std::array<int, 3> face;
			face[1] = tn[(j + 1) & 3];
			if (j & 1) {
				face[0] = tn[j];
				face[2] = tn[(j + 2) & 3];
			}
			else {
				face[2] = tn[j];
				face[0] = tn[(j + 2) & 3];
			}
			boundingNodes.insert(face.begin(), face.end());
			boundingTris.insert(std::move(face));
		}
	}
	for(auto bn : boundingNodes)
		deadNodes.erase(bn);
	boundingNodes.clear();
	int junk = 0;
	vnTets.clear();
	std::vector<int> vnTris;
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
		int* tr = _mt->triangleVertices(i);
		for (int j = 0; j < 3; ++j)
			if (vnVerts[tr[j]]) {
				vnTris.push_back(i);
				break;
			}
	}
	vnVerts.clear();
	std::vector<vnBccTetCutter_tbb::zIntrsct> zi_loc;
	_vbtc.zIntersectBoundingTris(boundingTris, zi_loc);
	junk = 0;
} */

void multiresBccTets::decimate(int level) {
	if (level < 2)
		throw(std::logic_error("Tet decimation below level 2 is requested.\n"));
	_decimatedNodes.clear();
	_decNodeConstraints.clear();
	int seed = 0, levelEnd = _tetNodes.size(), levelNow = 1, seedStart = 0;  // COURT nuke last if vec<bool> not used
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
	// binary tree structure for decimated nodes follows.  Only need closest parents.
	struct leafNode {
		int child;
		int parent0;
		int parent1;
		float depth;
	};
	std::list<leafNode> parents, leaves;
	auto processLeaf = [&]() {
		auto ln = leaves.begin();
		leafNode parent;
		parent.depth = ln->depth * 0.5f;
		parent.child = ln->parent0;
		if (nodeMap[parent.child] > -1)
			parents.push_back(parent);
		else {
			DNIT nit = _decimatedNodes.find(parent.child);
			parent.parent0 = nit->second.first;
			parent.parent1 = nit->second.second;
			leaves.push_back(parent);
		}
		parent.child = ln->parent1;
		if (nodeMap[parent.child] > -1)
			parents.push_back(parent);
		else {
			DNIT nit = _decimatedNodes.find(parent.child);
			parent.parent0 = nit->second.first;
			parent.parent1 = nit->second.second;
			leaves.push_back(parent);
		}
		leaves.erase(ln);
	};

	int ic2 = 0, ic3 = 0, ic4 = 0;

	for (auto dn = _decimatedNodes.begin(); dn != _decimatedNodes.end(); ++dn) {
		if (dn->first < 0)
			continue;
		int dNode = nodeMap[dn->first];
		if (dNode > -1) {  // used. A T-junction exists here
			parents.clear();
			leaves.clear();
			// collect leaves of binary tree of decimated nodes
			leafNode bud;
			bud.depth = 1.0f;
			bud.child = dn->first;
			bud.parent0 = dn->second.first;
			bud.parent1 = dn->second.second;
			leaves.push_back(bud);
			while (!leaves.empty())
				processLeaf();
			auto dnc = _decNodeConstraints.insert(std::make_pair(dNode, decimatedFaceNode()));
			auto& fn = dnc.first->second.faceNodes;
			auto& fp = dnc.first->second.faceBarys;
			if (parents.size() < 3) {
				fn.reserve(2);
				fn.push_back(nodeMap[parents.front().child]);
				fn.push_back(nodeMap[parents.back().child]);
				fp.assign(2, 0.5f);
			}
			else {
				std::map<int, float> dnVerts;
				for (auto& p : parents) {
					auto pr = dnVerts.insert(std::make_pair(p.child, p.depth));
					if (!pr.second)
						pr.first->second += p.depth;
				}

				if (dnVerts.size() == 2)
					++ic2;
				if (dnVerts.size() == 3)
					++ic3;
				if (dnVerts.size() > 3)
					++ic4;

				fn.reserve(dnVerts.size());
				fp.reserve(dnVerts.size());
				for (auto& dnv : dnVerts) {
					fn.push_back(nodeMap[dnv.first]);
					fp.push_back(dnv.second);
				}
			}
		}
	}

	std::cout << "Model has " << ic2 << " two macroNode, " << ic3 << " three macroNode, and " << ic4 << " four or greater macronode internode constraints.\n";

	_decimatedNodes.clear();
}

void multiresBccTets::getInterNodeConstraints(std::vector<int> &subNodes, std::vector<std::vector<int> > &macroNodes, std::vector<std::vector<float> > &macroBarycentrics) {
	size_t snSize = _decNodeConstraints.size();
	subNodes.clear();
	macroNodes.clear();
	macroBarycentrics.clear();
	subNodes.reserve(snSize);
	macroNodes.reserve(snSize);
	macroBarycentrics.reserve(snSize);
	for (auto& dn : _decNodeConstraints) {  //.begin(); dn != _vnTets.decNodeConstraints.end(); ++dn)
		subNodes.push_back(dn.first);
		macroNodes.push_back(std::move(dn.second.faceNodes));
		macroBarycentrics.push_back(std::move(dn.second.faceBarys));
	}
	_decNodeConstraints.clear();
}

void multiresBccTets::gridLocusToBarycentricWeightLevel(const Vec3f& gridLocus, const bccTetCentroid& tc, Vec3f& barycentricWeight)
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


void multiresBccTets::centroidPromoteOneLevel(const bccTetCentroid& tcLevelIn, bccTetCentroid& tc) {
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

void multiresBccTets::barycentricWeightToGridLocusLevels(const bccTetCentroid& tetCentroid, const Vec3f& barycentricWeight, Vec3f& gridLocus)
{
	short gridLoci[4][3];
	centroidToNodeLoci(tetCentroid, gridLoci);
	gridLocus = Vec3f((const short(&)[3])gridLoci[0]) * (1.0f - barycentricWeight.X - barycentricWeight.Y - barycentricWeight.Z);
	for (int i = 1; i < 4; ++i)
		gridLocus += Vec3f((const short(&)[3])gridLoci[i]) * barycentricWeight[i - 1];
}


