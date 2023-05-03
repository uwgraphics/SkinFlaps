#include <assert.h>
#include <limits>
#include <algorithm>
#include <cmath>
#include <functional>
#include "Mat2x2d.h"
#include "Vec3d.h"
#include "Mat2x2f.h"
#include "Mat3x3f.h"
#include "Mat3x3d.h"
#include "triTriIntersect_Shen.h"
#include "tri_tri_intersect_Guigue.h"

#include <chrono>  // for openMP timing
#include <ctime>  
#include <iostream>
#include <fstream>
#include <omp.h>

#include "vnBccTetCutter_omp.h"

std::vector<vnBccTetCutter_omp::tetTriangles> vnBccTetCutter_omp::_tetTris;
std::unordered_map<bccTetCentroid, int, vnBccTetCutter_omp::bccTetCentroidHasher> vnBccTetCutter_omp::_centroidIndices;


bool vnBccTetCutter_omp::makeFirstVnTets(materialTriangles *mt, vnBccTetrahedra *vbt, int maximumGridDimension)
{  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	if (maximumGridDimension > 0x8ffe)
		throw(std::logic_error("Maximum grid dimension requested must be less than 32K."));
	// WARNING - no complete tests are done to check for non-self-intersecting closed manifold triangulated surface input!!
	// This is essential. findAdjacentTriangles() is the closest test this routine provides.  Test externally.
	_mt = mt;
	_vbt = vbt;
	_vbt->_mt = mt;
	_vbt->_barycentricWeights.clear();
	_vbt->_barycentricWeights.assign(mt->numberOfVertices(), Vec3f());
	_vbt->_tetHash.clear();
	_vbt->_tetNodes.clear();
	evenXy.clear();
	oddXy.clear();
	if(_mt->findAdjacentTriangles(true))	return false;
	if (!setupBccIntersectionStructures(maximumGridDimension))
		return false;
	return tetCutCore();
}



bool vnBccTetCutter_omp::remakeVnTets(materialTriangles* mt)
{  // uses old vbt to get material coords of mt vertices and grid data, from which it makes new vbt.
	_mt = mt;
	_vbt->_mt = mt;


/*	_vbt->_barycentricWeights.clear();
	_vbt->_barycentricWeights.assign(mt->numberOfVertices(), Vec3f());
	_vbt->_tetHash.clear();
	_vbt->_tetNodes.clear();
	evenXy.clear();
	oddXy.clear();
	if (_mt->findAdjacentTriangles(true))	return false;
	if (!setupBccIntersectionStructures(90))
		return false; */




//	int newVertexStart = _vMatCoords.size();
	_vMatCoords.resize(_mt->numberOfVertices());
	// all oldVertices should have an unchanged grid locus
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i) {  // after incisions all new vertices should have been assigned a grid locus associated with the old _vbt
//	for (int n = _mt->numberOfVertices(), i = newVertexStart; i < n; ++i) {  // after incisions all new vertices should have been assigned a grid locus associated with the old _vbt
		if (_vbt->getVertexTetrahedron(i) < 0)  // deleted vertex
			continue;
		Vec3f Vf;
		_vbt->vertexGridLocus(i, Vf);
		_vMatCoords[i].set(Vf);
	}


/*	_vMatCoords.clear();
	_vMatCoords.assign(_mt->numberOfVertices(), Vec3f());
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i) {
		Vec3f Vf;
		Vf.set((const float(&)[3]) * _mt->vertexCoordinate(i));
		_vMatCoords[i].set((Vf - _vbt->_minCorner) * (float)_vbt->_unitSpacingInv);
	} */



	_vbt->_tetHash.clear();
	_vbt->_tetNodes.clear();
	_vbt->_nodeGridLoci.clear();
	_vbt->_tetCentroids.clear();
	_vbt->_vertexTets.clear();
	_vbt->_vertexTets.assign(mt->numberOfVertices(), -1);

	_vbt->_nodeSpatialCoords = nullptr;
	if (_mt->findAdjacentTriangles(true))
		throw(std::logic_error("topological error in vnBccTetCutter_omp::remakeVnTets()\n"));
	_vertexTetCentroids.clear();
	_vertexTetCentroids.assign(_mt->numberOfVertices(), bccTetCentroid());
//	_vbt->_barycentricWeights.clear();  // unnecessary as memory already allocated. Data will be changed here.
//	_vbt->_barycentricWeights.assign(_mt->numberOfVertices(), Vec3f());
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i) {
		_vbt->gridLocusToLowestTetCentroid(_vMatCoords[i], _vertexTetCentroids[i]);
		// set barycentric coordinate within that tet
		auto& bw = _vbt->_barycentricWeights[i];
		_vbt->gridLocusToBarycentricWeight(_vMatCoords[i], _vertexTetCentroids[i], bw);
		// any vertex on a tet face boundary should be nudged inside
		bool nudge = false;
		for (int j = 0; j < 3; ++j) {
			if (bw[j] < 1e-4f || bw[j] > 0.9999f)
				nudge = true;
		}
		if (bw[0] + bw[1] + bw[2] > 0.9999)
			nudge = true;
		if (nudge) {  // vertices on tet boundarys must be pulled inside to avoid dual identities
			Vec3f nV = (Vec3f(0.25f, 0.25f, 0.25f) - bw) * 0.0002f;
			bw += nV;
			_vbt->barycentricWeightToGridLocus(_vertexTetCentroids[i], bw, _vMatCoords[i]);
		}
	}
	// setup lines parallel with Z axis
	for (auto& exyRow : evenXy) {  // since calling this routine many times just emptying contents
		for(auto &exyCol : exyRow)
			exyCol.clear();
	}
	for (auto& oxyRow : oddXy) {  // since calling this routine many times just emptying contents
		for (auto& oxyCol : oxyRow)
			oxyCol.clear();
	}
	return tetCutCore();
}

bool vnBccTetCutter_omp::tetCutCore() {  // same cut core for first cut and recuts
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	int it, nOmpRange = _mt->numberOfTriangles();
	// private variables
	std::unordered_map<bccTetCentroid, std::vector<int>, bccTetCentroidHasher> centTris_local;
	std::vector<std::pair<bccTetCentroid, std::vector<int> > > cT_global;
	std::vector<zIntrsct> zInt_local;
	// COURT note the nowait.  Also #pragma statements cannot have a curly brace immediately after them, but must be on a new line.

//	_centroidIndices, , oddXy, evenXy, , _tetTris

// #pragma omp parallel shared(nOmpRange, cT_global) private(it, centTris_local, zInt_local)
	{
		centTris_local.clear();
		zInt_local.clear();
// #pragma omp for schedule(static, 1)
		for (it = 0; it < nOmpRange; ++it)
			inputTriangle(it, centTris_local, zInt_local);
// #pragma omp critical (centTris_local)
		for (auto ctl : centTris_local) {

//			cT_global.push_back(std::move(ctl));  // std::make_pair(ctl.first, std::move(ctl.second)));

			auto pr = _centroidIndices.insert(std::make_pair(ctl.first, -1));
			if (pr.second) {
				tetTriangles tt;
				tt.tc = ctl.first;
				tt.tetindx = -1;
				tt.tris = std::move(ctl.second);
				pr.first->second = _tetTris.size();
				_tetTris.push_back(tt);
			}
			else {
				auto& ttr = _tetTris[pr.first->second].tris;
				ttr.insert(ttr.end(), ctl.second.begin(), ctl.second.end());
			}
		}
		std::cout << "centTris_local.size = " << centTris_local.size() << " for thread id: " << omp_get_thread_num() << "\n";
// #pragma omp critical (zInt_local)
		for (auto& zi : zInt_local) {
			if (zi.odd)
				oddXy[zi.x][zi.y].insert(std::make_pair(zi.zInt, zi.solidBegin));
			else
				evenXy[zi.x][zi.y].insert(std::make_pair(zi.zInt, zi.solidBegin));
		}
	}

/*	for (auto& ctl : cT_global) {
		auto pr = _centroidIndices.insert(std::make_pair(ctl.first, -1));
		if (pr.second) {
			tetTriangles tt;
			tt.tc = ctl.first;
			tt.tetindx = -1;
			tt.tris = std::move(ctl.second);
			pr.first->second = _tetTris.size();
			_tetTris.push_back(tt);
		}
		else {
			auto& ttr = _tetTris[pr.first->second].tris;
			ttr.insert(ttr.end(), ctl.second.begin(), ctl.second.end());
		}
	} */

	std::cout << "After inputTriangles() _centroidIndices size is " << _centroidIndices.size() << " and _tetTris.size is " << _tetTris.size() << "\n";

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	std::string message("Inputting triangles took ");
	message += std::to_string(elapsed_seconds.count());
	message += " seconds for ";
	message += std::to_string(_vbt->_tetNodes.size());
	message += " tets.";
	std::cout << message << "\n";
	//	_surgAct->sendUserMessage(message.c_str(), "Timer");


	start = std::chrono::system_clock::now();

	// create and hash all interior nodes
	createInteriorNodes();
//	evenXy.clear();  // am reusing these structures for remakeVnTets()
//	oddXy.clear();

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	end_time = std::chrono::system_clock::to_time_t(end);
	message.assign("Interior node creation took ");
	message += std::to_string(elapsed_seconds.count());
	message += " seconds for ";
	message += std::to_string(_vbt->_tetNodes.size());
	message += " tets.";
	std::cout << message << "\n";


	_interiorNodes.clear();  // only nodes created thus far are interior nodes
	_interiorNodes.reserve(_vbt->_nodeGridLoci.size());
	for (int n = _vbt->_nodeGridLoci.size(), j = 0; j < n; ++j)
		_interiorNodes.insert(std::make_pair(_vbt->_nodeGridLoci[j], j));
	//	_exteriorNodeIndices.clear();
	//	_exteriorNodeIndices.reserve(_tetTris.size() >> 1);    // COURT later check this reserve size
	_vbt->_tetCentroids.clear();
	_vbt->_tetCentroids.reserve(_tetTris.size() >> 1);    // COURT perhaps assign() for multi threading

	start = std::chrono::system_clock::now();

	nOmpRange = _tetTris.size();
	_nSurfaceTets.store(0);
	nts_global.clear();
	nts_vec.clear();
	nts_locs.clear();
	nts_global.reserve(nOmpRange << 2);  // COURT check these reserve sizes for good guess
	nts_vec.reserve(nOmpRange << 2);
	nts_locs.reserve(nOmpRange << 2);
	_centroidTriangles.clear();
	_centroidTriangles.reserve(nOmpRange << 2);

#pragma omp parallel
	{
		// private variables
		std::vector<newTet> newTets_local;
		std::unordered_map<std::array<short, 3>, std::list<nodeTetSegment>, arrayShort3Hasher> nts_local;

#pragma omp for schedule(dynamic) nowait
		// COURT note the nowait.  Also #pragma statements cannot have a curly brace immediately after them, but must be on a new line.
		for (int i = 0; i < nOmpRange; ++i) {
			getConnectedComponents(_tetTris[i], newTets_local, nts_local);  // for this centroid split its triangles into solid connected components
		}
#pragma omp critical
		{
			int maxIdx = newTets_local.back().tetIdx;
			if (maxIdx >= _vbt->_tetCentroids.size()) {
				int incr = maxIdx - _vbt->_tetCentroids.size() + 1;
				_vbt->_tetCentroids.insert(_vbt->_tetCentroids.end(), incr, bccTetCentroid());
				_vbt->_tetNodes.insert(_vbt->_tetNodes.end(), incr, std::array<int, 4>());
			}
			for (auto& nt : newTets_local) {
				_vbt->_tetCentroids[nt.tetIdx] = nt.tc;  // COURT don't hash them here
				_vbt->_tetNodes[nt.tetIdx] = std::move(nt.tetNodes);
				auto pr = _centroidTriangles.insert(std::make_pair(nt.tc, std::list<tetTris>()));
				pr.first->second.push_back(tetTris());
				pr.first->second.back().tetIdx = nt.tetIdx;
				pr.first->second.back().tris = std::move(nt.tris);
			}
			for (auto& nts : nts_local) {
				auto pr = nts_global.insert(std::make_pair(nts.first, -1));
				if (pr.second) {
					pr.first->second = nts_vec.size();
					nts_vec.push_back(std::move(nts.second));
					nts_locs.push_back(pr.first->first);
				}
				else {
					auto& vr = nts_vec[pr.first->second];
					for (auto& nr : nts.second)
						vr.push_back(std::move(nr));
				}
			}
		}
	}

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	end_time = std::chrono::system_clock::to_time_t(end);
	message.assign("Getting connected components took ");
	message += std::to_string(elapsed_seconds.count());
	message += " seconds for ";
	message += std::to_string(_vbt->_tetNodes.size());
	message += " tets.";
	std::cout << message << "\n";

	// get tets where vertices reside
	_vbt->_vertexTets.clear();
	_vbt->_vertexTets.assign(_mt->numberOfVertices(), -1);
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i) {
		auto ci = _centroidTriangles.find(_vertexTetCentroids[i]);
		if (ci == _centroidTriangles.end())
			throw(std::logic_error("Couldn't find a tetrahedron containing a vertex.\n"));
		else if (ci->second.size() < 2)
			_vbt->_vertexTets[i] = ci->second.front().tetIdx;
		else {
			auto ltit = ci->second.begin();
			while (ltit != ci->second.end()) {
				auto tit = ltit->tris.begin();
				while (tit != ltit->tris.end()) {
					int* tr = _mt->triangleVertices(*tit);
					int j;
					for (j = 0; j < 3; ++j) {
						if (tr[j] == i) {
							_vbt->_vertexTets[i] = ltit->tetIdx;
							break;
						}
					}
					if (j < 3)
						break;
					++tit;
				}
				if (tit != ltit->tris.end())
					break;
				++ltit;
			}
			assert(_vbt->_vertexTets[i] > -1);
		}
	}

	start = std::chrono::system_clock::now();

	nOmpRange = nts_vec.size();
#pragma omp parallel
	{
		// private variables
		std::vector<extNode> eNodes_local;
#pragma omp for schedule(dynamic) nowait
		// COURT note the nowait.  Also #pragma statements cannot have a curly brace immediately after them, but must be on a new line.
		for (int i = 0; i < nOmpRange; ++i) {
			assignExteriorTetNodes(i, eNodes_local);
		}
#pragma omp critical
		{
			for (auto& en : eNodes_local) {
				int eNode = _vbt->_nodeGridLoci.size();
				_vbt->_nodeGridLoci.push_back(std::move(en.loc));
				for (auto& ti : en.tiPairs)
					_vbt->_tetNodes[ti.first][ti.second] = eNode;
			}
		}
	}

	_vbt->_firstInteriorTet = _vbt->_tetNodes.size();

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	end_time = std::chrono::system_clock::to_time_t(end);
	message.assign("Assigning exterior nodes took ");
	message += std::to_string(elapsed_seconds.count());
	message += " seconds for ";
	message += std::to_string(_vbt->_tetNodes.size());
	message += " tets.";
	std::cout << message << "\n";


	start = std::chrono::system_clock::now();

	fillNonVnTetCenter();

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	end_time = std::chrono::system_clock::to_time_t(end);
	message.assign("Filling tet center took ");
	message += std::to_string(elapsed_seconds.count());
	message += " seconds for ";
	message += std::to_string(_vbt->_tetNodes.size());
	message += " tets.";
	std::cout << message << "\n";

	_vbt->_tetNodes.shrink_to_fit();
	_vbt->_tetCentroids.shrink_to_fit();
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_vbt->_tetCentroids.size());
	for (int n = _vbt->_tetCentroids.size(), i = 0; i < n; ++i)
		_vbt->_tetHash.insert(std::make_pair(_vbt->_tetCentroids[i], i));
	return true;
}

void vnBccTetCutter_omp::createInteriorNodes() {
	short z;
	std::array<short, 3> s3;
	auto solidFilter = [&](std::multimap<double, bool>& mm) {  // isolate solid runs from surface edge and vertex hits and correct switch of coincident surfaces.
		auto mit = mm.begin();
		bool even = true;
		while (mit != mm.end()) {
			auto mNext = mit;
			++mNext;
			if (mNext == mm.end())
				break;
			while (mNext != mm.end() && mNext->first - mit->first < 1e-4f) {
				if (mit->second == mNext->second) {  // surface edge or vertex hit
					mm.erase(mNext);
					mNext = mit;
					++mNext;
				}
				else if (mit->second != even) {  // coincident surfaces
					mit->second = even;
					mNext->second = !even;
					break;
				}
				else
					break;
			}
			even = !even;
			++mit;
		}

//#ifdef _DEBUG
//		if (mm.size() & 1)
//			throw(std::logic_error("Program error in cutter in createInteriorNodes()\n");
//		bool inside = true;
//		for (auto mit = mm.begin(); mit != mm.end(); ++mit) {
//			if (mit->second != inside)
//				throw(std::logic_error("Program error in cutter in createInteriorNodes()\n");
//			inside = !inside;
//		}
// #endif

	};
	auto runInteriorNodes = [&](std::multimap<double, bool>& mm, bool evenLine) {
		auto mit = mm.begin();
		auto mend = mm.end();
		if (mit == mend)
			return;
		if (!mit->second)
			throw(std::runtime_error("Model is not a closed manifold surface.\n"));
		while (mit != mend) {
			// create runs of interior vertices
			z = (short)std::floor(mit->first);
			if (z < mit->first)
				++z;
			if ((z & 1) == evenLine)
				++z;
			++mit;
			if (mit->second) {  // correct any coincident pairs or micro collisions
				auto mit2 = mit;
				while (mit2 != mend && mit2->second && mit2->first - mit->first < 2e-5f)
					++mit2;
				if (mit2 == mend || mit2->second)
					throw(std::runtime_error("Model has a self collision.\n"));
				mit = mit2;
			}
			for (; z < mit->first; z += 2) {
				s3[2] = z;
				_vbt->_nodeGridLoci.push_back(s3);
			}
			++mit;
			if (mit != mend && !mit->second) {  // correct any coincident pairs or micro collisions
				auto mit2 = mit;
				while (mit2 != mend && !mit2->second && mit2->first - mit->first < 2e-5f)
					++mit2;
				if (mit2 == mend)
					return;
				if (!mit2->second)
					throw(std::runtime_error("Model has a self collision.\n"));
				mit = mit2;
			}
		}
	};
	for (short xi = 0; xi < oddXy.size(); ++xi) {
		for (short yi = 0; yi < oddXy[xi].size(); ++yi) {
			if (oddXy[xi][yi].empty())
				continue;
			s3[0] = xi * 2 + 1;
			s3[1] = yi * 2 + 1;
			solidFilter(oddXy[xi][yi]);
			runInteriorNodes(oddXy[xi][yi], false);
		}
	}
	for (short xi = 0; xi < evenXy.size(); ++xi) {
		for (short yi = 0; yi < evenXy[xi].size(); ++yi) {
			if (evenXy[xi][yi].empty())
				continue;
			s3[0] = (xi + 1) * 2;
			s3[1] = (yi + 1) * 2;
			solidFilter(evenXy[xi][yi]);
			runInteriorNodes(evenXy[xi][yi], true);
		}
	}
}

void vnBccTetCutter_omp::assignExteriorTetNodes(int exteriorNodeNum, std::vector<extNode>& eNodes) {
	auto &ntsList = nts_vec[exteriorNodeNum];
	auto &loc = nts_locs[exteriorNodeNum];
	struct tetIndex {
		int tet;
		int nodeIndex;
	};
	struct tetTris {
		std::set<int> tris;
		std::list<tetIndex> tetIndices;
		int node;
	};
	std::list<tetTris> triPools;
	auto tetsConnect = [&](const std::vector<int>& tris0, const std::set<int> &poolTris) ->bool{
		for (auto t : tris0) {
			if (poolTris.find(t) != poolTris.end())
				return true;
		}
		return false;
	};
	triPools.clear();
	for (auto nit = ntsList.begin(); nit != ntsList.end(); ++nit) {
		bool makeNewPool = true;
		auto firstConnectedPool = triPools.end();
		auto pit = triPools.begin();
		while ( pit != triPools.end()) {
			if (tetsConnect(nit->tetNodeTris, pit->tris)) {
				if (firstConnectedPool == triPools.end()) {
					firstConnectedPool = pit;
					tetIndex ti;
					ti.tet = nit->tetIdx;
					ti.nodeIndex = nit->tetNodeIndex;
					firstConnectedPool->tetIndices.push_back(ti);
					firstConnectedPool->tris.insert(nit->tetNodeTris.begin(), nit->tetNodeTris.end());
					++pit;
				}
				else {
					// merge pools
					firstConnectedPool->tris.insert(pit->tris.begin(), pit->tris.end());
					firstConnectedPool->tetIndices.splice(firstConnectedPool->tetIndices.end(), pit->tetIndices);
					pit = triPools.erase(pit);
				}
				makeNewPool = false;
			}
			else
				++pit;
		}
		if(makeNewPool){
			tetTris tt;
			tt.tris.insert(nit->tetNodeTris.begin(), nit->tetNodeTris.end());
			tetIndex ti;
			ti.tet = nit->tetIdx;
			ti.nodeIndex = nit->tetNodeIndex;
			tt.tetIndices.push_back(ti);
			triPools.push_back(tt);
		}
	}
	// any remaining triPools should get their own exterior, possibly virtual, node
	for (auto tp = triPools.begin(); tp != triPools.end(); ) {  // not multithreaded
		extNode en;
		en.loc = std::move(loc);
		for (auto& ti : tp->tetIndices)
			en.tiPairs.push_back(std::make_pair(ti.tet, ti.nodeIndex));
		eNodes.push_back(en);
		tp = triPools.erase(tp);
	}
}

int vnBccTetCutter_omp::nearestRayPatchHit(const Vec3f &rayBegin, Vec3f rayEnd, const std::vector<int>& tris, float& distanceSq) {  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
	distanceSq = FLT_MAX;
	Vec3d rE(rayBegin - rayEnd);
//	rayEnd -= rayBegin;
	int closeT = -1;
	Vec3d N;
	for (auto& t : tris) {
		int* tr = _mt->triangleVertices(t);
		Vec3d tv[3];
		for (int i = 0; i < 3; ++i)
			tv[i] = _vMatCoords[tr[i]];
		tv[1] -= tv[0];
		tv[2] -= tv[0];
		double dSq;
		Vec3d P, R;
		Mat3x3d M;
		M.Initialize_With_Column_Vectors(tv[1], tv[2], rE);  // Vec3d(- rayEnd)
		R = M.Robust_Solve_Linear_System(Vec3d(rayBegin) - tv[0]);
		if (R[0] < 0.0f || R[0] > 1.0f || R[1] < 0.0f || R[1] > 1.0f || R[2] < 0.0f || R[2] > 1.0f || R[0] + R[1] > 1.0f)
			continue;
//		P = rayEnd * R[2];
		P = rE * R[2];
		dSq = P * P;
		if (distanceSq > dSq) {
			distanceSq = dSq;
			closeT = t;
			N = tv[1] ^ tv[2];
		}
	}
	if (distanceSq < FLT_MAX) {
		if(distanceSq < 1e-3)  // coincident surfaces don't connect
			return 1;
		if (N * rE > 0.0f)
			return 1;
		else
			return -1;
	}
	else
		return 0;
}

bool vnBccTetCutter_omp::nearestPatchPoint(const short (&gl)[4][3], const int tetIdx, const std::vector<int>& tris, Vec3f& closeP, float& distanceSq) {
	// This routine assumes a patch intersection with an adjacent tet face.  It also hands off cases of opposing face intersections and closed patches entirely with a large tet.
	Vec3d closestP, P = {(double)gl[tetIdx][0], (double)gl[tetIdx][1], (double)gl[tetIdx][2]};
	double d, minD;
	int triEdge[4] = { -1, -1, -1, -1 }, teNext;
	for (int i = 0; i < 3; ++i) {  // search the 3 faces surrounding tetIdx
		minD = DBL_MAX;
		Vec3d tri[3];
		for (int j = 0; j < 3; ++j) {
			int nIdx = (tetIdx + 2 + j + i) & 3;
			tri[j].set((double)gl[nIdx][0], (double)gl[nIdx][1], (double)gl[nIdx][2]);  // each face of tet with last one opposite this tetIdx
		}
		teNext = 0;
		for (auto& t : tris) {
			Vec3d tt[3], source, target;
			int* tr = _mt->triangleVertices(t);
			for (int k = 0; k < 3; ++k)
				tt[k].set(_vMatCoords[tr[k]]);
			int coplanar = 0;
			if (tri_tri_intersection_test_3d(tri[0].xyz, tri[1].xyz, tri[2].xyz, tt[0].xyz, tt[1].xyz, tt[2].xyz, &coplanar, source.xyz, target.xyz)) {
				if (coplanar == 1 || source == target)
					continue;
				Vec3d N = target - source;
				double s = (P - source) * N / (N * N);
				if (s > 1.0)
					;
				else if (s < 0.0)
					target = source;
				else {
					target *= s;
					target += source * (1.0 - s);
				}
				d = (P - target).length2();
				if (d - 1e-5 < minD) {  // looking for a surface edge so admit double hit
					if (abs(d - minD) < 1e-5) {
						if (teNext > 3) {  // exceedingly rare vertex hit on a face
							throw(std::logic_error("Vertex intersection in nearestPatchPoint() not yet programmed.\n"));
						}
						triEdge[teNext++] = t;
					}
					else {
						triEdge[0] = t;
						teNext = 1;
						minD = d;
						closestP = target;
					}
				}
			}
		}
		if (triEdge[0] > -1)
			break;
	}
	if(triEdge[0] < 0)  // this patch is a closed lacuna completely inside this tet, or intersects with tet face opposite tetIdx
		return closestPatchPoint(Vec3f((float)gl[tetIdx][0], (float)gl[tetIdx][1], (float)gl[tetIdx][2]), tris, closeP, distanceSq);
	closeP.set((float)closestP[0], (float)closestP[1], (float)closestP[2]);
	distanceSq = (float)minD;
	if (teNext > 2)  // 2 coincident edges so always inside
		return true;
	auto triNorm = [&](const int triangle, Vec3d &N) {
		int* tr = _mt->triangleVertices(triangle);
		Vec3d U = _vMatCoords[tr[1]] - _vMatCoords[tr[0]], V = _vMatCoords[tr[2]] - _vMatCoords[tr[0]];
		N = U ^ V;
		N.q_normalize();
	};
	Vec3d N0;
	triNorm(triEdge[0], N0);  // intersection of edge of tet face with a single triangle if no triEdge[1]
	if (triEdge[1] > -1){  // face intersection with edge between two triangles
		Vec3d N1;
		triNorm(triEdge[1], N1);
		N0 += N1;
		if (abs(N0.X) < 0.02 && abs(N0.Y) < 0.02 && abs(N0.Z) < 0.02)  // at bottom of coincident edge.  Guaranteed inside. Generous due to Quake normalization
			return true;
	}
	return (P - closeP) * N0 < 0.0;
}

int vnBccTetCutter_omp::nearestPatchEdgePoint(const short(&gl)[4][3], const int tetIdx, const std::vector<int>& tris, Vec3f& closeP, float& distanceSq) {
	// This routine looks for a patch intersection with a tet edge surrounding tetIdx. Return 1 if tetIdx inside patch, -1 outside, and 0 is none of the 3 edge intersects.
	Vec3d closestP, P = { (double)gl[tetIdx][0], (double)gl[tetIdx][1], (double)gl[tetIdx][2] };
	double minD;
	bool inside, both = false;
	for (int i = 1; i < 4; ++i) {
		minD = DBL_MAX;
		Vec3d Q;
		int nIdx = (tetIdx + i) & 3;
		Q.set((double)gl[nIdx][0], (double)gl[nIdx][1], (double)gl[nIdx][2]);  // each face of tet with last one opposite this tetIdx
		Q -= P;
		for (auto& t : tris) {
			Vec3d tt[3], R;
			int* tr = _mt->triangleVertices(t);
			for (int k = 0; k < 3; ++k)
				tt[k].set(_vMatCoords[tr[k]]);
			tt[1] -= tt[0];
			tt[2] -= tt[0];
			Mat3x3d M(tt[1], tt[2], -Q);
			R = M.Robust_Solve_Linear_System(P - tt[0]);
			if (R[2] < 0.0 || R[2] > 1.0 || R[0] < 0.0 || R[0] > 1.0 || R[1] < 0.0 || R[1] > 1.0 || R[0] + R[1] > 1.0)
				continue;
			if (R[2] - 3e-3 < minD) {
				if (minD - R[2] < 3e-3) {
					if(inside != ((tt[1] ^ tt[2]) * Q >= 0))
						both = true;
				}
				else {
					minD = R[2];
					distanceSq = (Q * R[2]).length2();
					closestP = P + Q * R[2];
					inside = ((tt[1] ^ tt[2]) * Q >= 0);
					both = false;
				}
			}
		}
		if (minD < DBL_MAX)
			break;
	}
	closeP.set(closestP.xyz);
	if (minD > 1e38)
		return 0;
	if (both)  // coincident surface so point always inside
		return 1;
	return inside ? 1 : -1;
}

bool vnBccTetCutter_omp::closestPatchPoint(const Vec3f& P, const std::vector<int>& tris, Vec3f &closeP, float &distanceSq) {
	// finds and returns closest point on patch triangles to input point P. Return if P is inside patch.
	float dSq;
	distanceSq = FLT_MAX;
	bool repeatFind, coincidentSurface = false;
	Vec3f closeN;
	int closestT = -1, closeTI, closestTriIndex = -1;
	struct edgePoint {
		int tri0;
		int tri1;
	}closeE, closestEdge;
	closestEdge.tri0 = -1;
	closestEdge.tri1 = -1;
	auto setCloseE = [&](int tri, int edge) {
		int aT[3], aE[3];
		_mt->triangleAdjacencies(tri, aT, aE);
		closeE.tri0 = tri;
		closeE.tri1 = aT[edge];
		closeTI = -1;
	};
	for (auto &t : tris) {
		closeTI = -1;
		closeE.tri0 = -1;
		int *tr = _mt->triangleVertices(t);
		Vec3f T[3];
		for (int i = 0; i < 3; ++i)
			T[i] = _vMatCoords[tr[i]];
		T[1] -= T[0];
		T[2] -= T[0];
		Mat2x2f M = {T[1]*T[1], T[1]*T[2], 0.0f, T[2] * T[2] };
		M.x[2] = M.x[1];
		Vec3f Q = P - T[0];
		Vec2f S(Q * T[1], Q * T[2]);
		Vec2f R = M.Robust_Solve_Linear_System(S);
		// as tri barycentrics not orthogonal, must find closest outside edge
		if (R[0] < 0.0f) {
			if( R[1] < 0.0f) {
				R[0] = 0.0f;
				R[1] = 0.0f;
				closeTI = 0;
			}
			else {
				Vec3f X = T[0] + T[1] * R[0] + T[2] * R[1];
				R[0] = 0.0f;
				R[1] = ((X - T[0]) * T[2]) / M.x[3];
				if (R[1] < 0.0f) {
					R[1] = 0.0f;
					closeTI = 0;
				}
				else if (R[1] > 1.0f) {
					R[1] = 1.0f;
					closeTI = 2;
				}
				else
					setCloseE(t, 2);
			}
		}
		else if (R[1] < 0.0f) {
			Vec3f X = T[0] + T[1] * R[0] + T[2] * R[1];
			R[1] = 0.0f;
			R[0] = ((X - T[0]) * T[1]) / M.x[0];
			if (R[0] < 0.0f) {
				R[0] = 0.0f;
				closeTI = 0;
			}
			else if (R[0] > 1.0f) {
				R[0] = 1.0f;
				closeTI = 1;
			}
			else
				setCloseE(t, 0);
		}
		else if (R[0] + R[1] > 1.0f) {
			Vec3f K, L, X = T[0] + T[1] * R[0] + T[2] * R[1];
			L = T[2] - T[1];
			K = T[1] + T[0];
			R[1] = ((X - K) * L) / (L * L);
			if (R[1] < 0.0f) {
				R[1] = 0.0f;
				closeTI = 1;
			}
			else if (R[1] > 1.0f) {
				R[1] = 1.0f;
				closeTI = 2;
			}
			else
				setCloseE(t, 1);
			R[0] = 1.0f - R[1];
		}
		else
			;
		Q -= T[1] * R[0] + T[2] * R[1];
		dSq = Q * Q;
		if (dSq <= distanceSq){
			repeatFind = false;
			if (dSq == distanceSq)
				repeatFind = true;
			if (closeTI > -1) {
				if (repeatFind && closestTriIndex > -1) {
					if (_mt->triangleVertices(closestT)[closestTriIndex] != _mt->triangleVertices(t)[closeTI])
						coincidentSurface = true;
				}
				else {
					closestTriIndex = closeTI;
					closestT = t;
				}
				closestEdge.tri0 = -1;
			}
			else if (closeE.tri0 > -1) {
				if (repeatFind && closestEdge.tri0 > -1) {
					if (closeE.tri0 != closestEdge.tri0 && closeE.tri0 != closestEdge.tri1)  // opposite triangle not found
						coincidentSurface = true;
				}
				else {
					closestEdge.tri0 = closeE.tri0;
					closestEdge.tri1 = closeE.tri1;
				}
				closestTriIndex = -1;
			}
			else {
				if (repeatFind && closestT > -1)
					coincidentSurface = true;
				closestTriIndex = -1;
				closestEdge.tri0 = -1;
				closeN = T[1] ^ T[2];
				closestT = t;
			}
			distanceSq = dSq;
			closeP = P - Q;
		}
	}
	auto matTriNormal = [&](const int triangle, Vec3f &N) {
		int* tr = _mt->triangleVertices(triangle);
		Vec3f W0, W1;
		W0 = _vMatCoords[tr[2]] - _vMatCoords[tr[0]];
		W1 = _vMatCoords[tr[1]] - _vMatCoords[tr[0]];
		N = W0 ^ W1;
		N.q_normalize();
	};
	if (coincidentSurface)  // this connected patch always inside
		return true;
	else if (closestEdge.tri0 > -1) {
		Vec3f N[2];
		matTriNormal(closestEdge.tri0, N[0]);
		matTriNormal(closestEdge.tri1, N[1]);
		N[0] += N[1];
		if (fabs(N[0].X) < 1e-5f && fabs(N[0].Y) < 1e-5f && fabs(N[0].Z) < 1e-5f)  // coincident surface at a common edge
			return true;
		else {
			closeN = N[0];
			assert(!(closeN[0] == 0.0f && closeN[1] == 0.0f && closeN[2] == 0.0f));
		}
	}
	else if (closestTriIndex > -1) {
		return pointInsidePatchVertex(P, closestT, closestTriIndex);
	}
	else
		;
	return (closeN * (P - closeP)) < 0.0f;
}

bool vnBccTetCutter_omp::pointInsidePatchVertex(const Vec3f& P, const int triangle, const int tIndex) {
	// get closest edge neighbor to P
	std::vector<int> nt, nv;
	_mt->triangleVertexNeighbors(triangle, tIndex, nt, nv);
	assert(nt[0] > -1);
	Vec3f Q, R, S;
	Q = _vMatCoords[_mt->triangleVertices(triangle)[tIndex]];
	S = P - Q;
	S.q_normalize();
	float maxDot = -1e32f, dot;
	int edge = -1;
	for (int n = nv.size(), i = 0; i < n; ++i) {
		R = _vMatCoords[nv[i]];
		R -= Q;
		R.q_normalize();
		dot = S * R;
		if (maxDot < dot) {
			maxDot = dot;
			edge = i;
		}
	}
	auto tNorm = [&](const int triangle, Vec3f& N) {
		int *tr = _mt->triangleVertices(triangle);
		Vec3f N0 = _vMatCoords[tr[1]] - _vMatCoords[tr[0]], N1 = _vMatCoords[tr[2]] - _vMatCoords[tr[0]];
		N = N0 ^ N1;
		N.q_normalize();
	};
	tNorm(nt[edge], R);
	tNorm(nt[(edge + 1) % nt.size()], Q);
	R += Q;
	if (fabs(R.X) < 1e-8f && fabs(R.Y) < 1e-8f && fabs(R.Z) < 1e-8f)  // coincident surface at a common edge
		return true;
	return R * S < 0.0f;
}

void vnBccTetCutter_omp::getConnectedComponents(tetTriangles& tt, std::vector<newTet>& nt_vec, std::unordered_map<std::array<short, 3>, std::list<nodeTetSegment>, arrayShort3Hasher>& local_nts) {
	// for this centroid split its triangles into single solid connected components
	std::set<int> trSet, ts;
	trSet.insert(tt.tris.begin(), tt.tris.end());
	tt.tris.clear();
	struct patch {
		std::vector<int> lv;
		int tetIndex;
		std::array<int, 4> tetNodes;
		Vec3f insideP;  // point just inside this patch ?nuke
		float dsq;
	}p;
	std::list<patch> patches;
	p.tetNodes = { -1, -1, -1, -1 };
	p.tetIndex = -1;
	patches.push_back(p);
	auto trVec = &patches.back().lv;  // on entry only one triangle list for this centroid
	while (!trSet.empty()) { // find list of connected triangles
		ts.clear();
		int t, at[3], ae[3];
		t = *trSet.begin();
		trSet.erase(trSet.begin());
		ts.insert(t);
		std::forward_list<int> tList;  // unprocessed adjacencies
		tList.push_front(t);
		while (!tList.empty()) {
			t = tList.front();
			tList.pop_front();
			_mt->triangleAdjacencies(t, at, ae);  // COURT consider a local implementation of this
			for (int i = 0; i < 3; ++i) {
				auto ti = trSet.find(at[i]);
				if (ti == trSet.end())
					continue;
				trSet.erase(ti);
				ts.insert(at[i]);
				tList.push_front(at[i]);
			}
		}
		trVec->assign(ts.begin(), ts.end());
		if (!trSet.empty()) {
			patches.push_back(p);
			trVec = &patches.back().lv;
		}
	}

	if (patches.size() < 2) {  // only one tet. no inter patch contention to resolve

	}

	if (tt.tc[0] == 52 && tt.tc[1] == 26 && tt.tc[2] == 19)
		int junk = 0;

	// find interior nodes for each component, joining patches as appropriate
	short gl[4][3];
	_vbt->centroidToNodeLoci(tt.tc, gl);
	for (int i = 0; i < 4; ++i) {
		std::array<short, 3> loc = { gl[i][0], gl[i][1], gl[i][2] };
		auto iit = _interiorNodes.find(loc);


		if (patches.size() < 2) {  // only one tet.  No inter patch contention to resolve
			if (iit == _interiorNodes.end())
				patches.front().tetNodes[i] = -1;
			else
				patches.front().tetNodes[i] = iit->second;
		}
		else {  // put patch dissection block inside when tools are solid


			Vec3f node(loc[0], loc[1], loc[2]);  // COURT switch routine to Vec3d
			std::list<std::list<patch>::iterator> outside, inside;

			if (i == 0 && tt.tc[0] == 64 && tt.tc[1] == 9 && tt.tc[2] == 14)
				int junk = 0;

			for (auto pit = patches.begin(); pit != patches.end(); ++pit) {
				bool connected;
				int npep = nearestPatchEdgePoint(gl, i, pit->lv, pit->insideP, pit->dsq);  // fast and, even with coincident surfaces, never gets inside assignment wrong
				if (npep > 0)
					connected = true;
				else if (npep < 0)
					connected = false;
				else
					connected = nearestPatchPoint(gl, i, pit->lv, pit->insideP, pit->dsq);
				if (connected)
					inside.push_back(pit);
				else
					outside.push_back(pit);
			}
			if (inside.empty()) {  // No patch bound to an interior node. ? nuke
				if (iit != _interiorNodes.end())
					int junk = 0;
				//				throw(std::logic_error("An existing tet with an interior node has no patches attached to it in getConnectedComponents()\n"));
			}
			else if (!outside.empty()) {
				for (auto ip = inside.begin(); ip != inside.end(); ) {
					auto oMerge = outside.end();
					for (auto op = outside.begin(); op != outside.end(); ++op) {
						int connected = nearestRayPatchHit((*ip)->insideP, node, (*op)->lv, (*op)->dsq);  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
						if (connected < 0 && (*op)->dsq > 1e-5f) {  // rule out possible coincident surfaces
							if (oMerge == outside.end() || (*op)->dsq < (*oMerge)->dsq) {
								oMerge = op;
							}
						}
					}
					if (oMerge != outside.end()) {
						(*oMerge)->lv.insert((*oMerge)->lv.end(), (*ip)->lv.begin(), (*ip)->lv.end());
						patches.erase(*ip);
						ip = inside.erase(ip);
					}
					else
						++ip;
				}
			}
			// any remaining inside patches should be merged along with any interior nodes
			if (inside.size() > 1) {
				if (iit == _interiorNodes.end())
					int junk = 0;;
				auto ip0 = inside.begin();
				auto ip = ip0;
				++ip;
				while (ip != inside.end()) {
					(*ip0)->lv.insert((*ip0)->lv.end(), (*ip)->lv.begin(), (*ip)->lv.end());
					patches.erase(*ip);
					ip = inside.erase(ip);
				}
			}
			if (iit != _interiorNodes.end()) {
				if (inside.empty())
					int junk = 0;
				//				throw(std::logic_error("gCC() not working in cutter.\n"));
				else
					inside.front()->tetNodes[i] = iit->second;
			}
		}
	}

	for (auto& cTet : patches) {
		nt_vec.push_back(newTet());
		nt_vec.back().tetIdx = _nSurfaceTets.fetch_add(1);
		cTet.tetIndex = nt_vec.back().tetIdx;
		nt_vec.back().tc = tt.tc;
		nt_vec.back().tetNodes = cTet.tetNodes;
		nt_vec.back().tris.assign(cTet.lv.begin(), cTet.lv.end());
	}
	for (int i = 0; i < 4; ++i) {
		bool enNotEntered = true;
		std::pair< std::unordered_map<std::array<short, 3>, std::list<nodeTetSegment>, arrayShort3Hasher>::iterator, bool> pr;
		for (auto& p : patches) {
			if (p.tetNodes[i] < 0) {
				if (enNotEntered) {
					std::array<short, 3> loc = { gl[i][0], gl[i][1], gl[i][2] };
					pr = local_nts.insert(std::make_pair(loc, std::list<nodeTetSegment>()));
					enNotEntered = false;
				}
				nodeTetSegment nts;
				nts.tetNodeIndex = i;
				nts.tetIdx = p.tetIndex;
				nts.tetNodeTris.assign(p.lv.begin(), p.lv.end());
				pr.first->second.push_back(nts);
			}
		}
	}
}

bool vnBccTetCutter_omp::setupBccIntersectionStructures(int maximumGridDimension)
{
	boundingBox<float> bbf;
	bbf.Empty_Box();
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i)
		bbf.Enlarge_To_Include_Point((const float(&)[3])(*_mt->vertexCoordinate(i)));
	bbf.Minimum_Corner(_vbt->_minCorner.xyz);
	bbf.Maximum_Corner(_vbt->_maxCorner.xyz);
	// In this model all grid distances are 1 and all odd or even Cartesian distances are 2.
	_vbt->_unitSpacing = -1.0f;
	int bigDim;
	for (int i = 0; i < 3; ++i){
		float dimSize = bbf.val[(i << 1) + 1] - bbf.val[i << 1];
		if (dimSize > _vbt->_unitSpacing){
			_vbt->_unitSpacing = dimSize;
			bigDim = i;
		}
	}
	// object must be inside positive octant
	float offset = (_vbt->getMaximumCorner() - _vbt->getMinimumCorner()).length() * 0.0001f;
	_vbt->_minCorner -= Vec3f(offset, offset, offset);
	_vbt->_maxCorner += Vec3f(offset, offset, offset);
	double cs = _vbt->_maxCorner.xyz[bigDim] - _vbt->_minCorner.xyz[bigDim];
	_vbt->_unitSpacing = cs / maximumGridDimension;
	_vbt->_unitSpacingInv = 1.0 / _vbt->_unitSpacing;
	Vec3f maxMaterialCorner = (_vbt->_maxCorner - _vbt->_minCorner)*(float)_vbt->_unitSpacingInv;
	for (int i = 0; i < 3; ++i)
		_vbt->_gridSize[i] = _gridSize[i] = 1 + (int)std::floor(maxMaterialCorner.xyz[i]);
	_vMatCoords.clear();
	_vMatCoords.assign(_mt->numberOfVertices(), Vec3f());
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i){
		Vec3f Vf;
		Vf.set((const float(&)[3])*_mt->vertexCoordinate(i));
		_vMatCoords[i].set((Vf -_vbt->_minCorner)* (float)_vbt->_unitSpacingInv);
	}
	_vertexTetCentroids.clear();
	_vertexTetCentroids.assign(_mt->numberOfVertices(), bccTetCentroid());
	_vbt->_barycentricWeights.clear();
	_vbt->_barycentricWeights.assign(_mt->numberOfVertices(), Vec3f());
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i){
		_vbt->gridLocusToLowestTetCentroid(_vMatCoords[i], _vertexTetCentroids[i]);
		// set barycentric coordinate within that tet
		auto& bw = _vbt->_barycentricWeights[i];
		_vbt->gridLocusToBarycentricWeight(_vMatCoords[i], _vertexTetCentroids[i], bw);
		// any vertex on a tet face boundary should be nudged inside
		bool nudge = false;
		for (int j = 0; j < 3; ++j) {
			if (bw[j] < 1e-4f || bw[j] > 0.9999f)
				nudge = true;
		}
		if (bw[0] + bw[1] + bw[2] > 0.9999)
			nudge = true;
		if (nudge) {  // vertices on tet boundarys must be pulled inside to avoid dual identities
			Vec3f nV = (Vec3f(0.25f, 0.25f, 0.25f) - bw) * 0.0002f;
			bw += nV;
			_vbt->barycentricWeightToGridLocus(_vertexTetCentroids[i], bw, _vMatCoords[i]);
		}
	}
	// setup lines parallel with Z axis
	evenXy.assign(_gridSize[0] >> 1, std::vector<std::multimap<double, bool> >());  // 0th i always empty
	oddXy.assign(_gridSize[0] >> 1, std::vector<std::multimap<double, bool> >());
	int gsy = _gridSize[1] >> 1;
	for (int n = _gridSize[0]>> 1, i = 0; i < n; ++i) {
		evenXy[i].assign(gsy, std::multimap<double, bool>());  // 0th j always empty
		oddXy[i].assign(gsy, std::multimap<double, bool>());
	}
	return true;
}

void vnBccTetCutter_omp::inputTriangle(int tri, std::unordered_map<bccTetCentroid, std::vector<int>, bccTetCentroidHasher>& tc_loc, std::vector<zIntrsct>& zi_loc) {  // main routine for finding triangle intersections with Z grid lines and with tetCentroid locations
	int xy[4] = { INT_MAX, -2, INT_MAX, -2 };
	int* tr = _mt->triangleVertices(tri);
	Vec3f T[3];
	bccTetCentroid tc[3];
	for (int i = 0; i < 3; ++i) {
		T[i] = _vMatCoords[tr[i]];
		tc[i] = _vertexTetCentroids[tr[i]];
		int f = (int)std::floor(T[i][0]);
		if (f < xy[0])
			xy[0] = f;
		if (f > xy[1])
			xy[1] = f;
		f = (int)std::floor(T[i][1]);
		if (f < xy[2])
			xy[2] = f;
		if (f > xy[3])
			xy[3] = f;
	}
	// triangle-triangle intersection used to non-recursively grow tets who have a face intersecting triangle
	// Of 3 algorithms; Moller, Devillers-Guigue, and Shen-Heng-Tang this routine uses the last of these.
	// Routine modified to set one triangle fixed to evaluate against multiple others.
	std::set<bccTetCentroid> triTetsS;
	for (int i = 0; i < 3; ++i)
		triTetsS.insert(tc[i]);
	auto recurseTriangleTets = [&](bccTetCentroid& tetC) {
		std::map<bccTetCentroid, short> doTets;  // tet that need processing and its face which can be ignored
		doTets.insert(std::make_pair(tetC, -1));
		short gridLoci[4][3];
		auto dtit = doTets.begin();
		while (dtit != doTets.end()) {
			triTetsS.insert(dtit->first);
			_vbt->centroidToNodeLoci(dtit->first, gridLoci);
			bccTetCentroid tAdj;
			short adjFace;
			for (int i = 0; i < 4; ++i) {
				if (i == dtit->second)  // face already processed
					continue;
				adjFace = _vbt->faceAdjacentTet(dtit->first, i, tAdj);
				if (triTetsS.find(tAdj) != triTetsS.end()) // already found
					continue;
				Vec3f V0(gridLoci[i]), V1(gridLoci[(i + 1) & 3]), V2(gridLoci[(i + 2) & 3]);
				if (tri_fixedTri_intersect(V0.xyz, V1.xyz, V2.xyz))
					doTets.insert(std::make_pair(tAdj, adjFace));
			}
			doTets.erase(dtit);
			dtit = doTets.begin();
		}
	};
	if (triTetsS.size() > 2 || (triTetsS.size() > 1 && !_vbt->adjacentCentroids(*triTetsS.begin(), *triTetsS.rbegin()))) {
		triTetsS.clear();
		if (!setFixedTriangle_Shen(T[0].xyz, T[1].xyz, T[2].xyz)) {
			std::string str("Triangle ");
			str.append(std::to_string(tri));
			str.append(" is degenerate with zero area and can't be processed.");
			throw(std::runtime_error(str.c_str()));
		}
		recurseTriangleTets(tc[0]);
		for (int i = 1; i < 3; ++i){
			if (triTetsS.find(tc[i]) == triTetsS.end()) // should be there if recurse OK
				recurseTriangleTets(tc[i]);
		}
	}
//	tetTriangles ttr;
//	ttr.tetindx = -1;
//	ttr.tris.clear();
	for (auto& tt : triTetsS) {
		auto pr = tc_loc.insert(std::make_pair(tt, std::vector<int>()));
		pr.first->second.push_back(tri);
	}
	// now get any Z line intersects
	T[1] -= T[0];
	T[2] -= T[0];
	bool solidBegin = T[1][0] * T[2][1] - T[1][1] * T[2][0] < 0.0f;  // negative Z starts a solid
	auto triIntersectZ = [&](const int x, const int y, double& z) ->bool {
		Mat2x2d M;
		M.x[0] = T[1][0];  M.x[1] = T[1][1];  M.x[2] = T[2][0]; M.x[3] = T[2][1];
		Vec2d R = M.Robust_Solve_Linear_System(Vec2d(x - T[0][0], y - T[0][1]));
		if (R[0] < -1e-5f || R[0] >= 1.00001f || R[1] < -1e-5f || R[1] >= 1.00001f || R[0] + R[1] >= 1.00001f)
			return false;
		z = (T[0][2] + T[1][2] * R[0] + T[2][2] * R[1]);
		return true;
	};
	for (int i = xy[0] + 1; i <= xy[1]; ++i) {
		bool odd = i & 1;
		for (int j = xy[2] + 1; j <= xy[3]; ++j) {
			if (odd != (bool)(j & 1))
				continue;
			zIntrsct zi;
			if (odd) {
				if (triIntersectZ(i, j, zi.zInt)) {
					zi.odd = true;
					zi.x = (i - 1) >> 1;
					zi.y = (j - 1) >> 1;
					zi.solidBegin = solidBegin;
					zi_loc.push_back(zi);
				}
			}
			else {
				if (triIntersectZ(i, j, zi.zInt)){
					zi.odd = false;
					zi.x = (i - 2) >> 1;
					zi.y = (j - 2) >> 1;
					zi.solidBegin = solidBegin;
					zi_loc.push_back(zi);
				}
			}
		}
	}
	return;
}

void vnBccTetCutter_omp::fillNonVnTetCenter()
{
	_vbt->_firstInteriorTet = _vbt->_tetNodes.size();
	//  Use sequential z to get surrounding frame with less hash table searches.
	auto nextZnode = [&](const int prevNode, std::array<short, 3>& prevLocus) ->int {  // returns next node in z sequence from prevNode at locus. If former is -1 look for locus. If not found return -1
		prevLocus[2] += 2;
		if (prevNode > -1) {
			auto pn = _vbt->_nodeGridLoci[prevNode + 1].data();
			if (prevLocus[0] == pn[0] && prevLocus[1] == pn[1] && prevLocus[2] == pn[2])
				return prevNode + 1;
		}
		auto in = _interiorNodes.find(prevLocus);
		if (in != _interiorNodes.end())
			return in->second;
		else return -1;
	};
	auto findNode = [&](std::array<short, 3>& locus) ->int {
		auto in = _interiorNodes.find(locus);
		if (in != _interiorNodes.end())
			return in->second;
		else
			return -1;

	};
	auto zBlock = [&](int start, int stop) {
		std::array<int, 4> prevQuad, nextQuad;
		std::array<short, 3> locs[4], yLoc, xLoc = _vbt->_nodeGridLoci[start];
		yLoc = xLoc;
		yLoc[1] += 2;
		int nextY = findNode(yLoc);
		xLoc[0] += 2;
		int nextX = findNode(xLoc);
		locs[0] = xLoc;
		locs[0][0] -= 3;
		--locs[0][1];
		--locs[0][2];
		locs[1] = locs[0];
		locs[1][0] += 2;
		locs[2] = locs[0];
		locs[2][1] += 2;
		locs[3] = locs[2];
		locs[3][0] += 2;
		for (int i = 0; i < 4; ++i)
			prevQuad[i] = findNode(locs[i]);
		for (int i = start; i < stop; ++i) {
			for (int j = 0; j < 4; ++j)
				nextQuad[j] = nextZnode(prevQuad[j], locs[j]);
			// make tets
			bccTetCentroid center = { (unsigned short)(xLoc[0] - 2), (unsigned short)xLoc[1], (unsigned short)(xLoc[2] + 1) };
			for (int j = 0; j < 3; ++j)
				center[j] <<= 1;
			std::array<int, 4> tet;
			if (i < stop - 1) { // next even z ray
				tet[0] = i;
				tet[1] = i + 1;
				if (nextQuad[0] > -1 && nextQuad[1] > -1) {
					tet[2] = nextQuad[0];
					tet[3] = nextQuad[1];
					auto tc = center;
					--tc[1];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (nextQuad[2] > -1 && nextQuad[3] > -1) {
					tet[2] = nextQuad[3];
					tet[3] = nextQuad[2];
					auto tc = center;
					++tc[1];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (nextQuad[0] > -1 && nextQuad[2] > -1) {
					tet[0] = nextQuad[0];
					tet[1] = nextQuad[2];
					tet[2] = i + 1;
					tet[3] = i;
					auto tc = center;
					--tc[0];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (nextQuad[1] > -1 && nextQuad[3] > -1) {
					tet[0] = nextQuad[1];
					tet[1] = nextQuad[3];
					tet[2] = i;
					tet[3] = i + 1;
					auto tc = center;
					++tc[0];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
			}
			center[2] -= 2;  // next even x ray
			center[0] += 2;
			if (nextX > -1) {
				tet[0] = i;
				tet[1] = nextX;
				if (nextQuad[3] > -1 && nextQuad[1] > -1) {
					tet[2] = nextQuad[3];
					tet[3] = nextQuad[1];
					auto tc = center;
					++tc[2];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[1] > -1 && prevQuad[3] > -1) {
					tet[2] = prevQuad[1];
					tet[3] = prevQuad[3];
					auto tc = center;
					--tc[2];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[1] > -1 && nextQuad[1] > -1) {
					tet[0] = prevQuad[1];
					tet[1] = nextQuad[1];
					tet[2] = nextX;
					tet[3] = i;
					auto tc = center;
					--tc[1];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[3] > -1 && nextQuad[3] > -1) {
					tet[0] = prevQuad[3];
					tet[1] = nextQuad[3];
					tet[2] = i;
					tet[3] = nextX;
					auto tc = center;
					++tc[1];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
			}
			center[0] -= 2;
			center[1] += 2;
			if (nextY > -1) {
				tet[0] = i;
				tet[1] = nextY;
				if (prevQuad[2] > -1 && nextQuad[2] > -1) {
					tet[2] = prevQuad[2];
					tet[3] = nextQuad[2];
					auto tc = center;
					--tc[0];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[3] > -1 && nextQuad[3] > -1) {
					tet[2] = nextQuad[3];
					tet[3] = prevQuad[3];
					auto tc = center;
					++tc[0];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[2] > -1 && prevQuad[3] > -1) {
					tet[0] = prevQuad[2];
					tet[1] = prevQuad[3];
					tet[2] = nextY;
					tet[3] = i;
					auto tc = center;
					--tc[2];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (nextQuad[2] > -1 && nextQuad[3] > -1) {
					tet[0] = nextQuad[2];
					tet[1] = nextQuad[3];
					tet[2] = i;
					tet[3] = nextY;
					auto tc = center;
					++tc[2];
					if (_centroidIndices.find(tc) == _centroidIndices.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
			}
			// advance to next frame
			for (int j = 0; j < 4; ++j)
				prevQuad[j] = nextQuad[j];
			nextX = nextZnode(nextX, xLoc);
			nextY = nextZnode(nextX, yLoc);
		}
	};
	// interior tet nodes are created sequentially in z.  Do each z block as nucleus around even lattice
	for (int n = _interiorNodes.size(), i = 0; i < n; ) {
		auto lp0 = _vbt->_nodeGridLoci[i].data();
		if (lp0[0] & 1) {
			++i;
			continue;
		}
		int j = i + 1;
		auto lp1 = _vbt->_nodeGridLoci[j].data();
		while (lp0[0] == lp1[0] && lp0[1] == lp1[1] && lp0[2] + 2 == lp1[2]) { // this z segment exists
			++j;
			lp0 = lp1;
			lp1 = _vbt->_nodeGridLoci[j].data();
		}
		zBlock(i, j);
		i = j;
	}
}

