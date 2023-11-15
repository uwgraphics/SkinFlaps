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

#include "centroidBccTetCutter_tbb.h"

std::vector<centroidBccTetCutter_tbb::tetTriangles> centroidBccTetCutter_tbb::_tetTris;
std::unordered_map<bccTetCentroid, int, centroidBccTetCutter_tbb::bccTetCentroidHasher> centroidBccTetCutter_tbb::_centroidIndices;


bool centroidBccTetCutter_tbb::makeFirstVnTets(materialTriangles *mt, vnBccTetrahedra *vbt, int maximumGridDimension)
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



bool centroidBccTetCutter_tbb::remakeVnTets(materialTriangles* mt)
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
		throw(std::logic_error("topological error in centroidBccTetCutter_tbb::remakeVnTets()\n"));
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

bool centroidBccTetCutter_tbb::tetCutCore() {  // same cut core for first cut and recuts
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	int it, nOmpRange = _mt->numberOfTriangles();
	// private variables
	std::unordered_map<bccTetCentroid, std::vector<int>, bccTetCentroidHasher> centTris_local;
	std::vector<std::pair<bccTetCentroid, std::vector<int> > > cT_global;
	std::vector<zIntrsct> zInt_local;
	// COURT note the nowait.  Also #pragma statements cannot have a curly brace immediately after them, but must be on a new line.

//	_centroidIndices, , oddXy, evenXy, , _tetTris

//#pragma omp parallel shared(nOmpRange, cT_global) private(it, centTris_local, zInt_local)
	{
		centTris_local.clear();
		zInt_local.clear();
//#pragma omp for schedule(static, 1)
		for (it = 0; it < nOmpRange; ++it)
			inputTriangle(it, centTris_local, zInt_local);
//#pragma omp critical (centTris_local)
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
//#pragma omp critical (zInt_local)
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

//#pragma omp parallel
	{
		// private variables
		std::vector<newTet> newTets_local;
		std::unordered_map<std::array<short, 3>, std::list<nodeTetSegment>, arrayShort3Hasher> nts_local;

//#pragma omp for schedule(dynamic) nowait
		// COURT note the nowait.  Also #pragma statements cannot have a curly brace immediately after them, but must be on a new line.
		for (int i = 0; i < nOmpRange; ++i) {
			getConnectedComponents(_tetTris[i], newTets_local, nts_local);  // for this centroid split its triangles into solid connected components
		}
//#pragma omp critical
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

void centroidBccTetCutter_tbb::createInteriorNodes() {
	short z;
	std::array<short, 3> s3;
	double zTol = _gridSize[2] * 1e-3;
	auto solidFilter = [&](std::multimap<double, bool>& mm) {  // isolate solid runs from surface edge and vertex hits and correct switch of coincident surfaces.
		auto mit = mm.begin();
		bool inSolid = true;
		while (mit != mm.end()) {
			auto mNext = mit;
			++mNext;
			if (mNext == mm.end())
				break;
			while ( mNext->first - mit->first < zTol) {
				if (mit->second == mNext->second) {  // surface edge or vertex hit
					mm.erase(mNext);
					mNext = mit;
					++mNext;
					if (mNext == mm.end())
						break;
				}
				else if (mit->second != inSolid) {  // coincident surfaces
					mit->second = inSolid;
					mNext->second = !inSolid;
					break;
				}
				else
					break;
			}
			if (mit->second != inSolid)
				throw(std::logic_error("Solid ordering error in createInteriorNodes()\n"));
			inSolid = !inSolid;
			++mit;
		}
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

			if (xi == 20 && yi == 3)
				int junk = 0;

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

void centroidBccTetCutter_tbb::assignExteriorTetNodes(int exteriorNodeNum, std::vector<extNode>& eNodes) {
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

int centroidBccTetCutter_tbb::nearestRayPatchHit(const Vec3d &rayBegin, Vec3d rayEnd, const std::vector<int>& tris, Vec3d &hitP, double& distanceSq) {  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
	distanceSq = DBL_MAX;
	Vec3d rE(rayEnd - rayBegin);
	Vec3d N;
	for (auto& t : tris) {
		int* tr = _mt->triangleVertices(t);
		Vec3d tv[3];
		for (int i = 0; i < 3; ++i)
			tv[i].set(_vMatCoords[tr[i]]);
		tv[1] -= tv[0];
		tv[2] -= tv[0];
		double dSq;
		Vec3d P, R;
		Mat3x3d M;
		M.Initialize_With_Column_Vectors(tv[1], tv[2], -rE);
		R = M.Robust_Solve_Linear_System(Vec3d(rayBegin) - tv[0]);
		if (R[0] < -1e-8 || R[0] > 1.00000001 || R[1] < -1e-8 || R[1] > 1.00000001 || R[2] < 0.0 || R[2] > 1.0001 || R[0] + R[1] > 1.00000001)
			continue;
		P = rE * R[2];
		dSq = P * P;
		if (distanceSq > dSq) {
			distanceSq = dSq;
			hitP = rayBegin + rE * R[2];
			N = tv[1] ^ tv[2];
		}
	}
	if (distanceSq < DBL_MAX) {
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


bool centroidBccTetCutter_tbb::isInsidePatch(const Vec3d& P, const std::vector<int>& tris, Vec3d &closestP) {
	double d, minD = DBL_MAX;
	Vec3d N;
	for (auto t : tris) {
		int* tr = _mt->triangleVertices(t);
		Vec3d T[3], Q;
		for (int i = 0; i < 3; ++i)
			T[i].set(_vMatCoords[tr[i]]);
		T[1] -= T[0];
		T[2] -= T[0];
		Mat2x2d M(T[1] * T[1], T[1] * T[2], 0.0, T[2] * T[2]);
		M.x[2] = M.x[1];
		Q = P - T[0];
		Vec2d R = M.Robust_Solve_Linear_System(Vec2d(T[1] * Q, T[2] * Q));
		if (R[0] < 0.02)
			R[0] = 0.02;
		else if (R[0] > 0.98)
			R[0] = 0.98;
		else;
		if (R[1] < 0.02)
			R[1] = 0.02;
		else if (R[1] > 0.98)
			R[1] = 0.98;
		else;
		if (R[0] + R[1] > 0.98) {
			double d = 0.98/(R[0] + R[1]);
			R[0] *= d;
			R[1] *= d;
		}
		Q = T[0] + T[1] * R[0] + T[2] * R[1];
		d = (P - Q).length2();
		if (minD > d) {
			minD = d;
			closestP = Q;
			N = T[1] ^ T[2];
		}
	}
	return N * (closestP - P) >= 0.0;
}

void centroidBccTetCutter_tbb::getConnectedComponents(tetTriangles& tt, std::vector<newTet>& nt_vec, std::unordered_map<std::array<short, 3>, std::list<nodeTetSegment>, arrayShort3Hasher>& local_nts) {
	// for this centroid split its triangles into single solid connected components
	std::set<int> trSet, ts;
	trSet.insert(tt.tris.begin(), tt.tris.end());
	tt.tris.clear();
	struct patch {
		std::vector<int> tris;
		int tetIndex;
		std::array<int, 4> tetNodes;
	}p;
	std::list<patch> patches;
	typedef std::list<patch>::iterator PListIt;
	p.tetNodes = { -1, -1, -1, -1 };
	p.tetIndex = -1;
	patches.push_back(p);
	auto trVec = &patches.back().tris;  // on entry only one triangle list for this centroid
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
			trVec = &patches.back().tris;
		}
	}
	// find interior nodes for each component, joining patches as appropriate
	short gl[4][3];
	_vbt->centroidToNodeLoci(tt.tc, gl);
	int inodes[4] = { -1, -1, -1, -1 };
	for (int i = 0; i < 4; ++i) {
		std::array<short, 3> loc = { gl[i][0], gl[i][1], gl[i][2] };
		auto iit = _interiorNodes.find(loc);
		if (iit != _interiorNodes.end())
			inodes[i] = iit->second;
	}
	if (patches.size() < 2) {  // no interpatch computation required
		for (int i = 0; i < 4; ++i)
			patches.front().tetNodes[i] = inodes[i];
	}
	else {  // resolve which patches should be combined and which left separate, based on enclosure of solid section
		Vec3d tetNodes[4];
		for (int i = 0; i < 4; ++i)
			tetNodes[i].set(gl[i]);
		struct edgeIntrsct {
			bool nextOutside;
			patch* pptr;
		};
		std::list<patch*> unsolvedInsidePatches;
		std::multimap<double, edgeIntrsct> edges[6];
		for (auto pit = patches.begin(); pit != patches.end(); ++pit) {
			bool edgeCut = false;
			for (auto& t : pit->tris) {
				int* tr = _mt->triangleVertices(t);
				Vec3d T[3], P, Q;
				for (int i = 0; i < 3; ++i)
					T[i] = _vMatCoords[tr[i]];
				T[1] -= T[0];
				T[2] -= T[0];
				auto getEdgeIntersect = [&](int edge) {
					Mat3x3d M;
					M.Initialize_With_Column_Vectors(T[1], T[2], P - Q);
					Vec3d R = M.Robust_Solve_Linear_System(P - T[0]);
					if (R[0] < -1e-8 || R[0] > 1.00000001 || R[1] < -1e-8 || R[1] > 1.00000001 || R[2] <= 0.0 || R[2] >= 1.0 || R[0] + R[1] >= 1.00000001)  // don't want triangle vertex hits, only edges
						return;
					edgeIntrsct ei;
					ei.pptr = &(*pit);
					ei.nextOutside = (T[1] ^ T[2]) * (Q - P) >= 0;
					edges[edge].insert(std::make_pair(R[2], ei));
					edgeCut = true;
				};
				P = tetNodes[3];
				for (int i = 0; i < 4; ++i) {
					Q = tetNodes[i];
					getEdgeIntersect(i);
					P = Q;
				}
				P = tetNodes[0];
				Q = tetNodes[2];
				getEdgeIntersect(4);
				P = tetNodes[1];
				Q = tetNodes[3];
				getEdgeIntersect(5);
			}
			if (!edgeCut) {  // This patch is interior to all edges
				// if facing exterior to any tet node, must be independent of all patches with an edge intersection processed above and must be exterior to all tet nodes.
				// if interior could be interior to any existing solid
				Vec3d P;
				if (isInsidePatch(tetNodes[0], pit->tris, P))  // interior faces all nodes.
					unsolvedInsidePatches.push_back(&(*pit));
			}
		}
		std::list<std::set<patch*> > combineP;
		bool interiorLinks[2][4] = { false,false,false,false,false,false,false,false }, zeroEmpty = true, secondPair = false;
		auto combinePatches = [&](std::list<patch*>& pl) {
			auto cp = combineP.begin();
			while (cp != combineP.end()) {
				for (auto pp : pl) {
					if (cp->find(pp) != cp->end()) {
						for (auto pa : pl)
							cp->insert(pa);
						return;
					}
				}
				++cp;
			}
			if (cp == combineP.end()) {
				combineP.push_back(std::set<patch*>());
				for (auto pp : pl)
					combineP.back().insert(pp);
			}
		};
		auto getEdgeNodeIndices = [](const int edge, int& idx0, int& idx1) {
			if (edge < 1) { idx0 = 3; idx1 = 0; }
			else if (edge < 4) { idx0 = edge - 1; idx1 = edge; }
			else if (edge < 5) { idx0 = 0; idx1 = 2; }
			else { idx0 = 1; idx1 = 3; }
		};
		int nbegin, nend;
		double sameHit = std::max({ _gridSize[0], _gridSize[1], _gridSize[2] }) * 6e-4;
		for (int i = 0; i < 6; ++i) {
			getEdgeNodeIndices(i, nbegin, nend);
			if (!edges[i].empty()) {
				bool inSolid = inodes[nbegin] > -1;
				auto eit = edges[i].begin();
				auto enext = eit;
				++enext;
				while (enext != edges[i].end()) {
					if (enext->first - eit->first < sameHit) {  // possible coincident surface or minor collision hit.  Both sides guaranteed inside.
						if (eit->second.nextOutside == enext->second.nextOutside) {  // possible surface edge or vertex multi hit
							if (inSolid != eit->second.nextOutside) {
								edges[i].erase(enext);
								enext = eit;
								++enext;
								if (enext == edges[i].end())
									break;
							}
						}
						else if (inSolid != eit->second.nextOutside) {  // flip if necessary
							auto swap = enext->second;
							enext->second = eit->second;
							eit->second = swap;
						}
						else
							;
					}
					if (inSolid != eit->second.nextOutside)
						throw(std::logic_error("Solid ordering error in getConnectedComponents()\n"));
					inSolid = !inSolid;
					eit = enext;
					++enext;
				}
				eit = edges[i].begin();
				assert(eit->second.nextOutside == inodes[nbegin] > -1);
				enext = eit;
				++enext;
				while (enext != edges[i].end()) {
					if (!eit->second.nextOutside && enext->second.nextOutside) {  // solid interval. Combine patches
						if (eit->second.pptr != enext->second.pptr) {  // can be equal
							std::list<patch*> pl;
							pl.push_back(eit->second.pptr);
							pl.push_back(enext->second.pptr);
							combinePatches(pl);
						}
					}
					eit = enext;
					++enext;
				}
				assert(eit->second.nextOutside == inodes[nend] < 0);
				auto et = edges[i].begin();
				if (et->second.nextOutside) {
					assert(inodes[nbegin] > -1);
					et->second.pptr->tetNodes[nbegin] = inodes[nbegin];
				}
				else
					assert(inodes[nbegin] < 0);
				et = edges[i].end();
				--et;
				if (!et->second.nextOutside) {
					assert(inodes[nend] > -1);
					et->second.pptr->tetNodes[nend] = inodes[nend];
				}
				else
					assert(inodes[nend] < 0);
			}
			else {
				if (inodes[nbegin] > -1) {
					assert(inodes[nend] > -1);
					if (zeroEmpty) {
						interiorLinks[0][nbegin] = true;
						interiorLinks[0][nend] = true;
						zeroEmpty = false;
					}
					else if (interiorLinks[0][nbegin] || interiorLinks[0][nend]) { // 3rd or 4th member
						if (secondPair) {
							for (int k = 0; k < 4; ++k)
								interiorLinks[0][k] = true;
							secondPair = false;
						}
						else {
							interiorLinks[0][nbegin] = true;
							interiorLinks[0][nend] = true;
						}
					}
					else {  // can have 2 pairs
						interiorLinks[1][nbegin] = true;
						interiorLinks[1][nend] = true;
						secondPair = true;
					}
				}
			}
		}
		// add linked interior nodes to patches
		for (auto& p : patches) {
			for (int i = 0; i < 4; ++i) {
				if (p.tetNodes[i] > -1) {
					if (interiorLinks[0][i]) {
						for (int j = 0; j < 4; ++j) {
							if (i == j)
								continue;
							if (interiorLinks[0][j])
								p.tetNodes[j] = inodes[j];
						}
						if (!secondPair)
							break;
					}
					else if (secondPair && interiorLinks[1][i]) {
						for (int j = 0; j < 4; ++j) {
							if (i == j)
								continue;
							if (interiorLinks[1][j]) {
								p.tetNodes[j] = inodes[j];
								break;
							}
						}
					}
					else
						;
				}
			}
		}
		// look for patches that combine due to shared interior node (ie solid span bridges more than one edge)
		for (int i = 0; i < 4; ++i) {
			if (inodes[i] < 0)
				continue;
			std::list<patch*> pl;
			for (auto& p : patches) {
				if (p.tetNodes[i] > -1) {
					assert(p.tetNodes[i] == inodes[i]);
					pl.push_back(&p);
				}
			}
			if (pl.size() > 1)
				combinePatches(pl);
		}
		// combine patches that enclose a solid
		auto removePatch = [&](patch* p) {
			for (auto pi = patches.begin(); pi != patches.end(); ++pi)
				if (&(*pi) == p) {
					patches.erase(pi);
					break;
				}
		};
		for (auto& cp : combineP) {
			assert(cp.size() > 1);
			auto cpfirst = cp.begin();
			auto cpnext = cpfirst;
			++cpnext;
			while (cpnext != cp.end()) {
				(*cpfirst)->tris.insert((*cpfirst)->tris.end(), (*cpnext)->tris.begin(), (*cpnext)->tris.end());
				for (int i = 0; i < 4; ++i) {
					if ((*cpnext)->tetNodes[i] > -1)
						(*cpfirst)->tetNodes[i] = (*cpnext)->tetNodes[i];
				}
				removePatch(*cpnext);
				++cpnext;
			}
		}
		if (unsolvedInsidePatches.size() == patches.size()) {
			if (inodes[0] < 0 || inodes[1] < 0 || inodes[2] < 0 || inodes[3] < 0)
				throw(std::logic_error("An all inside node tet incorrectly processed in getConnectedComponents().\n"));
			if (unsolvedInsidePatches.size() > 1) {
				auto ui0 = unsolvedInsidePatches.begin();
				auto uiNext = ui0;
				++uiNext;
				while (uiNext != unsolvedInsidePatches.end()) {
					(*ui0)->tris.insert((*ui0)->tris.end(), (*uiNext)->tris.begin(), (*uiNext)->tris.end());
					removePatch(*uiNext);
				}
				for (int i = 0; i < 4; ++i)
					(*ui0)->tetNodes[i] = inodes[i];
			}
			unsolvedInsidePatches.clear();
		}
		for (auto ui : unsolvedInsidePatches) {  // these interior patches must be inside another solid
			// must find point on patch guaranteed inside tet closest to node 0
			Vec3d P((const unsigned short(&)[3])(*tt.tc.data())), centerP;
			P *= 0.5;  // closest point to tet center
			isInsidePatch(P, ui->tris, centerP);
			double dsq, minD = DBL_MAX;
			int connected = nearestRayPatchHit(tetNodes[0], centerP, ui->tris, P, dsq);
			patch* enclosingPatch = nullptr;
			for (auto pit = patches.begin(); pit != patches.end(); ++pit) {
				if (&(*pit) == ui)
					continue;
				Vec3d hitP;
				int connected = nearestRayPatchHit(P, tetNodes[0], pit->tris, hitP, dsq);  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
				if (connected < 0 && dsq < minD) {
					minD = dsq;
					enclosingPatch = &(*pit);
				}
			}
			if (enclosingPatch) {
				enclosingPatch->tris.insert(enclosingPatch->tris.end(), ui->tris.begin(), ui->tris.end());
				assert(ui->tetNodes[0] < 0 && ui->tetNodes[0] < 0 && ui->tetNodes[0] < 0 && ui->tetNodes[0] < 0);
				removePatch(ui);
			}
		}
		unsolvedInsidePatches.clear();
	}
	for (auto& cTet : patches) {
		nt_vec.push_back(newTet());
		nt_vec.back().tetIdx = _nSurfaceTets.fetch_add(1);
		cTet.tetIndex = nt_vec.back().tetIdx;
		nt_vec.back().tc = tt.tc;
		nt_vec.back().tetNodes = cTet.tetNodes;
		nt_vec.back().tris.assign(cTet.tris.begin(), cTet.tris.end());
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
				nts.tetNodeTris.assign(p.tris.begin(), p.tris.end());
				pr.first->second.push_back(nts);
			}
		}
	}
}

bool centroidBccTetCutter_tbb::setupBccIntersectionStructures(int maximumGridDimension)
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

void centroidBccTetCutter_tbb::inputTriangle(int tri, std::unordered_map<bccTetCentroid, std::vector<int>, bccTetCentroidHasher>& tc_loc, std::vector<zIntrsct>& zi_loc) {  // main routine for finding triangle intersections with Z grid lines and with tetCentroid locations

//	COURT - fill these instead for tbb implementation
//	tbb::concurrent_vector<zIntrsct> _zIntersects;
//	tbb::concurrent_unordered_map<bccTetCentroid, std::vector<int>, bccTetCentroidHasher> _centroidTriangles;  // careful does not support concurrent erase


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
	for (auto& tt : triTetsS) {
		auto pr = tc_loc.insert(std::make_pair(tt, std::vector<int>()));
		pr.first->second.push_back(tri);
	}
	// now get any Z line intersects
	Vec3d Td[3];
	for (int i = 0; i < 3; ++i)
		Td[i].set(T[i]);
	Td[1] -= Td[0];
	Td[2] -= Td[0];
	bool solidBegin = Td[1][0] * Td[2][1] - Td[1][1] * Td[2][0] < 0.0f;  // negative Z starts a solid
	auto triIntersectZ = [&](const int x, const int y, double& z) ->bool {
		Mat2x2d M;
		M.x[0] = Td[1][0];  M.x[1] = Td[1][1];  M.x[2] = Td[2][0]; M.x[3] = Td[2][1];
		Vec2d R = M.Robust_Solve_Linear_System(Vec2d(x - Td[0][0], y - Td[0][1]));
		if (R[0] < -1e-8f || R[0] >= 1.00000001 || R[1] < -1e-8f || R[1] >= 1.00000001 || R[0] + R[1] >= 1.00000001)  // intentionally permissive
			return false;
		z = (Td[0][2] + Td[1][2] * R[0] + Td[2][2] * R[1]);
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

void centroidBccTetCutter_tbb::fillNonVnTetCenter()
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

