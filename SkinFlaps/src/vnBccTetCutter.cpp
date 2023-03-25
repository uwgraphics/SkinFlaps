#include <assert.h>
#include <limits>
#include <algorithm>
#include <cmath>
#include <functional>
#include "Mat2x2d.h"
#include "Vec3d.h"
#include "Mat2x2f.h"
#include "Mat3x3f.h"
#include "triTriIntersect_Shen.h"
#include "tri_tri_intersect_Guigue.h"

#include <chrono>  // for openMP timing
#include <ctime>  // nuke after openMP debug
#include <fstream>

#include "vnBccTetCutter.h"

bool vnBccTetCutter::subcutMacroTets(std::vector<int>& macroTets ) {
	if (_mt->findAdjacentTriangles(true))	return false;
	std::vector<int> surfaceTris;
	std::vector<std::array<int, 3> > boundaryNodeTris;
	std::sort(macroTets.begin(), macroTets.end());  // COURT if lots should be hash table
	std::vector<int> subTris;
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
		int* tr = _mt->triangleVertices(i);
		for (int j = 0; j < 3; ++j) {  //  COURT good only once.  Fix later.
			if (std::binary_search(macroTets.begin(), macroTets.end(), _vbt->_vertexTets[tr[j]])) {
				subTris.push_back(i);
				break;
			}
		}
	}
	// delete old macrotets but keep all boundary triangles
	std::set<int> oldNodes;
	for (int& t : macroTets) {
		auto& tc = _vbt->_tetCentroids[t];
		if (!std::binary_search(_macroTetCentroidSubset.begin(), _macroTetCentroidSubset.end(), tc))
			_macroTetCentroidSubset.push_back(tc);
		for (int i = 0; i < 4; ++i)
			oldNodes.insert(_vbt->_tetNodes[t][i]);
		_vbt->_tetNodes[t][0] = -1;
	}
	std::sort(_macroTetCentroidSubset.begin(), _macroTetCentroidSubset.end());
	for (int n = _vbt->_tetNodes.size(), i = 0; i < n; ++i) {
		auto& tn = _vbt->_tetNodes[i];
		if (tn[0] < 0)
			continue;
		int on[4] = { -1, -1, -1, -1 }, nFound = 0;;
		for (int j = 0; j < 4; ++j) {
			if (oldNodes.find(tn[j]) != oldNodes.end()) {
				on[j] = tn[j];
				++nFound;
			}
		}
		if (nFound > 2) {
			for (int j = 0; j < 4; ++j) {
				if (on[j] > -1 && on[(j + 1) & 3] > -1 && on[(j + 2) & 3] > -1) {
					boundaryNodeTris.push_back(std::array<int, 3>());
					auto& an = boundaryNodeTris.back();
					an[0] = on[j];
					if (j & 1) {  // get correct clockwiseness
						an[1] = on[(j + 2) & 3];
						an[2] = on[(j + 1) & 3];
					}
					else {
						an[1] = on[(j + 1) & 3];
						an[2] = on[(j + 2) & 3];

					}
				}
			}
		}
	}
		
	// leave current intersection data structures except;
//	_vMatCoords.clear();
//	_vMatCoords.assign(_mt->numberOfVertices(), Vec3f());
//	_vertexTetCentroids.clear();
//	_vertexTetCentroids.assign(_mt->numberOfVertices(), bccTetCentroid());
//	_vbt->_barycentricWeights.clear();
//	_vbt->_barycentricWeights.assign(_mt->numberOfVertices(), Vec3f());
	std::set<int> surfaceVerts;
	for (auto& st : subTris) {
		int* tr = _mt->triangleVertices(st);
		for (int i = 0; i < 3; ++i)
			surfaceVerts.insert(tr[i]);
	}
	// leave old vertices in macrotets where they were barycentrically, but find new locs in microtet space.
	for (auto& sv : surfaceVerts) {
		if (!std::binary_search(macroTets.begin(), macroTets.end(), _vbt->_vertexTets[sv]))
			continue;
		_vbt->gridLocusToTetCentroid(_vMatCoords[sv], _vertexTetCentroids[sv]);
		// set barycentric coordinate within that tet
		auto& bw = _vbt->_barycentricWeights[sv];
		_vbt->gridLocusToBarycentricWeight(_vMatCoords[sv], _vertexTetCentroids[sv], bw);
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
			_vbt->barycentricWeightToGridLocus(_vertexTetCentroids[sv], bw, _vMatCoords[sv]);
		}
	}
	// We really don't need all this memory, but it's a pain to subcut it.
	// setup lines parallel with Z axis
	evenXy.assign(_gridSize[0] >> 1, std::vector<std::multimap<double, bool> >());  // 0th i always empty
	oddXy.assign(_gridSize[0] >> 1, std::vector<std::multimap<double, bool> >());
	int gsy = _gridSize[1] >> 1;
	for (int n = _gridSize[0] >> 1, i = 0; i < n; ++i) {
		evenXy[i].assign(gsy, std::multimap<double, bool>());  // 0th j always empty
		oddXy[i].assign(gsy, std::multimap<double, bool>());
	}
	_centroidTriangles.clear();  // COURT later reserve this hash table
	for (int n = subTris.size(), i = 0; i < n; ++i)
		inputTriangle(subTris[i]);
	// cull out data not in this tet subset
	for (auto ct = _centroidTriangles.begin(); ct != _centroidTriangles.end(); ){
		// all these are microtets so if its centroid not inside subset delete it
		Vec3f tc(ct->first[0], ct->first[1], ct->first[2]);
		tc *= 0.5f;
		bool found = false;
		for (auto& mtc : _macroTetCentroidSubset) {
			if (_vbt->insideTet(mtc, tc)) {
				found = true;
				break;
			}
		}
		if (found)
			++ct;
		else
			ct = _centroidTriangles.erase(ct);
	}
	for (int n = evenXy.size(), i = 0; i < n; ++i) {
		for (int m = evenXy[i].size(), j = 0; j < m; ++j) {
			auto& iLine = evenXy[i][j];
			if (iLine.empty())
				continue;
			Vec3f loc;
			loc.X = (i << 1);
			loc.Y = (j << 1);
			for (auto zit = iLine.begin(); zit != iLine.end(); ) {
				loc.Z = (float)zit->first;
				bool found = false;
				for (auto& mtc : _macroTetCentroidSubset) {
					if (_vbt->insideTet(mtc, loc)) {
						found = true;
						break;
					}
				}
				if (found)
					++zit;
				else
					zit = iLine.erase(zit);
			}
		}
	}
	for (int n = oddXy.size(), i = 0; i < n; ++i) {
		for (int m = oddXy[i].size(), j = 0; j < m; ++j) {
			auto &iLine = oddXy[i][j];
			if (iLine.empty())
				continue;
			Vec3f loc;
			loc.X = (i << 1) + 1;
			loc.Y = (j << 1) + 1;
			for (auto zit = iLine.begin(); zit != iLine.end(); ) {
				loc.Z = (float)zit->first;
				bool found = false;
				for (auto& mtc : _macroTetCentroidSubset) {
					if (_vbt->insideTet(mtc, loc)) {
						found = true;
						break;
					}
				}
				if (found)
					++zit;
				else
					zit = iLine.erase(zit);
			}
		}
	}
	// backstop even and oddXy with boundaryNodeTris for creation of new interior nodes
	for (auto& bt : boundaryNodeTris) {
		// get Z line intersects.  From inputTriangle()
		Vec3f T[3];
		for (int i = 0; i < 3; ++i)
			T[i].set((short (&)[3])*_vbt->_nodeGridLoci[bt[i]].data());
		int xy[4] = { INT_MAX, -2, INT_MAX, -2 };
		for (int i = 0; i < 3; ++i) {
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
				double z;
				// all boundaryTri intersects should make a legal bound to surfaceTri intersects.
				// surfaceTri intersects will overshoot bounds of macrotet.
				if (odd) {
					if (triIntersectZ(i, j, z))
						oddXy[(i - 1) >> 1][(j - 1) >> 1].insert(std::make_pair(z, solidBegin));
				}
				else {
					if (triIntersectZ(i, j, z))
						evenXy[(i - 2) >> 1][(j - 2) >> 1].insert(std::make_pair(z, solidBegin));
				}
			}
		}

	}
	// create and hash all interior nodes
	createInteriorNodes();

	return true;
}

bool vnBccTetCutter::makeFirstVnTets(materialTriangles *mt, vnBccTetrahedra *vbt, int maximumGridDimension)
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
	if(_mt->findAdjacentTriangles(true))	return false;
	if (!setupBccIntersectionStructures(maximumGridDimension))
		return false;
	_centroidTriangles.clear();  // COURT later reserve this hash table
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i)
		inputTriangle(i);
	// create and hash all interior nodes
	createInteriorNodes();
	_interiorNodes.clear();  // only nodes created thus far are interior nodes
	_interiorNodes.reserve(_vbt->_nodeGridLoci.size());
	for (int n = _vbt->_nodeGridLoci.size(), j = 0; j < n; ++j)
		_interiorNodes.insert(std::make_pair(_vbt->_nodeGridLoci[j], j));
	_exteriorNodeIndices.clear();
	_exteriorNodeIndices.reserve(_centroidTriangles.size() >> 1);    // COURT later check this reserve size
	_vbt->_tetCentroids.clear();
	_vbt->_tetCentroids.reserve(_centroidTriangles.size() >> 1);    // COURT perhaps assign() for multi threading
	for (auto& ct : _centroidTriangles)
		getConnectedComponents(ct.first, ct.second);  // for this centroid split its triangles into solid connected components
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_vbt->tetNumber());
	for (int n = _vbt->tetNumber(), i = 0; i < n; ++i)
		_vbt->_tetHash.insert(std::make_pair(_vbt->tetCentroid(i), i));
	// get tets where vertices reside
	_vbt->_vertexTets.clear();
	_vbt->_vertexTets.assign(mt->numberOfVertices(), -1);
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i) {
		auto ct = _centroidTriangles.find(_vertexTetCentroids[i]);
		if (ct == _centroidTriangles.end())
			throw(std::logic_error("Couldn't find a tetrahedron containing a vertex.\n"));
		else if (ct->second.size() < 2)
			_vbt->_vertexTets[i] = ct->second.front().tetindx;
		else {
			auto ltit = ct->second.begin();
			while (ltit != ct->second.end()) {
				auto tit = ltit->tris.begin();
				while (tit != ltit->tris.end()) {
					int *tr = _mt->triangleVertices(*tit);
					int j;
					for (j = 0; j < 3; ++j) {
						if (tr[j] == i) {
							_vbt->_vertexTets[i] = ltit->tetindx;
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
		}
	}
	assignExteriorTetNodes();
	_vbt->_firstInteriorTet = _vbt->_tetNodes.size();
	fillNonVnTetCenter();
	_vbt->_tetNodes.shrink_to_fit();
	_vbt->_tetCentroids.shrink_to_fit();
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_vbt->_tetCentroids.size());
	for (int n = _vbt->_tetCentroids.size(), i = 0; i < n; ++i)
		_vbt->_tetHash.insert(std::make_pair(_vbt->_tetCentroids[i], i));
	return true;
}

void vnBccTetCutter::createInteriorNodes() {
	short z;
	std::array<short, 3> s3;
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
			s3[0] = xi * 2 + 1;
			s3[1] = yi * 2 + 1;
			runInteriorNodes(oddXy[xi][yi], false);
		}
	}
	for (short xi = 0; xi < evenXy.size(); ++xi) {
		for (short yi = 0; yi < evenXy[xi].size(); ++yi) {
			s3[0] = (xi + 1) * 2;
			s3[1] = (yi + 1) * 2;
			runInteriorNodes(evenXy[xi][yi], true);
		}
	}
	evenXy.clear();
	oddXy.clear();
}

void vnBccTetCutter::assignExteriorTetNodes() {
	struct tetIndex {
		int tet;
		int nodeIndex;
	};
	struct tetTris {
		std::set<int> tris;
		std::list<tetIndex> tetIndices;
		int node;
		tetTriangles *seedTetPtr;
	};
	std::list<tetTris> triPools;
	auto tetsConnect = [&](const std::vector<int>& tris0, const std::set<int> &poolTris) ->bool{
		for (auto t : tris0) {
			if (poolTris.find(t) != poolTris.end())
				return true;
		}
		return false;
	};
	for (auto& tni : _exteriorNodeIndices) {
		triPools.clear();
		for (auto nit = tni.second.begin(); nit != tni.second.end(); ++nit) {
			bool makeNewPool = true;
			auto firstConnectedPool = triPools.end();
			auto pit = triPools.begin();
			while ( pit != triPools.end()) {
				if (tetsConnect(nit->ttPtr->tris, pit->tris)) {
					if (firstConnectedPool == triPools.end()) {
						firstConnectedPool = pit;
						tetIndex ti;
						ti.tet = nit->ttPtr->tetindx;
						ti.nodeIndex = nit->nodeIndex;
						firstConnectedPool->tetIndices.push_back(ti);
						firstConnectedPool->tris.insert(nit->ttPtr->tris.begin(), nit->ttPtr->tris.end());
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
				tt.tris.insert(nit->ttPtr->tris.begin(), nit->ttPtr->tris.end());
				tetIndex ti;
				ti.tet = nit->ttPtr->tetindx;
				ti.nodeIndex = nit->nodeIndex;
				tt.tetIndices.push_back(ti);
				tt.seedTetPtr = nit->ttPtr;
				triPools.push_back(tt);
			}
		}
		// any remaining triPools should get their own exterior, possibly virtual, node
		for (auto tp = triPools.begin(); tp != triPools.end(); ) {
			int eNode = _vbt->_nodeGridLoci.size();
			_vbt->_nodeGridLoci.push_back(tni.first);
			for (auto& ti : tp->tetIndices)
				_vbt->_tetNodes[ti.tet][ti.nodeIndex] = eNode;
			tp = triPools.erase(tp);
		}
	}
}

int vnBccTetCutter::nearestRayPatchHit(const Vec3f &rayBegin, Vec3f rayEnd, const std::vector<int>& tris, float& distanceSq) {  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
	distanceSq = FLT_MAX;
	rayEnd -= rayBegin;
	int closeT = -1;
	Vec3f N;
	for (auto& t : tris) {
		int* tr = _mt->triangleVertices(t);
		Vec3f tv[3];
		for (int i = 0; i < 3; ++i)
			tv[i] = _vMatCoords[tr[i]];
		tv[1] -= tv[0];
		tv[2] -= tv[0];
		float dSq;
		Vec3f P, R;
		Mat3x3f M;
		M.Initialize_With_Column_Vectors(tv[1], tv[2], -rayEnd);
		R = M.Robust_Solve_Linear_System(rayBegin - tv[0]);
		if (R[0] < 0.0f || R[0] > 1.0f || R[1] < 0.0f || R[1] > 1.0f || R[2] < 0.0f || R[2] > 1.0f || R[0] + R[1] > 1.0f)
			continue;
		P = rayEnd * R[2];
		dSq = P * P;
		if (distanceSq > dSq) {
			distanceSq = dSq;
			closeT = t;
			N = tv[1] ^ tv[2];
		}
	}
	if (distanceSq < FLT_MAX) {
		if (N * rayEnd > 0.0f)
			return -1;
		else
			return 1;
	}
	else
		return 0;
}

bool vnBccTetCutter::nearestPatchPoint(const short (&gl)[4][3], const int tetIdx, const std::vector<int>& tris, Vec3f& closeP, float& distanceSq) {
	// This routine assumes a patch intersection with a tet face.  It also handles possibility of closed patches entirely with a large tet.
	Vec3d closestP, P = {(double)gl[tetIdx][0], (double)gl[tetIdx][1], (double)gl[tetIdx][2]};
	double d, minD = DBL_MAX;
	int triEdge[2] = { -1, -1 };
	for (int i = 0; i < 4; ++i) {
		Vec3d tri[3];
		for (int j = 0; j < 3; ++j) {
			int nIdx = (tetIdx + 2 + j + i) & 3;
			tri[j].set((double)gl[nIdx][0], (double)gl[nIdx][1], (double)gl[nIdx][2]);  // each face of tet with last one opposite this tetIdx
		}
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
				if (d - minD < 1e-12) {
					if (d - minD > -1e-12)
						triEdge[1] = t;
					else {
						triEdge[0] = t;
						triEdge[1] = -1;
					}
					minD = d;
					closestP = target;
				}
			}
		}
		if (triEdge[0] > -1)
			break;
	}
	if(triEdge[0] < 0)  // this patch is a closed lacuna completely inside this tet with no tet face intersections
		return closestPatchPoint(Vec3f((float)gl[tetIdx][0], (float)gl[tetIdx][1], (float)gl[tetIdx][2]), tris, closeP, distanceSq);
	closeP.set((float)closestP[0], (float)closestP[1], (float)closestP[2]);
	distanceSq = (float)minD;
	auto triNorm = [&](const int triangle, Vec3f &N) {
		int* tr = _mt->triangleVertices(triangle);
		Vec3f U = _vMatCoords[tr[1]] - _vMatCoords[tr[0]], V = _vMatCoords[tr[2]] - _vMatCoords[tr[0]];
		N = U ^ V;
		N.q_normalize();
	};
	Vec3f N0;
	triNorm(triEdge[0], N0);  // intersection of edge of tet face with a single triangle if no triEdge[1]
	if (triEdge[1] > -1){  // face intersection with edge between two triangles
		Vec3f N1;
		triNorm(triEdge[1], N1);
		N0 += N1;
		if (fabs(N0.X) < 1e-5f && fabs(N0.Y) < 1e-5f && fabs(N0.Z) < 1e-5f)  // at bottom of coincident edge.  Guaranteed inside.
			return true;
	}
	return (Vec3f(P.xyz) - closeP) * N0 < 0.0f;
}

bool vnBccTetCutter::closestPatchPoint(const Vec3f& P, const std::vector<int>& tris, Vec3f &closeP, float &distanceSq) {
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

bool vnBccTetCutter::pointInsidePatchVertex(const Vec3f& P, const int triangle, const int tIndex) {
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

void vnBccTetCutter::getConnectedComponents(const bccTetCentroid& tc, std::list<tetTriangles>& ct) {  // for this centroid split its triangles into single solid connected components
	auto trVec = &ct.front().tris;  // on entry only one triangle list for this centroid
	ct.front().tetindx = -1;
	if (trVec->size() > 1) {
		std::set<int> trSet, ts;
		trSet.insert(trVec->begin(), trVec->end());
		trVec->clear();
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
				ct.push_back(tetTriangles());
				ct.back().tetindx = -1;
				trVec = &ct.back().tris;
			}
		}
	}
	struct patch {
		std::list<tetTriangles>::iterator pit;
		std::array<int, 4> tetNodes;
		Vec3f insideP;  // point just inside this patch
//		Vec3f N;  // normal pointing outside patch from insideP - COURT NUKE?
		float dsq;
	}p;
	std::list<patch> patches;
	p.tetNodes = { -1, -1, -1, -1 };
	auto pit = ct.begin();
	while (pit != ct.end()) {
		p.pit = pit;
		patches.push_back(p);
		++pit;
	}
	// find interior nodes for each component, joining patches as appropriate
	short gl[4][3];
	_vbt->centroidToNodeLoci(tc, gl);
	for (int i = 0; i < 4; ++i) {
		std::array<short, 3> loc = { gl[i][0], gl[i][1], gl[i][2] };
		auto iit = _interiorNodes.find(loc);
		Vec3f node(loc[0], loc[1], loc[2]);
		std::list<std::list<patch>::iterator> outside, inside;
		for (auto pit = patches.begin(); pit != patches.end(); ++pit) {
			bool connected = nearestPatchPoint(gl, i, pit->pit->tris, pit->insideP, pit->dsq);
			if (connected)
				inside.push_back(pit);
			else
				outside.push_back(pit);
		}
		if (inside.empty()) {  // No patch bound to an interior node.
//			if (iit != _interiorNodes.end())
//				throw(std::logic_error("An existing tet with an interior node has no patches attached to it in getConnectedComponents()\n"));
		}
		else if (!outside.empty()) {
			for (auto ip = inside.begin(); ip != inside.end(); ) {
				auto oMerge = outside.end();
				for (auto op = outside.begin(); op != outside.end(); ++op) {
					int connected = nearestRayPatchHit((*ip)->insideP, node, (*op)->pit->tris, (*op)->dsq);  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
					if (connected < 0 && (*op)->dsq > 1e-5f) {  // rule out possible coincident surfaces
						if (oMerge == outside.end() || (*op)->dsq < (*oMerge)->dsq) {
							oMerge = op;
						}
					}
				}
				if (oMerge != outside.end()) {
					(*oMerge)->pit->tris.insert((*oMerge)->pit->tris.end(), (*ip)->pit->tris.begin(), (*ip)->pit->tris.end());
					ct.erase((*ip)->pit);
					patches.erase(*ip);
					ip = inside.erase(ip);
				}
				else
					++ip;
			}
		}
		// any remaining inside patches should be merged along with any interior nodes
		if (inside.size() > 1) {
//			assert(iit != _interiorNodes.end());
			auto ip0 = inside.begin();
			auto ip = ip0;
			++ip;
			while (ip != inside.end()) {
				(*ip0)->pit->tris.insert((*ip0)->pit->tris.end(), (*ip)->pit->tris.begin(), (*ip)->pit->tris.end());
				ct.erase((*ip)->pit);
				patches.erase(*ip);
				ip = inside.erase(ip);
			}
		}
		if (iit != _interiorNodes.end()) {
			if (inside.empty())
				;
//				throw(std::logic_error("gCC() not working in cutter.\n"));
			else
				inside.front()->tetNodes[i] = iit->second;
		}
	}
	for (auto& cTet : patches) {
		cTet.pit->tetindx = _vbt->_tetNodes.size();
		_vbt->_tetNodes.push_back(cTet.tetNodes);
		_vbt->_tetCentroids.push_back(tc);  // careful later with multithreading
	}
	for (int i = 0; i < 4; ++i) {
		bool enNotEntered = true;
		std::pair< std::unordered_map<std::array<short, 3>, std::list<tetNodeIndex>, arrayShort3Hasher>::iterator, bool> pr;
		for (auto& p : patches) {
			if (p.tetNodes[i] < 0) {
				if (enNotEntered) {
					std::array<short, 3> loc = { gl[i][0], gl[i][1], gl[i][2] };
					pr = _exteriorNodeIndices.insert(std::make_pair(loc, std::list<tetNodeIndex>()));
					enNotEntered = false;
				}
				tetNodeIndex tni;
				tni.nodeIndex = i;
				tni.ttPtr = &(*p.pit);
				pr.first->second.push_back(tni);
			}
		}
	}
}

bool vnBccTetCutter::setupBccIntersectionStructures(int maximumGridDimension)
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
		_vbt->gridLocusToTetCentroid(_vMatCoords[i], _vertexTetCentroids[i]);
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

void vnBccTetCutter::inputTriangle(int tri) {  // main routine for finding triangle intersections with Z grid lines and with tetCentroid locations
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
	tetTriangles ttr;
	ttr.tetindx = -1;
	ttr.tris.clear();
	for (auto& tt : triTetsS) {
		auto pr = _centroidTriangles.insert(std::make_pair(tt, std::list<tetTriangles>()));
		if (pr.second)
			pr.first->second.push_back(ttr);
		pr.first->second.front().tris.push_back(tri);
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
			double z;
			if (odd) {
				if (triIntersectZ(i, j, z))
					oddXy[(i - 1) >> 1][(j - 1) >> 1].insert(std::make_pair(z, solidBegin));
			}
			else {
				if (triIntersectZ(i, j, z))
					evenXy[(i - 2) >> 1][(j - 2) >> 1].insert(std::make_pair(z, solidBegin));  // 0th in i and j always empty since object always inside box
			}
		}
	}
	return;
}

void vnBccTetCutter::fillNonVnTetCenter()
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
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (nextQuad[2] > -1 && nextQuad[3] > -1) {
					tet[2] = nextQuad[3];
					tet[3] = nextQuad[2];
					auto tc = center;
					++tc[1];
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
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
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
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
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
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
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[1] > -1 && prevQuad[3] > -1) {
					tet[2] = prevQuad[1];
					tet[3] = prevQuad[3];
					auto tc = center;
					--tc[2];
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
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
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
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
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
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
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[3] > -1 && nextQuad[3] > -1) {
					tet[2] = nextQuad[3];
					tet[3] = prevQuad[3];
					auto tc = center;
					++tc[0];
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
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
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
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
					if (_centroidTriangles.find(tc) == _centroidTriangles.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
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

