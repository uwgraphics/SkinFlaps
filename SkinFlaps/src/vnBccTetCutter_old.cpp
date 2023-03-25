#include <assert.h>
#include <limits>
#include <algorithm>
#include <cmath>
#include <functional>
#include "boundingBox.h"
#include "Mat2x2d.h"
#include "Mat3x3f.h"

#include <omp.h>
#include <chrono>  // for openMP timing
#include <ctime>  // nuke after openMP debug
#include <fstream>


#include "vnBccTetCutter.h"

vnBccTetCutter::vnBccTetCutter(void)
{
}

vnBccTetCutter::~vnBccTetCutter(void)
{
}

bool vnBccTetCutter::makeFirstVnTets(materialTriangles *mt, vnBccTetrahedra *vbt, int maximumGridDimension)
{  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	// WARNING - no complete tests are done to check for non-self-intersecting closed manifold triangulated surface input!!
	// This is essential. findAdjacentTriangles() is the closest test I provide.
	_mt = mt;
	_vbt = vbt;
	_vbt->_mt = mt;
	_vbt->_vertexTets.clear();
	_vbt->_vertexTets.assign(mt->numberOfVertices(), -1);
	_vbt->_barycentricWeights.clear();
	_vbt->_barycentricWeights.assign(mt->numberOfVertices(), Vec3f());
	_surfaceTetFaceNumber = 0;
	_vbt->_fixedNodes.clear();
	_vbt->_tetHash.clear();
	_vbt->_tetNodes.clear();
	if(!_mt->findAdjacentTriangles(true,false))	return false;
	if (!setupBccIntersectionStructures(maximumGridDimension))
		return false;
	for (int i = 0; i < 6; ++i){  // only 6 threads but less false sharing
		int n = _planeSets2[i].size();
		for (int j = 0; j < n; ++j)
			getPlanePolygons(i, j);
	}
	// first plane set gets all the interior tetrahedral nodes.
	_vbt->_nodeGridLoci.clear();
	_vbt->_nodeGridLoci.reserve((_planeSets2[0].size()*_vbt->_gridSize[1] * _vbt->_gridSize[2]) >> 2);  // reasonable guess
	for (int i = 0; i < 6; ++i){
		int j, n = _planeSets2[i].size();
		for (j = 0; j < n; ++j)
			processIntersectionPlane(i, j);
	}
	_interiorNodes.clear();  // only nodes created thus far are interior nodes
	_interiorNodes.reserve(_vbt->_nodeGridLoci.size());
	for (int n = _vbt->_nodeGridLoci.size(), j = 0; j < n; ++j)
		_interiorNodes.insert(std::make_pair(_vbt->_nodeGridLoci[j], j));
	_surfaceTetFaces.clear();
	_surfaceTetFaces.reserve(_surfaceTetFaceNumber << 1);
	collectSurfaceTetCentersFaces();
	createVirtualNodedSurfaceTets();
	_surfaceTetFaces.clear();
	createSurfaceTetNodes();
	int nSurfaceTets = _vbt->_tetNodes.size();
	fillNonVnTetCenter();
	// In some complex solids createVirtualNodedSurfaceTets() will create multiple tets for the same tet locus as they will have no shared connected components,
	// but they can have the same nodes since the tets surrounding them may have connected components.  While these are not strictly an error, they do the user no good
	// since the tets are incapable of independent movement.  Further these are tets that are sparsely filled with solid.  Combining them provides somewhat more accurate
	// physics behavior. The next section removes these identical multiplicities.
	int k = 0;
	long *nk = _vbt->_tetNodes[0].data();
	std::vector<long> tetDispl;
	tetDispl.assign(nSurfaceTets, -1);
	tetDispl[0] = 0;
	bool dupFound = false;
	for (int i = 1; i < nSurfaceTets; ++i){
		tetDispl[i] = k + 1;
		long *ni = _vbt->_tetNodes[i].data();
		if (ni[0] == nk[0] && ni[1] == nk[1] && ni[2] == nk[2] && ni[3] == nk[3]){  // due to createVirtualNodedSurfaceTets() these identicals will all be sequential
			--tetDispl[i];
			dupFound = true;
			continue;
		}
		++k;
		nk = _vbt->_tetNodes[k].data();;
		if (dupFound){
			_vbt->_tetNodes[k] = _vbt->_tetNodes[i];
			_vbt->_tetCentroids[k] = _vbt->_tetCentroids[i];
		}
	}
	if (dupFound){
		++k;
		_vbt->_tetNodes.erase(_vbt->_tetNodes.begin() + k, _vbt->_tetNodes.begin() + nSurfaceTets);
		_vbt->_tetCentroids.erase(_vbt->_tetCentroids.begin() + k, _vbt->_tetCentroids.begin() + nSurfaceTets);
		for (int n = _vbt->_vertexTets.size(), i = 0; i < n; ++i){
			if (_vbt->_vertexTets[i] > -1)
				_vbt->_vertexTets[i] = tetDispl[_vbt->_vertexTets[i]];
		}
	}
	_vbt->_tetNodes.shrink_to_fit();
	_vbt->_tetCentroids.shrink_to_fit();
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_vbt->_tetCentroids.size());
	for (int n = _vbt->_tetCentroids.size(), i = 0; i < n; ++i)
		_vbt->_tetHash.insert(std::make_pair(_vbt->_tetCentroids[i].ll, i));
	// for debug - nuke later
 /* #ifdef _DEBUG
	for (auto &th : _vbt->_tetHash){
		Vec3f tV[4];
		std::array<short, 3> mean, tmp;
		mean.assign(0);
		bccTetCentroid tc;
		tc.ll = th.first;
		for (int i = 0; i < 4; ++i){
			tmp = _vbt->_nodeGridLoci[_vbt->_tetNodes[th.second][i]];
			if (i < 1){
				if ((tc.xyz[tc.halfCoordAxis] + tc.xyz[(tc.halfCoordAxis + 1) % 3]) & 1){
					assert(tmp[tc.halfCoordAxis] == tc.xyz[tc.halfCoordAxis] && tmp[(tc.halfCoordAxis + 2) % 3] == tc.xyz[(tc.halfCoordAxis + 2) % 3] && tmp[(tc.halfCoordAxis + 1) % 3] + 1 == tc.xyz[(tc.halfCoordAxis + 1) % 3]);
				}
				else{
					assert(tmp[tc.halfCoordAxis] - 1 == tc.xyz[tc.halfCoordAxis] && tmp[(tc.halfCoordAxis + 2) % 3] == tc.xyz[(tc.halfCoordAxis + 2) % 3] && tmp[(tc.halfCoordAxis + 1) % 3] + 1 == tc.xyz[(tc.halfCoordAxis + 1) % 3]);
				}
			}
			else if (i < 2){
				if ((tc.xyz[tc.halfCoordAxis] + tc.xyz[(tc.halfCoordAxis + 1) % 3]) & 1){
					assert(tmp[tc.halfCoordAxis] == tc.xyz[tc.halfCoordAxis] && tmp[(tc.halfCoordAxis + 2) % 3] == tc.xyz[(tc.halfCoordAxis + 2) % 3] && tmp[(tc.halfCoordAxis + 1) % 3] - 1 == tc.xyz[(tc.halfCoordAxis + 1) % 3]);
				}
				else{
					assert(tmp[tc.halfCoordAxis] - 1 == tc.xyz[tc.halfCoordAxis] && tmp[(tc.halfCoordAxis + 2) % 3] == tc.xyz[(tc.halfCoordAxis + 2) % 3] && tmp[(tc.halfCoordAxis + 1) % 3] - 1 == tc.xyz[(tc.halfCoordAxis + 1) % 3]);
				}
			}
			else if (i < 3){
				if ((tc.xyz[tc.halfCoordAxis] + tc.xyz[(tc.halfCoordAxis + 1) % 3]) & 1){
					assert(tmp[tc.halfCoordAxis] - 1 == tc.xyz[tc.halfCoordAxis] && tmp[(tc.halfCoordAxis + 2) % 3] - 1 == tc.xyz[(tc.halfCoordAxis + 2) % 3] && tmp[(tc.halfCoordAxis + 1) % 3] == tc.xyz[(tc.halfCoordAxis + 1) % 3]);
				}
				else{
					assert(tmp[tc.halfCoordAxis] == tc.xyz[tc.halfCoordAxis] && tmp[(tc.halfCoordAxis + 2) % 3] + 1 == tc.xyz[(tc.halfCoordAxis + 2) % 3] && tmp[(tc.halfCoordAxis + 1) % 3] == tc.xyz[(tc.halfCoordAxis + 1) % 3]);
				}
			}
			else{
				if ((tc.xyz[tc.halfCoordAxis] + tc.xyz[(tc.halfCoordAxis + 1) % 3]) & 1){
					assert(tmp[tc.halfCoordAxis] - 1 == tc.xyz[tc.halfCoordAxis] && tmp[(tc.halfCoordAxis + 2) % 3] + 1 == tc.xyz[(tc.halfCoordAxis + 2) % 3] && tmp[(tc.halfCoordAxis + 1) % 3] == tc.xyz[(tc.halfCoordAxis + 1) % 3]);
				}
				else{
					assert(tmp[tc.halfCoordAxis] == tc.xyz[tc.halfCoordAxis] && tmp[(tc.halfCoordAxis + 2) % 3] - 1 == tc.xyz[(tc.halfCoordAxis + 2) % 3] && tmp[(tc.halfCoordAxis + 1) % 3] == tc.xyz[(tc.halfCoordAxis + 1) % 3]);
				}
			}
			tV[i].set(tmp[0], tmp[1], tmp[2]);
			for (int j = 0; j < 3; ++j)
				mean[j] += tmp[j];
		}
		for (int i = 0; i < 3; ++i){
			if (mean[i] & 2)
				assert(i == tc.halfCoordAxis);
			mean[i] >>= 2;
			tV[i + 1] -= tV[0];
		}
		assert(mean == tc.xyz);
		// test for positiveness in tet ordering
		float vol = tV[1] * (tV[2] ^ tV[3]);
		assert(vol > 0.0f);
	}
#endif */

	// all vertices in duplicated, virtual noded tets have their tet index already assigned.  Vertices in unique tets still unassigned.
	for (int n = _vbt->_vertexTets.size(), i = 0; i < n; ++i){
		if (_vbt->_vertexTets[i] < 0 && *_mt->vertexFaceTriangle(i) < 0x80000000) {  // second test ignores excised vertices
			auto rng = _vbt->_tetHash.equal_range(_vertexTetLoci[i].ll);
			assert(std::distance(rng.first, rng.second) == 1);
			_vbt->_vertexTets[i] = rng.first->second;
		}
	}

	// COURT - no longer doing this here
/*	_vbt->_nodeSpatialCoords.clear();
	_vbt->_nodeSpatialCoords.assign(_vbt->_nodeGridLoci.size(), Vec3f());
	for (int n = _vbt->_nodeGridLoci.size(), i = 0; i < n; ++i){
		const short *np = _vbt->_nodeGridLoci[i].data();
		Vec3f *vp = &_vbt->_nodeSpatialCoords[i];
		vp->set((float)np[0], (float)np[1], (float)np[2]);
		*vp *= (float)_vbt->_unitSpacing;
		*vp += _vbt->_minCorner;
	} */

	return true;
}

void vnBccTetCutter::createSurfaceTetNodes()
{
	struct sharedTetNode{
		bool internal;
		std::set<unsigned long> nSet;
	}stn;
	stn.internal = false;
	stn.nSet.clear();
	std::unordered_map<std::array<short, 3>, std::list<sharedTetNode>, arrayShort3Hasher> surfaceTetNodes;
	surfaceTetNodes.reserve(_surfaceTetFaceNumber << 2);
	std::array<short, 3> nodeLoc;
	auto addTetNodePair = [&](int triVert, vnTetFace *tfp){
		unsigned long v0 = tfp->tetNodes[0][triVert], v1 = tfp->tetNodes[1][triVert];
		assert(v0 < 0xffffffff && v1 < 0xffffffff);  // Should never happen.
		bool internalNode;
		if (tfp->interiorNodes & (1 << triVert))
			internalNode = true;
		else
			internalNode = false;
		auto tnSets = &surfaceTetNodes.insert(std::make_pair(nodeLoc, std::list<sharedTetNode>())).first->second;
		auto tsit = tnSets->begin();
		while (tsit != tnSets->end()){
			if (tsit->nSet.find(v0) != tsit->nSet.end()){
				tsit->internal |= internalNode;  // if any face says it is an internal node, they all are
				if (tsit->nSet.insert(v1).second){
					auto tsit2 = tnSets->begin();
					while (tsit2 != tnSets->end()){
						if (tsit == tsit2){
							++tsit2;
							continue;
						}
						assert(tsit2->nSet.find(v0) == tsit2->nSet.end());
						if (tsit2->nSet.find(v1) != tsit2->nSet.end()){
							tsit->internal |= tsit2->internal;
							tsit->nSet.insert(tsit2->nSet.begin(), tsit2->nSet.end());
							tnSets->erase(tsit2);
							tsit2 = tnSets->end();
							break;
						}
						++tsit2;
					}
				}
				break;
			}
			if (tsit->nSet.find(v1) != tsit->nSet.end()){
				tsit->internal |= internalNode;  // if any face says it is an internal node, they all are
				tsit->nSet.insert(v0);  // must happen
				auto tsit2 = tnSets->begin();
				while (tsit2 != tnSets->end()){
					if (tsit == tsit2){
						++tsit2;
						continue;
					}
					assert(tsit2->nSet.find(v1) == tsit2->nSet.end());
					if (tsit2->nSet.find(v0) != tsit2->nSet.end()){
						tsit->internal |= tsit2->internal;
						tsit->nSet.insert(tsit2->nSet.begin(), tsit2->nSet.end());
						tnSets->erase(tsit2);
						tsit2 = tnSets->end();
						break;
					}
					++tsit2;
				}
				break;
			}
			++tsit;
		}
		if (tsit == tnSets->end()){
			tnSets->push_back(stn);
			tnSets->back().internal = internalNode;
			tnSets->back().nSet.insert(v0);
			tnSets->back().nSet.insert(v1);
		}
	};
	for (int i = 0; i < 6; ++i){
		for (int n = _planeSets2[i].size(), j = 0; j < n; ++j){
			int c0 = i >> 1, oddPlane = i&1;
			for (auto &tf : _planeSets2[i][j]->vnTetFaces){
				unsigned short faceOdd = (tf.first.first + tf.first.second) & 1;
				for (int k = 0; k < 3; ++k){
					nodeLoc[c0] = oddPlane ? _planeSets2[i][j]->D + tf.first.second + faceOdd : _planeSets2[i][j]->D - tf.first.second - faceOdd;
					nodeLoc[(c0 + 1) % 3] = tf.first.second + faceOdd;
					nodeLoc[(c0 + 2) % 3] = tf.first.first;
					if (k == 1)
						nodeLoc[(c0 + 2) % 3] += 2;
					if (k == 2){
						++nodeLoc[(c0 + 2) % 3];
						int vOffset = faceOdd ? -1 : 1;
						nodeLoc[(c0 + 1) % 3] += vOffset;
						nodeLoc[c0] += oddPlane ? vOffset : -vOffset;
					}
					addTetNodePair(k, &tf.second);
				}
			}
		}
	}
	// create and assign all shared surface tet nodes just found
	std::array<long, 4> emptyTet;
	emptyTet.fill(-1);
	_vbt->_tetNodes.assign(_surfaceTetNumber, emptyTet);
	for (auto &stn : surfaceTetNodes){
		for (auto &sharedNode : stn.second){
			long node;
			if (sharedNode.internal){
				auto in = _interiorNodes.find(stn.first);
				assert(in != _interiorNodes.end());
				node = in->second;
			}
			else{
				node = _vbt->_nodeGridLoci.size();
				_vbt->_nodeGridLoci.push_back(stn.first);
			}
			for (auto &tetNode : sharedNode.nSet){

				if ((tetNode >> 2) > 199 && (tetNode >> 2) < 203)
					int junk = 0;

				_vbt->_tetNodes[tetNode >> 2][tetNode & 3] = node;
			}
		}
	}
	_vbt->_tetHash.reserve(_vbt->_tetNodes.size());  // COURT - fix by reserving expected # of interior cubes to avoid rehashing
	for (int n = _vbt->_tetNodes.size(), i = 0; i < n; ++i){
		_vbt->_tetHash.insert(std::make_pair(_vbt->_tetCentroids[i].ll, i));
		for (int j = 0; j < 4; ++j){
			long *tn = &_vbt->_tetNodes[i][j];
			// single interior polygon face tet can have 3 faces unpenetrated leaving one vertex isolated.
			if (*tn < 0){
				std::array<short, 3> locus = _vbt->nodeGridLocation(_vbt->_tetCentroids[i], j);
				*tn = _vbt->_nodeGridLoci.size();
				_vbt->_nodeGridLoci.push_back(locus);
			}
		}
	}
}

bool vnBccTetCutter::getPlanePolygons(const int planeSet, const int plane)
{
	bccPlane *pp = _planeSets2[planeSet][plane].get();
	std::vector<double> planeDist;
	planeDist.reserve(_mt->numberOfVertices());
	int c2, c1, c0 = planeSet >> 1;
	c1 = (c0 + 1) % 3;
	c2 = (c0 + 2) % 3;
	bool odd = (planeSet & 1), somePos = false, someNeg = false;
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i){
		if (odd)
			planeDist.push_back(-_vMatCoords[i]._v[c0] + _vMatCoords[i]._v[c1] + pp->D);
		else
			planeDist.push_back(-_vMatCoords[i]._v[c0] - _vMatCoords[i]._v[c1] + pp->D);
		std::signbit(planeDist.back()) ? somePos = true : someNeg = true;
	}
	if (!(somePos & someNeg))
		return false;
	long *tr;
	auto getIntersect = [&](int edge, Vec2d &intrsct){
		double denom = abs(planeDist[tr[edge]] - planeDist[tr[(edge + 1) % 3]]);
		if (denom < 1e-16f){
			intrsct.X = _vMatCoords[tr[edge]][c2];
			intrsct.Y = _vMatCoords[tr[edge]][c1];
		}
		else{
			denom = 1.0 / denom;
			double *dp = _vMatCoords[tr[(edge + 1) % 3]]._v;
			Vec2d I1(dp[c2], dp[c1]);
			dp = _vMatCoords[tr[edge]]._v;
			intrsct.set(dp[c2], dp[c1]);
			intrsct *= abs(planeDist[tr[(edge + 1) % 3]]) * denom;
			intrsct += I1 * abs(planeDist[tr[edge]]) * denom;
		}
	};
	std::vector<char> trisDone;
	trisDone.assign(_mt->numberOfTriangles(),0);
	for (int n = _mt->numberOfTriangles(), j, i = 0; i < n; ++i){
		if (trisDone[i] || _mt->triangleMaterial(i) < 0)
			continue;
		tr = _mt->triangleVertices(i);
		for (j = 0; j < 2; ++j){
			if (std::signbit(planeDist[tr[j]]) != std::signbit(planeDist[tr[(j + 1) % 3]]))  // polygon start
				break;
		}
		if (j > 1)
			continue;
		// list polygon counterclockwise in the plane as the triangle edge polygons will be listed that way as well
		// counterclockwise going edge of triangle will always have signbit going from false to true
		// triSegment2 elements in polygon list triangle and its EXITING uv.  Entering uv is not listed in a ts2.
		if (j > 0){
			if (!std::signbit(planeDist[tr[1]]))
				j = 2;
		}
		else{
			if (!std::signbit(planeDist[tr[0]])){  // not edge 0. Will be remaining edge crossing
				if (std::signbit(planeDist[tr[2]]) != std::signbit(planeDist[tr[0]]))
					j = 2;
				else
					j = 1;
			}
		}
		// this polygon starts at triangle i edge j
		pp->polygons2.push_back(std::list<triSegment2>());
		triSegment2 ts;
		ts.triangle = i;
		do{
			getIntersect(j, ts.uv);
			pp->polygons2.back().push_back(ts);
			unsigned long adj = _mt->triAdjs(ts.triangle)[j];
			if (adj == 3)
				return false;
			ts.triangle = adj >> 2;
			tr = _mt->triangleVertices(ts.triangle);
			trisDone[ts.triangle] = 1;
			j = adj & 3;
			for (int k = 1; k < 3; ++k){
				if (std::signbit(planeDist[tr[(j + k) % 3]]) != std::signbit(planeDist[tr[(j + k + 1) % 3]])){  // next edge
					j = (j + k) % 3;
					break;
				}
			}
		} while (ts.triangle != i);
	}
	return true;
}

bool vnBccTetCutter::setupBccIntersectionStructures(int maximumGridDimension)
{
	boundingBox<float> bbf;
	bbf.Empty_Box();
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i)
		bbf.Enlarge_To_Include_Point((const float(&)[3])(*_mt->vertexCoordinate(i)));
	bbf.Minimum_Corner(_vbt->_minCorner._v);
	bbf.Maximum_Corner(_vbt->_maxCorner._v);
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
	float offset = FLT_EPSILON * (float)_vbt->_unitSpacing;
	_vbt->_minCorner -= Vec3f(offset, offset, offset);
	_vbt->_maxCorner += Vec3f(offset, offset, offset);
	double cs = _vbt->_maxCorner._v[bigDim] - _vbt->_minCorner._v[bigDim];
	_vbt->_unitSpacing = cs / maximumGridDimension;
	_vbt->_unitSpacingInv = 1.0 / _vbt->_unitSpacing;
	{
		Vec3f box = (_vbt->_maxCorner - _vbt->_minCorner)*(float)_vbt->_unitSpacingInv * 0.5f;
		for (int i = 0; i < 3; ++i){
			_vbt->_gridSize[i] = 1 + (int)std::floor(box._v[i]);
		}
	}
	_vMatCoords.clear();
	_vMatCoords.assign(_mt->numberOfVertices(), Vec3d());
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i){
		Vec3d Vf;
		Vf.set((const float(&)[3])*_mt->vertexCoordinate(i));
		_vMatCoords[i].set((Vf -_vbt->_minCorner)* _vbt->_unitSpacingInv);
	}
	for (int j, i = 0; i < 6; ++i)	{
		// while there will be physics nodes at the boundaries of the grid, there will be no intersections there as the object is inside the boundary
		int nf = _vbt->_gridSize[i >> 1], ns = _vbt->_gridSize[((i >> 1) + 1) % 3];
		_planeSets2[i].reserve(nf + ns - 1);
		for (j = 0; j < nf + ns - 1; ++j)
			_planeSets2[i].push_back(std::unique_ptr<bccPlane>(new bccPlane));
		if (i & 1){
			for (j = 0; j < ns - 1; ++j)
				_planeSets2[i][j]->D = -((ns - j - 1) << 1);
			for (j = 0; j < nf; ++j)
				_planeSets2[i][ns - 1 + j]->D = j << 1;
		}
		else{
			for (j = 0; j < nf - 1; ++j)
				_planeSets2[i][j]->D = (j + 1) << 1;
			for (j = 0; j < ns; ++j)
				_planeSets2[i][j + nf - 1]->D = (j + nf) << 1;
		}
	}
	_vertexTetLoci.clear();
	_vertexTetLoci.assign(_mt->numberOfVertices(), bccTetCentroid());
	_vbt->_barycentricWeights.clear();
	_vbt->_barycentricWeights.assign(_mt->numberOfVertices(), Vec3f());
	int n = _mt->numberOfVertices();
	for (int i = 0; i < n; ++i){
		Vec3f B;
		B.set(_vMatCoords[i]._v);
		_vbt->gridLocusToTetCentroid(B, _vertexTetLoci[i]);
		// set barycentric coordinate within that tet
		_vbt->gridLocusToBarycentricWeight(B, _vertexTetLoci[i], _vbt->_barycentricWeights[i]);
	}
	return true;
}

bool vnBccTetCutter::remakeVnTets(materialTriangles *mt)
{  // uses old vbt to get material coords of mt vertices and grid data, from which it makes new vbt.
	int newVertexStart = _vMatCoords.size();
	_vMatCoords.resize(_mt->numberOfVertices());
#ifdef _DEBUG
	for (int i = 0; i < newVertexStart; ++i){
		if (_vbt->getVertexTetrahedron(i) < 0)
			continue;
		Vec3f Vf;
		_vbt->vertexGridLocus(i, Vf);
		assert((Vf - Vec3f(_vMatCoords[i]._v)).length2() < 1e-10f);  // old vertices should retain their gridLocus
	}
#endif
	// all oldVertices should have an unchanged grid locus
	for (int n = _mt->numberOfVertices(), i = newVertexStart; i < n; ++i){  // after incisions all new vertices should have been assigned a grid locus associated with the old _vbt
		if (_vbt->getVertexTetrahedron(i) < 0)
			continue;
		Vec3f Vf;
		_vbt->vertexGridLocus(i, Vf);
		_vMatCoords[i].set(Vf);
	}
	_surfaceTetFaceNumber = 0;
	_vbt->_fixedNodes.clear();
	_vbt->_tetHash.clear();
	_vbt->_tetNodes.clear();
	if (!_mt->findAdjacentTriangles(true, false))	return false;
	for (int i = 0; i < 6; ++i)	{
		// while there will be physics nodes at the boundaries of the grid, there will be no intersections there as the object is inside the boundary
		for (int n = _planeSets2[i].size(), j = 0; j < n; ++j){
			_planeSets2[i][j]->polygons2.clear();
			_planeSets2[i][j]->vnTetFaces.clear();
		}
	}
	_vertexTetLoci.clear();
	bccTetCentroid emptyTC;
	emptyTC.ll = LLONG_MAX;
	_vertexTetLoci.assign(_mt->numberOfVertices(), emptyTC);
	_vbt->_barycentricWeights.clear();
	_vbt->_barycentricWeights.assign(_mt->numberOfVertices(), Vec3f());
	// set tet centroids and barycentric coordinates within that tet
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i){
		if (_vbt->_vertexTets[i] < 0)
			continue;
		Vec3f Vf;
		Vf.set((double(&)[3])_vMatCoords[i]._v);
		_vbt->gridLocusToTetCentroid(Vf, _vertexTetLoci[i]);
		_vbt->gridLocusToBarycentricWeight(Vf, _vertexTetLoci[i], _vbt->_barycentricWeights[i]);
	}
	_vbt->_vertexTets.clear();
	_vbt->_vertexTets.assign(mt->numberOfVertices(), -1);
	for (int i = 0; i < 6; ++i){
		for (int n = _planeSets2[i].size(), j = 0; j < n; ++j)
			getPlanePolygons(i, j);
	}
	// first plane set gets all the interior tetrahedral nodes.  To avoid need for mutex, don't multithread set 0.
	_vbt->_nodeGridLoci.clear();
	_vbt->_nodeGridLoci.reserve((_planeSets2[0].size()*_vbt->_gridSize[1] * _vbt->_gridSize[2]) >> 2);  // reasonable guess
	for (int i = 0; i < 6; ++i){
		int n = _planeSets2[i].size();
		for (int j = 0; j < n; ++j)
			processIntersectionPlane(i, j);
	}
	_interiorNodes.clear();
	_interiorNodes.reserve(_vbt->_nodeGridLoci.size());
	for (int n = _vbt->_nodeGridLoci.size(), j = 0; j < n; ++j)
		_interiorNodes.insert(std::make_pair(_vbt->_nodeGridLoci[j], j));
	_surfaceTetFaces.clear();
	_surfaceTetFaces.reserve(_surfaceTetFaceNumber << 1);
	collectSurfaceTetCentersFaces();
	createVirtualNodedSurfaceTets();
	_surfaceTetFaces.clear();
	createSurfaceTetNodes();
	int nSurfaceTets = _vbt->_tetNodes.size();
	fillNonVnTetCenter();
	// In some complex solids createVirtualNodedSurfaceTets() will create multiple tets for the same tet locus as they will have no shared connected components,
	// but they can have the same nodes since the tets surrounding them may have connected components.  While these are not strictly an error, they do the user no good
	// since the tets are incapable of independent movement.  Further these are tets that are sparsely filled with solid.  Combining them provides somewhat more accurate
	// physics behavior. The next section removes these identical multiplicities.
	int k = 0;
	long *nk = _vbt->_tetNodes[0].data();
	std::vector<long> tetDispl;
	tetDispl.assign(nSurfaceTets, -1);
	tetDispl[0] = 0;
	bool dupFound = false;
	for (int i = 1; i < nSurfaceTets; ++i){
		tetDispl[i] = k + 1;
		long *ni = _vbt->_tetNodes[i].data();
		if (ni[0] == nk[0] && ni[1] == nk[1] && ni[2] == nk[2] && ni[3] == nk[3]){  // due to createVirtualNodedSurfaceTets() these identicals will all be sequential
			--tetDispl[i];
			dupFound = true;
			continue;
		}
		++k;
		nk = _vbt->_tetNodes[k].data();;
		if (dupFound){
			_vbt->_tetNodes[k] = _vbt->_tetNodes[i];
			_vbt->_tetCentroids[k] = _vbt->_tetCentroids[i];
		}
	}
	if (dupFound){
		++k;
		_vbt->_tetNodes.erase(_vbt->_tetNodes.begin() + k, _vbt->_tetNodes.begin() + nSurfaceTets);
		_vbt->_tetCentroids.erase(_vbt->_tetCentroids.begin() + k, _vbt->_tetCentroids.begin() + nSurfaceTets);
		for (int n = _vbt->_vertexTets.size(), i = 0; i < n; ++i){
			if (_vbt->_vertexTets[i] > -1)
				_vbt->_vertexTets[i] = tetDispl[_vbt->_vertexTets[i]];
		}
	}
	_vbt->_tetNodes.shrink_to_fit();
	_vbt->_tetCentroids.shrink_to_fit();
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_vbt->_tetCentroids.size());
	for (int n = _vbt->_tetCentroids.size(), i = 0; i < n; ++i)
		_vbt->_tetHash.insert(std::make_pair(_vbt->_tetCentroids[i].ll, i));
	// all vertices in duplicated, virtual noded tets have their tet index already assigned.  Vertices in unique tets still unassigned.
	for (int n = _vbt->_vertexTets.size(), i = 0; i < n; ++i){
		if (_vbt->_vertexTets[i] < 0){
			if (_vertexTetLoci[i].ll < LLONG_MAX) {
				auto rng = _vbt->_tetHash.equal_range(_vertexTetLoci[i].ll);
				assert(std::distance(rng.first, rng.second) == 1);
				_vbt->_vertexTets[i] = rng.first->second;
			}
			else {
				assert(*_mt->vertexFaceTriangle(i) == 0x80000000);  // a deleted vertex from an excision
			}
		}
	}

	// COURT no longer doing this here
/*	_vbt->_nodeSpatialCoords.clear();
	_vbt->_nodeSpatialCoords.assign(_vbt->_nodeGridLoci.size(), Vec3f());
	for (int n = _vbt->_nodeGridLoci.size(), i = 0; i < n; ++i){
		const short *np = _vbt->_nodeGridLoci[i].data();
		Vec3f *vp = &_vbt->_nodeSpatialCoords[i];
		vp->set((float)np[0], (float)np[1], (float)np[2]);
		*vp *= (float)_vbt->_unitSpacing;
		*vp += _vbt->_minCorner;
	} */

	return true;
}

/* bool vnBccTetCutter::getTriangleAdjacencies(int nTris, const int (*triangles)[3])
{
	assert(false);
	typedef std::map<longPair,int,longPairTest> edgeMap;
	typedef edgeMap::iterator edgeIt;
	std::pair <edgeIt,bool> P;
	edgeIt ei;
	edgeMap M;
	M.clear();
	_adjTris.clear();
	_adjTris.assign(nTris*3,0x00000003);
	long maxVert=0;
	for(int j,oldCode,newCode,i=0; i<nTris; ++i)	{
		for(j=0; j<3; j++)	{
			longPair E(triangles[i][j],triangles[i][(j+1)%3]);
			if(E.lMax>maxVert)	maxVert=E.lMax;
			newCode = (i<<2)+j;
			P = M.insert(std::make_pair(E,newCode));
			if(P.second==false)	// edge match found
			{
				oldCode = P.first->second;
				_adjTris[i*3+j] = oldCode;
				_adjTris[(oldCode>>2)*3+(oldCode&0x00000003)] = newCode;
				M.erase(P.first);
			}
		}
	}
	_vertexTris.assign(++maxVert,-1);
	for(int j,i=0; i<nTris; ++i)	{
		for(j=0; j<3; ++j)
			_vertexTris[triangles[i][j]]=i;	}
	if(M.empty()) return true;
	return false;
} */

/*void vnBccTetCutter::getTriangleVertexNeighbors(int triVert, std::vector<int> &triangles, std::vector<int> &vertices)
{	// this one only works for a closed manifold surface.  Use other version if not closed.
	assert(false); 
	triangles.clear();	vertices.clear();
	int i,tr=-1,tStart=_vertexTris[triVert];
	for(i=0; i<3; ++i)
		if(tri[tStart][i]==triVert)	break;
	unsigned int adj=_adjTris[tStart*3+i];
	while(tr!=tStart)	{
		i = adj&0x00000003;
		tr = adj>>2;
		triangles.push_back(tr);
		vertices.push_back(tri[tr][i]);
		adj=_adjTris[tr*3+((i+1)%3)];	}
} */


void vnBccTetCutter::processIntersectionPlane(const int planeSet, const int plane)
{
	std::vector<PLANE_LINE> planeHorizLines, planeDiagonals[2];
	int first = planeSet >> 1;
	int second = (first + 1) % 3;
	int third = (first + 2) % 3;
	int diagDim = _vbt->_gridSize[second] + _vbt->_gridSize[third] + 1;
	planeDiagonals[0].resize(diagDim);
	planeDiagonals[1].resize(diagDim);
	int phlDim = (_vbt->_gridSize[second] << 1) + 1;
	planeHorizLines.resize(phlDim);
	_planeSets2[planeSet][plane]->vnTetFaces.clear();
	std::list<std::list<triSegment2> > *L = &_planeSets2[planeSet][plane]->polygons2;
	if (L->empty())
		return;
	std::multimap<std::pair<short, short>, std::list<triSegment2> > planeFaces;
	for (auto lit = L->begin(); lit != L->end(); ++lit)	{
		Vec2d lastUv = lit->back().uv;
		std::pair<short, short> lp, lpLast = uvToPlaneFace(lastUv);
		std::list<triSegment2> segmentLast, segmentNow;
		std::list<std::pair<short, short> > triFaces;
		bool segmentLastDone = false;
		for (auto pit = lit->begin(); pit != lit->end(); ++pit){  // pit contains triangle and its exiting uv.  Entering uv not listed.
			if (segmentLastDone)
				segmentNow.push_back(*pit);
			else
				segmentLast.push_back(*pit);
			lp = uvToPlaneFace(pit->uv);
			if (lp != lpLast){
				if (segmentLastDone){
					segmentNow.back().uv[0] = DBL_MAX;  // last uv not inside previous face.  Mark as an unclosed polygon segment
					planeFaces.insert(std::make_pair(lpLast, segmentNow));
				}
				else{
					segmentLast.back().uv[0] = DBL_MAX;  // same logic as above
					segmentLastDone = true;
				}
				// next routine collects face path of this polygon segment and logs any triangle edge intersections
				getFaceEdgePath(&lastUv, lpLast, &pit->uv, lp, pit->triangle, second, triFaces, planeHorizLines, planeDiagonals);
				segmentNow.clear();
				// this uv is in the new face
				segmentNow.push_back(*pit);
				if (!triFaces.empty()){  // complete triangle pass through these faces with no uv inside
					double firstU = pit->uv[0];
					segmentNow.back().uv[0] = DBL_MAX;
					for (auto &tf : triFaces)
						planeFaces.insert(std::make_pair(tf, segmentNow));
					segmentNow.back().uv[0] = firstU;
				}
				lpLast = lp;
			}
			lastUv = pit->uv;
		}
		if (segmentLastDone){
			segmentNow.splice(segmentNow.end(), segmentLast);
			planeFaces.insert(std::make_pair(lpLast, segmentNow));
		}
		else  // closed poplygon inside one triangular face
			planeFaces.insert(std::make_pair(lpLast, segmentLast));
	}
	L->clear();
	PLANE_LINE::iterator hit;
	for (int k, j = 0; j<diagDim; ++j)	{
		if (!planeDiagonals[0][j].empty())	{
			k = 1;
			for (hit = planeDiagonals[0][j].begin(); hit != planeDiagonals[0][j].end(); ++hit)	{
				if (hit->second.solidRight != k){  // set from polygon in getFaceEdgePath()
					// may be from floating point roundoff of coincident surface
					planeLineCrossing plc = hit->second;
					auto hit2 = hit;
					++hit2;
					assert(hit2 != planeDiagonals[0][j].end() && hit2->first - hit->first < 1e-6);
					hit->second = hit2->second;
					hit2->second = plc;
				}
				k ^= 1;
			}
			assert(k>0);
		}
		if (planeDiagonals[1][j].empty())	continue;
		k = 1;
		for (hit = planeDiagonals[1][j].begin(); hit != planeDiagonals[1][j].end(); ++hit)	{
			if (hit->second.solidRight != k){  // set from polygon in getFaceEdgePath()
				// may be from floating point roundoff of coincident surface
				planeLineCrossing plc = hit->second;
				auto hit2 = hit;
				++hit2;
				assert(hit2 != planeDiagonals[1][j].end() && hit2->first - hit->first < 1e-6);
				hit->second = hit2->second;
				hit2->second = plc;
			}
			k ^= 1;
		}
		assert(k>0);
	}
	for (int k, j = 0; j<phlDim; ++j)	{
		if (planeHorizLines[j].empty())	continue;
		k = 1;
		for (hit = planeHorizLines[j].begin(); hit != planeHorizLines[j].end(); ++hit)	{
			if (hit->second.solidRight != k){  // set from polygon in getFaceEdgePath()
				// may be from floating point roundoff of coincident surface
				planeLineCrossing plc = hit->second;
				auto hit2 = hit;
				++hit2;
				assert(hit2 != planeHorizLines[j].end() && hit2->first - hit->first < 1e-6);
				hit->second = hit2->second;
				hit2->second = plc;
			}
			k ^= 1;
		}
		assert(k>0);
	}
	if (planeSet < 1)	{	// save interior tet nodes.  Careful here with openMP. Would need to create mutex for interiorNodeLatticeLoci. ?not worth the trouble.
		assert (planeHorizLines[0].empty());
		assert(planeHorizLines[phlDim - 1].empty());
		for (int j = 1; j < phlDim - 1; ++j){
			if (planeHorizLines[j].empty())	continue;
			std::array<short, 3> ll;
			ll[0] = _planeSets2[0][plane]->D - j;
			ll[1] = j;
			for (auto hit = planeHorizLines[j].begin(); hit != planeHorizLines[j].end(); ++hit)	{
				short z = (short)std::floor(hit->first);
				assert(hit->second.solidRight);
				++z;	++hit;
				if (hit->second.solidRight){  // this is floating point indecision at collision of coincident surfaces
					assert(false);
					// write me
				}
				// COURT - currently making them, but not putting into hash table.
				while (z < hit->first)	{
					if ((j & 1) == (z & 1)){  // all nodes in lattice have either all odd or all even coordinates
						ll[2] = z;
						_vbt->_nodeGridLoci.push_back(ll);  // no mutex at present
					}
					++z;
				}
			}
		}
	}
	auto fit = planeFaces.begin();
	while (fit != planeFaces.end()){
		auto pr = planeFaces.equal_range(fit->first);
		makeConnectedComponentTetFaces(planeSet, plane, pr, planeHorizLines, planeDiagonals);
		fit = pr.second;
	}
}

void vnBccTetCutter::makeConnectedComponentTetFaces(int planeSet, int planeNumber, std::pair<PFIT, PFIT> &face,
	std::vector<PLANE_LINE> &planeHorizLines, std::vector<PLANE_LINE>(&planeDiagonals)[2])
{
	std::pair<short, short> fac = face.first->first;
	bool oddFace = ((fac.first + fac.second) & 1);
	PLANE_LINE *h, *d0, *d1;
	if (oddFace)
		h = &planeHorizLines[fac.second + 1];
	else
		h = &planeHorizLines[fac.second];
	d0 = &planeDiagonals[0][_vbt->_gridSize[((planeSet >> 1) + 1) % 3] + ((fac.first - fac.second + 1) >> 1)];
	d1 = &planeDiagonals[1][(fac.second + fac.first + 2) >> 1];
	std::list<std::list<triSegment2> > edgePolygons;
	// first get all open edgePolygons
	triSegment2 ts;
	auto getEdgeSolids = [&](PLANE_LINE *pl, float low, float high, bool fillX, bool reverse) {
		auto lb = pl->lower_bound(low);
		auto ub = pl->upper_bound(high);
		if (lb == ub){
			if (lb != pl->end() && !lb->second.solidRight){  // entire edge is interior
				std::list<triSegment2> lt;
				if (fillX)
					ts.uv.X = low;
				else
					ts.uv.Y = low;
				ts.triangle = -1;
				lt.push_back(ts);
				if (fillX)
					ts.uv.X = high;
				else
					ts.uv.Y = high;
				lt.push_back(ts);
				if (reverse){
					lt.reverse();
					edgePolygons.push_front(lt);
				}
				else
					edgePolygons.push_back(lt);
			}
			return;
		}
		if (!lb->second.solidRight){
			std::list<triSegment2> lt;
			if (fillX)
				ts.uv.X = low;
			else
				ts.uv.Y = low;
			ts.triangle = -2;
			lt.push_back(ts);
			if (fillX)
				ts.uv.X = lb->first;
			else
				ts.uv.Y = lb->first;
			ts.triangle = lb->second.triangle;
			lt.push_back(ts);
			if (reverse){
				lt.reverse();
				edgePolygons.push_front(lt);
			}
			else
				edgePolygons.push_back(lt);
			++lb;
		}
		while (lb != ub){
			assert(lb->second.solidRight);
			std::list<triSegment2> lt;
			if (fillX)
				ts.uv.X = lb->first;
			else
				ts.uv.Y = lb->first;
			ts.triangle = lb->second.triangle;
			lt.push_back(ts);
			++lb;
			if (lb == ub){
				if (fillX)
					ts.uv.X = high;
				else
					ts.uv.Y = high;
				ts.triangle = -2;
			}
			else{
				if (fillX)
					ts.uv.X = lb->first;
				else
					ts.uv.Y = lb->first;
				ts.triangle = lb->second.triangle;
			}
			lt.push_back(ts);
			if (reverse){
				lt.reverse();
				edgePolygons.push_front(lt);
			}
			else
				edgePolygons.push_back(lt);
			if (lb != ub)
				++lb;
		}
	};
	if (oddFace){
		ts.uv.X = -3.0f;
		getEdgeSolids(d1, (float)fac.second, (float)fac.second + 1.0f, false, true);
		ts.uv.X = -1.0f;
		getEdgeSolids(d0, (float)fac.second, (float)fac.second + 1.0f, false, false);
		ts.uv.Y = fac.second + 1.0f;
		getEdgeSolids(h, (float)fac.first, (float)fac.first + 2.0f, true, true);
		// compute diagonal x values
		for (auto &tp : edgePolygons){
			for (auto &p : tp){
				if (p.uv.X < -2.0f)
					p.uv.X = fac.first + 1.0f - (p.uv.Y - fac.second);
				else if (p.uv.X < -0.5f)
					p.uv.X = fac.first + 1.0f + (p.uv.Y - fac.second);
				else ;
			}
		}
	}
	else{
		ts.uv.Y = (float)fac.second;
		getEdgeSolids(h, (float)fac.first, (float)fac.first + 2.0f, true, false);
		ts.uv.X = -1.0f;
		getEdgeSolids(d0, (float)fac.second, (float)fac.second + 1.0f, false, true);
		ts.uv.X = -3.0f;
		getEdgeSolids(d1, (float)fac.second, (float)fac.second + 1.0f, false, false);
		// compute diagonal x values
		for (auto &tp : edgePolygons){
			for (auto &p : tp){
				if (p.uv.X < -2.0f)
					p.uv.X = fac.first + 2.0f - (p.uv.Y - fac.second);
				else if (p.uv.X < -0.5f)
					p.uv.X = fac.first + (p.uv.Y - fac.second);
				else;
			}
		}
	}
	// remove any doubled corner points
	auto epit = edgePolygons.begin();
	while (epit != edgePolygons.end()){
		if (epit->back().triangle < 0){  // all corners points should be doubled
			auto pNext = epit;
			++pNext;
			if (pNext == edgePolygons.end())
				pNext = edgePolygons.begin();
			assert(pNext->front().triangle < 0 && epit->back().uv == pNext->front().uv);
			epit->pop_back();
			if (epit != pNext){
				epit->splice(epit->end(), *pNext);
				edgePolygons.erase(pNext);
			}
			else
				break;
		}
		else
			++epit;
	}
	std::list<std::list<triSegment2> > facePolygons;
	// now splice or add (if closed interior polygon) triangle-edge strings traversing this triangle+
	while (face.first != face.second){
		if (face.first->second.back().uv.X < FLT_MAX){  // a closed interior polygon
			facePolygons.push_back(std::list<triSegment2>());
			facePolygons.back().splice(facePolygons.back().end(), face.first->second);
		}
		else{
			auto lp = &face.first->second;
			assert(lp->back().uv[0] > 1e32f);  // unclosed chain with this uv not inside face
			epit = edgePolygons.begin();
			while (epit != edgePolygons.end()){
				if (epit->front().triangle == lp->back().triangle){
					lp->pop_back();
					epit->splice(epit->begin(), *lp);
					if (epit->front().triangle == epit->back().triangle){  // polygon closed - done
						facePolygons.splice(facePolygons.end(), edgePolygons, epit);
						epit = edgePolygons.end();
					}
					break;
				}
				++epit;
			}
			if (epit == edgePolygons.end()){
				assert(face.first->second.empty());
				++face.first;
				continue;
			}
			auto epit2 = edgePolygons.begin();
			while (epit2 != edgePolygons.end()){
				if (epit2 == epit){
					++epit2;
					continue;
				}
				if (epit2->back().triangle == epit->front().triangle){
					epit->splice(epit->begin(), *epit2);
					edgePolygons.erase(epit2);
					break;
				}
				++epit2;
			}
		}
		assert(face.first->second.empty());
		++face.first;
	}
	if (!edgePolygons.empty()){  // single polygon of only interior corner points surrounding an interior polygon(s)
		assert(!facePolygons.empty());
		assert(edgePolygons.size() < 2);
		facePolygons.splice(facePolygons.end(), edgePolygons);
	}
	// All polygons now are closed polygons with inside edge vertices labelled as < 0
	createSurfaceTetFaces(fac, facePolygons, _planeSets2[planeSet][planeNumber].get());
}

void vnBccTetCutter::createSurfaceTetFaces(std::pair<short, short> &face, std::list<std::list<triSegment2> > &facePolygons, bccPlane *planeData)
{
	std::list<std::list<triSegment2> > exteriorPolygons, interiorPolygons;
	std::list<vnTetFace> exteriorTetFaces;
	bool oddFace = ((face.first + face.second) & 1);
	auto counterClockwise = [](std::list<triSegment2> &poly) ->bool {
		std::list<triSegment2>::iterator minIt, pit;
		double minX = DBL_MAX;
		for (pit = poly.begin(); pit != poly.end(); ++pit) {
			if (pit->uv[0] < minX) {
				minX = pit->uv[0];
				minIt = pit;
			}
		}
		Vec2d uv0, uv1;
		pit = minIt;
		do {
			if (pit == poly.begin())
				pit = poly.end();
			--pit;
		} while (pit->uv == minIt->uv && pit != minIt);
		if (pit == minIt)  // zero area polygon
			return true;
		uv0 = pit->uv - minIt->uv;
		pit = minIt;
		do {
			++pit;
			if (pit == poly.end())
				pit = poly.begin();
			uv1 = pit->uv - minIt->uv;
			double a = uv0[0] * uv1[1] - uv0[1] * uv1[0];
			if (a > 0.0)
				return false;
			else if (a < 0.0)
				return true;
			else;
		} while (pit != minIt);
		return true;  // zero area polygon
	};
	auto windingNumber = [](Vec2d &uv, std::list<triSegment2> &poly) ->int{
		// Winding number test for a point in a polygon.  Assumes closed polygon without repeat of first point.
		auto lit = poly.end();
		--lit;
		int wn = 0;
		for (auto pit = poly.begin(); pit != poly.end(); ++pit){
			if (lit->uv.Y <= uv.Y) {  // start
				if (pit->uv.Y  > uv.Y)  // an upward crossing
					if (((pit->uv.X - lit->uv.X) * (uv.Y - lit->uv.Y) - (uv.X - lit->uv.X) * (pit->uv.Y - lit->uv.Y)) > 0.0)  // uv left of  edge
						++wn;            // have  a valid up intersect
			}
			else {  // lit->uv.Y > uv.y
				if (pit->uv.Y <= uv.Y)     // a downward crossing
					if (((pit->uv.X - lit->uv.X) * (uv.Y - lit->uv.Y) - (uv.X - lit->uv.X) * (pit->uv.Y - lit->uv.Y)) < 0.0)  // uv right of  edge
						--wn;            // have  a valid down intersect
			}
			lit = pit;
		}
		return wn;  // winding number.  0 only when P is outside. Sign gives clockwiseness.
	};
	while (!facePolygons.empty()){
		auto fpit = facePolygons.begin();
		vnTetFace tf;
		tf.interiorNodes = 0;
		auto lastIt = fpit->begin();
		bool exterior = false;
		for (auto tsit = fpit->end(); tsit != fpit->begin();){
			--tsit;
			if (tsit->triangle < 0){  // interior corner point
				exterior = true;
				// COURT - decided to have plane horizontal be vertex0 to vertex1 for both odd and even faces. This means even is counterclockwise and odd clockwise
				if (tsit->uv[0] == face.first){
					tf.interiorNodes |= 1;
					assert(tsit->uv[1] == face.second + ((face.first + face.second) & 1) ? 1 : 0);
				}
				else if (tsit->uv[0] == face.first + 2.0){
					tf.interiorNodes |= 2;
					assert(tsit->uv[1] == face.second + ((face.first + face.second) & 1) ? 1 : 0);
				}
				else if (tsit->uv[0] == face.first + 1.0){  // doesn't matter if up or down
					tf.interiorNodes |= 4;
					assert(tsit->uv[1] == face.second + ((face.first + face.second) & 1) ? 0 : 1);
				}
				else
					assert(false);
				if (tsit != fpit->begin()){
					--tsit;
					if (tsit->triangle > -1)
						tf.edgeTriangles.push_back(tsit->triangle);
					else
						++tsit;
				}
			}
			else{
				if (tsit->triangle == lastIt->triangle){  // contains a border segment so is an exterior solid polygon
					tf.edgeTriangles.push_back(tsit->triangle);
					--tsit;
					if (tsit->triangle < 0){
						++tsit;
						continue;
					}
					assert(tsit->triangle != tf.edgeTriangles.back());
					tf.edgeTriangles.push_back(tsit->triangle);
					exterior = true;
				}
				else{
					if (tf.edgeTriangles.empty() || tsit->triangle != tf.edgeTriangles.front())
						tf.interiorTriangles.push_back(tsit->triangle);
				}

			}
			lastIt = tsit;
		}
		if (exterior){
			exteriorTetFaces.push_back(std::move(tf));
			assert(counterClockwise(*fpit));
			exteriorPolygons.splice(exteriorPolygons.end(), facePolygons, fpit);
		}
		else{  // test this interior polygon for clockwiseness
			if (counterClockwise(*fpit)){
				exteriorTetFaces.push_back(std::move(tf));
				exteriorPolygons.splice(exteriorPolygons.end(), facePolygons, fpit);
			}
			else
				interiorPolygons.splice(interiorPolygons.end(), facePolygons, fpit);
		}
	}
	if (!interiorPolygons.empty()){  // any interior must be inside an exterior one, so add its triangles to its interior tf
		for (auto &ip : interiorPolygons){
			auto efit = exteriorTetFaces.begin();
			auto epit = exteriorPolygons.begin();
			while (epit != exteriorPolygons.end()){
				if (windingNumber(ip.front().uv, *epit) != 0){
					for (auto &ts2 : ip)
						efit->interiorTriangles.push_back(ts2.triangle);
					break;
				}
				++epit;
				++efit;
			}
			assert(epit != exteriorPolygons.end());
		}
	}
	// create a new tet face for each exteriorPolygon
	_surfaceTetFaceNumber += exteriorTetFaces.size();
	for (auto &etf : exteriorTetFaces)
		planeData->vnTetFaces.insert(std::make_pair(face, std::move(etf)));
}

void vnBccTetCutter::collectSurfaceTetCentersFaces()
{
	// Every tet face intersected by the surface will have at one surface tet on either side.
	// This routine gets tet centers on either side of each intersected face, hashes them, and adds the face to that hash.
	// Also fills in any interior nodes in the face.
	for (int i = 0; i < 6; ++i){
		int j, n = _planeSets2[i].size();
		for (j = 0; j < n; ++j){
			int planeD = _planeSets2[i][j]->D;
			auto tfit = _planeSets2[i][j]->vnTetFaces.begin();
			auto tfEnd = _planeSets2[i][j]->vnTetFaces.end();
			while (tfit != tfEnd){
				// get tet centers on either side of this face
				bccTetCentroid tc0, tc1;
				int axis = i >> 1, oddFace = (tfit->first.first + tfit->first.second) & 1;
				// midpoint of face horizontal
				tc0.xyz[axis] = (i & 1) ? planeD + tfit->first.second + oddFace : planeD - tfit->first.second - oddFace;
				tc0.xyz[(axis + 1) % 3] = tfit->first.second + oddFace;
				tc0.xyz[(axis + 2) % 3] = tfit->first.first + 1;
				tc0.halfCoordAxis = axis;
				assert((tc0.xyz[axis] & 1) == (tc0.xyz[(axis + 1) % 3] & 1) && (tc0.xyz[axis] & 1) != (tc0.xyz[(axis + 2) % 3] & 1));
				tc1.xyz = tc0.xyz;
				tc1.halfCoordAxis = ((axis + 1) % 3);
				if (oddFace){
					--tc1.xyz[(axis + 1) % 3];
					if (i & 1)
						--tc0.xyz[axis];
				}
				else{
					if ((i & 1) < 1)
						--tc0.xyz[axis];
				}
				auto fl0 = &_surfaceTetFaces.insert(std::make_pair(tc0.ll, std::list<vnTetFace*>())).first->second;
				auto fl1 = &_surfaceTetFaces.insert(std::make_pair(tc1.ll, std::list<vnTetFace*>())).first->second;
				auto tfLast = tfit;
				while (tfit->first == tfLast->first){
					tfit->second.set = i;
					fl0->push_back(&tfit->second);
					fl1->push_back(&tfit->second);
					++tfit;
					if (tfit == tfEnd)
						break;
				}
			}
		}
	}
}

std::pair<short, short> vnBccTetCutter::uvToPlaneFace(const Vec2d &v)
{  // plane face indexing has V component in the horizontal span minimum
	// The U component alternates with even permutation U+V has even U as upward pointing triangles, odd U down pointing.
	// Similarly odd permutation U+V has odd U as upward pointing triangles, even U down pointing.
	std::pair<short, short> sp;
	sp.first = (short)std::floor(v.X);
	sp.second = (short)std::floor(v.Y);
	if ((sp.first + sp.second) & 1){
		if (v.Y - sp.second < sp.first + 1.0 - v.X)
			--sp.first;
	}
	else{
		if (v.Y - sp.second > v.X - sp.first)
			--sp.first;
	}
	return sp;
}

void vnBccTetCutter::getFaceEdgePath(const Vec2d *startUv, std::pair<short, short> startFace, const Vec2d *endUv, std::pair<short, short> endFace, int triangle, int secondAxis,
	std::list<std::pair<short, short> > &triFaces, std::vector<PLANE_LINE> &horizLines, std::vector<PLANE_LINE>(&diagonals)[2])
{
	triFaces.clear();
	int previousEdge = -1;
	planeLineCrossing plc;
	plc.solidRight = 0;
	plc.triangle = triangle;
	Vec2d N;
	N = *endUv - *startUv;  // due to counterclockwiseness around solid, its solid surface normal is y, -x
	double s, t;  // s is parameter along line from startUv to endUv. t is parameter along a candidate triangle edge.
	std::pair<short, short> faceNow = startFace;
	while (faceNow != endFace)	{
		triFaces.push_back(faceNow);
		bool oddTri = ((faceNow.first + faceNow.second) & 1);
		int i;
		for (i = 0; i < 3; ++i)	{
			if (i == previousEdge)
				continue;
			if (i < 1)	{  // horizontal axis crossing?
				if (faceNow.second == endFace.second)
					continue;
				int horizIndx = faceNow.second + (oddTri? 1 : 0);
				s = (horizIndx - startUv->Y) / N[1];
				if (s < 0.0 || s > 1.0)
					continue;
				if (!(s > 0.0) && N.Y >= 0.0){  // On horizLine but can't be a starting edge in this direction
					assert(previousEdge < 0);
					continue;
				}
				t = s * N[0] + startUv->X;
				if (t < (double)faceNow.first || t >= (double)faceNow.first + 2.0)
					continue;
				plc.solidRight = N[1] < 0.0 ? 1 : 0;
				horizLines[horizIndx].insert(std::make_pair(t, plc));
				if (oddTri)
					++faceNow.second;
				else
					--faceNow.second;
				previousEdge = i;
				break;
			}
			else if (i < 2)	{  // diag0 axis crossing - increasing direction in uv is [1, 1].
				if (N.X == N.Y) // crossing not possible
					continue;
				double dy = startUv->Y - faceNow.second + (oddTri ? 1.0 : 0.0), dx = startUv->X - faceNow.first;
				if (dy == dx && ((oddTri && -N.X + N.Y >= 0.0) || (!oddTri && N.X - N.Y >= 0.0))){    // use same test for being exactly on diagonal line as uvToPlaneFace() or possible roundoff error
					assert(previousEdge < 0);
					continue;
				}
				s = dx - dy;
				s /= (N.Y - N.X);
				t = startUv->X - faceNow.first + N.X * s + (oddTri ? -1.0 : 0.0);
				if (s < 0.0 || s > 1.000001 || t < 0.0 || t > 1.0)  // allow a little destination roundoff error
					continue;
				int diagIndx = _vbt->_gridSize[secondAxis] + ((faceNow.first - faceNow.second + 1) >> 1);
				plc.solidRight = (N[1] - N[0]) < 0.0 ? 1 : 0;
				diagonals[0][diagIndx].insert(std::make_pair(t + (double)faceNow.second, plc));
				if (oddTri)
					++faceNow.first;
				else
					--faceNow.first;
				previousEdge = i;
				break;
			}
			else	{	// diag1 axis crossing - increasing direction in uv is [-1, 1].
				if (N.X == -N.Y) // crossing not possible
					continue;
				double dy = startUv->Y - faceNow.second, dx = faceNow.first + (oddTri ? 1.0 : 2.0) - startUv->X;
				if (dx == dy && ((oddTri && N.X + N.Y >= 0.0) || (!oddTri && -N.X + -N.Y >= 0.0))){    // use same test for being exactly on diagonal line as uvToPlaneFace() or possible roundoff error
					assert(previousEdge < 0);
					continue;
				}
				s = dx - dy;
				s /= N.Y + N.X;
				t = startUv->Y - faceNow.second + N.Y * s;
				if (s < 0.0 || s > 1.000001 || t < 0.0 || t > 1.0)  // allow a little destination roundoff error
					continue;
				int diagIndx = (faceNow.second + faceNow.first + 2) >> 1;
				plc.solidRight = (-N[1] - N[0]) < 0.0 ? 1 : 0;
				diagonals[1][diagIndx].insert(std::make_pair(t + (double)faceNow.second, plc));
				if (oddTri)
					--faceNow.first;
				else
					++faceNow.first;
				previousEdge = i;
				break;
			}
		}
		assert(i < 3);
	}
	assert(triFaces.front() == startFace);
	triFaces.pop_front();
}

void vnBccTetCutter::tetConnectedSurface(bccTetCentroid tc, std::set<long> &triangles, std::vector<long> &vertices)
{  // Inputs tc and surface triangles seed.  Returns all the surface triangles and vertices inside the tet.
	std::set<long> verts;
	std::forward_list<long> edgeTris;
	edgeTris.assign(triangles.begin(), triangles.end());
	triangles.clear();
	// input triangles should already contain triangles intersecting tet faces. Must find all new interior triangles using their interior vertices.
	std::function<void(long)> recurseInteriorTriangles = [&](long triangle){
		if (!triangles.insert(triangle).second)
			return;
		long *tr = _mt->triangleVertices(triangle);
		for (int i = 0; i < 3; ++i){
			if (_vertexTetLoci[tr[i]] == tc){
				if (!verts.insert(tr[i]).second)
					continue;
				unsigned long adj = _mt->triAdjs(triangle)[i];
				while (adj >> 2 != triangle){
					recurseInteriorTriangles(adj >> 2);
					assert(_mt->triangleVertices(adj >> 2)[((adj & 3) + 1) % 3] == tr[i]);
					adj = _mt->triAdjs(adj >> 2)[((adj&3)+ 1) % 3];
				}
			}
		}
	};
	for (auto &et : edgeTris)
		recurseInteriorTriangles(et);
	vertices.assign(verts.begin(), verts.end());
}

void vnBccTetCutter::createVirtualNodedSurfaceTets()
{
	bccTetCentroid tc;
	int firstAxis;  // secondAxis will be halfCoordAxis
	bool firstAxisUp;
	auto createTetAssignNodes = [&](std::list<vnTetFace*> &faces) ->long{
		long newTet = _surfaceTetNumber.fetch_add(1);
		_vbt->_tetCentroids.push_back(tc);
		newTet <<= 2;
		for (auto &f : faces){
			assert(f->set >> 1 != ((tc.halfCoordAxis + 1) % 3));
			if ((f->set >> 1) == firstAxis){
				int faceSide = f->set & 1 ? 0 : 1;
				f->tetNodes[faceSide][0] = newTet;
				f->tetNodes[faceSide][1] = newTet | 1;
				if (f->set&1)
					f->tetNodes[faceSide][2] = newTet | 2;
				else
					f->tetNodes[faceSide][2] = newTet | 3;
			}
			else{  // second axis face
				int faceSide = f->set & 1;
				if (firstAxisUp){
					f->tetNodes[faceSide][0] = newTet | 2;
					f->tetNodes[faceSide][1] = newTet | 3;
					if (f->set & 1)
						f->tetNodes[faceSide][2] = newTet | 1;
					else
						f->tetNodes[faceSide][2] = newTet;
				}
				else{
					f->tetNodes[faceSide][0] = newTet | 3;
					f->tetNodes[faceSide][1] = newTet | 2;
					if (f->set & 1)
						f->tetNodes[faceSide][2] = newTet;
					else
						f->tetNodes[faceSide][2] = newTet | 1;
				}
			}
		}
		return newTet >> 2;
	};
	_surfaceTetNumber.store(0);
	_vbt->_tetCentroids.clear();  //  _surfaceTetCenters.clear();
	_vbt->_tetCentroids.reserve(_surfaceTetFaceNumber);
	for (auto &st : _surfaceTetFaces){
		tc.ll = st.first;
		firstAxis = (tc.halfCoordAxis + 2) % 3;
		firstAxisUp = (tc.xyz[tc.halfCoordAxis] + tc.xyz[firstAxis]) & 1;
		// get groups of faces sharing common edge triangle intersects - most common surface tets
		struct edgeIntersectGroup{
			std::set<long> edgeTriangles;
			std::list<vnTetFace*> faces;
		};
		std::list<edgeIntersectGroup> eig;
		auto tfit = st.second.begin();
		while (tfit != st.second.end()){
			if ((*tfit)->edgeTriangles.empty()){
				++tfit;
				continue;
			}
			edgeIntersectGroup *prevFused = NULL;
			auto egit = eig.begin();
			while (egit != eig.end()){
				if (prevFused == NULL){
					for (auto &fet : (*tfit)->edgeTriangles){
						if (egit->edgeTriangles.find(fet) != egit->edgeTriangles.end()){
							egit->edgeTriangles.insert((*tfit)->edgeTriangles.begin(), (*tfit)->edgeTriangles.end());
							egit->faces.push_back(*tfit);
							prevFused = &(*egit);
							break;
						}
					}
					++egit;
				}
				else{
					auto pfit = (*tfit)->edgeTriangles.begin();
					while (pfit != (*tfit)->edgeTriangles.end()){
						if (egit->edgeTriangles.find(*pfit) != egit->edgeTriangles.end()){  // fuse these 2 groups
							prevFused->edgeTriangles.insert(egit->edgeTriangles.begin(), egit->edgeTriangles.end());
							prevFused->faces.splice(prevFused->faces.end(), egit->faces);
							egit = eig.erase(egit);
							break;
						}
						++pfit;
					}
					if (pfit == (*tfit)->edgeTriangles.end())
						++egit;
				}
			}
			if (prevFused == NULL){
				eig.push_back(edgeIntersectGroup());
				eig.back().edgeTriangles.insert((*tfit)->edgeTriangles.begin(), (*tfit)->edgeTriangles.end());
				eig.back().faces.push_back(*tfit);
			}
			tfit = st.second.erase(tfit);
		}
		// have all independent edge groups and interior faces. Do any interior tests for all except single tet centroid locations.
		if (eig.size() < 2 && st.second.empty())
			createTetAssignNodes(eig.front().faces);
		else if (eig.empty() && st.second.size() < 2){
			std::list<vnTetFace*> oneFace;
			oneFace.splice(oneFace.end(), st.second, st.second.begin());
			createTetAssignNodes(oneFace);
		}
		else{  // have gone as far as possible using simple common tet edge intersections.  Now must use all triangles and do
			// interior tests for possible further agglomeration or multiple tets at this location
			for (auto &eg : eig){
				for (auto &fac : eg.faces)
					eg.edgeTriangles.insert(fac->interiorTriangles.begin(), fac->interiorTriangles.end());
			}
			for (auto &fac : st.second){
				eig.push_back(edgeIntersectGroup());
				assert(fac->edgeTriangles.empty());
				eig.back().edgeTriangles.insert(fac->interiorTriangles.begin(), fac->interiorTriangles.end());
				eig.back().faces.push_back(fac);
			}
			st.second.clear();
			auto eit = eig.begin();
			while (eit != eig.end()){
				std::vector<long> vertices;
				tetConnectedSurface(tc, eit->edgeTriangles, vertices);
				auto eit2 = eit;
				++eit2;
				bool recollectInterior = false;  // "dents" into a solid can lead to omissions in gathering interior
				while (eit2 != eig.end()){
					bool notFound = true;
					for (auto tri : eit2->edgeTriangles){
						if (eit->edgeTriangles.find(tri) != eit->edgeTriangles.end()){
							for (auto tri2 : eit2->edgeTriangles){  // the full interior connected search should have gotten everything
								if (eit->edgeTriangles.find(tri2) == eit->edgeTriangles.end()){
									eit->edgeTriangles.insert(eit2->edgeTriangles.begin(), eit2->edgeTriangles.end());
									recollectInterior = true;
								}
							}
							eit->faces.splice(eit->faces.end(), eit2->faces);
							eit2 = eig.erase(eit2);
							notFound = false;
							break;
						}
					}
					if (notFound)
						++eit2;
				}
				// now have agglomerated connected faces in eit
				long thisTet = createTetAssignNodes(eit->faces);
				// since there are likely multiple tets for this BCC locus, assign this tet to its internal vertices
				if (recollectInterior){
					vertices.clear();
					tetConnectedSurface(tc, eit->edgeTriangles, vertices);
				}
				for (auto &v : vertices){
					assert(_vbt->_vertexTets[v] < 0);  // should only happen in one surface tet
					_vbt->_vertexTets[v] = thisTet;
				}
				eit = eig.erase(eit);
			}
		}
	}
}

void vnBccTetCutter::fillNonVnTetCenter()
{
	_vbt->_firstInteriorTet = _vbt->_tetNodes.size();
	// interior tet nodes are created sequentially in z.  Use this to get all z-x and z-y tets with less hash table searches per tet.
	for (int n = _interiorNodes.size() - 1, i = 0; i < n; ++i){
		auto lpi = _vbt->_nodeGridLoci[i].data();
		auto lpj = _vbt->_nodeGridLoci[i + 1].data();
		if (lpi[0] == lpj[0] && lpi[1] == lpj[1] && lpi[2] + 2 == lpj[2]){ // this z segment exists.  Look for z-x and z-y tets
			long nn[4];
			for (int j = 0; j < 4; ++j){
				auto ll = _vbt->_nodeGridLoci[i];
				ll[0] += j & 1 ? 1 : -1;
				ll[1] += j & 2 ? 1 : -1;
				++ll[2];
				auto in = _interiorNodes.find(ll);
				if (in != _interiorNodes.end()){
					nn[j] = in->second;
				}
				else
					nn[j] = -1;
			}
			// See vnBccTetrahedra.h for how tetrahedral nodes are ordered.
			for (int j = 0; j < 3; j += 2){
				if (nn[j] > -1 && nn[j + 1] > -1){
					std::array<long, 4> nodes;
					nodes[0] = i;
					nodes[1] = i + 1;
					bccTetCentroid tc;
					tc.halfCoordAxis = 1;
					tc.xyz = _vbt->_nodeGridLoci[i];  // *lpi
					tc.xyz[2] += 1;
					if (j){
						nodes[2] = nn[j + 1];
						nodes[3] = nn[j];
					}
					else{
						--tc.xyz[tc.halfCoordAxis];
						nodes[2] = nn[j];
						nodes[3] = nn[j + 1];
					}
					if (_vbt->_tetHash.find(tc.ll) == _vbt->_tetHash.end()){  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(nodes);
					}
				}
			}
			for (int j = 0; j < 2; ++j){
				if (nn[j] > -1 && nn[j + 2] > -1){
					std::array<long, 4> nodes;
					nodes[0] = nn[j];
					nodes[1] = nn[j + 2];
					bccTetCentroid tc;
					tc.halfCoordAxis = 0;
					tc.xyz = _vbt->_nodeGridLoci[i];  // *lpi
					tc.xyz[2] += 1;
					if (j){
						nodes[2] = i;
						nodes[3] = i + 1;
					}
					else{
						--tc.xyz[tc.halfCoordAxis];
						nodes[2] = i + 1;
						nodes[3] = i;
					}
					if (_vbt->_tetHash.find(tc.ll) == _vbt->_tetHash.end()){  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(nodes);
					}
				}
			}
		}
	}
	// now create all x-y tets
	for (int n = _interiorNodes.size(), i = 0; i < n; ++i){
		auto ll = _vbt->_nodeGridLoci[i];
		ll[0] += 2;
		auto in = _interiorNodes.find(ll);
		if (in == _interiorNodes.end())
			continue;
		long x2 = in->second;
		bccTetCentroid tc;
		tc.xyz = _vbt->_nodeGridLoci[i];
		tc.xyz[0] += 1;
		for (int j = 1; j > -1; --j){
			ll = tc.xyz;
			ll[2] += j ? 1 : -1;
			--ll[1];
			in = _interiorNodes.find(ll);
			if (in == _interiorNodes.end())
				continue;
			long y0 = in->second;
			ll[1] += 2;
			in = _interiorNodes.find(ll);
			if (in == _interiorNodes.end())
				continue;
			tc.halfCoordAxis = 2;
			if (j < 1)
				--tc.xyz[2];
			std::array<long, 4> nodes;
			nodes[0] = i;
			nodes[1] = x2;
			nodes[2] = j ? in->second : y0;
			nodes[3] = j ? y0 : in->second;
			if (_vbt->_tetHash.find(tc.ll) == _vbt->_tetHash.end()){  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
				_vbt->_tetCentroids.push_back(tc);  // perhaps just push back if not multithreading
				_vbt->_tetNodes.push_back(nodes);
			}
		}
	}
}

