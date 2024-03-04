#include <assert.h>
#include "Vec2d.h"
#include "Vec3f.h"
#include "Mat3x3f.h"
#include "Mat2x2d.h"
#include "Mat2x2f.h"
#include "insidePolygon.h"
#include "boundingBox.h"
#include "vnBccTetrahedra.h"
#include <functional>
#include <exception>
#include "closestPointOnTriangle.h"
#include "fence.h"
#include "FacialFlapsGui.h"
#include "clockwise.h"
#include "deepCut.h"

// below for threading building blocks
#include "oneapi/tbb.h"
//#include "tbb/concurrent_vector.h"
//#include "tbb/parallel_for.h"
//#include "tbb/blocked_range.h"
//#include "tbb/tick_count.h"

// this version replaced my old hand written constrained Delaunay triangulation with an excellent library
// which can be found at: https://github.com/artem-ogre/CDT
// Please visit that site regarding MPL licensing and further implementation details.
#include "CDT.h"

float deepCut::_cutSpacingInv = 15.0f;  // inverse of deep cut interior point spacing. COURT - Model dependent. Should be moved into scene file later.

bool deepCut::cutDeep()  // interpost connection data already loaded in _deepPosts
{
	_previousSkinTopEnd = -1;
	_loopSkinTopBegin = -1;
	_endPlanes[0].P.X = DBL_MAX;
	_endPlanes[1].P.X = DBL_MAX;
	_holePolyLines.clear();
	if (_deepPosts.size() < 2)
		return false;
	std::map<std::pair<int, int>, int > rtiHits;  // COURT - get rtiHits as separate pass
	auto getRtiHits = [&]() {
		for (auto& dp : _deepPosts) {
			for (auto& ti : dp.triIntersects) {
				if (ti.scl.rtiIndexTo > -1) {
					auto pr = rtiHits.insert(std::make_pair(std::make_pair(ti.postNum, ti.rayIndex), 1));
					if (!pr.second)
						++pr.first->second;
					pr = rtiHits.insert(std::make_pair(std::make_pair(ti.scl.rtiPostTo, ti.scl.rtiIndexTo), 1));
					if (!pr.second)
						++pr.first->second;
				}
			}
		}
	};
	rtiHits.clear();
	getRtiHits();
	for (auto rit = rtiHits.begin(); rit != rtiHits.end(); ) {  // these must be per post pairs. Searching for interior intersected holes in interpost quads only.
		int post = rit->first.first;
		auto& ti = _deepPosts[post].triIntersects;
		int i0 = rit->first.second + 1;
		assert(i0 == 1);
		++rit;
		if (rit->first.first != post)
			throw(std::logic_error("Program error. Topological connection problem 0 in cutDeep().\n"));
		int i1 = rit->first.second;
		while (i0 < i1) {  // COURT - check that open end gaps handled correctly
			bool interpostPath = false;
			for(int k = i0 + 1; k< i1; k += 2){
				double lowV;
				if (post > 0 && ti[i0].scl.rtiIndexTo < 0) {
					if (surfacePath(ti[i0], ti[k], false, lowV) < DBL_MAX) {
						ti[i0].scl.rtiPostTo = post;
						ti[i0].scl.rtiIndexTo = k;
						ti[i0].scl.lowestV = lowV;
						interpostPath = true;
					}
				}
				if (post < _deepPosts.size() - 1 && ti[k].scl.rtiIndexTo < 0) {
					if (surfacePath(ti[k], ti[i0], false, lowV) < DBL_MAX) {
						ti[k].scl.rtiPostTo = post;
						ti[k].scl.rtiIndexTo = i0;
						ti[k].scl.lowestV = lowV;
						interpostPath = true;
					}
				}
				if (interpostPath) {
					i0 = k + 1;
					break;
				}
			}
			if(!interpostPath)
				throw(std::logic_error("Program error. Topological connection problem 1 in cutDeep().\n"));
		}
		++rit;
		if (rit != rtiHits.end() && rit->first.first == post)
			throw(std::logic_error("Program error. Topological connection problem 2 in cutDeep().\n"));
	}
	rtiHits.clear();
	getRtiHits();
	if (!_deepPosts.front().closedEnd)
		if (!connectOpenEnd(0))
			return false;
	if (!_deepPosts.back().closedEnd)
		if (!connectOpenEnd(_deepPosts.size() - 1))
			return false;
	rtiHits.clear();
	getRtiHits();
	// trim unused back end intersects from posts. Leave any open end gap intersects so scl indexing isn't broken.
	for (auto rit = rtiHits.begin(); rit != rtiHits.end(); ) {
		int post = rit->first.first;
		int i0;
		do {
			i0 = rit->first.second + 1;
			++rit;
		} while (rit != rtiHits.end() && rit->first.first == post);
		_deepPosts[post].triIntersects.erase(_deepPosts[post].triIntersects.begin() + i0, _deepPosts[post].triIntersects.end());
	}
	_preDeepCutVerts = _mt->numberOfVertices() - 1;
	auto punchDeepVert = [&](rayTriangleIntersect& ri) ->bool{
		if (ri.deepVert > -1)
			return true;
		int mat = _mt->triangleMaterial(ri.triangle);
		float uv[2] = { (float)ri.uv[0], (float)ri.uv[1] };
		if (mat == 2) {
			int topVertex, bottomVertex;
			createFlapTopBottomVertices(ri.triangle, uv, topVertex, bottomVertex);
			ri.mat2Vert = topVertex;
			ri.deepVert = bottomVertex;
		}
		else if (mat == 1 || (mat > 4 && mat < 10)) {
			Vec3f gridLocus, bw;
			int tet = _vbt->parametricTriangleTet(ri.triangle, uv, gridLocus);
			_vbt->gridLocusToBarycentricWeight(gridLocus, _vbt->_tetCentroids[tet], bw);
			ri.deepVert = _mt->addNewVertexInMidTriangle(ri.triangle, uv);
			assert(ri.deepVert == _vbt->_vertexTets.size());
			_vbt->_vertexTets.push_back(tet);
			_vbt->_barycentricWeights.push_back(bw);
		}
		else {
			std::cout << "trying to create a deep cut point on material " << mat << "\n";
			return false;
		}
		ri.triangle = -1;  // triangle not valid after punch
		return true;
	};

	// punch rtiHit triangles with an interior vertex
	for (auto &r : rtiHits)
		if (!punchDeepVert(_deepPosts[r.first.first].triIntersects[r.first.second]))
			return false;
	_mt->findAdjacentTriangles(true);
	getDeepSpatialCoordinates();  // after punches need to redo
	// create bilinear surface precomputations for each patch per nVidia subroutine
	for (int i = 1; i < _deepPosts.size(); ++i) {
		auto &ti = _deepPosts[i].triIntersects;
		int i0, i1;
		for (int n = ti.size(), j = 1; j < n; j += 2) {
			if (ti[j].scl.rtiPostTo == i - 1) {  // interpost connection levels for this quad
				i1 = j;
				i0 = ti[j].scl.rtiIndexTo;
				break;
			}
		}
		Vec3d P00 =_deepXyz[_deepPosts[i - 1].triIntersects[i0].deepVert];
		Vec3d P10 = _deepXyz[_deepPosts[i].triIntersects[i1].deepVert];
		Vec3d P01 = _deepXyz[_deepPosts[i - 1].triIntersects.front().deepVert];
		Vec3d P11 = _deepXyz[_deepPosts[i].triIntersects.front().deepVert];
		makeBilinearPatch(P00, P10, P01, P11, _deepPosts[i].bl);
	}
	// now get all interior solid vertices on posts and collect vertices which should not be duplicated
	std::set<int> nonDupedVertices;
	for (auto rit = rtiHits.begin(); rit != rtiHits.end(); ++rit) {  // these must be pairs searching for interior intersected holes
		int post = rit->first.first;
		int i0 = rit->first.second;
		int hits = rit->second;
		++rit;
		if (post != rit->first.first || rit->first.second != i0 + 1 || rit->second != hits)
			throw(std::logic_error("Program error making solid vertices in cutDeep().\n"));
		auto& dti = _deepPosts[post].triIntersects;
		getDeepCutLine(dti[i0], dti[i0 + 1]);
		if (hits == 1) {  // blind end cut into a solid
			nonDupedVertices.insert(dti[i0].deepVert);
			for (auto& dv : dti[i0].dcl)
				nonDupedVertices.insert(dv);
			nonDupedVertices.insert(dti[i0 + 1].deepVert);
		}
	}


	// make skin incisions and collect all connected deep surface lines on incised surface for split later
	std::list<std::list<int> > surfacePolyLines;  // deepVert lines
	surfacePolyLines.push_back(std::list <int>());

	bool setTopLoopBegin;
	int postTo, idxTo;
	auto inputScl = [&](std::map<std::pair<int, int>, int >::iterator &rit) ->bool{
		double stub;
		auto& ti = _deepPosts[rit->first.first].triIntersects[rit->first.second];
		postTo = ti.scl.rtiPostTo;
		idxTo = ti.scl.rtiIndexTo;
		if (surfacePath(ti, _deepPosts[postTo].triIntersects[idxTo], true, stub) == DBL_MAX)
			return false;
		if (setTopLoopBegin && ti.mat2Vert > -1)
			_loopSkinTopBegin = ti.mat2Vert;
		setTopLoopBegin = false;
		_previousSkinTopEnd = _deepPosts[postTo].triIntersects[idxTo].mat2Vert;
		_mt->findAdjacentTriangles(true);
		updateDeepSpatialCoordinates();  // after punches need to redo
		auto dit = ti.scl.deepVertsTris.begin();
		if (!surfacePolyLines.back().empty()) // skip duplicated entry
			++dit;
		while (dit != ti.scl.deepVertsTris.end()) {
			surfacePolyLines.back().push_back(*dit);
			++dit;
		}
		rtiHits.erase(rit);
		rit = rtiHits.find(std::make_pair(postTo, idxTo));
		return true;
	};
	setTopLoopBegin = true;
	auto rhit = rtiHits.begin();
	while (rhit != rtiHits.end()) {
		if (rhit->second > 1 || _deepPosts[rhit->first.first].triIntersects[rhit->first.second].scl.rtiIndexTo < 0) {
			++rhit;
			continue;
		}
		// at beginning of non-closed loop line
		do {
			if (!inputScl(rhit))
				return false;
		} while (_deepPosts[rhit->first.first].triIntersects[rhit->first.second].scl.rtiIndexTo > -1);
		rtiHits.erase(rhit);
		rhit = rtiHits.begin();
		if (rhit != rtiHits.end())
			surfacePolyLines.push_back(std::list <int>());
	}
	while(!rtiHits.empty()){ // now do any closed loops that remain
		assert(surfacePolyLines.back().empty());
		setTopLoopBegin = true;
		rhit = rtiHits.begin();
		int firstPost = rhit->first.first, firstIdx = rhit->first.second;
		do {
			if (!inputScl(rhit))
				return false;
		}while (postTo != firstPost || idxTo != firstIdx);
		rhit = rtiHits.begin();
		if (rhit != rtiHits.end())
			surfacePolyLines.push_back(std::list <int>());
	}
	// now cut the quads
	updateDeepSpatialCoordinates();  // necessary for hole finding in quads and end planes
	for (int np = _deepPosts.size(), i = 1; i < np; i++)
		deepCutQuad(i);
	if (_endPlanes[0].P.X < DBL_MAX)
		deepCutEndPlane(0);
	if (_endPlanes[1].P.X < DBL_MAX)
		deepCutEndPlane(1);
	struct posTex {
		int pos;
		int tex;
	} ptex;
	ptex.pos = -1;
	ptex.tex = -1;
	std::map<int, posTex> oppositeSideVertices;  // relates vertices on first side of cut to opposite side which are created here.
	auto cloneVertex = [&](const int sourceVertex) ->int {
		// makes a new vertex which is a copy of sourceVertex
		int retval = _mt->addVertices(1);
		_mt->setVertexCoordinate(retval, (const float(&)[3])* _mt->vertexCoordinate(sourceVertex));
		return retval;
	};
	auto dupVertex = [&](int vert) {
		auto pr = oppositeSideVertices.insert(std::make_pair(vert, ptex));
		if (pr.second) {
			if (nonDupedVertices.find(vert) != nonDupedVertices.end())  // no dup
				pr.first->second.pos = vert;
			else {  // all deep vertices so don't add to _deepBed
				pr.first->second.pos = cloneVertex(vert);
				assert(pr.first->second.pos == _vbt->_vertexTets.size());
				_vbt->_vertexTets.push_back(_vbt->_vertexTets[vert]);
				_vbt->_barycentricWeights.push_back(_vbt->_barycentricWeights[vert]);
			}
		}
	};
	for (int n = (int)_deepPosts.size(), i = 1; i < n; ++i) {
		for(auto &t : _deepPosts[i].quadTriangles){
			for(int j=0; j<3; ++j)
				dupVertex(t.v[j]);
		}
	}
	if (_endPlanes[0].P.X < DBL_MAX) {
		for (auto& t : _endPlanes[0].quadTriangles) {
			for (int j = 0; j < 3; ++j)
				dupVertex(t.v[j]);
		}
	}
	if (_endPlanes[1].P.X < DBL_MAX) {
		for (auto& t : _endPlanes[1].quadTriangles) {
			for (int j = 0; j < 3; ++j)
				dupVertex(t.v[j]);
		}
	}
	// do surfaceCutLine splits of original surface
	auto splitSurface = [&](std::list<std::list<int> >& spl) {
		for (auto& sp : spl) {
			auto vit = sp.begin();
			++vit;
			if (sp.front() == sp.back())  // closed polygon
				sp.push_back(oppositeSideVertices[*vit].pos);  // new end point after first replace
			unsigned int adj = 0xffffffff;
			int vNow = *vit, *tr;
			for (int nt = _mt->numberOfTriangles(), k, j = 0; j < nt; ++j) {
				if (_mt->triangleMaterial(j) == 2)  // must be a deep triangle
					continue;
				tr = _mt->triangleVertices(j);
				for (k = 0; k < 3; ++k) {
					if (tr[k] == sp.front() && tr[(k + 1) % 3] == vNow) {
						adj = (j << 2) + k;
						break;
					}
				}
				if (k < 3)
					break;
			}
			assert(adj != 0xffffffff);
			++vit;
			while (vit != sp.end()) {
				int loopV = -1;
				unsigned int topAdj, botAdj = _mt->triAdjs(adj >> 2)[adj & 3];
				if (_mt->triangleMaterial(botAdj >> 2) == 3) {  // correct the _deepBed entry for the top node if an incision edge
					// from incision convention
					topAdj = _mt->triAdjs(botAdj >> 2)[1];
					if ((topAdj >> 2) + 1 != (botAdj >> 2)) {
						topAdj = _mt->triAdjs(botAdj >> 2)[2];
						assert((topAdj >> 2) + 1 == (botAdj >> 2));
					}
					tr = _mt->triangleVertices(topAdj >> 2);
					int k;
					for (k = 0; k < 2; ++k) {
						auto dbit = _deepBed.find(tr[k]);
						if (dbit != _deepBed.end() && dbit->second.deepMtVertex == vNow) {
							dbit->second.deepMtVertex = oppositeSideVertices[vNow].pos;
							break;
						}
					}
				}
				int newTx = -1, lastTx = -1;
				do {
					adj = _mt->triAdjs(adj >> 2)[adj & 3];
					tr = _mt->triangleVertices(adj >> 2);
					tr[adj & 3] = oppositeSideVertices[vNow].pos;
					if (_mt->triangleMaterial(adj >> 2) != 3) {
						int* txp = _mt->triangleTextures(adj >> 2);
						if (lastTx != txp[adj & 3]) {  // remember possible texture seam crossing
							lastTx = txp[adj & 3];
							newTx = _mt->addTexture();
							_mt->setTexture(newTx, (const float(&)[2]) * _mt->getTexture(txp[adj & 3]));
						}
						txp[adj & 3] = newTx;
					}
					loopV = tr[((adj & 3) + 2) % 3];
					adj = (adj & 0xfffffffc) + (((adj & 3) + 2) % 3);
				} while (loopV != *vit);
				vNow = *vit;
				++vit;
				adj = _mt->triAdjs(adj >> 2)[adj & 3];
			}
		}
	};
	splitSurface(surfacePolyLines);
	// counter clockwise listing of holes must be reversed to split original surface correctly
	for (auto& hpl : _holePolyLines)
		hpl.reverse();
	splitSurface(_holePolyLines);
	// fill in new cut surface triangles
	auto addDeepTriangles = [&](std::vector<matTriangle>& tris) {
		for(auto &t : tris){
			int V[3], T[3];
			// remove any degenrate triangles along straight line of non-duped vertices
			V[0] = t.v[0];
			V[1] = t.v[2];
			V[2] = t.v[1];
			T[0] = t.tex[0];
			T[1] = t.tex[2];
			T[2] = t.tex[1];
			_mt->addTriangle(V, 6, T);
			for (int i = 0; i < 3; ++i) {
				auto ovit = oppositeSideVertices.find(t.v[i]);
				V[i] = ovit->second.pos;
				if (ovit->second.tex < 0) {
					ovit->second.tex = T[i] = _mt->addTexture();
					_mt->setTexture(T[i], (const float (&)[2])*_mt->getTexture(t.tex[i]));
				}
				else
					T[i] = ovit->second.tex;
			}
			_mt->addTriangle(V, 6, T);
		}
	};
	for (int n = (int)_deepPosts.size(), i = 1; i < n; ++i)
		addDeepTriangles(_deepPosts[i].quadTriangles);
	if (_endPlanes[0].P.X < DBL_MAX)
		addDeepTriangles(_endPlanes[0].quadTriangles);
	if (_endPlanes[1].P.X < DBL_MAX)
		addDeepTriangles(_endPlanes[1].quadTriangles);
	if (_mt->findAdjacentTriangles(true))
		return false;  // should throw
	return true;
}

void deepCut::getDeepCutLine(rayTriangleIntersect &top, rayTriangleIntersect &bot) {
	Vec3f V, C = Vec3f(top.intersect.xyz), D = Vec3f(bot.intersect.xyz) - C;
	float len = D.length();
	int n = (int)(D.length() * _vbt->getTetUnitSizeInv());
	D /= (float)n;
	for (int i = 1; i < n - 1; ++i) {
		V = C + D * (float)i;
		int tet;
		Vec3f baryWeight;
		if (!uniqueSpatialTet(V, tet, baryWeight))
			continue;
		top.dcl.push_back(_mt->addVertices(1));
		_mt->setVertexCoordinate(top.dcl.back(), V.xyz);  // texture not set

		// COURT this was commented out.
		//  
		assert(_vbt->_vertexTets.size() == top.dcl.back());
		_vbt->_vertexTets.push_back(tet);
		_vbt->_barycentricWeights.push_back(baryWeight);
	}
}

bool deepCut::uniqueSpatialTet(const Vec3f pos, int& tet, Vec3f& baryWeight) {

//	std::chrono::time_point<std::chrono::system_clock> start, end;
//	start = std::chrono::system_clock::now();

	tet = -1;  // only found once so no write contention.
	auto tetInside = [&](int tetid, Vec3f& bw) ->bool {
		boundingBox<float> bb;
		bb.Empty_Box();
		Vec3f vp[4];
		const int* tn = _vbt->tetNodes(tetid);  //  > _tetNodes[i].data();
		for (int j = 0; j < 4; ++j) {
			vp[j] = _vbt->_nodeSpatialCoords[tn[j]];
			bb.Enlarge_To_Include_Point(vp[j].xyz);
		}
		if (bb.Outside(pos.xyz))
			return false;
		Mat3x3f M;
		M.Initialize_With_Column_Vectors(vp[1] - vp[0], vp[2] - vp[0], vp[3] - vp[0]);
		Vec3f R = M.Robust_Solve_Linear_System(pos - vp[0]);
		if (R[0] <= 0.0f || R[1] <= 0.0f || R[2] <= 0.0f || R[0] >= 1.0f || R[1] >= 1.0f || R[2] >= 1.0f || R[0] + R[1] + R[2] >= 1.0f)
			return false;
		bw = R;
		return true;
	};
	for (int i = 0; i < _vbt->_nMegatets; ++i) {  // in this new megatet environment check these largest unique tets first
		if (tetInside(i, baryWeight)){
			tet = i;
			return true;
		}
	}
//	for (int i = _vbt->firstInteriorTet(); i < _vbt->tetNumber(); ++i) {  // Now check the unique interior tets 
	tbb::parallel_for(tbb::blocked_range<std::size_t>(_vbt->firstInteriorTet(), _vbt->tetNumber()), [&](tbb::blocked_range<size_t> r) {
		for (int i = r.begin(); i != r.end(); ++i) {
			Vec3f bw;
			if (tetInside(i, bw)) {
				tet = i;
				baryWeight = bw;
				break;
			}
		}
		if (tet > -1)
			oneapi::tbb::task_group_context().cancel_group_execution();
	});

//	end = std::chrono::system_clock::now();
//	std::chrono::duration<double> elapsed_seconds = end - start;
//	std::cout << "interiorSpatialTet() took " << elapsed_seconds.count() << "\n";

	// worst case seems to be 0.008 seconds without tbb. tbb version worst 0.001 but some 15x faster on my desktop
	if (tet < 0)
		return false;
	return true;
}

bool deepCut::connectOpenEnd(int postNum) {
	auto &post = _deepPosts[postNum].triIntersects;
	endPlane *ep = postNum < 1 ? &_endPlanes[0] : &_endPlanes[1];
	auto makePlane = [&](int topIndex) {
		ep->P = post[topIndex].intersect;
		ep->V = post.front().intersect - ep->P;
		ep->VlengthSqInv = 1.0/ep->V.length2();
		if(postNum > 0)
			ep->U = post.front().intersect - _deepPosts[postNum - 1].triIntersects.front().intersect;
		else
			ep->U = _deepPosts[1].triIntersects.front().intersect - post.front().intersect;
		ep->N = ep->V ^ ep->U;
		ep->N.normalize();
		ep->d = ep->N * ep->P;
		ep->U = ep->N ^ ep->V;
	};
	makePlane(post.size() - 1);
	int interPostIdx, maxIdx;  // start with maximum index from interpost connection
	if (postNum < 1) {
		for (maxIdx = 1; maxIdx < _deepPosts[1].triIntersects.size(); maxIdx += 2) {
			if (_deepPosts[1].triIntersects[maxIdx].scl.rtiPostTo < 1) {
				maxIdx = _deepPosts[1].triIntersects[maxIdx].scl.rtiIndexTo;
				break;
			}
		}
	}
	else {
		for (maxIdx = 1; maxIdx < post.size(); maxIdx += 2) {
			if (post[maxIdx].scl.rtiPostTo < _deepPosts.size() - 1)
				break;
		}
	}
	interPostIdx = maxIdx;
	if (postNum < 1) {  // front
		int k, j, i;
		for (k = 0; k < post.size(); k += 2) {
			for (i = 1; i < post.size(); i += 2) {
				double lowV;
				if (surfacePath(post[i], post[k], false, lowV) < DBL_MAX) {
					if (k < interPostIdx || i <= interPostIdx) {
						post[i].scl.rtiPostTo = postNum;
						post[i].scl.rtiIndexTo = k;
						post[i].scl.lowestV = lowV;
						if (k < 1 && i > maxIdx)
							maxIdx = i;
						break;
					}
				}
			}
			if (k < 1) {
				if (i > post.size())
					return false;
			}
		}
		makePlane(maxIdx);
	}
	else {  // back end
		int k, j, i;
		for (k = 0; k < (k < 1 ? post.size() : maxIdx); k += 2) {
			for (i = 1; i < (k < 1 ? post.size() : maxIdx + 1); i += 2) {
				double lowV;
				if (surfacePath(post[k], post[i], false, lowV) < DBL_MAX) {
					if (k < interPostIdx || i <= interPostIdx) {
						post[k].scl.rtiPostTo = postNum;
						post[k].scl.rtiIndexTo = i;
						post[k].scl.lowestV = lowV;
						if (k<1 && i > maxIdx)
							maxIdx = i;
						break;

					}
				}
			}
			if (k < 1 && i > post.size())
					return false;
		}
		makePlane(maxIdx);
	}
	return true;
}

bool deepCut::deepCutQuad(int postNum) {
	// get CW quad border
	auto& pti0 = _deepPosts[postNum - 1].triIntersects;
	auto& pti1 = _deepPosts[postNum].triIntersects;
	std::list<int> poly, tmp;
	// keep poly and its UV in sync;
	std::list<Vec2d> polyUV, tmpUV;
	poly.assign(pti0[0].scl.deepVertsTris.begin(), pti0[0].scl.deepVertsTris.end());
	polyUV.assign(pti0[0].scl.deepUVs.begin(), pti0[0].scl.deepUVs.end());
	const bilinearPatch& bl = _deepPosts[postNum].bl;
	Vec2d postTx(0.0, 1.0f);
	polyUV.front() = postTx;
	postTx.X = 1.0;
	polyUV.back() = postTx;
	bool lastDeep = false;  // openEnd = false, 
	int interPostIdx = -1;
//	if (!pti1[0].scl.deepVertsTris.empty() && postNum == _deepPosts.size() - 1) // open end
//		openEnd = true;
	// find first interpost return
	for (int i = 1; i < pti1.size(); i += 2) {
		if (pti1[i].scl.rtiPostTo == postNum - 1) {
			interPostIdx = i;
			break;
		}
	}
	if (interPostIdx < 0 || pti1[interPostIdx].scl.rtiIndexTo < 0)
		throw(std::logic_error("Topological connection error in deepQuadCut()."));
	auto getPostV = [&](int vertex, bool post0) ->double {
		float* vp = _mt->vertexCoordinate(vertex);
		Vec3d V((double)vp[0], (double)vp[1], (double)vp[2]);
		if (post0)
			return (V - bl.P00).length() / bl.e00.length();
		else
			return (V - bl.P10).length() / bl.e11.length();
	};
	for (int i = 0; i != interPostIdx; ) {
		if (i & 1){
			assert(!pti1[i].scl.deepVertsTris.empty());
			tmp.assign(pti1[i].scl.deepVertsTris.begin(), pti1[i].scl.deepVertsTris.end());
			tmpUV.assign(pti1[i].scl.deepUVs.begin(), pti1[i].scl.deepUVs.end());
			postTx.set(1.0, getPostV(tmp.front(), false));
			tmpUV.front() = postTx;
			postTx.Y = getPostV(tmp.back(), false);
			tmpUV.back() = postTx;
			if (!lastDeep) {
				assert(poly.back() == tmp.front());
				poly.pop_back();
				polyUV.pop_back();
			}
			poly.splice(poly.end(), tmp);
			polyUV.splice(polyUV.end(), tmpUV);
			lastDeep = false;
			i = pti1[i].scl.rtiIndexTo;
		}
		else {
			if (lastDeep) {
				poly.push_back(pti1[i].deepVert);
				postTx.set(1.0, getPostV(pti1[i].deepVert, false));
				polyUV.push_back(postTx);
			}
			if (i > interPostIdx) {
				--i;
				tmp.assign(pti1[i].dcl.begin(), pti1[i].dcl.end());
				tmp.reverse();
			}
			else {
				tmp.assign(pti1[i].dcl.begin(), pti1[i].dcl.end());
				++i;
			}
			for (auto v : tmp) {
				postTx.set(1.0, getPostV(v, false));
				polyUV.push_back(postTx);
			}
			poly.splice(poly.end(), tmp);
			lastDeep = true;
		}
	}
	tmp.assign(pti1[interPostIdx].scl.deepVertsTris.begin(), pti1[interPostIdx].scl.deepVertsTris.end());
	tmpUV.assign(pti1[interPostIdx].scl.deepUVs.begin(), pti1[interPostIdx].scl.deepUVs.end());
	postTx.set(1.0, getPostV(tmp.front(), false));
	tmpUV.front() = postTx;
	postTx.set(0.0, getPostV(tmp.back(), true));
	tmpUV.back() = postTx;
	interPostIdx = pti1[interPostIdx].scl.rtiIndexTo;
	if (!lastDeep) {
		assert(poly.back() == tmp.front());
		poly.pop_back();
		polyUV.pop_back();
	}
	poly.splice(poly.end(), tmp);
	polyUV.splice(polyUV.end(), tmpUV);
	lastDeep = false;
	for (int i = interPostIdx; i != 0; ) {
		if ((i & 1) < 1) {
			assert(!pti0[i].scl.deepVertsTris.empty());
			tmp.assign(pti0[i].scl.deepVertsTris.begin(), pti0[i].scl.deepVertsTris.end());
			tmpUV.assign(pti0[i].scl.deepUVs.begin(), pti0[i].scl.deepUVs.end());
			postTx.X = 0.0;
			postTx.Y = getPostV(tmp.front(), true);
			tmpUV.front() = postTx;
			postTx.Y = getPostV(tmp.back(), true);
			tmpUV.back() = postTx;
			if (!lastDeep) {
				assert(poly.back() == tmp.front());
				poly.pop_back();
				polyUV.pop_back();
			}
			poly.splice(poly.end(), tmp);
			polyUV.splice(polyUV.end(), tmpUV);
			lastDeep = false;
			i = pti0[i].scl.rtiIndexTo;
		}
		else {
			if (lastDeep) {
				poly.push_back(pti0[i].deepVert);
				postTx.X = 0.0;
				postTx.Y = getPostV(pti0[i].deepVert, true);
				polyUV.push_back(postTx);
			}
			if (i <= interPostIdx) {
				--i;
				tmp.assign(pti0[i].dcl.begin(), pti0[i].dcl.end());
				tmp.reverse();
			}
			else {
				tmp.assign(pti0[i].dcl.begin(), pti0[i].dcl.end());
				++i;
			}
			for (auto v : tmp) {
				postTx.set(0.0, getPostV(v, true));
				polyUV.push_back(postTx);
			}
			poly.splice(poly.end(), tmp);
			lastDeep = true;
		}
	}
	if (poly.back() == poly.front()) {
		poly.pop_back();
		polyUV.pop_back();
	}
	poly.reverse(); // to make CCW
	polyUV.reverse();
	makePolygonTriangles(poly, polyUV, &bl, nullptr, _deepPosts[postNum].quadTriangles);
	return true;
}

bool deepCut::deepCutEndPlane(int endPlane) {
	_endPlanes[endPlane].quadTriangles.clear();
	std::list<int> polyVerts, tmp;
	std::list<Vec2d> polyUV, tmpUV;
	auto vCoord = [&](int vertex) ->Vec3d {
		float* fp = _mt->vertexCoordinate(vertex);
		Vec3d ret(fp[0], fp[1], fp[2]);
		return ret;
	};
	// unlike an interpost cut which generates only one polygon, open ends can generate more than one
	std::vector<rayTriangleIntersect>* tip;
	if (endPlane)  // open end on last point
		tip = &_deepPosts.back().triIntersects;
	else
		tip = &_deepPosts.front().triIntersects;
	std::vector<bool> iPoints(tip->size(), true);
	double postLen = (tip->back().intersect - tip->front().intersect).length();
	auto dclPath = [&](int idx, int start) ->bool {
		if ((*tip)[idx].scl.rtiIndexTo < 0 || (*tip)[idx].scl.rtiPostTo != (*tip)[idx].postNum)  // end into solid edge
			return true;
		if (endPlane < 1) {
			if ((*tip)[idx].scl.rtiIndexTo < idx)
				return idx != start;
		}
		else {
			if ((*tip)[idx].scl.rtiIndexTo > idx)
				return idx != start;
		}
		return false;
	};
	auto cutPolygon = [&](int start) {
		iPoints[start] = false;
		polyVerts.assign((*tip)[start].scl.deepVertsTris.begin(), (*tip)[start].scl.deepVertsTris.end());
		polyUV.assign((*tip)[start].scl.deepUVs.begin(), (*tip)[start].scl.deepUVs.end());
		Vec2d postTx(0.0, 0.0);
		postTx.Y = (vCoord(polyVerts.front()) - _endPlanes[endPlane].P).length() / postLen;
		polyUV.front() = postTx;
		postTx.Y = (vCoord(polyVerts.back()) - _endPlanes[endPlane].P).length() / postLen;
		polyUV.back() = postTx;
		if ((*tip)[start].scl.rtiPostTo != (*tip)[start].postNum)
			throw(std::logic_error("Can't make polygon connection in end plane deep cut."));
		int idx = (*tip)[start].scl.rtiIndexTo;
		iPoints[idx] = false;
		while (idx != start) { // follow polygon
			if (dclPath(idx, start)) {  // end into solid edge
				if (endPlane) {
					--idx;
					tmp.assign((*tip)[idx].dcl.begin(), (*tip)[idx].dcl.end());
				}
				else {
					tmp.assign((*tip)[idx].dcl.begin(), (*tip)[idx].dcl.end());
					++idx;
				}
				iPoints[idx] = false;
				if (endPlane)
					tmp.reverse();
				for (auto& v : tmp) {
					postTx.Y = (vCoord(v) - _endPlanes[endPlane].P).length() / postLen;
					polyUV.push_back(postTx);
				}
				polyVerts.splice(polyVerts.end(), tmp);
			}
			else {
				tmp.assign((*tip)[idx].scl.deepVertsTris.begin(), (*tip)[idx].scl.deepVertsTris.end());
				tmpUV.assign((*tip)[idx].scl.deepUVs.begin(), (*tip)[idx].scl.deepUVs.end());
				postTx.Y = ((*tip)[idx].intersect - _endPlanes[endPlane].P).length() / postLen;
				tmpUV.front() = postTx;
				if ((*tip)[idx].scl.rtiPostTo != (*tip)[idx].postNum)
					throw(std::logic_error("Can't make polygon connection in end plane deep cut."));
				idx = (*tip)[idx].scl.rtiIndexTo;
				iPoints[idx] = false;
				postTx.Y = ((*tip)[idx].intersect - _endPlanes[endPlane].P).length() / postLen;
				tmpUV.back() = postTx;
				polyVerts.splice(polyVerts.end(), tmp);
				polyUV.splice(polyUV.end(), tmpUV);
			}
		}
		assert(polyVerts.back() != polyVerts.front());
		polyVerts.reverse(); // to make CCW
		polyUV.reverse();
		makePolygonTriangles(polyVerts, polyUV, nullptr, &_endPlanes[endPlane], _endPlanes[endPlane].quadTriangles);
	};
	if (endPlane) {  // open end on last point
		for (int i = 0; i < tip->size(); i += 2) {
			if (!iPoints[i] || (*tip)[i].scl.rtiIndexTo < 0)  // not previously touched and a valid start
				continue;
			cutPolygon(i);
		}
	}
	else{  // open first point
		for (int i = tip->size() - 1; i > -1;  i -= 2) {
			if (!iPoints[i] || (*tip)[i].scl.rtiIndexTo < 0)
				continue;
			cutPolygon(i);
		}
	}
	return true;
}

void deepCut::makePolygonTriangles(const std::list<int> &polyVerts, const std::list<Vec2d> &polyUV, const bilinearPatch *blp, const endPlane *ep, std::vector<matTriangle> &polyTriangles){
	// polyVerts contains polygonVertices and polyUV their patch parametric locations
	std::vector<std::pair<int, Vec2d> > deepOuterPolygon;
	deepOuterPolygon.reserve(polyVerts.size());
	// put in missing uv coords for this quad
	Vec2d minC(DBL_MAX, DBL_MAX), maxC(-DBL_MAX, -DBL_MAX);
	std::vector<Vec2d> bUv;
	bUv.assign(polyUV.begin(), polyUV.end());
	auto uvit = bUv.begin();
	for (auto pit = polyVerts.begin(); pit != polyVerts.end(); ++pit) {
		deepOuterPolygon.push_back(std::make_pair(*pit, *uvit));
		if (minC.X > uvit->X)
			minC.X = uvit->X;
		if (maxC.X < uvit->X)
			maxC.X = uvit->X;
		if (minC.Y > uvit->Y)
			minC.Y = uvit->Y;
		if (maxC.Y < uvit->Y)
			maxC.Y = uvit->Y;
		++uvit;
	}
	std::vector<int> bTx, holeTx;
	bTx.reserve(bUv.size());
	float tex[2];
	for (int n = bUv.size(), i = 0; i < n; ++i) {
		bTx.push_back(_mt->addTexture());
		tex[0] = (float)bUv[i].X;
		tex[1] = (float)bUv[i].Y;
		_mt->setTexture(bTx.back(), tex);
	}
	int holesSize = 0;
	std::list< vertUvPolygon > holes;
	findCutInteriorHoles(blp, ep, bUv, deepOuterPolygon, holes);
	for (auto hole : holes)
		holesSize += hole.vertices.size();
	holeTx.reserve(holesSize);
	for (auto hole : holes) {
		for(auto &huv : hole.uvs){
			holeTx.push_back(_mt->addTexture());
			float tx[2] = { (float)huv.X, (float)huv.Y };
			_mt->setTexture(holeTx.back(), tx);
		}
	}
	// do Delaunay triangulation of border with interior points
	int nOuterPolygon = (int)bUv.size();
	// recursively add internal Delaunay vertices to mtFace
	double incr = (maxC.X - minC.X) + (maxC.Y - minC.Y);
	incr /= (polyVerts.size() * 0.33);
	int n = (int)((maxC.X - minC.X) / incr), m = (int)((maxC.Y - minC.Y) / incr);
	Vec2d V2;
	insidePolygon ip;
	struct interiorPoint {
		int newV;
		int tx;
		Vec2d uv;
	}intPt;
	std::vector<interiorPoint> ceiVec;
	ceiVec.reserve(m*n);
	for (int i = 1; i < n; ++i) {
		for (int j = 1; j < m; ++j) {
			V2.X = minC.X + (incr + 5.0e-18) * i;  // Delaunay routines don't like equal edge lengths
			V2.Y = minC.Y + (incr + 3.0e-18) * j;
			if (ip.insidePolygon2d(V2, bUv)) {
				auto hit = holes.begin();
				while (hit != holes.end()) {
					if (ip.insidePolygon2d(V2, hit->uvs))
						break;
					++hit;
				}
				if (hit != holes.end())
					continue;
				Vec3d Nd;
				if (blp != nullptr)
					Nd = blp->P00 + blp->e10 * V2.X + blp->e00 * V2.Y + (blp->e11 - blp->e00) * (V2.X * V2.Y);
				else 
					Nd = ep->P + ep->U * V2.X + ep->V * V2.Y;
				int tet;
				Vec3f N(Nd.xyz), baryWeight;
				if (!uniqueSpatialTet(Vec3f(N.xyz), tet, baryWeight))
					continue;
				intPt.newV = _mt->addVertices(1);
				_mt->setVertexCoordinate(intPt.newV, N.xyz);  // texture not set
				assert(_vbt->_vertexTets.size() == intPt.newV);
				_vbt->_vertexTets.push_back(tet);
				_vbt->_barycentricWeights.push_back(baryWeight);
				intPt.uv = V2;
				intPt.tx = _mt->addTexture();
				tex[0] = (float)V2.X;
				tex[1] = (float)V2.Y;
				_mt->setTexture(intPt.tx, tex);
				ceiVec.push_back(intPt);
			}
		}
	}
	// Using new constrained Delaunay triangulation library. See top of file
	typedef CDT::V2d<double> V2d;
	typedef CDT::Triangulation<double> Triangulation;
	Triangulation cdt = Triangulation(CDT::FindingClosestPoint::ClosestRandom, 10);
	std::vector<V2d> DelPts;
	DelPts.reserve(nOuterPolygon + holesSize + ceiVec.size());
	for (int i = 0; i < nOuterPolygon; ++i)
		DelPts.push_back(V2d().make(bUv[i].X, bUv[i].Y));
	for (auto& h : holes) {
		for(auto &uv : h.uvs)
			DelPts.push_back(V2d().make(uv.X, uv.Y));
	}
	for (int r = ceiVec.size(), i = 0; i < r; ++i)
		DelPts.push_back(V2d().make(ceiVec[i].uv.X, ceiVec[i].uv.Y));
	cdt.insertVertices(DelPts);
	std::vector<CDT::Edge> edges;
	edges.reserve(nOuterPolygon + holesSize);
	for(int i=1; i< nOuterPolygon; ++i)
		edges.push_back(CDT::Edge(i-1, i));
	edges.push_back(CDT::Edge(nOuterPolygon - 1, 0));
	int edgeSize = nOuterPolygon;
	for (auto& h : holes) {
		for (int i = 1; i < h.uvs.size(); ++i)
			edges.push_back(CDT::Edge(i - 1 + edgeSize, i + edgeSize));
		edges.push_back(CDT::Edge(h.uvs.size() - 1 + edgeSize, edgeSize));
		edgeSize += h.uvs.size();
	}
//	for (auto hole : holeUvs) {
//		for (int i = 1; i < hole.size(); ++i)
//			edges.push_back(CDT::Edge(i - 1 + edgeSize, i + edgeSize));
//		edges.push_back(CDT::Edge(hole.size() - 1 + edgeSize, edgeSize));
//		edgeSize += hole.size();
//	}
	cdt.insertEdges(edges);
	cdt.eraseOuterTrianglesAndHoles();
	polyTriangles.reserve(cdt.triangles.size());
	for (size_t n = cdt.triangles.size(), i = 0; i < n; ++i) {
//		std::vector<int> tri;
		matTriangle tri;
		tri.material = 6;
		for (int j = 0; j < 3; ++j) {
			int v;
			if ((v = cdt.triangles[i].vertices[j]) < nOuterPolygon) {
				tri.v[2 - j] = deepOuterPolygon[v].first;
				tri.tex[2 - j] = bTx[v];
			}
			else if (v < nOuterPolygon + holesSize) {
				v -= nOuterPolygon;
				tri.tex[2 - j] = holeTx[v];
				for (auto hole : holes) {
					if (v < hole.uvs.size()) {
						tri.v[2 - j] = hole.vertices[v];
						break;
					}
					else
						v -= hole.vertices.size();
				}
			}
			else {
				tri.v[2 - j] = ceiVec[v - nOuterPolygon - holesSize].newV;
				tri.tex[2 - j] = ceiVec[v - nOuterPolygon - holesSize].tx;
			}
		}
		polyTriangles.push_back(tri);
	}
}

bool deepCut::updateDeepSpatialCoordinates()
{
	// adds any new vertices created
//	_deepXyz.reserve(_mt->numberOfVertices()); // , Vec3d(DBL_MAX, 0.0f, 0.0f)
	for (int n = _mt->numberOfVertices(), i = _deepXyz.size(); i < n; ++i) {
		float v[3];
		_mt->getVertexCoordinate(i, v);
		_deepXyz.push_back(Vec3d(v));
		auto dbit = _deepBed.find(i);
		if(dbit == _deepBed.end())
			continue;
		if (dbit->second.deepMtVertex > -1) {  // mark material 3-4 vertices invalid
			std::vector<materialTriangles::neighborNode> nei;
			_mt->getNeighbors(dbit->second.deepMtVertex, nei);
			auto nit = nei.begin();
			bool has4 = false, has5 = false;
			while (nit != nei.end()) {
				if (_mt->triangleMaterial(nit->triangle) == 5)
					has5 = true;
				else if (_mt->triangleMaterial(nit->triangle) == 4)
					has4 = true;
				else
					;
				++nit;
			}
			if (has4 && !has5) {  // invalid vertex for deep cut. Must be on a flap bottom.
				_deepXyz[dbit->first].X = DBL_MAX;
				continue;
			}
		}
		Vec3f bw;
		int botTet = deepPointTetWeight(dbit, bw);
		if (botTet < 0) {
			_deepXyz.clear();  // COURT - this is bad.  Input deepBed needs to be recomputed.
			return false;
		}
		else {
			Vec3f v;
			_vbt->getBarycentricTetPosition(botTet, bw, v);
			_deepXyz[dbit->first] = Vec3d(v.xyz);
		}
	}
	return true;
}

bool deepCut::getDeepSpatialCoordinates()
{
	_deepXyz.clear();
	_deepXyz.reserve(_mt->numberOfVertices()); // , Vec3d(DBL_MAX, 0.0f, 0.0f)
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i) {
		float v[3];
		_mt->getVertexCoordinate(i, v);
		_deepXyz.push_back(Vec3d(v));
	}
	for (auto dbit = _deepBed.begin(); dbit != _deepBed.end(); ++dbit){  // avoid logN hash searches
		if (dbit->second.deepMtVertex > -1){  // mark material 3-4 vertices invalid
			std::vector<materialTriangles::neighborNode> nei;
			_mt->getNeighbors(dbit->second.deepMtVertex, nei);
			auto nit = nei.begin();
			bool has4 = false, has5 = false;
			while (nit != nei.end()) {
				if (_mt->triangleMaterial(nit->triangle) == 5)
					has5 = true;
				else if (_mt->triangleMaterial(nit->triangle) == 4)
					has4 = true;
				else
					;
				++nit;
			}
			if (has4 && !has5) {  // invalid vertex for deep cut. Must be on a flap bottom.
				_deepXyz[dbit->first].X = DBL_MAX;
				continue;
			}
			Vec3f v;
			_vbt->getBarycentricTetPosition(_vbt->getVertexTetrahedron(dbit->second.deepMtVertex), *_vbt->getVertexWeight(dbit->second.deepMtVertex), v);
			_deepXyz[dbit->first] = Vec3d(v.xyz);
		}
		else {
			Vec3f bw;
			int botTet = deepPointTetWeight(dbit, bw);
			if (botTet < 0) {  // COURt if there are a lot of these or any occur over an important area of the model should recompute deep bed.
				std::cout << "Deep bed point at vertex " << dbit->first << " with deep vertex " << dbit->second.deepMtVertex << " not connected to a tet.\n";
				_deepXyz[dbit->first].X = DBL_MAX;  // mark invalid for deep cut.
			}
			else {
				Vec3f v;
				_vbt->getBarycentricTetPosition(botTet, bw, v);
				_deepXyz[dbit->first] = Vec3d(v.xyz);
			}
		}
	}
	boundingBox<double> bb;
	bb.Empty_Box();
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i){
		if (_deepXyz[i].X < DBL_MAX)  // don't include invalid vertex
			bb.Enlarge_To_Include_Point(_deepXyz[i].xyz);
	}
	Vec3d vn, vx;
	bb.Maximum_Corner(vx.xyz);
	bb.Minimum_Corner(vn.xyz);
	_maxSceneSize = (float)((vx - vn).length());
	return true;
}

void deepCut::getDeepPosts(std::vector<Vec3f>& xyz, std::vector<Vec3f>& nrm) {
	xyz.clear();
	nrm.clear();
	xyz.reserve(_deepPosts.size());
	nrm.reserve(_deepPosts.size());
	for (int n = _deepPosts.size(), i = 0; i < n; ++i) {
		xyz.push_back(Vec3f(_deepPosts[i].triIntersects[0].intersect.xyz));
		nrm.push_back(Vec3f(_deepPosts[i].triIntersects.back().intersect.xyz));
		nrm.back() -= xyz.back();
	}
}

int deepCut::preventPreviousCrossover(const int postNum) {
	// assumes all posts have previous surface connections
	assert(postNum > 0);  // already checked in calling routine
	deepPost* dp = &_deepPosts[postNum];
	Vec3d P = dp->triIntersects[0].intersect, N;;
	int k;
	for (k = 1; k < dp->triIntersects.size(); k += 2) {
		if (dp->triIntersects[k].scl.rtiPostTo == postNum - 1)
			break;
	}
	if (k > dp->triIntersects.size()) {
		return postNum;
	}
	N = dp->triIntersects[k].intersect;
	N -= P;
	bool normalChanged = false;
	for (int i = 0; i < postNum; ++i) {
		if (postNum < 2) {
			int j = dp->triIntersects[k].scl.rtiIndexTo;
			Vec3d Q = _deepPosts[i].triIntersects[0].intersect, M = _deepPosts[i].triIntersects[j].intersect;
			M -= Q;
			double MdN = M * N, mn = M.length() * N.length();
			if (MdN/mn > 0.998) // lines parallel
				continue;
			Mat2x2d Mat(M*M, MdN, -MdN, -(N * N));
			P -= Q;
			Vec2d R = Mat.Robust_Solve_Linear_System(Vec2d(P * M, P * N));
			double sep = (Q + M * R.X - P - N * R.Y).length2();
			if (sep/mn > 0.001) // no crossover
				continue;
			return postNum;
		}
		else {
			if (i < 1)
				continue;
			// get this active patch
			int j;
			for (j = 1; j < _deepPosts[i].triIntersects.size(); ++j) {
				if (_deepPosts[i].triIntersects[j].scl.rtiPostTo == i - 1)
					break;
			}
			if (j > _deepPosts[i].triIntersects.size())
				throw(std::logic_error("Previous processing should have removed this case."));
			bilinearPatch bl;
			int idx = _deepPosts[i].triIntersects[j].scl.rtiIndexTo;
			makeBilinearPatch(_deepPosts[i - 1].triIntersects[idx].intersect, _deepPosts[i].triIntersects[j].intersect, _deepPosts[i - 1].triIntersects.front().intersect, _deepPosts[i].triIntersects.front().intersect, bl);
			double rP[2];
			Vec2d fP[2];
			if (bilinearRayIntersection(P, N, bl, rP, fP)) {
				int k;
				for (k = 0; k < 2; ++k) {
					if (rP[k] < 0.0 || rP[k] > 1.0)
						continue;
					if (fP[k].Y > _deepPosts[i].triIntersects[j].scl.lowestV)
						return postNum;
				}
			}
		}
	}
	return 0;
}

bool deepCut::inputCorrectFence(fence* fp, FacialFlapsGui* ffg) {
	// interactive deep cut builder.  On successful exit deep posts have interpost connections set up.
	_deepPosts.clear();
	if (fp->numberOfPosts() < 2) {
		char str[200];
		sprintf(str, "You need at least 2 posts to create a deep cut.");
		ffg->sendUserMessage(str, "Invalid deep cut-");
		return false;
	}
	std::vector<Vec3f> positions, normals;
	std::vector<int> triangles;
	std::vector<float> uv;
	bool edgeStart, edgeEnd, startOpen, endOpen;
	int n = fp->getPostData(positions, normals, triangles, uv, edgeStart, edgeEnd, startOpen, endOpen);
	for (int i = 0; i < n; ++i) {
		bool closedEnd = false;
		if (i < 1)
			closedEnd = !startOpen;
		if (i == n-1)
			closedEnd = !endOpen;
		if (addDeepPost(triangles[i], (const float(&)[2])uv[i * 2], -Vec3d(normals[i]), closedEnd) < 0) {
			char str[200];
			sprintf(str, "Post number %d needs direction adjustment.", i + 1);
			ffg->sendUserMessage(str, "Please correct deep cut-");
			return false;
		}
	}
	for (int i = 1; i < n; ++i) {
		if (!topConnectToPreviousPost(i)) {
			char str[200];
			sprintf(str, "Post number %d has no top side connection to previous post.\nDelete it and try again", i + 1);
			ffg->sendUserMessage(str, "Please correct deep cut-");
			return false;
		}
	}
	for (int i = 1; i < n; ++i) {
		if (!deepConnectToPreviousPost(i)) {
			char str[200];
			sprintf(str, "Post number %d has no bottom connection to previous post.\nAdjust its post direction or previous post direction.", i + 1);
			ffg->sendUserMessage(str, "Please correct deep cut-");
			return false;
		}
	}
	if (n > 1) {  // no corrections necessary til 2
		for (int i = 1; i < n; ++i) {
			if (preventPreviousCrossover(i) > 0) {
				char str[200];
				sprintf(str, "Solid post line %d intersects a previous cut.\nPlease adjust it's direction.", i + 1);
				ffg->sendUserMessage(str, "Please correct deep cut-");
				return false;
			}
		}
	}
	return true;
}

int deepCut::addDeepPost(const int triangle, const float (&uv)[2], const Vec3d& rayDirection, bool closedEnd) {
	if (_deepPosts.empty()) {
		if (!getDeepSpatialCoordinates()) {
			return -1;
		}
	}
	int ret = (int)_deepPosts.size();
	_deepPosts.push_back(deepPost());
	_deepPosts.back().closedEnd = closedEnd;
	_deepPosts.back().rayDirection = rayDirection;
	int* tr = _mt->triangleVertices(triangle);
	Vec3d rayStart;
	if (_deepXyz[tr[0]].X == DBL_MAX || _deepXyz[tr[1]].X == DBL_MAX || _deepXyz[tr[2]].X == DBL_MAX) {  // invalid top triangle
		rayStart = Vec3d((const float (&)[3])*_mt->vertexCoordinate(tr[0])) * (1.0 - uv[0] - uv[1]);
		rayStart += Vec3d((const float(&)[3])*_mt->vertexCoordinate(tr[1])) * uv[0];
		rayStart += Vec3d((const float(&)[3])*_mt->vertexCoordinate(tr[2])) * uv[1];
	}
	else {
		rayStart = _deepXyz[tr[0]] * (1.0 - uv[0] - uv[1]);
		rayStart += _deepXyz[tr[1]] * uv[0];
		rayStart += _deepXyz[tr[2]] * uv[1];
	}
	if (!rayIntersectMaterialTriangles(rayStart, rayDirection, _deepPosts.back().triIntersects))
		return -1;
	for (auto& ti : _deepPosts.back().triIntersects)
		ti.postNum = ret;
	return ret;
}

bool deepCut::topConnectToPreviousPost(int postNum) {
	assert(postNum > 0);
	double lowV;
	if (surfacePath(_deepPosts[postNum - 1].triIntersects[0], _deepPosts[postNum].triIntersects[0], false, lowV) == DBL_MAX)
		return false;
	auto fp = &_deepPosts[postNum - 1].triIntersects[0].scl;
	fp->rtiPostTo = postNum;
	fp->rtiIndexTo = 0;
	fp->lowestV = lowV;
//	fp->deepUVs.clear();  // clear saved topo search data
//	fp->deepVertsTris.clear();
	return true;
}

bool deepCut::deepConnectToPreviousPost(int postNum) {
	// find shortest return path
	double pathLen = DBL_MAX, minPath = DBL_MAX, lowV;
	int minJ = 10000;
	int i, j, n = _deepPosts[postNum].triIntersects.size();
	for (i = 1; i < n; i += 2) {
		for (j = 1; j < _deepPosts[postNum - 1].triIntersects.size(); j += 2) {
			pathLen = surfacePath(_deepPosts[postNum].triIntersects[i], _deepPosts[postNum - 1].triIntersects[j], false, lowV);
//			_deepPosts[postNum].triIntersects[i].scl.deepUVs.clear();  // clear topo search data
//			_deepPosts[postNum].triIntersects[i].scl.deepVertsTris.clear();
			if (pathLen < DBL_MAX)
				break;
		}
		if (pathLen < DBL_MAX)
			break;
	}
	if (pathLen < DBL_MAX) {
		auto fp = &_deepPosts[postNum].triIntersects[i].scl;
		fp->rtiPostTo = postNum - 1;
		fp->rtiIndexTo = j;  //  minJ;
		fp->lowestV = lowV;  //  minLowV;
//		fp->deepUVs.clear();  // clear saved topo search data
//		fp->deepVertsTris.clear();
		return true;
	}
	return false;
}

/* bool deepCut::connectToPreviousPost(int postNum) {

	assert(false);  // no longer used

	assert(postNum > 0);
	double lowV;
	if (surfacePath(_deepPosts[postNum - 1].triIntersects[0], _deepPosts[postNum].triIntersects[0], false, lowV) == DBL_MAX)
		return false;
	auto fp = &_deepPosts[postNum - 1].triIntersects[0].scl;
	fp->rtiPostTo = postNum;
	fp->rtiIndexTo = 0;
	fp->lowestV = lowV;
	fp->deepUVs.clear();
	fp->deepVertsTris.clear();
	// find shortest return path
	double pathLen = DBL_MAX, minPath = DBL_MAX;
	int minJ = 10000;
	int i, j, n = _deepPosts[postNum].triIntersects.size();
	for (i = 1; i < n; i += 2) {
		for (j = 1; j < _deepPosts[postNum - 1].triIntersects.size(); j += 2) {
			pathLen = surfacePath(_deepPosts[postNum].triIntersects[i], _deepPosts[postNum - 1].triIntersects[j], false, lowV);
			if (pathLen < DBL_MAX)
				break;
		}
		if (pathLen < DBL_MAX)
			break;
	}
	if (pathLen < DBL_MAX){
		auto fp = &_deepPosts[postNum].triIntersects[i].scl;
		fp->rtiPostTo = postNum - 1;
		fp->rtiIndexTo = j;  //  minJ;
		fp->lowestV = lowV;  //  minLowV;
		fp->deepUVs.clear();
		fp->deepVertsTris.clear();
		return true;
	}
	return false;
} */

bool deepCut::rayIntersectMaterialTriangles(const Vec3d& rayStart, const Vec3d& rayDirection, std::vector<rayTriangleIntersect>& intersects) {
	std::multimap<double, rayTriangleIntersect> rtiMap;
	boundingBox<double> bb, tbb;
	bb.Empty_Box();
	bb.Enlarge_To_Include_Point(rayStart.xyz);
	bb.Enlarge_To_Include_Point((rayStart + rayDirection * _maxSceneSize).xyz);
	Vec3d P, T[3], N;
	// do slightly permissive find
	for (int n = _mt->numberOfTriangles(), j, i = 0; i < n; ++i) {
		int tm = _mt->triangleMaterial(i);
		if (tm == 3 || tm == 4 || tm < 0)  // only look for permissible deep cut triangles.
			continue;
		int* tr = _mt->triangleVertices(i);
		tbb.Empty_Box();
		for (j = 0; j < 3; ++j) {
			if (_deepXyz[tr[j]].X > 1e22)  // invalid vertex, thus so is this triangle
				break;
			T[j].set(_deepXyz[tr[j]]);
			tbb.Enlarge_To_Include_Point(T[j].xyz);
		}
		if (j < 3)
			continue;
		if (!bb.Intersection(tbb))
			continue;
		Mat3x3d M;
		M.Initialize_With_Column_Vectors(T[1] - T[0], T[2] - T[0], -rayDirection);
		P = M.Robust_Solve_Linear_System(rayStart - T[0]);
		if (P.X < -1e-16 || P.Y < -1e-16 || P.X > 1.00001 || P.Y > 1.00001 || P.X + P.Y >= 1.00001)
			continue;
		if (P.Z <= -1e-8)  // look only in correct direction
			continue;
		rayTriangleIntersect rti;
		rti.scl.rtiIndexTo = -1;
		rti.scl.rtiPostTo = 1000;  // code for no connection.  diagonals are negative
		rti.deepVert = -1;
		rti.mat2Vert = -1;
		N = (T[1] - T[0]) ^ (T[2] - T[0]);
		if (N * rayDirection < 0.0f)
			rti.solidDown = 1;
		else
			rti.solidDown = 0;
		rti.triangle = i;
		rti.uv[0] = P.X;
		rti.uv[1] = P.Y;
		rti.rayParam = P.Z;
		rti.intersect = T[0] * (1.0 - P.X - P.Y) + T[1] * P.X + T[2] * P.Y;
		rtiMap.insert(std::make_pair(P.Z, rti));
	}
	// in cases where an undermine has been done there can be a double hit both at the deep material 5 triangle and its top material 2 counterpart.
	// When this occurs remove the material 2 entry and only leave its corresponding material 5 triangle.
	std::multimap<double, rayTriangleIntersect>::iterator rit, ritLast;
	rit = ritLast = rtiMap.begin();
	while (rit != rtiMap.end()) {
		++rit;
		if (rit == rtiMap.end())
			break;
		if (rit->second.solidDown == ritLast->second.solidDown && abs(rit->first - ritLast->first) < 1e-16) {  // double hit
			if (_mt->triangleMaterial(rit->second.triangle) == 2)
				rit = rtiMap.erase(rit);
			else if (_mt->triangleMaterial(ritLast->second.triangle) == 2)
				rtiMap.erase(ritLast);
			else  // doesn't matter which one you choose
				rtiMap.erase(ritLast);
		}
		ritLast = rit;
	}
	// consistency check
	rit = rtiMap.begin();
	if (rit == rtiMap.end())
		return false;
	if (rit->second.solidDown == false)
		rit = rtiMap.erase(rit);
	if (rit == rtiMap.end())
		return false;
	// initially demanded solid pairs, but in spatial coords can have minor solid collisions. Now these are allowed if minor.
	if (rtiMap.size() & 1) // insist on pairs
		return false;
	intersects.clear();
	intersects.reserve(rtiMap.size());
	int intsctIndex = 0;
	while (rit != rtiMap.end()) {  // create pairs enclosing solids, but allow minor soft-soft collisions
		assert (rit->second.solidDown == true);
		rit->second.rayIndex = intsctIndex++;
		intersects.push_back(rit->second);
		rit = rtiMap.erase(rit);
		if (rit == rtiMap.end())
			return false;
		if (rit->second.solidDown == true) {
			++rit;
			if (rit == rtiMap.end() || rit->second.solidDown == true)
				return false;
			rit->second.rayIndex = intsctIndex++;
			intersects.push_back(rit->second);
			rtiMap.erase(rit);
		}
		else {
			rit->second.rayIndex = intsctIndex++;
			intersects.push_back(rit->second);
			rtiMap.erase(rit);
		}
		rit = rtiMap.begin();
	}
	return true;
}

int deepCut::addPeriostealUndermineTriangle(const int triangle, const Vec3f &linePickDirection, bool incisionConnect)
{
	int perioTri = 0x7fffffff, mat = _mt->triangleMaterial(triangle);
	if (mat == 7 || mat == 8 || mat == 10)  // a periosteal triangle possibly in the middle of an undermine
		perioTri = triangle;
	else{
		std::vector<Vec3f> positions;
		std::vector<float> rayParams;
		std::vector<int> rayTris;
		Vec3f startPos;
		_mt->getVertexCoordinate(_mt->triangleVertices(triangle)[0], startPos.xyz);
		_mt->linePick(startPos, linePickDirection, positions, rayTris, rayParams, 7);
		if (rayTris.empty())
			return -1;
		for (int n = rayTris.size(), i = 0; i < n; ++i) {
			int mat = _mt->triangleMaterial(rayTris[i]);
			if(mat > 6){  // periosteal triangle or one already undermined
				perioTri = rayTris[i];
				break;
			}
		}
	}
	if (perioTri == 0x7fffffff)
		return 0x7fffffff;
	addUndermineTriangle(perioTri, 7, incisionConnect);
	return perioTri;
}

void deepCut::rayBilinearPatchIntersection(Vec3d& rayStart, Vec3d& rayDir, const Vec3d& P00, const Vec3d& P10, const Vec3d& P01, const Vec3d& P11,
	double(&rayParam)[2], double(&faceParam)[2][2]) {
	// slightly modified code from https://research.nvidia.com/sites/default/files/pubs/2019-03_Cool-Patches%3A-A//Chapter_08.pdf
	Vec3d e10 = P10 - P00;
	Vec3d e11 = P11 - P10;
	Vec3d e00 = P01 - P00;
	Vec3d qn = (P10 - P00) ^ (P01 - P11);
	Vec3d q00 = P00 - rayStart;
	Vec3d q10 = P10 - rayStart;
	double a = (q00 ^ rayDir) * e00;
	double c = qn * rayDir;
	double b = (q10 ^ rayDir) * e11;
	b -= a + c;
	double det = b * b - 4.0 * a * c;
	if (det < 0.0) {
		rayParam[0] = -DBL_MAX;
		rayParam[1] = -DBL_MAX;
		return;
	}
	det = sqrt(det);
	double uVec[2];
	if (c == 0.0) {  // trapezoid so one root
		uVec[0] = -a / b;
		uVec[1] = -DBL_MAX;
	}
	else {
		uVec[0] = (-b - copysign(det, b)) * 0.5;
		uVec[1] = a / uVec[0];
		uVec[0] /= c;
	}
	for (int i = 0; i < 2; ++i) {
		if (uVec[i] < -1e32) {
			rayParam[i] = -DBL_MAX;
			break;
		}
		Vec3d pa = q00 + (q10 - q00) * uVec[i];
		Vec3d pb = e00 + (e11 - e00) * uVec[i];
		Vec3d n = rayDir ^ pb;
		det = n * n;
		n = n ^ pa;
		double t1 = n * pb;
		double v1 = n * rayDir;
		if (t1 > 0.0) {
			rayParam[i] = t1 / det;
			faceParam[0][i] = uVec[i];
			faceParam[1][i] = v1 / det;
		}
		else
			rayParam[i] = -DBL_MAX;
	}
}

void deepCut::makeBilinearPatch(const Vec3d& P00, const Vec3d& P10, const Vec3d& P01, const Vec3d& P11, bilinearPatch& bl) {
	bl.P00 = P00;  bl.P10 = P10;
	bl.e10 = P10 - P00;
	bl.e01 = P11 - P01;
	bl.e11 = P11 - P10;
	bl.e00 = P01 - P00;
	bl.qn = bl.e10 ^ (P01 - P11);
	bl.P01 = P01;
	bl.P11 = P11;
}

int deepCut::bilinearRayIntersection(const Vec3d& rayStart, const Vec3d& rayDir, const bilinearPatch &bl, double (&rayParam)[2], Vec2d (&faceParams)[2]) {
	// slightly modified code from https://research.nvidia.com/sites/default/files/pubs/2019-03_Cool-Patches%3A-A//Chapter_08.pdf
	// preliminary patch data precomputed
	// new version 7/11/22 only finds patch intersections with u in the range of 0 to 1
	rayParam[0] = rayParam[1] = -DBL_MAX;
	Vec3d q00 = bl.P00 - rayStart;
	Vec3d q10 = bl.P10 - rayStart;
	double a = (q00 ^ rayDir) * bl.e00;
	double c = bl.qn * rayDir;
	double b = (q10 ^ rayDir) * bl.e11;
	b -= a + c;
	double det = b * b - 4.0 * a * c;
	if (det < 0.0)
		return 0;
	det = sqrt(det);
	double uVec[2];
	if (c == 0.0) {  // trapezoid so one root
		uVec[0] = -a / b;
		uVec[1] = -DBL_MAX;
	}
	else {
		uVec[0] = (-b - copysign(det, b)) * 0.5;
		uVec[1] = a / uVec[0];
		uVec[0] /= c;
	}
	int nSols = 0;
	// nVidia routine for finding v and t not sufficiently accurate in extrapolation regions.  Using more precise form instead.
	// Also due to extrapolation need both answers since a valid value of t may have extreme values of u and v.
	// As of 7/11/2022 will no longer extrapolate. u must be in range of 0 to 1
	for (int i = 0; i < 2; ++i) {
		if (uVec[i] < 0.0 || uVec[i] > 1.0)
			continue;
		Vec3d pa = bl.P00 + bl.e10 * uVec[i];
		Vec3d T, pb = bl.P01 * (1.0 - uVec[i]) + bl.P11 * uVec[i];
		pb -= pa;
		T = rayStart - pa;
		Mat2x2d M(pb * pb, rayDir * pb, 0, -rayDir * rayDir);
		M.x[2] = -M.x[1];
		Vec2d B(T * pb, T*rayDir), R;
		R = M.Robust_Solve_Linear_System(B);
		if (R.Y > 0.0 && R.Y <= 1.0) {
			// not allowing extrapolation any longer
			rayParam[nSols] = R.Y;
			faceParams[nSols].X = uVec[i];
			faceParams[nSols].Y = R.X;
			++nSols;
		}
	}
	return nSols;
}

bool deepCut::planeRayIntersection(const Vec3d P, const Vec3d R, const endPlane *ep, double& rayParam, Vec2d& faceParam) {
	// Same behavior as bilinearRayIntersection() but intersection surface is an end plane.
	double d0 = P * ep->N - ep->d, d1 = (P + R) * ep->N - ep->d;
	if (signbit(d0) == signbit(d1))
		return false;
	rayParam = d0 / (d0 - d1);
	Vec3d I = P + R * rayParam - ep->P;
	faceParam.X = I * ep->U;
	if (ep != &_endPlanes[0]) {
		if (faceParam.X <= 0.0)
			return false;
	}
	else {
		if (faceParam.X >= 0.0)
			return false;
	}
	faceParam.Y = I * ep->V;
	faceParam *= ep->VlengthSqInv;
	return true;
}

void deepCut::findCutInteriorHoles(const bilinearPatch* blp, const endPlane* ep, const std::vector<Vec2d>& bUv, const std::vector<std::pair<int, Vec2d> >& deepOuterPolygon, std::list< vertUvPolygon >& holes) {
	holes.clear();
	boundingBox<double> bb;
	bb.Empty_Box();
	for (auto& bp : deepOuterPolygon)
		bb.Enlarge_To_Include_Point(_deepXyz[bp.first].xyz);
	Vec3d E, nE;
	double rayParams[2];
	Vec2d faceParams[2];
	auto edgeIntersect = [&]() ->bool {
		if (blp == nullptr) {
			if (!planeRayIntersection(E, nE - E, ep, rayParams[0], faceParams[0]))
				return false;
		}
		else {
			if (bilinearRayIntersection(E, nE - E, *blp, rayParams, faceParams) != 1)
				return false;
		}
		return true;
	};
	std::vector<bool> trisUsed(_mt->numberOfTriangles(), false);
	struct interiorHole {
		std::vector<Vec2d> uvs;
		Vec2d uv2[2];
		unsigned int te2[2];
		double param2[2];
	};
	std::list< interiorHole > holeList;
	for (int n = _mt->numberOfTriangles(), j, i = 0; i < n; ++i) {
		int mat = _mt->triangleMaterial(i);
		if (mat < 0 || (mat >2 && mat < 5))
			continue;
		if (trisUsed[i])
			continue;
		int* tr = _mt->triangleVertices(i);
		boundingBox<double> bt;
		bt.Empty_Box();
		for (j = 0; j < 3; ++j) {
			if (tr[j] > _preDeepCutVerts)  // new triangle already processed by this deep cut.
				break;
			bt.Enlarge_To_Include_Point(_deepXyz[tr[j]].xyz);
		}
		if (j < 3 || !bb.Intersection(bt))
			continue;

		Vec2d uv2[2];
		unsigned int te2[2] = {0xffffffff, 0xffffffff};
		double param2[2];
		nE = _deepXyz[tr[0]];
		for (j = 2; j > -1; --j) {
			E = _deepXyz[tr[j]];
			if(j<1 && te2[1] > 0xfffffffe)
				continue;
			if (edgeIntersect()) {
				if (te2[1] < 0xffffffff) {
					te2[0] = (i << 2) + j;  // reference original triangle
					param2[0] = rayParams[0];
					uv2[0] = faceParams[0];
				}
				else {
					insidePolygon ip;
					if (!ip.insidePolygon2d(faceParams[0], bUv))
						break;
					te2[1] = _mt->triAdjs(i)[j];
					param2[1] = 1.0f - rayParams[0];
					uv2[1] = faceParams[0];
				}
			}
			nE = E;
		}
		if (te2[0] > 0xfffffffe)
			continue;
		// te3[0] is original triangle data
		// get all possible hole paths
		holeList.push_back(interiorHole());
		surfaceCutLine scl;
		double minV;
		if (surfacePathSub(-1, -1, -1, -1, te2[1], param2[1], uv2[1], te2[0] >> 2, blp, ep, false, scl, minV) < DBL_MAX) {  // only take in connected holes. May have read in a line of coincident surface which should be ignored.
			auto& h = holeList.back();
			h.uvs.assign(scl.deepUVs.begin(), scl.deepUVs.end());
			h.uvs.push_back(uv2[1]);
			h.uv2[0] = uv2[0];
			h.uv2[1] = uv2[1];
			h.param2[0] = param2[0];
			h.param2[1] = param2[1];
			h.te2[0] = te2[0];
			h.te2[1] = te2[1];
		}
		else
			holeList.pop_back();
		for (auto& t : scl.deepVertsTris)  // don't look along this line again
			trisUsed[t] = true;
	}
	if (holeList.empty())
		return;
	// only want first level holes. Solids within holes are not cut.
	auto hit = holeList.begin();
	while (hit != holeList.end()) {
		auto hit2 = hit;
		++hit2;
		while (hit2 != holeList.end()) {
			insidePolygon ip;
			if (ip.insidePolygon2d(hit2->uvs.front(), hit->uvs)) {
				hit2 = holeList.erase(hit2);
			}
			else if (ip.insidePolygon2d(hit->uvs.front(), hit2->uvs)) {
				hit = holeList.erase(hit);
				break;
			}
			else
				++hit2;
		}
		if (hit2 == holeList.end())
			++hit;
	}
	auto uvFromEdge = [&](int edge, double param, Vec2f& uv) {
		uv.set(0.0, 0.0);
		if (edge < 1)
			uv[0] = param;
		else if(edge > 1)
			uv[1] = 1.0 - param;
		else {
			uv[0] = 1.0 - param;
			uv[1] = param;
		}
	};
	auto deepSurfacePoint = [&](unsigned int &te, Vec2f &uv) ->int {
		Vec3f gridLocus, bw;
		int tet = _vbt->parametricTriangleTet(te >> 2, uv.xy, gridLocus);
		if (tet < 0)
			throw(std::logic_error("Program error in deepCut hole finder."));
		_vbt->gridLocusToBarycentricWeight(gridLocus, _vbt->_tetCentroids[tet], bw);
		int ret = _mt->addNewVertexInMidTriangle(te >> 2, uv.xy);
		assert(ret == _vbt->_vertexTets.size());
		_vbt->_vertexTets.push_back(tet);
		_vbt->_barycentricWeights.push_back(bw);
		return ret;
	};
	for (auto& h : holeList) {
		int j, newP, endTriangle, cutP[2][2], mat = _mt->triangleMaterial(h.te2[0] >> 2);
		auto adj0 = _mt->triAdjs(h.te2[0] >> 2)[h.te2[0] & 3];
		Vec2f uv[2], uvt;
		uvFromEdge(h.te2[0] & 3, h.param2[0], uv[0]);
		auto te1 = _mt->triAdjs(h.te2[1] >> 2)[h.te2[1] & 3];
		uvFromEdge(te1 & 3, 1.0 - h.param2[1], uv[1]);
		uvt = uv[0] * 0.6667 + uv[1] * 0.3333;
		if (mat == 2) {
			createFlapTopBottomVertices(h.te2[0] >> 2, uvt.xy, cutP[0][0], cutP[0][1]);
			newP = cutP[0][0];
		}
		else {
			cutP[0][0] = -1;
			cutP[0][1] = deepSurfacePoint(h.te2[0], uvt);
			newP = cutP[0][1];
		}
		te1 = _mt->triAdjs(h.te2[1] >> 2)[h.te2[1] & 3];
		auto tr = _mt->triangleVertices(te1>>2);
		for (j = 0; j < 3; ++j)
			if (tr[j] == newP)
				break;
		if (j < 1)
			uvt.set(0.0f, 0.0f);
		else if (j > 1)
			uvt.set(0.0f, 1.0f);
		else
			uvt.set(1.0f, 0.0f);
		uvFromEdge(te1 & 3, 1.0 - h.param2[1], uv[1]);
		uvt *= 0.5f;
		uvt += uv[1] * 0.5f;
		if (mat == 2)
			createFlapTopBottomVertices(te1 >> 2, uvt.xy, cutP[1][0], cutP[1][1]);
		else {
			cutP[1][0] = -1;
			cutP[1][1] = deepSurfacePoint(te1, uvt);
		}
		endTriangle = _mt->triAdjs(adj0 >> 2)[adj0 & 3] >> 2;
		updateDeepSpatialCoordinates();
		surfaceCutLine scl;
		double minimumBilinearV;
		if (surfacePathSub(cutP[1][0], cutP[1][1], cutP[0][0], cutP[0][1], h.te2[1], h.param2[1], h.uv2[1], endTriangle, blp, ep, true, scl, minimumBilinearV) == DBL_MAX)
			throw(std::logic_error("Couldn't track a valid hole path in hole finder."));
		if (mat == 2) {
			std::list<int> topVerts, botVerts;
			topVerts.push_back(cutP[0][0]);
			topVerts.push_back(cutP[1][0]);
			botVerts.push_back(cutP[0][1]);
			botVerts.push_back(cutP[1][1]);
			_loopSkinTopBegin = cutP[1][0];
			_previousSkinTopEnd = cutP[0][0];
			if (!topDeepSplit_Sub(topVerts, botVerts, true, true))
				throw(std::logic_error("Program error in skinCutLine()."));
		}
		vertUvPolygon vup;
		vup.vertices.assign(scl.deepVertsTris.begin(), scl.deepVertsTris.end());
		vup.uvs.assign(scl.deepUVs.begin(), scl.deepUVs.end());
		vup.uvs[vup.uvs.size() - 1] = vup.uvs[1] * 0.3333f + vup.uvs[vup.uvs.size() - 2] * 0.6667f;
		vup.uvs[0] = vup.uvs[1] * 0.6667f + vup.uvs[vup.uvs.size() - 2] * 0.3333f;
		// holes must be listed counter clockwise
		if (clockwise(vup.uvs) < 2) {
			std::reverse(vup.uvs.begin(), vup.uvs.end());
			std::reverse(vup.vertices.begin(), vup.vertices.end());
			scl.deepVertsTris.reverse();
		}
		holes.push_back(vup);
		// save for surface splitter
		_holePolyLines.push_back(std::move(scl.deepVertsTris));
		_holePolyLines.back().push_back(_holePolyLines.back().front());  // closed loop code
	}
}

double deepCut::surfacePath(rayTriangleIntersect& from, const rayTriangleIntersect& to, const bool cutPath, double& minimumBilinearV) {
	const bilinearPatch* blp = nullptr;
	const endPlane* ep = nullptr;
	if (from.postNum == to.postNum) {
		if (from.postNum == 0) {
			if ((from.rayIndex & 1) == 0 && from.rayIndex > to.rayIndex)
				ep = nullptr;
			else {
				ep = &_endPlanes[0];  // endPlaneCut = 1;
				assert(_endPlanes[0].P.X != DBL_MAX);
			}
		}
		if (from.postNum == _deepPosts.size() - 1) {
			if (from.rayIndex & 1 && from.rayIndex < to.rayIndex)
				ep = nullptr;
			else {
				ep = &_endPlanes[1];  // endPlaneCut = 2;
				assert(ep->P.X != DBL_MAX);
			}
		}
	}
	minimumBilinearV = DBL_MAX;
//	if (cutPath) {
//		if (!from.scl.deepVertsTris.empty())
//			throw(std::logic_error("surfacePath() called on a path that is not empty;"));
//	}
	bilinearPatch bl;
	Vec3d P00, P10, P01, P11;  //  , lastI = from.intersect;
	if (cutPath && ep == nullptr) {  // use final bilinear surface definition
		if (from.postNum < to.postNum)
			blp = &_deepPosts[to.postNum].bl;
		else if (from.postNum > to.postNum)
			blp = &_deepPosts[from.postNum].bl;
		else if (from.rayIndex > to.rayIndex) {  // COURT - no longer allowed
			if (from.postNum < 1) {  // open end case
				if (to.rayIndex & 1) // hole
					blp = &_deepPosts[from.postNum + 1].bl;
				else  // other open end case
					throw(std::logic_error("Program me."));
			}
			else if (from.postNum == _deepPosts.size() - 1) // open end case in hole
				throw(std::logic_error("Program me."));
			else
				blp = &_deepPosts[from.postNum + 1].bl;
		}
		else {  // (from.rayIndex < to.rayIndex)
			if (from.postNum < 1)  // open end case
				throw(std::logic_error("Program me."));
			else if (from.postNum == _deepPosts.size() - 1) { // open end case in hole
				if (from.rayIndex & 1) // hole
					blp = &_deepPosts[from.postNum].bl;
				else  // other open end case
					throw(std::logic_error("Program me."));
			}
			else {
				blp = &_deepPosts[from.postNum].bl;
			}
		}
	}
	else if (ep == nullptr) {  // get a provisional bilinear surface
		int intersectLevel = _deepPosts[from.postNum].triIntersects.size() - 1;  // keep bilinear surface as square as possible
		if (intersectLevel > _deepPosts[to.postNum].triIntersects.size() - 1)
			intersectLevel = _deepPosts[to.postNum].triIntersects.size() - 1;
		if (from.postNum < to.postNum) {
			assert((from.rayIndex & 1) < 1);
			P01 = _deepPosts[from.postNum].triIntersects.front().intersect;
			P11 = _deepPosts[to.postNum].triIntersects.front().intersect;
			P00 = _deepPosts[from.postNum].triIntersects[intersectLevel].intersect;
			P10 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
		}
		else if (from.postNum > to.postNum) {
			assert(from.rayIndex & 1);
			P10 = _deepPosts[from.postNum].triIntersects[intersectLevel].intersect;
			P00 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
			P11 = _deepPosts[from.postNum].triIntersects.front().intersect;
			P01 = _deepPosts[to.postNum].triIntersects.front().intersect;
		}
		else if (from.rayIndex > to.rayIndex) {  // COURT - no longer allowed
			if (from.postNum == _deepPosts.size() - 1) { // open end case
				if (intersectLevel > _deepPosts[from.postNum - 1].triIntersects.size() - 1)
					intersectLevel = _deepPosts[from.postNum - 1].triIntersects.size() - 1;
				P11 = _deepPosts[from.postNum].triIntersects.front().intersect;
				P10 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
				P01 = _deepPosts[from.postNum - 1].triIntersects.front().intersect;
				P00 = _deepPosts[to.postNum - 1].triIntersects[intersectLevel].intersect;
			}
			else {
				if (from.rayIndex & 1) {  // other open end case.  COURT - again no longer allowed
					assert(from.postNum < 1);
				}
				if (intersectLevel > _deepPosts[from.postNum + 1].triIntersects.size() - 1)
					intersectLevel = _deepPosts[from.postNum + 1].triIntersects.size() - 1;
				P01 = _deepPosts[from.postNum].triIntersects.front().intersect;
				P00 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
				P11 = _deepPosts[from.postNum + 1].triIntersects.front().intersect;
				P10 = _deepPosts[to.postNum + 1].triIntersects[intersectLevel].intersect;
			}
		}
		else {
			if (from.postNum < 1) {  // open end case
				if (intersectLevel > _deepPosts[from.postNum + 1].triIntersects.size() - 1)
					intersectLevel = _deepPosts[from.postNum + 1].triIntersects.size() - 1;
				P01 = _deepPosts[from.postNum].triIntersects.front().intersect;
				P00 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
				P11 = _deepPosts[from.postNum + 1].triIntersects.front().intersect;
				P10 = _deepPosts[to.postNum + 1].triIntersects[intersectLevel].intersect;
			}
			else {
				if (to.rayIndex & 1) {  // other open end case
					assert(from.postNum == _deepPosts.size() - 1);
				}
				if (intersectLevel > _deepPosts[from.postNum - 1].triIntersects.size() - 1)
					intersectLevel = _deepPosts[from.postNum - 1].triIntersects.size() - 1;
				P11 = _deepPosts[from.postNum].triIntersects.front().intersect;
				P10 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
				P01 = _deepPosts[from.postNum - 1].triIntersects.front().intersect;
				P00 = _deepPosts[to.postNum - 1].triIntersects[intersectLevel].intersect;
			}
		}
		makeBilinearPatch(P00, P10, P01, P11, bl);
		blp = &bl;
	}
	else  // is an end plane cut
		assert(ep != nullptr);
	Vec3d E, nE; // , I;
	double rayParams[2];
	Vec2d faceParams[2];
	auto edgeIntersect = [&]() ->bool {
		if (ep != nullptr) {
			if (!planeRayIntersection(E, nE - E, ep, rayParams[0], faceParams[0]))
				return false;
		}
		else {
			if (bilinearRayIntersection(E, nE - E, *blp, rayParams, faceParams) != 1)
				return false;
		}
		if (minimumBilinearV > faceParams[0].Y)
			minimumBilinearV = faceParams[0].Y;
		return true;
	};
	unsigned int te = 0xffffffff;
	int* tr;
	int nEdges = 1, prevMat;
	// next lambda is for difficult triangles whose edges don't cut cleanly with bilinearRayIntersection()
	auto getStartEndTE = [&](const rayTriangleIntersect& rti) ->unsigned int {
		std::vector<materialTriangles::neighborNode> nei;
		_mt->getNeighbors(rti.mat2Vert < 0 ? rti.deepVert : rti.mat2Vert, nei);
		E = _deepXyz[nei.back().vertex];
		int mat;
		for (auto& n : nei) {
			nE = _deepXyz[n.vertex];
			if ((mat = _mt->triangleMaterial(n.triangle)) < 3 || mat > 4) {
				assert(E.X != DBL_MAX);
				if (edgeIntersect()) {
					tr = _mt->triangleVertices(n.triangle);
					for (int i = 0; i < 3; ++i) {
						if (tr[i] == n.vertex) {
							prevMat = mat;
							return _mt->triAdjs(n.triangle)[(i + 2) % 3];
						}
					}
					break;
				}
			}
			E = nE;
		}
		return 3UL;  // empty triEdge
	};
	int endTriangle = to.triangle;
	if (from.deepVert < 0) {
		prevMat = _mt->triangleMaterial(from.triangle) == 2 ? 2 : 5;  // 5 or 6, 7, or 8
		tr = _mt->triangleVertices(from.triangle);
		int i;
		nE = _deepXyz[tr[0]];
		for (i = 2; i > -1; --i) {
			E = _deepXyz[tr[i]];
			assert(E.X != DBL_MAX);
			if (edgeIntersect()) {
				te = _mt->triAdjs(from.triangle)[i];
				break;
			}
			nE = E;
		}
		if (i < 0)
			throw(std::logic_error("This surfacePath() call does not have a valid starting triangle."));
	}
	else {
		if ((te = getStartEndTE(to)) == 3)  // COURT use triangleEdgeCuts()
			throw(std::logic_error("Program error start/ending a deepCut."));
		endTriangle = _mt->triAdjs(te >> 2)[te & 3] >> 2;
		if ((te = getStartEndTE(from)) == 3)
			throw(std::logic_error("Program error start/ending a deepCut."));
	}
	int topStartV = from.mat2Vert, deepStartV = from.deepVert, topEndV = to.mat2Vert, deepEndV = to.deepVert;
	return surfacePathSub(topStartV, deepStartV, topEndV, deepEndV, te, 1.0 - rayParams[0], faceParams[0], endTriangle, blp, ep, cutPath, from.scl, minimumBilinearV);
}

double deepCut::surfacePathSub(int topStartV, int deepStartV, int topEndV, int deepEndV, const unsigned int startTE, const double &startParam, const Vec2d &startUV, const int endTriangle,
	const bilinearPatch* bl, const endPlane* ep, const bool cutPath, surfaceCutLine& scl, double& minimumBilinearV){
	// broken into a subroutine for use in finding holes.  We are still using minimumBilinearV to prevent deepPost crossover.
	// To process holes will fill scl.deepUVs even when not cutPath looking for topological connections.  Delete these before calling cutPath == true subsequently.
	if (cutPath) {
		scl.deepUVs.clear();  // delete previous intermediate results used in path searching. Intentionally kept after a search for hole finding.
		scl.deepVertsTris.clear();
	}
	double len = 0.0;
	int nEdges = 0;
	Vec3d E, nE, I, lastI;
	double rayParams[2];
	Vec2d faceParams[2];
	scl.deepUVs.clear();  // next lines clear from a previous topological path check
	scl.deepVertsTris.clear();
	if (!cutPath)  // used for hole finding
		scl.deepVertsTris.push_back(endTriangle);
	auto edgeIntersect = [&]() ->bool {
		if (ep != nullptr) {
			if (!planeRayIntersection(E, nE - E, ep, rayParams[0], faceParams[0]))
				return false;
		}
		else {
			if (bilinearRayIntersection(E, nE - E, *bl, rayParams, faceParams) != 1)
				return false;
		}
		// above routine no longer allows u outside patch results
		I = E * (1.0 - rayParams[0]) + nE * rayParams[0];
		if(nEdges)
			len += (I - lastI).length();
		lastI = I;
		if (minimumBilinearV > faceParams[0].Y)
			minimumBilinearV = faceParams[0].Y;
		return true;
	};
	unsigned int te = startTE;
	std::vector<unsigned int> topTe;
	std::vector<float> topParams;
	std::vector<Vec2d> topUVs;
	if (cutPath) {
		topTe.push_back(te);
		topParams.push_back((float)startParam);
		topUVs.push_back(startUV);
	}
	int prevMat = _mt->triangleMaterial(_mt->triAdjs(te>>2)[te&3] >> 2);  //  , deepStart = from.deepVert, topStart = from.mat2Vert;
	Vec2d TinUV(DBL_MAX, 0.0), UinUV(DBL_MAX, 0.0), deepUvStart(DBL_MAX, 0.0), topUvStart(DBL_MAX, 0.0), UoutUV(-DBL_MAX, 0.0), UintUV(-DBL_MAX, 0.0);  // This routine will always start and finish on a deep post point.
	int TinTri = -1, ToutTri = -1, mat2BorderVertex = -1, mat2BorderTexture = -1;
	unsigned int UinTe = 3, UoutTe = 3;
	float TinParam, ToutParam, UinParam, UoutParam;
	bool Tin = false;
	do{
		int splitTri, splitEdge, splitMat;
		splitTri = te >> 2, splitEdge = te & 3;
		splitMat = _mt->triangleMaterial(splitTri);
		if (splitMat < 2) {
			if (splitMat > 0)
				splitMat = 5;
			else if(splitMat < 0)
				throw(std::logic_error("Trying to deep cut through a deleted triangle.|n"));
			else
				throw(std::logic_error("Model error with triangle material = 0.|n"));
		}
		if (splitTri == endTriangle)
			break;
		auto tr = _mt->triangleVertices(splitTri);
		if (splitMat == 4) {  // mat 5 march to an undermine border edge. Pop back on top
			for (auto& db : _deepBed) {
				if (db.second.deepMtVertex == tr[(splitEdge + 1) % 3]) {
					std::vector<materialTriangles::neighborNode> nei;
					_mt->getNeighbors(db.first, nei);
					auto nit = nei.begin();
					while (nit != nei.end()) {
						if (_deepBed[nit->vertex].deepMtVertex == tr[splitEdge]) {
							tr = _mt->triangleVertices(nit->triangle);
							int i;
							for (i = 0; i < 3; ++i) {
								if (tr[i] == nit->vertex) {
									te = (nit->triangle << 2) + i;
									if (cutPath) {
										topTe.pop_back();  // don't cut this one. It will be cut by this skin split.
										topParams.pop_back();
										UinUV = topUVs.back();
										topUVs.pop_back();
										cutDeepSurface(deepStartV, deepUvStart, -1, deepUvStart, topTe, topParams, topUVs, scl);  // COURT - watch for deepUvStart
										UinTe = _mt->triAdjs(nit->triangle)[i];
										UinParam = 1.0f - (float)rayParams[0];
										Tin = false;
										prevMat = 2;
									}
									break;
								}
							}
							if (i > 2)
								return false;
							break;
						}
						++nit;
					}
					assert(nit != nei.end());
					++nEdges;
					break;
				}
			}
		}
		else if (splitMat == 3) {  // march over an incision edge until a mat 5 or 2 triangle is found
			if (prevMat == 2) {  // this is a T out
				auto adjs = _mt->triAdjs(splitTri + 1);  // incision convention
				if (_mt->triangleMaterial(adjs[0] >> 2) == 3) {  // non undermined incision edge

					// COURT - write me

					adjs = _mt->triAdjs((adjs[0] >> 2) - 1);  // incision convention again
					assert(_mt->triangleMaterial(adjs[0] >> 2) == 2);
					te = adjs[0];
					prevMat = 2;
				}
				else {
					assert(_mt->triangleMaterial(adjs[0] >> 2) > 4);
					if (cutPath) {
						ToutTri = topTe.back() >> 2;
						ToutParam = topParams.back();
						deepStartV = -1;
						deepUvStart.set(DBL_MAX, 0.0);
					}
					te = adjs[0];
					prevMat = 5;
				}
			}
			else {  // Tin
				assert(prevMat == 5);
				if (cutPath) {
					TinTri = (te >> 2) - 1;
					TinParam = (float)rayParams[0];
					TinUV = faceParams[0];
					topTe.pop_back();
					topParams.pop_back();
					topUVs.pop_back();
				}
				auto adjs = _mt->triAdjs((te >> 2) - 1);  // incision convention
				te = adjs[0];
				assert(_mt->triangleMaterial(te >> 2) == 2);
				prevMat = 2;
			}
			++nEdges;
		}
		else {
			E = _deepXyz[tr[(splitEdge + 1) % 3]];
			int i;
			for (i = 2; i < 4; ++i) {  // always find the next edge before cutting previous triangle
				nE = _deepXyz[tr[(splitEdge + i) % 3]];
				if (nE.X == DBL_MAX) {  // signals invalid triangle beyond an undermined edge
					UoutTe = te;
					UoutParam = 1.0f - (float)rayParams[0];
					int v50 = _deepBed[tr[splitEdge]].deepMtVertex;
					int v51 = _deepBed[tr[(splitEdge + 1) % 3]].deepMtVertex;
					assert(v50 > -1 && v51 > -1);
					std::vector<materialTriangles::neighborNode> nei;
					_mt->getNeighbors(v51, nei);
					auto nit = nei.begin();
					while (nit != nei.end()) {
						if (nit->vertex == v50) {
							assert(_mt->triangleMaterial(nit->triangle) == 5);
							int* trp = _mt->triangleVertices(nit->triangle);
							int k;
							for (k = 0; k < 3; ++k) {
								if (trp[k] == v50)
									break;
							}
							assert(k < 3);
							te = (nit->triangle << 2) + k;  // deep mat 5 tris always listed in same order as the top tri that generated it.
							prevMat = 5;
							break;
						}
						++nit;
					}
					assert(nit != nei.end());
					break;
				}
				if (edgeIntersect()) {  //  || splitTri == endTriangle condition for open ends
					te = _mt->triAdjs(splitTri)[(splitEdge + i - 1) % 3];  // these two lines are for the next cut
					++nEdges;
					if (cutPath) {
						// do previous Tin or Tout only after next te secured as the T op will split the top edge triangle
						if (TinTri > -1) {
							cutDeepSurface(deepStartV, deepUvStart, -1, deepUvStart, topTe, topParams, topUVs, scl);
							topStartV = TinSub(TinTri, TinParam);
							topUvStart = TinUV;
							topTe.clear();
							topParams.clear();
							topUVs.clear();
							Tin = true;
							TinTri = -1;
						}
						if (ToutTri > -1) {
							topTe.pop_back();
							topParams.pop_back();
							auto ToutUV = topUVs.back();
							topUVs.pop_back();
							unsigned int ate;
							if(!topTe.empty())
								ate = _mt->triAdjs(topTe.back() >> 2)[topTe.back() & 3];
							int tv = TinSub(ToutTri, ToutParam);
							if (!topTe.empty())
								topTe.back() = _mt->triAdjs(ate >> 2)[ate & 3];
							cutSkinLine(topStartV, topUvStart, tv, ToutUV, topTe, topParams, topUVs, Tin, true, scl);
							if (mat2BorderVertex > -1) {
								auto vit = scl.deepVertsTris.begin();
								do {
									++vit;
								} while (*vit != mat2BorderVertex);
								++vit;
								int topV = -1;
								for (auto dp : _deepBed) {
									if (dp.second.deepMtVertex == *vit)
										topV = dp.first;
								}
								if (topV < 0)
									throw(std::logic_error("Program error at material 2-border vertex.\n"));
								mat2BorderSplit(mat2BorderVertex, mat2BorderTexture, topV);
								mat2BorderVertex = -1;
							}
							topTe.clear();
							topParams.clear();
							topUVs.clear();
							deepStartV = -1;
							deepUvStart.set(-DBL_MAX, 0.0);
							ToutTri = -1;
							Tin = false;
						}
						if (UinTe != 3) {
							unsigned int ate = _mt->triAdjs(UinTe >> 2)[UinTe & 3];
							assert(_deepBed[_mt->triangleVertices(ate >> 2)[((ate & 3) + 2) % 3]].deepMtVertex < 0);
							std::list<int> topVerts, deepVerts;
							topVerts.push_back( _mt->triangleVertices(UinTe >> 2)[((UinTe&3)+2)%3]);
							deepVerts.push_back(_deepBed[topVerts.back()].deepMtVertex);
							float uv[2] = {0.333f, 0.333f};
							int bottomVertex;
							createFlapTopBottomVertices(UinTe >> 2, uv, topStartV, bottomVertex);
							topVerts.push_back(topStartV);
							deepVerts.push_back(bottomVertex);
							if ((ate & 3) < 1) {
								uv[0] = 1.0f - UinParam;
								uv[1] = 0.0f;
							}
							else if ((ate & 3) > 1) {
								uv[1] = UinParam;
								uv[0] = 0.0f;
							}
							else {
								uv[0] = 1.0f - UinParam;
								uv[1] = UinParam;
							}
							createFlapTopBottomVertices(ate >> 2, uv, topStartV, bottomVertex);
							topVerts.push_back(topStartV);
							deepVerts.push_back(bottomVertex);
							topDeepSplit_Sub(topVerts, deepVerts, false, false);
							topStartV = topVerts.back();
							topUvStart = UinUV;
							topTe.clear();
							topParams.clear();
							topUVs.clear();
							Tin = true;
							UinTe = 3;
						}
						if (UoutTe != 3) {  // already have a skin top cut in progress
							unsigned int ate = _mt->triAdjs(UoutTe >> 2)[UoutTe & 3];
							std::list<int> topVerts, deepVerts;
							topVerts.push_back(_mt->triangleVertices(UoutTe >> 2)[((UoutTe & 3) + 2) % 3]);
							deepVerts.push_back(_deepBed[topVerts.back()].deepMtVertex);
							float uv[2] = { 0.333f, 0.333f };
							int topVertex, bottomVertex;
							createFlapTopBottomVertices(UoutTe >> 2, uv, topVertex, bottomVertex);
							topVerts.push_back(topVertex);
							deepVerts.push_back(bottomVertex);
							UoutTe = _mt->triAdjs(ate >> 2)[ate & 3];
							if ((UoutTe & 3) < 1) {
								uv[0] = 1.0f - UoutParam;
								uv[1] = 0.0f;
							}
							else if ((UoutTe & 3) > 1) {
								uv[1] = UoutParam;
								uv[0] = 0.0f;
							}
							else {
								uv[0] = 1.0f - UoutParam;
								uv[1] = UoutParam;
							}
							topTe.pop_back();
							topParams.pop_back();
							auto UoutUV = topUVs.back();
							topUVs.pop_back();
							unsigned int teEnd = _mt->triAdjs(topTe.back() >> 2)[topTe.back() & 3];
							createFlapTopBottomVertices(UoutTe >> 2, uv, topVertex, bottomVertex);
							topVerts.push_back(topVertex);
							deepVerts.push_back(bottomVertex);
							topDeepSplit_Sub(topVerts, deepVerts, false, false);
							topTe.back() = _mt->triAdjs(teEnd >> 2)[teEnd & 3];
							cutSkinLine(topStartV, topUvStart, topVertex, UoutUV, topTe, topParams, topUVs, Tin, true, scl);
							if (mat2BorderVertex > -1) {
								auto vit = scl.deepVertsTris.begin();
								do {
									++vit;
								} while (*vit != mat2BorderVertex);
								++vit;
								int topV = -1;
								for (auto dp : _deepBed) {
									if (dp.second.deepMtVertex == *vit)
										topV = dp.first;
								}
								if (topV < 0)
									throw(std::logic_error("Program error at material 2-border vertex.\n"));
								mat2BorderSplit(mat2BorderVertex, mat2BorderTexture, topV);
								mat2BorderVertex = -1;
							}
							topStartV = -1;
							deepStartV = bottomVertex;
							topTe.clear();
							topParams.clear();
							topUVs.clear();
							Tin = false;
							UoutTe = 3;
						}
						if (splitMat == 2 && prevMat > 4) {  // skin start from a boundary or periosteal edge
							topTe.pop_back();
							float lastParam = topParams.back();
							topParams.pop_back();
							auto lastUV = topUVs.back();
							topUVs.pop_back();
							cutDeepSurface(deepStartV, deepUvStart, -1, deepUvStart, topTe, topParams, topUVs, scl);
							Vec3f gridLocus, bw;
							int* st = _mt->triangleVertices(splitTri);
							int tet = _vbt->parametricEdgeTet(st[splitEdge], st[(splitEdge + 1) % 3], lastParam, gridLocus);
							_vbt->gridLocusToBarycentricWeight(gridLocus, _vbt->tetCentroid(tet), bw);
							mat2BorderVertex = _mt->splitTriangleEdge(splitTri, splitEdge, lastParam);
							mat2BorderTexture = _mt->numberOfTextures() - 1;  // added in above call
							assert(mat2BorderVertex == _vbt->_vertexTets.size());
							_vbt->_vertexTets.push_back(tet);
							_vbt->_barycentricWeights.push_back(bw);
							scl.deepVertsTris.push_back(mat2BorderVertex);
							scl.deepUVs.push_back(lastUV);
							topStartV = -1;
							topUvStart = faceParams[0];
							topTe.clear();
							topParams.clear();
							topUVs.clear();
							Tin = false;
							TinTri = -1;
							prevMat = 2;
						}
						if (splitMat > 4 && prevMat == 2) {  // skin end from a boundary or periosteal edge
							unsigned int lastTe = topTe.back();
							topTe.pop_back();
							float lastParam = topParams.back();
							topParams.pop_back();
							auto lastUV = topUVs.back();
							topUVs.pop_back();
							unsigned int lastTe2 = topTe.back();
							topTe.pop_back();
							float lastParam2 = topParams.back();
							topParams.pop_back();
							auto lastUV2 = topUVs.back();
							topUVs.pop_back();
							unsigned int ate = _mt->triAdjs(topTe.back() >> 2)[topTe.back() & 3];
							float uv[2] = { 0.0f, 0.0f };
							if ((lastTe2&3) < 1)
								uv[0] = lastParam2;
							else if ((lastTe2 & 3) > 1)
								uv[1] = 1.0f - lastParam2;
							else {
								uv[0] = 1.0f - lastParam2;
								uv[1] = lastParam2;
							}
							int topVertex, bottomVertex;
							createFlapTopBottomVertices(lastTe2>>2, uv, topVertex, bottomVertex);
							topTe.back() = _mt->triAdjs(ate >> 2)[ate & 3];
							cutSkinLine(topStartV, topUvStart, topVertex, lastUV2, topTe, topParams, topUVs, Tin, false, scl);
							Vec3f gridLocus, bw;
							int* st = _mt->triangleVertices(lastTe >> 2);
							int tet = _vbt->parametricEdgeTet(st[lastTe & 3], st[((lastTe & 3) + 1) % 3], lastParam, gridLocus);
//							int tet = parametricMTedgeTet(lastTe >> 2, lastTe & 3, lastParam, gridLocus);
							_vbt->gridLocusToBarycentricWeight(gridLocus, _vbt->tetCentroid(tet), bw);
							int newV = _mt->splitTriangleEdge(lastTe >> 2, lastTe & 3, lastParam);
							int newTx = _mt->numberOfTextures() - 1;  // added in above call
							assert(newV == _vbt->_vertexTets.size());
							_vbt->_vertexTets.push_back(tet);
							_vbt->_barycentricWeights.push_back(bw);
							mat2BorderSplit(newV, newTx, topVertex);
							topTe.clear();
							topParams.clear();
							topUVs.clear();
							deepStartV = newV;
							deepUvStart = lastUV;
							ToutTri = -1;
							TinTri = -1;
							Tin = false;
						}
						topTe.push_back(te);
						topParams.push_back(1.0f - (float)rayParams[0]);
						topUVs.push_back(faceParams[0]);
					}
					else {  // Use intersect UV and triangle path in hole finder
						scl.deepUVs.push_back(faceParams[0]);
						scl.deepVertsTris.push_back(splitTri);
					}
					if (splitMat == 2)
						prevMat = 2;
					else
						prevMat = 5;  // includes 1 & 5-9
					break;
				}
				E = nE;
			}
			if (i > 3)
				return DBL_MAX;
		}
	}while ((te >> 2) != endTriangle && nEdges < 500);
	if (cutPath) {
		Vec2d toUV(DBL_MAX, 0.0);
		if (prevMat == 2) {
			cutSkinLine(topStartV, topUvStart, topEndV, toUV, topTe, topParams, topUVs, Tin, false, scl);
			if (mat2BorderVertex > -1) {
				auto vit = scl.deepVertsTris.begin();
				do {
					++vit;

				} while (*vit != mat2BorderVertex);
				++vit;
				int topV = -1;
				for(auto dp : _deepBed){
					if (dp.second.deepMtVertex == *vit)
						topV = dp.first;
				}
				if (topV < 0)
					throw(std::logic_error("Program error at material 2-border vertex.\n"));
				mat2BorderSplit(mat2BorderVertex, mat2BorderTexture, topV);
			}
		}
		else
			cutDeepSurface(deepStartV, deepUvStart, deepEndV, toUV, topTe, topParams, topUVs, scl);
	}
	if (nEdges > 499)
		return DBL_MAX;
	return len;
}

void deepCut::mat2BorderSplit(int borderV, int borderTx, int incisionTopV) {  // At material 2 triangle edge with mat 1, or > 4 border, this split must be made.
	// The deep vertex of an incision end must directly connect to this borderV for a topologically correct deep split.  This could have been prevented
	// if model had material 3 edges around all material 2 triangles at the borders, usually with boundary or bone triangles.
	// This convention was not adopted prompting the need for this routine.
	auto dbit = _deepBed.find(incisionTopV);
	if (dbit == _deepBed.end())
		throw(std::logic_error("Program error: mat2BorderSplit() called with improper setup.\n"));
	int oppVert = _mt->addVertices(1), oppTex = _mt->addTexture(), bottomVertex = dbit->second.deepMtVertex, topTx = -1;
	assert(oppVert == _vbt->_vertexTets.size());
	_vbt->_vertexTets.push_back(_vbt->getVertexTetrahedron(incisionTopV));
	_vbt->_barycentricWeights.push_back(*_vbt->getVertexWeight(incisionTopV));
	Vec3f V;
	_mt->getVertexCoordinate(incisionTopV, V.xyz);
	_mt->setVertexCoordinate(oppVert, V.xyz);
	_mt->findAdjacentTriangles(true);
	std::vector<materialTriangles::neighborNode> nei;
	_mt->getNeighbors(incisionTopV, nei);
	auto nit = nei.begin();
	for (; nit != nei.end(); ++nit) {
		if (nit->vertex == bottomVertex)
			break;
	}
	while (nit->vertex != borderV) {
		++nit;
		if (nit == nei.end())
			nit = nei.begin();
		auto tr = _mt->triangleVertices(nit->triangle);
		auto triTx = _mt->triangleTextures(nit->triangle);
		int i;
		for (i = 0; i < 3; ++i) {
			if (tr[i] == incisionTopV) {
				if (_mt->triangleMaterial(nit->triangle) == 2) {
					topTx = triTx[i];
					Vec2f tx;
					auto fp = _mt->getTexture(topTx);
					tx.set(fp[0], fp[1]);
					_mt->setTexture(oppTex, tx.xy);
				}
				tr[i] = oppVert;
				triTx[i] = oppTex;
				break;
			}
		}
		if (i > 2)
			throw(std::logic_error("Program error in cutting surface path.\n"));
	}
	// COURT - could make graphics better.
	int bottomTx = 0;  // stub
	int v3[3] = { borderV, incisionTopV, bottomVertex }, tex[3] = { borderTx, topTx, bottomTx };
	int firstTri = _mt->addTriangle(v3, 6, tex);  // splitMat
	v3[1] = v3[2];
	v3[2] = oppVert;
	tex[1] = tex[2];
	tex[2] = oppTex;
	int oppTri = _mt->addTriangle(v3, 6, tex);
	_mt->findAdjacentTriangles(true);
	unsigned int te = _mt->triAdjs(firstTri)[1];
	int firstTex = _mt->triangleTextures(te >> 2)[te & 3];
	_mt->triangleTextures(firstTri)[2] = firstTex;
	te = _mt->triAdjs(oppTri)[1];
	bottomTx = _mt->triangleTextures(te >> 2)[((te & 3) + 1) % 3];
	_mt->triangleTextures(oppTri)[1] = bottomTx;
}

void deepCut::cutDeepSurface(int startV, Vec2d& startUV, int endV, Vec2d& endUV, std::vector<unsigned int>& te, std::vector<float>& params, std::vector<Vec2d>& UVs, surfaceCutLine& scl) {
	if ((scl.deepVertsTris.empty() || startV != scl.deepVertsTris.back()) && startV > -1) {
		scl.deepVertsTris.push_back(startV);
		scl.deepUVs.push_back(startUV);
	}
	auto pit = params.begin();
	auto uvit = UVs.begin();
	for (auto& t : te) {
		Vec3f gridLocus;
		int tri = t >> 2, edge = t & 3;
		int* st = _mt->triangleVertices(tri);
		int tet = _vbt->parametricEdgeTet(st[edge], st[(edge + 1) % 3], *pit, gridLocus);
//		int tet = parametricMTedgeTet(tri, edge, *pit, gridLocus);
		int bedVertex = _mt->splitTriangleEdge(tri, edge, *pit);
		assert(bedVertex == _vbt->_vertexTets.size());
		_vbt->_vertexTets.push_back(tet);
		_vbt->_barycentricWeights.push_back(Vec3f());
		_vbt->gridLocusToBarycentricWeight(gridLocus, _vbt->_tetCentroids[tet], _vbt->_barycentricWeights.back());
		scl.deepVertsTris.push_back(bedVertex);
		scl.deepUVs.push_back(*uvit);
		++pit;
		++uvit;
	}
	if (endV > -1) {
		scl.deepVertsTris.push_back(endV);
		scl.deepUVs.push_back(endUV);
	}
}

void deepCut::cutSkinLine(int startV, Vec2d &startUV, int endV, Vec2d& endUV, std::vector<unsigned int>& te, std::vector<float>& params, std::vector<Vec2d>& UVs, bool Tin, bool Tout, surfaceCutLine& scl) {
	assert(UVs.size() == te.size());
	if (startV == _previousSkinTopEnd)
		Tin = true;
	if(endV == _loopSkinTopBegin)
		Tout = true;
	auto pit = params.begin();
	std::list<int> topVerts, botVerts;
	std::list<Vec2d> botUVs;
	if (startV > -1) {  // usual situation unless startV was a material 2 border vertex, which won't have a deep bed vertex.
		topVerts.push_back(startV);
		auto dbit = _deepBed.find(startV);
		if (dbit == _deepBed.end() || dbit->second.deepMtVertex < 0)
			throw(std::logic_error("Program error in skinCutLine()."));
		botVerts.push_back(dbit->second.deepMtVertex);
		botUVs.push_back(startUV);
	}
	int nUV = 0;
	for (auto teit = te.begin(); teit != te.end(); ++teit) {
		int tri = *teit >> 2, edge = *teit & 3;
		float uv[2] = { 0.0f, 0.0f };
		if (edge < 1)
			uv[0] = *pit;
		else if (edge > 1)
			uv[1] = 1.0f - *pit;
		else {
			uv[1] = *pit;
			uv[0] = 1.0f - *pit;
		}
		int topVertex, bottomVertex;
		createFlapTopBottomVertices(tri, uv, topVertex, bottomVertex);
		topVerts.push_back(topVertex);
		botVerts.push_back(bottomVertex);
		botUVs.push_back(UVs[nUV++]);
		++pit;
	}
	if (endV > -1) {
		topVerts.push_back(endV);
		auto dbit = _deepBed.find(endV);
		if (dbit == _deepBed.end() || dbit->second.deepMtVertex < 0)
			throw(std::logic_error("Program error in skinCutLine()."));
		botVerts.push_back(dbit->second.deepMtVertex);
		botUVs.push_back(endUV);
	}
	_mt->findAdjacentTriangles(true);
	updateDeepSpatialCoordinates();
	if(!topDeepSplit_Sub(topVerts, botVerts, Tin, Tout))
		throw(std::logic_error("Program error in skinCutLine()."));
	auto bvit = botVerts.begin();
	auto uvit = botUVs.begin();
	for (; bvit != botVerts.end(); ++bvit) {
		scl.deepVertsTris.push_back(*bvit);
		scl.deepUVs.push_back(*uvit);
		++uvit;
	}
}

