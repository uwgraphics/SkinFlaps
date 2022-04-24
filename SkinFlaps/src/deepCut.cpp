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

bool deepCut::cutDeep()  // data already loaded in _deepPosts
{
	_previousSkinTopEnd = -1;
	_loopSkinTopBegin = -1;
	_holePolyLines.clear();
	if (_deepPosts.size() < 2)
		return false;
	// trim unused intersects from posts
	int prevMax = -1;
	std::vector<int> trimTi;
	trimTi.assign(_deepPosts.size(), -1);
	for (int j, i = _deepPosts.size() - 1; i > -1; --i) {
		for (j = _deepPosts[i].triIntersects.size() - 1; j > -1; --j) {
			if (_deepPosts[i].triIntersects[j].scl.rtiIndexTo > -1) {
				if (prevMax > j)
					j = prevMax;
				prevMax = _deepPosts[i].triIntersects[j].scl.rtiIndexTo;
				break;
			}
		}
		if (i == _deepPosts.size() - 1 && _deepPosts[i].triIntersects[0].scl.rtiIndexTo > j) {  // open back end
			assert(_deepPosts[i].triIntersects[0].scl.rtiPostTo == i);
			j = _deepPosts[i].triIntersects[0].scl.rtiIndexTo;
		}
		if(j< 0){
			if (prevMax > 0)
				j = prevMax;
			else
				return false;
		}
		auto ti = &_deepPosts[i].triIntersects;
		trimTi[i] = j;
	}
	// made all inter-ray surface connections
	if (!_deepPosts.front().closedEnd)
		if (!connectOpenEnd(0, trimTi[0]))
			return false;
	if (!_deepPosts.back().closedEnd)
		if (!connectOpenEnd(_deepPosts.size() - 1, trimTi[_deepPosts.size() - 1]))
			return false;
	// now make all single ray open scl connections
	for (int i = 0; i < _deepPosts.size(); ++i) {
		auto& dti = _deepPosts[i].triIntersects;
		for (int j = 1; j < dti.size() - 1; j += 2) {
			if (j >= trimTi[i])
				break;
			if (dti[j].scl.rtiIndexTo < 0 && i > 0) {
				for (int k = j + 1; k < dti.size(); k += 2) {
					if (k > trimTi[i])
						break;
					double lowV;
					if (surfacePath(dti[j], dti[k], nullptr, false, lowV) < DBL_MAX) {
						dti[j].scl.rtiPostTo = i;
						dti[j].scl.rtiIndexTo = k;
						dti[j].scl.lowestV = lowV;
						break;
					}
				}
			}
		}
		for (int j = 2; j < dti.size() - 1; j += 2) {
			if (j > trimTi[i])
				break;
			if (dti[j].scl.rtiIndexTo < 0 && i < _deepPosts.size() - 1) {
				for (int k = j - 1; k > -1; k -= 2) {
					double lowV;
					if (surfacePath(dti[j], dti[k], nullptr, false, lowV) < DBL_MAX) {
						dti[j].scl.rtiPostTo = i;
						dti[j].scl.rtiIndexTo = k;
						dti[j].scl.lowestV = lowV;
						break;
					}
				}
			}
		}
	}
	// now organize all skin cut sequences
	std::map<std::pair<int, int>, int > rtiHits;
	rayTriangleIntersect* last, * rip = &_deepPosts[0].triIntersects[0];
	last = rip;
	std::list< rayTriangleIntersect* > outerLoop;
	bool closedLoop = true;
	do {
		outerLoop.push_back(rip);
		auto pr = rtiHits.insert(std::make_pair(std::make_pair(rip->postNum, rip->rayIndex), 1));
		if (!pr.second)
			++pr.first->second;
		pr = rtiHits.insert(std::make_pair(std::make_pair(rip->scl.rtiPostTo, rip->scl.rtiIndexTo), 1));
		if (!pr.second)
			++pr.first->second;
		last = rip;
		rip = &_deepPosts[last->scl.rtiPostTo].triIntersects[last->scl.rtiIndexTo];
		if (rip->scl.rtiIndexTo < 0) {  // solid run on a ray possibly closed end. Always alternates with a scl.
			closedLoop = false;
			outerLoop.push_back(nullptr);
			if (last->scl.rtiIndexTo & 1) {
				if (last->scl.rtiIndexTo == 1 && last->scl.rtiPostTo < 1)  // closed front end
					break;
				else  // start back on the bottom
					rip = &_deepPosts[last->scl.rtiPostTo].triIntersects[last->scl.rtiIndexTo - 1];
			}
			else
				rip = &_deepPosts[last->scl.rtiPostTo].triIntersects[last->scl.rtiIndexTo + 1];
		}
	} while (rip != outerLoop.front());
	if (!closedLoop && outerLoop.back() != nullptr && outerLoop.back()->scl.rtiPostTo == 0 && outerLoop.back()->scl.rtiIndexTo == 0) {
		auto lit = outerLoop.end(); --lit;
		while (*lit != nullptr) {
			outerLoop.push_front(outerLoop.back());
			--lit;
			outerLoop.pop_back();
		}
	}
	// now look for interior skin rings which need to be cut
	for (int i = 1; i < _deepPosts.size() - 1; ++i) {
		for (int j = 1; j < _deepPosts[i].triIntersects.size() - 2; j += 2) {
			if (j > trimTi[i])
				break;
			if (rtiHits.find(std::make_pair(i, j)) != rtiHits.end())
				continue;
			// closed loop inside outer loop
			outerLoop.push_back(nullptr);
			rip = &_deepPosts[i].triIntersects[j];
			outerLoop.push_back(rip);
			if(rip->scl.rtiPostTo != rip->postNum || rip->scl.rtiIndexTo != j + 1)
				throw(std::logic_error("Error finding closed interior skin loop within a deep cut-\n"));
			auto pr = rtiHits.insert(std::make_pair(std::make_pair(i, j), 2));
			rip = &_deepPosts[rip->scl.rtiPostTo].triIntersects[rip->scl.rtiIndexTo];
			outerLoop.push_back(rip);
			pr = rtiHits.insert(std::make_pair(std::make_pair(rip->postNum, rip->rayIndex), 2));
			if (!pr.second || rip->scl.rtiPostTo != i || rip->scl.rtiIndexTo != j)
				throw(std::logic_error("Error finding closed interior skin loop within a deep cut-\n"));
		}
	}
	// get all rti with an scl entering and leaving, otherwise the rti blind ends into a solid
	for (auto rit = rtiHits.begin(); rit != rtiHits.end(); ) {
		if (rit->second != 2)
			rit = rtiHits.erase(rit);
		else
			++rit;
	}
	_preDeepCutVerts = _mt->numberOfVertices();
	auto punchDeepVert = [&](rayTriangleIntersect& ri) ->bool{
		if (ri.deepVert > -1)
			return true;
		int mat = _mt->triangleMaterial(ri.triangle);
		float uv[2] = { (float)ri.uv[0], (float)ri.uv[1] };
		if (mat == 2) {
			long topVertex, bottomVertex;
			createFlapTopBottomVertices(ri.triangle, uv, topVertex, bottomVertex);
			ri.mat2Vert = topVertex;
			ri.deepVert = bottomVertex;
		}
		else if (mat == 1 || (mat > 4 && mat < 10)) {
			Vec3f gridLocus, bw;
			long tet = parametricMTtriangleTet(ri.triangle, uv, gridLocus);
			_vbt->gridLocusToBarycentricWeight(gridLocus, *_vbt->tetCentroid(tet), bw);
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
	for (int n = _deepPosts.size(), i = 0; i < n; ++i) {
		for (int j = 0; j <= trimTi[i]; ++j) {
			if (!punchDeepVert(_deepPosts[i].triIntersects[j]))
				return false;
		}
	}
	_mt->findAdjacentTriangles(true, false);
	getDeepSpatialCoordinates();  // after punches need to redo
	// create bilinear surface precomputations for each patch per nVidia subroutine
	for (int i = 1; i < _deepPosts.size(); ++i) {
		int intersectLevel = trimTi[i - 1];	//  _deepPosts[i - 1].triIntersects.size() - 1;  // keep bilinear surface as square as possible
		if (intersectLevel > trimTi[i])  // _deepPosts[i].triIntersects.size() - 1)
			intersectLevel = trimTi[i];  //  _deepPosts[i].triIntersects.size() - 1;

		Vec3d P00 =_deepXyz[_deepPosts[i - 1].triIntersects[intersectLevel].deepVert];
		Vec3d P10 = _deepXyz[_deepPosts[i].triIntersects[intersectLevel].deepVert];
		Vec3d P01 = _deepXyz[_deepPosts[i - 1].triIntersects.front().deepVert];
		Vec3d P11 = _deepXyz[_deepPosts[i].triIntersects.front().deepVert];
		makeBilinearPatch(P00, P10, P01, P11, _deepPosts[i].bl);
	}
	// make skin incisions and collect all connected deep surface lines on incised surface for split later
	std::list<std::list<long> > surfacePolyLines;  // deepVert lines
	surfacePolyLines.push_back(std::list <long>());
	bool setTopLoopBegin = true;
	for (auto oit = outerLoop.begin(); oit != outerLoop.end(); ++oit) {
		if (*oit == nullptr) {
			if (surfacePolyLines.empty() || !surfacePolyLines.back().empty())
				surfacePolyLines.push_back(std::list<long>());
			setTopLoopBegin = true;
			continue;
		}
		double stub;
		if (surfacePath(**oit, _deepPosts[(*oit)->scl.rtiPostTo].triIntersects[(*oit)->scl.rtiIndexTo], nullptr, true, stub) == DBL_MAX)
			return false;
		if (setTopLoopBegin && (*oit)->mat2Vert > -1)
			_loopSkinTopBegin = (*oit)->mat2Vert;
		setTopLoopBegin = false;
		_previousSkinTopEnd = _deepPosts[(*oit)->scl.rtiPostTo].triIntersects[(*oit)->scl.rtiIndexTo].mat2Vert;

		_mt->findAdjacentTriangles(true, false);
		updateDeepSpatialCoordinates();  // after punches need to redo
		auto dit = (*oit)->scl.deepVerts.begin();
		if (!surfacePolyLines.back().empty() && surfacePolyLines.back().back() == *dit)
			++dit;
		for (; dit != (*oit)->scl.deepVerts.end(); ++dit) {
			surfacePolyLines.back().push_back(*dit);
		}
	}
	if (surfacePolyLines.front().empty())
		surfacePolyLines.pop_front();
	if (surfacePolyLines.back().empty())
		surfacePolyLines.pop_back();
	// splice any surfacePolyLines that need it
	if(surfacePolyLines.size()> 1){
		auto dpit0 = surfacePolyLines.begin();
		auto dpit1 = dpit0;
		while (dpit0 != surfacePolyLines.end()) {
			++dpit1;
			if (dpit1 == surfacePolyLines.end()) {
				++dpit0;
				dpit1 = dpit0;
				continue;
			}
			if (dpit0->front() == dpit1->back()) {
				dpit0->pop_front();
				dpit0->splice(dpit0->begin(), *dpit1);
				surfacePolyLines.erase(dpit1);
				dpit1 = dpit0;
				continue;
			}
			if (dpit0->back() == dpit1->front()) {
				dpit0->pop_back();
				dpit0->splice(dpit0->end(), *dpit1);
				surfacePolyLines.erase(dpit1);
				dpit1 = dpit0;
			}
		}
	}
	// now get all interior solid vertices on posts
	std::set<long> nonDupedVertices;
	for (int np = _deepPosts.size(), i = 0; i < np; i++) {
		auto &dti = _deepPosts[i].triIntersects;
		bool openEnd, openSolid;
		if (i < 1 && !dti.back().scl.deepVerts.empty())
			openEnd = true;
		else if (i == _deepPosts.size() - 1 && !dti.front().scl.deepVerts.empty())
			openEnd = true;
		else
			openEnd = false;
		for (int j = 0; j < dti.size() - 1; j += 2) {
			if (rtiHits.find(std::make_pair(i, j)) == rtiHits.end() || rtiHits.find(std::make_pair(i, j + 1)) == rtiHits.end())
				openSolid = true;
			else
				openSolid = false;
			if (openEnd && !openSolid)  // don't bother making this deep cut line as just one large polygon
				continue;
			getDeepCutLine(dti[j], dti[j + 1]);
			if(openSolid){  // this segment blind ends into a solid
				nonDupedVertices.insert(dti[j].deepVert);
				for (auto& dv : dti[j].dcl)
					nonDupedVertices.insert(dv);
				nonDupedVertices.insert(dti[j + 1].deepVert);
			}
		}
	}
	// now cut the quads
	updateDeepSpatialCoordinates();  // necessary for hole finding in quads
	for (int np = _deepPosts.size(), i = 1; i < np; i++)
		deepCutQuad(i);
	struct posTex {
		int pos;
		int tex;
	} ptex;
	ptex.pos = -1;
	ptex.tex = -1;
	std::map<int, posTex> oppositeSideVertices;  // relates vertices on first side of cut to opposite side which are created here.
	auto dupVertex = [&](long vert) {
		auto pr = oppositeSideVertices.insert(std::make_pair(vert, ptex));
		if (pr.second) {
			if (nonDupedVertices.find(vert) != nonDupedVertices.end())  // no dup
				pr.first->second.pos = vert;
			else {  // all deep vertices so don't add to _deepBed
				pr.first->second.pos = _mt->cloneVertex(vert);
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
	// do surfaceCutLine splits of original surface
	auto splitSurface = [&](std::list<std::list<long> >& spl) {
		for (auto& sp : spl) {
			auto vit = sp.begin();
			++vit;
			if (sp.front() == sp.back())  // closed polygon
				sp.push_back(oppositeSideVertices[*vit].pos);  // new end point after first replace
			unsigned long adj = 0xffffffff;
			long vNow = *vit, * tr;
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
				long loopV = -1;
				unsigned long topAdj, botAdj = _mt->triAdjs(adj >> 2)[adj & 3];
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
						long* txp = _mt->triangleTextures(adj >> 2);
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
	// clockwise listing of holes must be reversed to split original surface correctly
	for (auto& hpl : _holePolyLines)
		hpl.reverse();
	splitSurface(_holePolyLines);
	// fill in new cut surface triangles
	auto addDeepTriangles = [&](std::vector<materialTriangles::matTriangle>& tris) {
		for(auto &t : tris){
			long V[3], T[3];
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
	if (_mt->findAdjacentTriangles(true, false))
		return false;  // should throw
	return true;
}

void deepCut::getDeepCutLine(rayTriangleIntersect &top, rayTriangleIntersect &bot) {
	Vec3f V, C = Vec3f(top.intersect._v), D = Vec3f(bot.intersect._v) - C;
	float len = D.length();
	int n = (int)(D.length() * _vbt->getTetUnitSizeInv());
	D /= (float)n;
	for (int i = 1; i < n - 1; ++i) {
		V = C + D * (float)i;
		int tet;
		Vec3f baryWeight;
		if (!interiorSpatialTet(V, tet, baryWeight))
			continue;
		top.dcl.push_back(_mt->addVertices(1));
		_mt->setVertexCoordinate(top.dcl.back(), V._v);  // texture not set
		assert(_vbt->_vertexTets.size() == top.dcl.back());
		_vbt->_vertexTets.push_back(tet);
		_vbt->_barycentricWeights.push_back(baryWeight);
	}
}

bool deepCut::interiorSpatialTet(const Vec3f pos, int& tet, Vec3f& baryWeight) {

//	std::chrono::time_point<std::chrono::system_clock> start, end;
//	start = std::chrono::system_clock::now();

	tet = -1;  // only found once so no write contention. Could make atomic.
	int n = _vbt->tetNumber();  // , i
//	for (i = _vbt->firstInteriorTet(); i < n; ++i) {
	tbb::parallel_for(tbb::blocked_range<std::size_t>(_vbt->firstInteriorTet(), _vbt->tetNumber()), [&](tbb::blocked_range<size_t> r) {
		for (int i = r.begin(); i != r.end(); ++i) {
			boundingBox<float> bb;
			bb.Empty_Box();
			Vec3f* vp[4];
			long* tn = _vbt->_tetNodes[i].data();
			for (int j = 0; j < 4; ++j) {
				vp[j] = &_vbt->_nodeSpatialCoords[tn[j]];
				bb.Enlarge_To_Include_Point(vp[j]->_v);
			}
			if (bb.Outside(pos._v))
				continue;
			Mat3x3f M;
			M.Initialize_With_Column_Vectors(*vp[1] - *vp[0], *vp[2] - *vp[0], *vp[3] - *vp[0]);
			Vec3f R = M.Robust_Solve_Linear_System(pos - *vp[0]);
			if (R[0] <= 0.0f || R[1] <= 0.0f || R[2] <= 0.0f || R[0] >= 1.0f || R[1] >= 1.0f || R[2] >= 1.0f || R[0] + R[1] + R[2] >= 1.0f)
				continue;
			tet = i;
			baryWeight = R;
			break;  // internal tets are guaranteed unique barring rare inversions
		}
		if (tet > -1)
			oneapi::tbb::task_group_context().cancel_group_execution();
	});
//	}

//	end = std::chrono::system_clock::now();
//	std::chrono::duration<double> elapsed_seconds = end - start;
//	std::cout << "interiorSpatialTet() took " << elapsed_seconds.count() << "\n";

	// worst case seems to be 0.008 seconds without tbb. tbb version worst 0.001 but some 15x faster on my desktop
	if (tet < 0)
		return false;
	return true;
}

bool deepCut::connectOpenEnd(int postNum, int &interpostEnd) {
	auto post = &_deepPosts[postNum].triIntersects;
	if (postNum < 1) {  // front
		int j, i;
		std::vector<int> pathBack;
		for (i = 0; i < post->size(); ) {
			for (j = i + 1; j < post->size(); j += 2) {
				double lowV;
				if (surfacePath(_deepPosts[postNum].triIntersects[j], _deepPosts[postNum].triIntersects[i], nullptr, false, lowV) < DBL_MAX) {
						(*post)[j].scl.rtiPostTo = postNum;
					(*post)[j].scl.rtiIndexTo = i;
					(*post)[j].scl.lowestV = lowV;
					int k = j - 1;
					while (k > i) {
						pathBack.push_back(k);
						k -= 2;
					}
					if (j >= interpostEnd) {
						interpostEnd = j;
						i = post->size();  // stop looking
					}
					else
						i = j + 1;
					break;
				}
			}
			if (j > post->size())
				return false;
		}
		// now do any back segments as well
		for (auto& k : pathBack) {
			double lowV;
			if (surfacePath(_deepPosts[postNum].triIntersects[k - 1], _deepPosts[postNum].triIntersects[k], nullptr, false, lowV) < DBL_MAX) {
				(*post)[k - 1].scl.rtiPostTo = postNum;
				(*post)[k - 1].scl.rtiIndexTo = k;
				(*post)[k - 1].scl.lowestV = lowV;
			}
			else
				return false;
		}
	}
	else {  // back end
		int j, i;
		std::vector<int> pathBack;
		for (i = 0; i < post->size(); ) {
			for (j = i + 1; j < post->size(); j += 2) {
				double lowV;
				if (surfacePath(_deepPosts[postNum].triIntersects[i], _deepPosts[postNum].triIntersects[j], nullptr, false, lowV) < DBL_MAX) {
					(*post)[i].scl.rtiPostTo = postNum;
					(*post)[i].scl.rtiIndexTo = j;
					(*post)[i].scl.lowestV = lowV;
					int k = j - 1;
					while (k > i) {
						pathBack.push_back(k);
						k -= 2;
					}
					if (j >= interpostEnd) {
						interpostEnd = j;
						i = post->size();  // stop looking
					}
					else
						i = j + 1;
					break;
				}
			}
			if (j > post->size())
				return false;
		}
		// now do any back segments as well
		for (auto &k : pathBack) {
			double lowV;
			if (surfacePath(_deepPosts[postNum].triIntersects[k], _deepPosts[postNum].triIntersects[k-1], nullptr, false, lowV) < DBL_MAX) {
				(*post)[k].scl.rtiPostTo = postNum;
				(*post)[k].scl.rtiIndexTo = k - 1;
				(*post)[k].scl.lowestV = lowV;
			}
			else
				return false;
		}
	}
	return true;
}

bool deepCut::deepCutQuad(int postNum) {
	// get CW quad border
	auto& pti0 = _deepPosts[postNum - 1].triIntersects;
	auto& pti1 = _deepPosts[postNum].triIntersects;
	std::list<long> poly, tmp;
	poly.assign(pti0[0].scl.deepVerts.begin(), pti0[0].scl.deepVerts.end());
	bool openEnd = false, lastDeep = false;
	int interPostIdx = -1;
	if (!pti1[0].scl.deepVerts.empty() && postNum == _deepPosts.size() - 1) // open end
		openEnd = true;
	// find first interpost return
	for (int i = 1; i < pti1.size(); i += 2) {
		if (pti1[i].scl.rtiPostTo == postNum - 1) {
			interPostIdx = i;
			break;
		}
	}
	if(interPostIdx < 0 || pti1[interPostIdx].scl.rtiIndexTo < 0)
		throw(std::logic_error("Topological connection error in deepQuadCut()."));
	for (int i = 0; i != interPostIdx; ) {
		if ((openEnd || i&1) && !pti1[i].scl.deepVerts.empty()) {
			tmp.assign(pti1[i].scl.deepVerts.begin(), pti1[i].scl.deepVerts.end());
			if (!lastDeep) {
				assert(poly.back() == tmp.front());
				poly.pop_back();
			}
			poly.splice(poly.end(), tmp);
			lastDeep = false;
			i = pti1[i].scl.rtiIndexTo;
		}
		else {
			if (lastDeep)
				poly.push_back(pti1[i].deepVert);
			if (i > interPostIdx) {
				--i;
				tmp.assign(pti1[i].dcl.begin(), pti1[i].dcl.end());
				tmp.reverse();
			}
			else {
				tmp.assign(pti1[i].dcl.begin(), pti1[i].dcl.end());
				++i;
			}
			poly.splice(poly.end(), tmp);
			lastDeep = true;
		}
	}
	tmp.assign(pti1[interPostIdx].scl.deepVerts.begin(), pti1[interPostIdx].scl.deepVerts.end());
	if (!lastDeep) {
		assert(poly.back() == tmp.front());
		poly.pop_back();
	}
	poly.splice(poly.end(), tmp);
	lastDeep = false;
	interPostIdx = pti1[interPostIdx].scl.rtiIndexTo;
	if (!pti0[interPostIdx].scl.deepVerts.empty() && postNum == 1) // open front
		openEnd = true;
	else
		openEnd = false;
	for (int i = interPostIdx; i != 0; ) {
		if ((openEnd || (i&1) < 1) && !pti0[i].scl.deepVerts.empty()) {
			tmp.assign(pti0[i].scl.deepVerts.begin(), pti0[i].scl.deepVerts.end());
			if (!lastDeep) {
				assert(poly.back() == tmp.front());
				poly.pop_back();
			}
			poly.splice(poly.end(), tmp);
			lastDeep = false;
			i = pti0[i].scl.rtiIndexTo;
		}
		else {
			if (lastDeep)
				poly.push_back(pti0[i].deepVert);
			if (i <= interPostIdx) {
				--i;
				tmp.assign(pti0[i].dcl.begin(), pti0[i].dcl.end());
				tmp.reverse();
			}
			else {
				tmp.assign(pti0[i].dcl.begin(), pti0[i].dcl.end());
				++i;
			}
			poly.splice(poly.end(), tmp);
			lastDeep = true;
		}
	}
	if (poly.back() == poly.front())
		poly.pop_back();
	poly.reverse(); // to make CCW
	std::vector<std::pair<long, Vec2d> > deepOuterPolygon;
	deepOuterPolygon.reserve(poly.size());
	// put in missing uv coords for this quad
	Vec2d minC(DBL_MAX, DBL_MAX), maxC(-DBL_MAX, -DBL_MAX);
	const bilinearPatch& bl = _deepPosts[postNum].bl;
	std::vector<Vec2d> bUv;
	bUv.reserve(poly.size());
	Vec3d lastN;
	bilinearNormal(0.0, 1.0, bl, lastN);  // first poly point
	for (auto pit = poly.begin(); pit != poly.end(); ++pit) {
		Vec2d R;
		// get bilinear uv for deep non-surface points
		Vec3d E((const float(&)[3]) * _mt->vertexCoordinate(*pit));
		E -= lastN;
		double rp[2];
		Vec2d fp[2];
		int nSols = bilinearRayIntersection(E, lastN * 2.0, bl, rp, fp);
		if (nSols < 1)
			throw(std::logic_error("Bilinear intersection error in deepQuadCut()."));
		if (nSols > 1) {
			if (abs(rp[1] - 0.5) < abs(rp[0] - 0.5))
				fp[0] = fp[1];
		}
		bUv.push_back(fp[0]);
		deepOuterPolygon.push_back(std::make_pair(*pit, fp[0]));
		bilinearNormal(fp[0].X, fp[0].Y, bl, lastN);
		if (minC.X > fp[0].X)
			minC.X = fp[0].X;
		if (maxC.X < fp[0].X)
			maxC.X = fp[0].X;
		if (minC.Y > fp[0].Y)
			minC.Y = fp[0].Y;
		if (maxC.Y < fp[0].Y)
			maxC.Y = fp[0].Y;
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
	std::list< std::vector<std::pair<long, Vec2d> > > holes;
	findCutInteriorHoles(bl, bUv, deepOuterPolygon, holes);
	for (auto hole : holes)
		holesSize += hole.size();
	holeTx.reserve(holesSize);
	std::list< std::vector<Vec2d> > holeUvs;
	for (auto hole : holes) {
		holeUvs.push_back(std::vector<Vec2d>());
		for (auto hp = hole.begin(); hp != hole.end(); ++hp) {
			holeUvs.back().push_back(hp->second);
			holeTx.push_back(_mt->addTexture());
			float tx[2] = { (float)hp->second.X, (float)hp->second.Y };
			_mt->setTexture(holeTx.back(), tx);
		}
	}
	// do Delaunay triangulation of border with interior points
	int nOuterPolygon = (int)bUv.size();
	// recursively add internal Delaunay vertices to mtFace
	double incr = 2.0/ pti0.front().scl.deepVerts.size();
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
				if (!holeUvs.empty()) {
					auto hole = holeUvs.begin();
					while (hole != holeUvs.end()) {
						if (ip.insidePolygon2d(V2, *hole))
							break;
						++hole;
					}
					if (hole != holeUvs.end())
						continue;
				}
				Vec3d Nd = bl.P00 + bl.e10 * V2.X + bl.e00 * V2.Y + (bl.e11 - bl.e00) * (V2.X * V2.Y);
				int tet;
				Vec3f N(Nd._v), baryWeight;
				if (!interiorSpatialTet(Vec3f(N._v), tet, baryWeight))
					continue;
				intPt.newV = _mt->addVertices(1);
				_mt->setVertexCoordinate(intPt.newV, N._v);  // texture not set
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
	for (auto hole : holeUvs) {
		for (int r = hole.size(), i = 0; i < r; ++i)
			DelPts.push_back(V2d().make(hole[i].X, hole[i].Y));
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
	for (auto hole : holeUvs) {
		for (int i = 1; i < hole.size(); ++i)
			edges.push_back(CDT::Edge(i - 1 + edgeSize, i + edgeSize));
		edges.push_back(CDT::Edge(hole.size() - 1 + edgeSize, edgeSize));
		edgeSize += hole.size();
	}
	cdt.insertEdges(edges);
	cdt.eraseOuterTrianglesAndHoles();
	_deepPosts[postNum].quadTriangles.clear();
	_deepPosts[postNum].quadTriangles.reserve(cdt.triangles.size());
	for (size_t n = cdt.triangles.size(), i = 0; i < n; ++i) {
		materialTriangles::matTriangle tri;
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
					if (v < hole.size()) {
						tri.v[2 - j] = hole[v].first;
						break;
					}
					else
						v -= hole.size();
				}
			}
			else {
				tri.v[2 - j] = ceiVec[v - nOuterPolygon - holesSize].newV;
				tri.tex[2 - j] = ceiVec[v - nOuterPolygon - holesSize].tx;
			}
		}
		_deepPosts[postNum].quadTriangles.push_back(tri);
	}
	return true;
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
		long botTet = deepPointTetWeight(dbit, bw);
		if (botTet < 0) {
			_deepXyz.clear();  // COURT - this is bad.  Input deepBed needs to be recomputed.
			return false;
		}
		else {
			Vec3f v;
			_vbt->getBarycentricTetPosition(botTet, bw, v);
			_deepXyz[dbit->first] = Vec3d(v._v);
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
		}
		Vec3f bw;
		long botTet = deepPointTetWeight(dbit, bw);
		if (botTet < 0) {  // COURt if there are a lot of these or any occur over an important area of the model should recompute deep bed.
			std::cout << "Deep bed point at vertex " << dbit->first << " with deep vertex " << dbit->second.deepMtVertex << " not in a tet.\n";
		}
		else {
			Vec3f v;
			_vbt->getBarycentricTetPosition(botTet, bw, v);
			_deepXyz[dbit->first] = Vec3d(v._v);
		}
	}
	boundingBox<double> bb;
	bb.Empty_Box();
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i){
		if (_deepXyz[i].X < DBL_MAX)  // don't include invalid vertex
			bb.Enlarge_To_Include_Point(_deepXyz[i]._v);
	}
	Vec3d vn, vx;
	bb.Maximum_Corner(vx._v);
	bb.Minimum_Corner(vn._v);
	_maxSceneSize = (float)((vx - vn).length());
	return true;
}

void deepCut::getDeepPosts(std::vector<Vec3f>& xyz, std::vector<Vec3f>& nrm) {
	xyz.clear();
	nrm.clear();
	xyz.reserve(_deepPosts.size());
	nrm.reserve(_deepPosts.size());
	for (int n = _deepPosts.size(), i = 0; i < n; ++i) {
		xyz.push_back(Vec3f(_deepPosts[i].triIntersects[0].intersect._v));
		nrm.push_back(Vec3f(_deepPosts[i].triIntersects.back().intersect._v));
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
	char str[200];
	if (fp->numberOfPosts() < 2) {
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
			sprintf(str, "Post number %d needs direction adjustment.", i + 1);
			ffg->sendUserMessage(str, "Please correct deep cut-");
			return false;
		}
	}
	for (int i = 1; i < n; ++i) {
		if (!topConnectToPreviousPost(i)) {
			sprintf(str, "Post number %d has no top side connection to previous post.\nDelete it and try again", i + 1);
			ffg->sendUserMessage(str, "Please correct deep cut-");
			return false;
		}
	}
	for (int i = 1; i < n; ++i) {
		if (!deepConnectToPreviousPost(i)) {
			sprintf(str, "Post number %d has no bottom connection to previous post.\nAdjust its post direction or previous post direction.", i + 1);
			ffg->sendUserMessage(str, "Please correct deep cut-");
			return false;
		}
	}
	if (n > 1) {  // no corrections necessary til 2
		for (int i = 1; i < n; ++i) {
			if (preventPreviousCrossover(i) > 0) {
				sprintf(str, "Solid post line %d intersects a previous cut.\nPlease adjust it's direction.", i + 1);
				ffg->sendUserMessage(str, "Please correct deep cut-");
				return false;
			}
		}
	}
//	positions.clear();  // section no longer necessary as fence must be corrected on entry to pass this far
//	normals.clear();
//	positions.reserve(n);
//	normals.reserve(n);
//	for (int i = 0; i < n; ++i) {
//		positions.push_back(Vec3f(_deepPosts[i].triIntersects.front().intersect._v));
//		normals.push_back(-Vec3f(_deepPosts[i].rayDirection._v));
//	}
//	fp->updatePosts(positions, normals);
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
	long* tr = _mt->triangleVertices(triangle);
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
	if (surfacePath(_deepPosts[postNum - 1].triIntersects[0], _deepPosts[postNum].triIntersects[0], nullptr, false, lowV) == DBL_MAX)
		return false;
	auto fp = &_deepPosts[postNum - 1].triIntersects[0].scl;
	fp->rtiPostTo = postNum;
	fp->rtiIndexTo = 0;
	fp->lowestV = lowV;
	return true;
}

bool deepCut::deepConnectToPreviousPost(int postNum) {
	// find shortest return path
	double pathLen = DBL_MAX, minPath = DBL_MAX, lowV;
	int minJ = 10000;
	int i, j, n = _deepPosts[postNum].triIntersects.size();
	for (i = 1; i < n; i += 2) {
		for (j = 1; j < _deepPosts[postNum - 1].triIntersects.size(); j += 2) {
			pathLen = surfacePath(_deepPosts[postNum].triIntersects[i], _deepPosts[postNum - 1].triIntersects[j], nullptr, false, lowV);
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
		return true;
	}
	return false;
}

bool deepCut::connectToPreviousPost(int postNum) {
	assert(postNum > 0);
	double lowV;
	if (surfacePath(_deepPosts[postNum - 1].triIntersects[0], _deepPosts[postNum].triIntersects[0], nullptr, false, lowV) == DBL_MAX)
		return false;
	auto fp = &_deepPosts[postNum - 1].triIntersects[0].scl;
	fp->rtiPostTo = postNum;
	fp->rtiIndexTo = 0;
	fp->lowestV = lowV;
	// find shortest return path
	double pathLen = DBL_MAX, minPath = DBL_MAX;
	int minJ = 10000;
	int i, j, n = _deepPosts[postNum].triIntersects.size();
	for (i = 1; i < n; i += 2) {
		for (j = 1; j < _deepPosts[postNum - 1].triIntersects.size(); j += 2) {
			pathLen = surfacePath(_deepPosts[postNum].triIntersects[i], _deepPosts[postNum - 1].triIntersects[j], nullptr, false, lowV);
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
		return true;
	}
	return false;
}

bool deepCut::rayIntersectMaterialTriangles(const Vec3d& rayStart, const Vec3d& rayDirection, std::vector<rayTriangleIntersect>& intersects) {
	std::multimap<double, rayTriangleIntersect> rtiMap;
	boundingBox<double> bb, tbb;
	bb.Empty_Box();
	bb.Enlarge_To_Include_Point(rayStart._v);
	bb.Enlarge_To_Include_Point((rayStart + rayDirection * _maxSceneSize)._v);
	Vec3d P, T[3], N;
	// do slightly permissive find
	for (int n = _mt->numberOfTriangles(), j, i = 0; i < n; ++i) {
		long tm = _mt->triangleMaterial(i);
		if (tm == 3 || tm == 4 || tm < 0)  // only look for permissible deep cut triangles.
			continue;
		long* tr = _mt->triangleVertices(i);
		tbb.Empty_Box();
		for (j = 0; j < 3; ++j) {
			if (_deepXyz[tr[j]].X > 1e22)  // invalid vertex, thus so is this triangle
				break;
			T[j].set(_deepXyz[tr[j]]);
			tbb.Enlarge_To_Include_Point(T[j]._v);
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

long deepCut::addPeriostealUndermineTriangle(const int triangle, const Vec3f &linePickDirection, bool incisionConnect)
{
	long perioTri = 0x7fffffff, mat = _mt->triangleMaterial(triangle);
	if (mat == 7 || mat == 8 || mat == 10)  // a periosteal triangle possibly in the middle of an undermine
		perioTri = triangle;
	else{
		std::vector<float> positions, rayParams;
		std::vector<int> rayTris;
		_mt->linePick(_mt->vertexCoordinate(_mt->triangleVertices(triangle)[0]), linePickDirection._v, positions, rayTris, rayParams, 7);
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
	for (int i = 0; i < 2; ++i) {
		if (uVec[i] < -1e32 || uVec[i] > 8.0 || uVec[i] < -8.0)  // last 2 throw out wild extrapolation results for open end processing
			continue;
		Vec3d pa = bl.P00 + bl.e10 * uVec[i];
		Vec3d T, pb = bl.P01 * (1.0 - uVec[i]) + bl.P11 * uVec[i];
		pb -= pa;
		T = rayStart - pa;
		Mat2x2d M(pb * pb, rayDir * pb, 0, -rayDir * rayDir);
		M.x[2] = -M.x[1];
		Vec2d B(T * pb, T*rayDir), R;
		R = M.Robust_Solve_Linear_System(B);
		double t1 = R.Y;
		double v1 = R.X;
		if (t1 > 0.0) {
			if (t1 <= 1.0)
				++nSols;
			else
				continue;
			// again throw out wild extrapolation results for open end processing
			if ((uVec[i] > 1.0 || uVec[i] < 0.0) && (v1 < -8.0 || v1 > 2.0)) {
				--nSols;
				continue;
			}
			rayParam[nSols-1] = t1;
			faceParams[nSols-1].X = uVec[i];
			faceParams[nSols-1].Y = v1;
		}
	}
	return nSols;
}

void deepCut::findCutInteriorHoles(const bilinearPatch& bl, const std::vector<Vec2d>& bUv, const std::vector<std::pair<long, Vec2d> >& deepOuterPolygon, std::list< std::vector<std::pair<long, Vec2d> > >& holes) {
	holes.clear();
	boundingBox<double> bb;
	bb.Empty_Box();
	for (auto& bp : deepOuterPolygon)
		bb.Enlarge_To_Include_Point(_deepXyz[bp.first]._v);
	std::vector<Vec2d> uvs;
	std::vector<unsigned long> tes;
	std::vector<double> params;
	std::map<unsigned long, int> teMap;
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
		if (_mt->triangleMaterial(i) < 0)
			continue;
		long* tr = _mt->triangleVertices(i);
		for (int j = 0; j < 3; ++j) {
			int last = tr[j], next = j < 2 ? tr[j + 1] : tr[0];
			if (next < last && last < _preDeepCutVerts) {  // only do edge once and don't do new deep cut edges
				Vec3d lV = _deepXyz[last], V = _deepXyz[next];
				boundingBox<double> be;
				be.Empty_Box();  // is this necessary?
				be.Enlarge_To_Include_Point(lV._v);
				be.Enlarge_To_Include_Point(V._v);
				if (bb.Intersection(be)) {
					double rayParams[2];
					Vec2d faceParams[2];
					if (bilinearRayIntersection(lV, V - lV, bl, rayParams, faceParams) == 1) {
						insidePolygon ip;
						if (ip.insidePolygon2d(faceParams[0], bUv)) {
							tes.push_back((i << 2) + j);
							teMap.insert(std::make_pair(tes.back(), uvs.size()));
							teMap.insert(std::make_pair(_mt->triAdjs(i)[j], uvs.size()));
							params.push_back(rayParams[0]);
							uvs.push_back(faceParams[0]);
						}
					}
				}
			}
		}
	}
	if (tes.empty())
		return;
	struct holePoint {
		unsigned long te;
		double param;
		Vec2d uv;
	}hp;
	std::vector<std::list<holePoint> > holesIn;
	auto growLoop = [&](int idx) {
		holesIn.push_back(std::list<holePoint>());
		hp.te = tes[idx];
		hp.param = params[idx];
		hp.uv = uvs[idx];
		holesIn.back().push_back(hp);
		auto lowestV = holesIn.back().begin();
		unsigned long lastTe = tes[idx], te;
		tes[idx] = 3;
		while (true) {
			// get 2 possible opposite te and check for path
			int j;
			for (j = 1; j < 3; ++j) {
				te = _mt->triAdjs(lastTe >> 2)[((lastTe & 3) + j) % 3];
				auto tit = teMap.find(te);
				if (tit != teMap.end()) {
					idx = tit->second;
					if (tes[idx] == 3)
						continue;
					lastTe = hp.te = te;
					if (tes[idx] != te)
						hp.param = 1.0 - params[idx];
					else
						hp.param = params[idx];
					hp.uv = uvs[idx];
					holesIn.back().push_back(hp);
					if (lowestV->uv.Y > hp.uv.Y) {
						lowestV = holesIn.back().end();
						--lowestV;
					}
					tes[idx] = 3;  // mark empty
					break;
				}
			}
			if (j > 2)
				break;
		}
	};
	for (int n = tes.size(), i = 0; i < n; ++i) {
		if (tes[i] != 3)
			growLoop(i);
	}
	// occasionally get spurious edge capture from collisions which are not a valid hole
	auto hiit = holesIn.begin();
	while (hiit != holesIn.end()) {
		if (_mt->triAdjs(hiit->front().te >> 2)[hiit->front().te & 3] >> 2 != hiit->back().te >> 2) {
			hiit = holesIn.erase(hiit);
			continue;
		}
		int mat;
		while ((mat = _mt->triangleMaterial(hiit->front().te >> 2)) == 3 || mat == 4) {  // put a valid start in the front
			hiit->push_back(hiit->front());
			hiit->pop_front();
		}
		++hiit;
	}
	std::vector<std::vector<Vec2d> > holeBorders;
	for (auto& hole : holesIn) {
		holeBorders.push_back(std::vector<Vec2d>());
		holeBorders.back().reserve(hole.size());
		for (auto& hp : hole)
			holeBorders.back().push_back(hp.uv);
	}
	// only want first level holes so prune any holes interior to others
	auto hbit = holeBorders.begin();
	for (int i = 0; i < holesIn.size(); ++i) {
		auto hbit2 = holeBorders.begin();
		while(hbit2 != holeBorders.end()){
			if (hbit == hbit2) {
				++hbit2;
				continue;
			}
			insidePolygon ip;
			if (ip.insidePolygon2d(hbit->front(), *hbit2)) {
				holesIn[i].clear();
				hbit->clear();
				break;
			}
			++hbit2;
		}
		++hbit;
	}
	for (auto& hole : holesIn) {
		if (hole.empty())
			continue;
		std::vector<unsigned long> topTe;
		std::vector<float> topParams;
		std::list<Vec2d> topUvs;
		long deepStart = -1, topStart = -1;
		surfaceCutLine scl;
		bool Tin = false, closedLoop = true;
		unsigned long backTeRev = _mt->triAdjs(hole.back().te >> 2)[hole.back().te & 3];
		int splitTri, splitEdge;
		float splitParam;
		int prevMat = _mt->triangleMaterial(hole.front().te >> 2), Ttriangle = -1;
		for (auto pt = hole.begin(); pt != hole.end(); ++pt) {
			splitTri = pt->te >> 2, splitEdge = pt->te & 3;
			splitParam = (float)pt->param;
			long* tr = _mt->triangleVertices(splitTri);
			if (_mt->triangleMaterial(splitTri) == 4) {  // mat 5 march to an undermine border edge. Pop back on top
				closedLoop = false;
				for (auto& db : _deepBed) {
					if (db.second.deepMtVertex == tr[((pt->te & 3) + 1) % 3]) {
						std::vector<materialTriangles::neighborNode> nei;
						_mt->getNeighbors(db.first, nei);
						auto nit = nei.begin();
						while (nit != nei.end()) {
							if (_deepBed[nit->vertex].deepMtVertex == tr[pt->te & 3]) {
								tr = _mt->triangleVertices(nit->triangle);
								int i;
								for (i = 0; i < 3; ++i) {
									if (tr[i] == nit->vertex) {
//										pt->te = (nit->triangle << 2) + i;  // deep mat 5 tris always listed in same order as the top tri that generated it.
										break;
									}
								}
								if (i > 2)
									int junk = 0;
								//									return false;
								break;
							}
							++nit;
						}
						assert(nit != nei.end());
						break;
					}
				}
			}
			else if (_mt->triangleMaterial(splitTri) == 3) {  // march over an incision edge until a mat 5 or 2 triangle is found
				closedLoop = false;
				if (prevMat == 2) {  // this is a T out
					auto adjs = _mt->triAdjs(splitTri + 1);  // incision convention
					if (_mt->triangleMaterial(adjs[0] >> 2) == 3) {  // non undermined incision edge

						// COURT - write me

						prevMat = 2;
					}
					else {
						assert(_mt->triangleMaterial(adjs[0] >> 2) == 5 && _mt->triangleMaterial(_mt->triAdjs(splitTri)[0]>>2) == 2);
						unsigned long backTe = _mt->triAdjs(topTe.back() >> 2)[topTe.back() & 3];
						int tv = TinSub(splitTri, splitParam);
						topUvs.push_back(pt->uv);
						topTe.back() = _mt->triAdjs(backTe >> 2)[backTe & 3];
						cutSkinLine(topStart, tv, topTe, topParams, Tin, true, scl);
						topTe.clear();
						topParams.clear();
						deepStart = -1;
//						ToutTri = -1;
						Tin = false;
						++pt;
						++pt;
						deepStart = -1;
						prevMat = 5;
					}
				}
				else {  // Tin
					assert(prevMat == 5);
					cutDeepSurface(deepStart, -1, topTe, topParams, scl);
					++pt;
					topStart = TinSub(pt->te >> 2, (float)pt->param);
					topUvs.push_back(pt->uv);
					topTe.clear();
					topParams.clear();
					Tin = true;
					++pt;
					deepStart = -1;
					prevMat = 2;
				}
			}
			else {
				// do previous Tin or Tout only after next te secured as the T op will split the top edge triangle
				if (pt == hole.begin()) {
					if (prevMat == 2) {
						float uv[2] = { 0.0f, 0.0f };
						if (splitEdge < 1)
							uv[0] = splitParam;
						else if (splitEdge > 1)
							uv[1] = 1.0f - splitParam;
						else {
							uv[1] = splitParam;
							uv[0] = 1.0f - splitParam;
						}
						long topVertex, bottomVertex;
						createFlapTopBottomVertices(splitTri, uv, topVertex, bottomVertex);
						topStart = topVertex;
						_loopSkinTopBegin = topStart;
					}
					else {
						Vec3f gridLocus;
						long tet = parametricMTedgeTet(splitTri, splitEdge, splitParam, gridLocus);
						int bedVertex = _mt->splitTriangleEdge(splitTri, splitEdge, splitParam);
						assert(bedVertex == _vbt->_vertexTets.size());
						_vbt->_vertexTets.push_back(tet);
						_vbt->_barycentricWeights.push_back(Vec3f());
						_vbt->gridLocusToBarycentricWeight(gridLocus, *_vbt->tetCentroid(tet), _vbt->_barycentricWeights.back());
						deepStart = bedVertex;
					}
					topUvs.push_back(pt->uv);
					backTeRev = _mt->triAdjs(backTeRev >> 2)[backTeRev & 3];  // necessary because above split may have changed the last te.
					hole.back().te = backTeRev;
				}
				else {
					topTe.push_back(pt->te);
					topParams.push_back(splitParam);
					topUvs.push_back(pt->uv);
				}
			}
		}
		if (prevMat == 2) {
			if (closedLoop) {
				cutSkinLine(topStart, -1, topTe, topParams, false, false, scl);
				std::list<long> topVerts, botVerts;
				topVerts.push_back(topStart);
				botVerts.push_back(scl.deepVerts.front());
				auto dbit = _deepBed.begin();
				topStart = scl.deepVerts.back();
				while (dbit != _deepBed.end()) {
					if (dbit->second.deepMtVertex == topStart)
						break;
					++dbit;
				}
				if (dbit == _deepBed.end())
					throw(std::logic_error("Program error in interior hole surface cutter."));
				topVerts.push_back(dbit->first);
				botVerts.push_back(scl.deepVerts.back());
				if (!topDeepSplit_Sub(topVerts, botVerts, true, true))
					throw(std::logic_error("Program error in interior hole surface cutter."));
			}
			else {  // this segment will always connect to its first segment
				cutSkinLine(topStart, _loopSkinTopBegin, topTe, topParams, Tin, true, scl);  // 
				scl.deepVerts.pop_back();
			}
			Tin = false;  // COURT reinit all variables for next hole
		}
		else {
			if (closedLoop) {
				cutDeepSurface(deepStart, -1, topTe, topParams, scl);
			}
		}
		//		std::list< std::vector<std::pair<long, Vec2d> > >& holes
		holes.push_back(std::vector<std::pair<long, Vec2d> >());
		if(topUvs.size() != scl.deepVerts.size())
			throw(std::logic_error("Program error in interior hole surface cutter."));
		auto uit = topUvs.begin();
		for (auto dv : scl.deepVerts) {
			holes.back().push_back(std::make_pair(dv, *uit));
			++uit;
		}
		_holePolyLines.push_back(std::list<long>());
		_holePolyLines.back().splice(_holePolyLines.back().begin(), scl.deepVerts);
		_holePolyLines.back().push_back(_holePolyLines.back().front());  // closed loop code
	}
	// holes must be listed clockwise
	auto hlit = holes.begin();
	auto hpit = _holePolyLines.begin();
	hbit = holeBorders.begin();
	while (hbit != holeBorders.end()) {
		if (clockwise(*hbit) < 2) {
			hpit->reverse();
			std::reverse(hlit->begin(), hlit->end());
		}
		++hbit;
		++hlit;
	}
}

double deepCut::surfacePath(rayTriangleIntersect& from, const rayTriangleIntersect& to, const bilinearPatch *holeBl, const bool cutPath, double& minimumBilinearV) {
	minimumBilinearV = DBL_MAX;
	if (cutPath) {
		if(!from.scl.deepVerts.empty())
			throw(std::logic_error("surfacePath() called on a path that is not empty;"));
	}
	bilinearPatch bl;
	bool outsidePatch = false, firstUzero;
	Vec3d P00, P10, P01, P11, lastI = from.intersect;
	if (holeBl != nullptr) {
		P00 = holeBl->P00;
		P10 = holeBl->P10;
		P01 = holeBl->P01;
		P11 = holeBl->P11;
		outsidePatch = false;
		firstUzero = true;
	}
	else {
		int intersectLevel = _deepPosts[from.postNum].triIntersects.size() - 1;  // keep bilinear surface as square as possible
		if (intersectLevel > _deepPosts[to.postNum].triIntersects.size() - 1)
			intersectLevel = _deepPosts[to.postNum].triIntersects.size() - 1;
		if (from.postNum < to.postNum) {
			assert((from.rayIndex & 1) < 1);
			P01 = _deepPosts[from.postNum].triIntersects.front().intersect;
			P11 = _deepPosts[to.postNum].triIntersects.front().intersect;
			P00 = _deepPosts[from.postNum].triIntersects[intersectLevel].intersect;
			P10 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
			firstUzero = true;
		}
		else if (from.postNum > to.postNum) {
			assert(from.rayIndex & 1);
			P10 = _deepPosts[from.postNum].triIntersects[intersectLevel].intersect;
			P00 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
			P11 = _deepPosts[from.postNum].triIntersects.front().intersect;
			P01 = _deepPosts[to.postNum].triIntersects.front().intersect;
			firstUzero = false;
		}
		else if (from.rayIndex > to.rayIndex) {
			if (from.postNum == _deepPosts.size() - 1) { // open end case
				if (intersectLevel > _deepPosts[from.postNum - 1].triIntersects.size() - 1)
					intersectLevel = _deepPosts[from.postNum - 1].triIntersects.size() - 1;
				P11 = _deepPosts[from.postNum].triIntersects.front().intersect;
				P10 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
				P01 = _deepPosts[from.postNum - 1].triIntersects.front().intersect;
				P00 = _deepPosts[to.postNum - 1].triIntersects[intersectLevel].intersect;
				outsidePatch = true;
				firstUzero = false;
			}
			else {
				if (from.rayIndex & 1) {  // other open end case
					assert(from.postNum < 1);
					outsidePatch = true;
				}
				if (intersectLevel > _deepPosts[from.postNum + 1].triIntersects.size() - 1)
					intersectLevel = _deepPosts[from.postNum + 1].triIntersects.size() - 1;
				P01 = _deepPosts[from.postNum].triIntersects.front().intersect;
				P00 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
				P11 = _deepPosts[from.postNum + 1].triIntersects.front().intersect;
				P10 = _deepPosts[to.postNum + 1].triIntersects[intersectLevel].intersect;
				firstUzero = true;
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
				outsidePatch = true;
				firstUzero = true;
			}
			else {
				if (to.rayIndex & 1) {  // other open end case
					assert(from.postNum == _deepPosts.size() - 1);
					outsidePatch = true;
				}
				if (intersectLevel > _deepPosts[from.postNum - 1].triIntersects.size() - 1)
					intersectLevel = _deepPosts[from.postNum - 1].triIntersects.size() - 1;
				P11 = _deepPosts[from.postNum].triIntersects.front().intersect;
				P10 = _deepPosts[to.postNum].triIntersects[intersectLevel].intersect;
				P01 = _deepPosts[from.postNum - 1].triIntersects.front().intersect;
				P00 = _deepPosts[to.postNum - 1].triIntersects[intersectLevel].intersect;
				firstUzero = false;
			}
		}
	}
	makeBilinearPatch(P00, P10, P01, P11, bl);
	double len = 0.0;
	Vec3d lastN, E, nE, I;
	// get first surface normal
	if (from.rayIndex < 1) {
		if (!cutPath && to.triangle == from.triangle) {
			from.scl.deepVerts.clear();
			minimumBilinearV = 1.0;
			return 0.0;
		}
		bilinearNormal(firstUzero ? 0.0 : 1.0, 1.0, bl, lastN);
	}
	else {
		double v;
		if(firstUzero)
			v = (from.intersect - bl.P00).length2() / (bl.e00 * bl.e00);
		else
			v = (from.intersect - bl.P10).length2() / (bl.e11 * bl.e11);
		bilinearNormal(firstUzero ? 0.0 : 1.0, sqrt(v), bl, lastN);
		if (!cutPath && to.triangle == from.triangle) {
			from.scl.deepVerts.clear();
			minimumBilinearV = sqrt(v);
			return 0.0;
		}
	}
	double rayParams[2];
	Vec2d faceParams[2];
	auto edgeIntersect = [&]() ->bool {
		if (bilinearRayIntersection(E, nE - E, bl, rayParams, faceParams) != 1)
			return false;
		bool ret;
		if (faceParams[0].X < 0.0 || faceParams[0].X >= 1.0)
			ret = outsidePatch ? true : false;
		else
			ret = outsidePatch ? false : true;
		if (!ret)
			return false;
		I = E * (1.0 - rayParams[0]) + nE * rayParams[0];
		len += (I - lastI).length();
		lastI = I;
		bilinearNormal(faceParams[0].X, faceParams[0].Y, bl, lastN);
		if (minimumBilinearV > faceParams[0].Y)
			minimumBilinearV = faceParams[0].Y;
		return true;
	};
	unsigned long te = 0xffffffff;
	long *tr;
	int nEdges = 1, prevMat;
	// next lambda is for difficult triangles whose edges don't cut cleanly with bilinearRayIntersection()
	auto triangleEdgeCuts = [&](int tri, int(&edges)[2], double(&eParams)[2], Vec2d (&fParams)[2]) ->int{
		tr = _mt->triangleVertices(tri);
		double p[3], last;
		Vec2d fp[3];
		for (int i = 0; i < 3; ++i) {
			E = _deepXyz[tr[i]] - lastN;
			double rP[2];
			Vec2d fP[2];
			int nSols = bilinearRayIntersection(E, lastN * 2.0, bl, rP, fP);
			if(nSols < 1)
				return -1;
			if (nSols > 1) {
				if (abs(rP[1] - 0.5) < abs(rP[0] - 0.5)) {
					rP[0] = rP[1];
					fP[0] = fP[1];
				}
			}
			p[i] = rP[0] - 0.5;
			fp[i] = fP[0];
		}
		last = p[0];
		int eNum = 0;
		for (int i = 2; i > -1; --i) {
			if (signbit(last) != signbit(p[i])) {
				eParams[eNum] = p[i] / (p[i] - last);
				fParams[eNum] = fp[i] * (1.0 - eParams[eNum]) + fp[(i+1)%3] * eParams[eNum];
				if (fParams[eNum][0] < 0.0 || fParams[eNum][0] >= 1.0)
					if(outsidePatch)
						edges[eNum++] = i;
				else
					if (!outsidePatch)
						edges[eNum++] = i;
			}
			last = p[i];
		}
		if (eNum < 1)
			return -1;
		return eNum;
	};
	auto getStartEndTE = [&](const rayTriangleIntersect &rti) ->unsigned long{
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
							return _mt->triAdjs(n.triangle)[(i+2)%3];
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
	if(from.deepVert < 0) {
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
		if (i < 0) {
			int edgesCut[2];
			double edgeParams[2];
			Vec2d bilinParams[2];
			if (triangleEdgeCuts(from.triangle, edgesCut, edgeParams, bilinParams) != 1)
				return DBL_MAX;
//				throw(std::logic_error("This surfacePath() call does not have a valid starting triangle."));
			te = _mt->triAdjs(from.triangle)[edgesCut[0]];
			rayParams[0] = edgeParams[0];
			bilinearNormal(bilinParams[0].X, bilinParams[0].Y, bl, lastN);
		}
	}
	else {
		if ((te = getStartEndTE(to)) == 3)  // COURT use triangleEdgeCuts()
			throw(std::logic_error("Program error start/ending a deepCut."));
		endTriangle = _mt->triAdjs(te>>2)[te&3] >> 2;
		if ((te = getStartEndTE(from)) == 3)
			throw(std::logic_error("Program error start/ending a deepCut."));
	}
	std::vector<unsigned long> topTe;
	std::vector<float> topParams;
	if (cutPath) {
		topTe.push_back(te);
		topParams.push_back(1.0f - (float)rayParams[0]);
	}
	long deepStart = from.deepVert, topStart = from.mat2Vert;
	long TinTri = -1, ToutTri = -1;
	unsigned long UinTe = 3, UoutTe = 3;
	float TinParam, ToutParam, UinParam, UoutParam;
	bool Tin = false;
	do{
		int splitTri, splitEdge, splitMat;
		splitTri = te >> 2, splitEdge = te & 3;
		splitMat = _mt->triangleMaterial(splitTri);
		if (splitTri == endTriangle)
			break;
		tr = _mt->triangleVertices(splitTri);
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
										cutDeepSurface(deepStart, -1, topTe, topParams, from.scl);
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
						deepStart = -1;
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
					topTe.pop_back();
					topParams.pop_back();
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
					long v50 = _deepBed[tr[splitEdge]].deepMtVertex;
					long v51 = _deepBed[tr[(splitEdge + 1) % 3]].deepMtVertex;
					assert(v50 > -1 && v51 > -1);
					std::vector<materialTriangles::neighborNode> nei;
					_mt->getNeighbors(v51, nei);
					auto nit = nei.begin();
					while (nit != nei.end()) {
						if (nit->vertex == v50) {
							assert(_mt->triangleMaterial(nit->triangle) == 5);
							long* trp = _mt->triangleVertices(nit->triangle);
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
							cutDeepSurface(deepStart, -1, topTe, topParams, from.scl);
							topStart = TinSub(TinTri, TinParam);
							topTe.clear();
							topParams.clear();
							Tin = true;
							TinTri = -1;
						}
						if (ToutTri > -1) {
							topTe.pop_back();
							topParams.pop_back();
							unsigned long ate;
							if(!topTe.empty())
								ate = _mt->triAdjs(topTe.back() >> 2)[topTe.back() & 3];
							int tv = TinSub(ToutTri, ToutParam);
							if (!topTe.empty())
								topTe.back() = _mt->triAdjs(ate >> 2)[ate & 3];
							cutSkinLine(topStart, tv, topTe, topParams, Tin, true, from.scl);
							topTe.clear();
							topParams.clear();
							deepStart = -1;
							ToutTri = -1;
							Tin = false;
						}
						if (UinTe != 3) {
							unsigned long ate = _mt->triAdjs(UinTe >> 2)[UinTe & 3];
							assert(_deepBed[_mt->triangleVertices(ate >> 2)[((ate & 3) + 2) % 3]].deepMtVertex < 0);
							std::list<long> topVerts, deepVerts;
							topVerts.push_back( _mt->triangleVertices(UinTe >> 2)[((UinTe&3)+2)%3]);
							deepVerts.push_back(_deepBed[topVerts.back()].deepMtVertex);
							float uv[2] = {0.333f, 0.333f};
							long bottomVertex;
							createFlapTopBottomVertices(UinTe >> 2, uv, topStart, bottomVertex);
							topVerts.push_back(topStart);
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
							createFlapTopBottomVertices(ate >> 2, uv, topStart, bottomVertex);
							topVerts.push_back(topStart);
							deepVerts.push_back(bottomVertex);
							topDeepSplit_Sub(topVerts, deepVerts, false, false);
							topStart = topVerts.back();
							topTe.clear();
							topParams.clear();
							Tin = true;
							UinTe = 3;
						}
						if (UoutTe != 3) {
							unsigned long ate = _mt->triAdjs(UoutTe >> 2)[UoutTe & 3];
							std::list<long> topVerts, deepVerts;
							topVerts.push_back(_mt->triangleVertices(UoutTe >> 2)[((UoutTe & 3) + 2) % 3]);
							deepVerts.push_back(_deepBed[topVerts.back()].deepMtVertex);
							float uv[2] = { 0.333f, 0.333f };
							long topVertex, bottomVertex;
							createFlapTopBottomVertices(UoutTe >> 2, uv, topVertex, bottomVertex);
							topVerts.push_back(topVertex);
							deepVerts.push_back(bottomVertex);
							UoutTe = _mt->triAdjs(ate >> 2)[ate & 3];
							if ((ate & 3) < 1) {
								uv[0] = 1.0f - UoutParam;
								uv[1] = 0.0f;
							}
							else if ((ate & 3) > 1) {
								uv[1] = UoutParam;
								uv[0] = 0.0f;
							}
							else {
								uv[0] = 1.0f - UoutParam;
								uv[1] = UoutParam;
							}
							topTe.pop_back();
							topParams.pop_back();
							unsigned long teEnd = _mt->triAdjs(topTe.back() >> 2)[topTe.back() & 3];
							createFlapTopBottomVertices(ate >> 2, uv, topVertex, bottomVertex);
							topVerts.push_back(topVertex);
							deepVerts.push_back(bottomVertex);
							topDeepSplit_Sub(topVerts, deepVerts, false, false);
							topTe.back() = _mt->triAdjs(teEnd >> 2)[teEnd & 3];
							cutSkinLine(topStart, topVertex, topTe, topParams, Tin, true, from.scl);
							topStart = -1;
							deepStart = bottomVertex;
							topTe.clear();
							topParams.clear();
							Tin = false;
							UoutTe = 3;
						}
						if (splitMat == 2 && prevMat > 4) {  // skin start from a boundary or periosteal edge
							float lastParam = topParams.back();
							topTe.pop_back();
							topParams.pop_back();
							cutDeepSurface(deepStart, -1, topTe, topParams, from.scl);
							float uv[2] = { 0.0f, 0.0f };
							long topVertex, bottomVertex;
							if (splitEdge < 1)
								uv[0] = lastParam;
							else if (splitEdge > 1)
								uv[1] = 1.0f - lastParam;
							else {
								uv[0] = 1.0f - lastParam;
								uv[1] = lastParam;
							}
							createFlapTopBottomVertices(splitTri, uv, topVertex, bottomVertex);
							from.scl.deepVerts.push_back(topVertex);
							topStart = topVertex;
							topTe.clear();
							topParams.clear();
							Tin = false;
							TinTri = -1;
							prevMat = 2;
						}
						if (splitMat > 4 && prevMat == 2) {  // skin end from a boundary or periosteal edge
							float lastParam = topParams.back();
							topTe.pop_back();
							topParams.pop_back();
							float uv[2] = { 0.0f, 0.0f };
							long topVertex, bottomVertex;
							unsigned long lastTe = _mt->triAdjs(topTe.back()>>2)[topTe.back()&3], adj = _mt->triAdjs(splitTri)[splitEdge];
							if ((adj&3) < 1)
								uv[0] = 1.0f - lastParam;
							else if ((adj & 3) > 1)
								uv[1] = lastParam;
							else {
								uv[1] = 1.0f - lastParam;
								uv[0] = lastParam;
							}
							createFlapTopBottomVertices(adj>>2, uv, topVertex, bottomVertex);
							topTe.back() = _mt->triAdjs(lastTe >> 2)[lastTe & 3];
							cutSkinLine(topStart, topVertex, topTe, topParams, Tin, false, from.scl);
							from.scl.deepVerts.push_back(topVertex);
							topTe.clear();
							topParams.clear();
							deepStart = -1;
							ToutTri = -1;
							TinTri = -1;
							Tin = false;
						}
						topTe.push_back(te);
						topParams.push_back(1.0f - (float)rayParams[0]);
					}
					if (splitMat == 2)
						prevMat = 2;
					else
						prevMat = 5;  // includes 1 & 5-9
					break;
				}
				E = nE;
			}
			if (i > 3) {

				return DBL_MAX;  // COURT - reincorporate later

				int edgesCut[2];
				double edgeParams[2];
				Vec2d bilinParams[2];
				int nCuts = triangleEdgeCuts(splitTri, edgesCut, edgeParams, bilinParams);
				if (nCuts < 1)
					return DBL_MAX;
				if (nCuts < 2)
					return DBL_MAX;  // surface path returns to negative half space
				for (i = 0; i < 2; ++i) {
					if (edgesCut[i] == splitEdge)
						continue;
					else
						break;
				}
				if (i > 1)
					throw(std::logic_error("Error in deep cut surface path finder"));
				te = _mt->triAdjs(splitTri)[edgesCut[i]];  // these two lines are for the next cut
				++nEdges;
				if (cutPath) {
					topTe.push_back(te);
					topParams.push_back(1.0f - (float)rayParams[0]);
				}
				if (_mt->triangleMaterial(splitTri) == 2)
					prevMat = 2;
				else
					prevMat = 5;  // includes 1 & 5-9
				bilinearNormal(bilinParams[i].X, bilinParams[i].Y, bl, lastN);
			}
		}
	}while ((te >> 2) != endTriangle && nEdges < 500);
	if (cutPath) {
		if (prevMat == 2)
			cutSkinLine(topStart, to.mat2Vert, topTe, topParams, Tin, false, from.scl);
		else
			cutDeepSurface(deepStart, to.deepVert, topTe, topParams, from.scl);
	}
	if (nEdges > 499)
		return DBL_MAX;
	return len;
}

void deepCut::cutDeepSurface(int startV, int endV, std::vector<unsigned long>& te, std::vector<float>& params, surfaceCutLine& scl) {
	if ((scl.deepVerts.empty() || startV != scl.deepVerts.back()) && startV > -1)
		scl.deepVerts.push_back(startV);
	auto pit = params.begin();
	for (auto& t : te) {
		Vec3f gridLocus;
		long tri = t >> 2, edge = t & 3;
		long tet = parametricMTedgeTet(tri, edge, *pit, gridLocus);
		int bedVertex = _mt->splitTriangleEdge(tri, edge, *pit);
		assert(bedVertex == _vbt->_vertexTets.size());
		_vbt->_vertexTets.push_back(tet);
		_vbt->_barycentricWeights.push_back(Vec3f());
		_vbt->gridLocusToBarycentricWeight(gridLocus, *_vbt->tetCentroid(tet), _vbt->_barycentricWeights.back());
		scl.deepVerts.push_back(bedVertex);
		++pit;
	}
	if(endV > -1)
		scl.deepVerts.push_back(endV);
}

void deepCut::cutSkinLine(int startV, int endV, std::vector<unsigned long>& te, std::vector<float>& params, bool Tin, bool Tout, surfaceCutLine& scl) {
	if (startV == _previousSkinTopEnd)
		Tin = true;
	if(endV == _loopSkinTopBegin)
		Tout = true;
	auto pit = params.begin();
	std::list<long> topVerts, botVerts;
	topVerts.push_back(startV);
	auto dbit = _deepBed.find(startV);
	if (dbit == _deepBed.end() || dbit->second.deepMtVertex < 0)
		throw(std::logic_error("Program error in skinCutLine()."));
	botVerts.push_back(dbit->second.deepMtVertex);
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
		long topVertex, bottomVertex;
		createFlapTopBottomVertices(tri, uv, topVertex, bottomVertex);
		topVerts.push_back(topVertex);
		botVerts.push_back(bottomVertex);
		++pit;
	}
	if (endV > -1) {
		topVerts.push_back(endV);
		dbit = _deepBed.find(endV);
		if (dbit == _deepBed.end() || dbit->second.deepMtVertex < 0)
			throw(std::logic_error("Program error in skinCutLine()."));
		botVerts.push_back(dbit->second.deepMtVertex);
	}
	_mt->findAdjacentTriangles(true, false);
	updateDeepSpatialCoordinates();
	if(!topDeepSplit_Sub(topVerts, botVerts, Tin, Tout))
		throw(std::logic_error("Program error in skinCutLine()."));
	auto bvit = botVerts.begin();
	for ( ; bvit !=botVerts.end(); ++bvit)
		scl.deepVerts.push_back(*bvit);
}

