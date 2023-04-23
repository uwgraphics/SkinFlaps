////////////////////////////////////////////////////////////////////////////
// File: skinCutUndermineTets.cpp
// Author: Court Cutting
// Date: 6/2/2019
// Purpose: Interactive skin-mucosal incision class which operates in spatial coordinates. Requires setup with a triangulated subdermal surface corresponding to skin-mucosal surface of
//    closed manifold surface surrounding elastic solid modelled as virtual noded cubes.  Both superficial and deep surfaces in material coordinates on setup.
//    Allows user to make depth and normal agnostic incisions on skin-mucosa in spatial coordinates.  Also allows undermining of flaps aint an incision line
//    away from the deep bed.  Updates both the virtual noded cubes and the surface model.
////////////////////////////////////////////////////////////////////////////

#include <tuple>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <deque>
#include <fstream>
#include <unordered_set>
#include "Mat3x3f.h"
#include "Mat2x2f.h"
#include "boundingBox.h"
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"
#include "insidePolygon.h"
#include "skinCutUndermineTets.h"
#include "surgicalActions.h"
#include "gl3wGraphics.h"

// could also have done a singleton class, but more cumbersome
std::unordered_map<int, skinCutUndermineTets::deepPoint> skinCutUndermineTets::_deepBed;
gl3wGraphics* skinCutUndermineTets::_gl3w;
materialTriangles* skinCutUndermineTets::_mt;
vnBccTetrahedra* skinCutUndermineTets::_vbt;

bool skinCutUndermineTets::skinCut(std::vector<Vec3f> &topCutPoints, std::vector<Vec3f> &topNormals, bool startOpen, bool endOpen)
{  //  Only cuts material 2 triangles down to deep bed and creates single vertex deep cut line along with material 3 side triangles.
	// Does not cut cubes.  This doesn't happen until undermining is done.  Input is an array of cutter vertices.
	// Inputs start/endOpen indicate whether or not the first or last point "open" which means to freely cut the material 3 surface
	// to the nearest open edge.  cutBottom==true indicates that the straight line cuts between points is done on the bottom side.
	// These straight bottom cuts are necessary when deep cutting from the skin surface so that large deep planes are created.
	_startOpen = startOpen;
	_endOpen = endOpen;
	_solidRecutRequired = false;
	_firstTopVertex = -1;
	if (topCutPoints.size() < 2)
		return false;
	std::vector<int> topMtVertices, deepVertexLine;
	topMtVertices.assign(topCutPoints.size(), -1);
	deepVertexLine.assign(topCutPoints.size(), -1);
	float uv[2];
	for (int n=topCutPoints.size(),i = 0; i < n; ++i){
		int tri = _mt->numberOfTriangles();
		_mt->closestPoint(topCutPoints[i].xyz, tri, uv);
		Vec3f v3;
		_mt->getBarycentricPosition(tri, uv, v3.xyz);
		if ((i<1 && _startOpen) || (i>n - 2 && _endOpen)){  // COURT - fix me. Must allow re-entrant endOpen not yet programmed.
			topCutPoints[i].set(v3);
			if (i<1)
				topMtVertices[i] = addTinEdgeVertex(topCutPoints[0], topCutPoints[1]);
			else
				topMtVertices[i] = addTinEdgeVertex(topCutPoints[n-1], topCutPoints[n-2]);
			auto dbit = _deepBed.find(topMtVertices[i]);
			assert(dbit != _deepBed.end());
			deepVertexLine[i] = dbit->second.deepMtVertex;
			continue;
		}
		assert((topCutPoints[i] - v3).length2() < 0.001f);
		if (uv[0] < -0.0001 || uv[1]<-0.0001 || uv[0] + uv[1]>1.001f)
			assert(false);
		else
			createFlapTopBottomVertices(tri, uv, topMtVertices[i], deepVertexLine[i]);
		if (i < 1)
			_firstTopVertex = topMtVertices[i];
	}
	_mt->findAdjacentTriangles(true);
	if (!topDeepSplit(topMtVertices, deepVertexLine, startOpen, endOpen))
		return false;
	return true;
}

int skinCutUndermineTets::parametricMTedgeTet(const int triangle, const int edge, const float param, Vec3f &gridLocus)
{

	assert(false);  // NUKE THIS ROUTINE  use new on in vnt for multires

	int *tr = _mt->triangleVertices(triangle);
	Vec3f tV[2];
	_vbt->vertexGridLocus(tr[edge], tV[0]);
	_vbt->vertexGridLocus(tr[(edge + 1) % 3], tV[1]);
	gridLocus = tV[0] * (1.0f - param) + tV[1] * param;
	bccTetCentroid tC;
	_vbt->gridLocusToLowestTetCentroid(gridLocus, tC);
	if (tC == _vbt->tetCentroid(_vbt->getVertexTetrahedron(tr[edge])))
		return _vbt->getVertexTetrahedron(tr[edge]);
	if (tC == _vbt->tetCentroid(_vbt->getVertexTetrahedron(tr[(edge + 1) % 3])))
		return _vbt->getVertexTetrahedron(tr[(edge + 1) % 3]);
	// find candidate cubes
	std::list<int> cc, tp;
	auto pr = _vbt->_tetHash.equal_range(tC);
	while (pr.first != pr.second){
		cc.push_back(pr.first->second);
		++pr.first;
	}
	if (cc.size() < 1){
		assert(false);
		return -1;
	}
	if (cc.size() < 2)
		return cc.front();
	for (auto c : cc){
		if (_vbt->decreasingCentroidPath(c, _vbt->_vertexTets[tr[edge]], tp))
			return c;
		if (_vbt->decreasingCentroidPath(c, _vbt->_vertexTets[tr[(edge + 1) % 3]], tp))
			return c;
	}
	assert(false);
	return -1;
}

int skinCutUndermineTets::parametricMTtriangleTet(const int mtTriangle, const float(&uv)[2], Vec3f &gridLocus)
{  // in material coords

	assert(false);  // NUKE THIS ROUTINE  use new on in vnt for multires

	int *tr = _mt->triangleVertices(mtTriangle);
	Vec3f tV[3];
	for (int i = 0; i < 3; ++i)
		_vbt->vertexGridLocus(tr[i], tV[i]);
	gridLocus = tV[0] * (1.0f - uv[0] - uv[1]) + tV[1] * uv[0] + tV[2] * uv[1];
	bccTetCentroid tC;
	_vbt->gridLocusToLowestTetCentroid(gridLocus, tC);
	for (int i = 0; i < 3; ++i){
		if (tC == _vbt->tetCentroid(_vbt->_vertexTets[tr[i]]))
			return _vbt->_vertexTets[tr[i]];
	}
	// find candidate cubes
	std::list<int> cc, tp;
	auto pr = _vbt->_tetHash.equal_range(tC);
	while (pr.first != pr.second){
		cc.push_back(pr.first->second);
		++pr.first;
	}
	if (cc.size() < 1){
		assert(false);
		return -1;
	}
	if (cc.size() < 2)
		return cc.front();
	for (auto c : cc){
		for (int i = 0; i < 3; ++i){
			if (_vbt->decreasingCentroidPath(c, _vbt->_vertexTets[tr[i]], tp))
				return c;
		}
	}
	assert(false);
	return -1;
}

int skinCutUndermineTets::createDeepBedVertex(std::unordered_map<int, deepPoint>::iterator &dit)
{
	if (dit->second.deepMtVertex > -1)
		return dit->second.deepMtVertex;
	Vec3f bw;
	int tet = deepPointTetWeight(dit, bw);
	if (tet < 0)
		return -1;
	dit->second.deepMtVertex = _mt->numberOfVertices();
	_vbt->_vertexTets.push_back(tet);
	_vbt->_barycentricWeights.push_back(bw);
	_mt->addVertices(1);
	Vec3f pos;
	_vbt->vertexBarycentricPosition(dit->second.deepMtVertex, pos);
	_mt->setVertexCoordinate(dit->second.deepMtVertex, pos.xyz);
	// Don't add texture here anymore as not unique so make them when using vert to make triangle
	return dit->second.deepMtVertex;
}

int skinCutUndermineTets::addSurfaceVertex(const int tet, const Vec3f &gridLocus)
{  // adds both here and to materialTriangles. Return -1 means gridLocus not in tet.
	Vec3f bw;
	_vbt->gridLocusToBarycentricWeight(gridLocus, _vbt->tetCentroid(tet), bw);
	assert(bw[0] >= 0.0f && bw[1] >= 0.0f && bw[2] >= 0.0f && bw[0] <= 1.0f && bw[1] <= 1.0f && bw[2] <= 1.0f && bw[0] + bw[1] + bw[2] <= 1.0f);
	int newVert = _vbt->vertexNumber();
	_vbt->_vertexTets.push_back(tet);
	_vbt->_barycentricWeights.push_back(bw);
//	assert(_mt->numberOfVertices() == newVert);  No longer necessarily true with skin=border deep cut
	_mt->addVertices(1);
	Vec3f pos;
	_vbt->vertexBarycentricPosition(newVert, pos);
	_mt->setVertexCoordinate(newVert, pos.xyz);
	// texture too often unknown at time of call. Leave for later.
	return newVert;
}

void skinCutUndermineTets::createFlapTopBottomVertices(const int topTriangle, float (&uv)[2], int &topVertex, int &bottomVertex)
{  // If uv == 0, 0 || uv == 1, 0 || uv == 0, 1 existing topVertex returned, otherwise new one created. Creates/gets a corresponding bottomVertex.
	// If topTriangle already part of a flap a bottom vertex will be added. If no flap present an unconnected bottom edge vertex will be created.
	// If a flap bottom does not already exists on input, bottomVertex will be negated on output.
	assert(_mt->triangleMaterial(topTriangle) == 2);
	// look for a flap bottom replicant of topTriangle if it exists.
	int deepVerts[3], *tr = _mt->triangleVertices(topTriangle);
	int bottomTriangle = 0, n = _mt->numberOfTriangles(), j = -1;
	Vec3f deepLocus, gridLocus;
	float mult;
	deepLocus.set(0.0f, 0.0f, 0.0f);
	for (int i = 0; i < 3; ++i){
		auto dbit = _deepBed.find(tr[i]);
		assert(dbit != _deepBed.end());
		if (i < 1)
			mult = 1.0f - uv[0] - uv[1];
		else if (i < 2)
			mult = uv[0];
		else
			mult = uv[1];
		deepLocus += dbit->second.gridLocus*mult;
		if (dbit->second.deepMtVertex < 0){  // signal topTriangle without an undermine. Still possible that a bottom undermined doesn't exist.
			j = 3;
			bottomTriangle = n;
		}
		deepVerts[i] = dbit->second.deepMtVertex;
	}
	if (bottomTriangle < n){
		std::vector<materialTriangles::neighborNode> nei;
		_mt->getNeighbors(deepVerts[0], nei);
		auto nit = nei.begin();
		while ( nit != nei.end()){
			if (nit->vertex == deepVerts[2]){
				++nit;
				if (nit == nei.end()){
					if (nei.front().vertex == deepVerts[1]) {
						nit = nei.begin();
						break;
					}
					else
						throw(std::logic_error("Topological eror found in createFlapTopBottomVertices()."));
				}
				else if (nit->vertex == deepVerts[1])
					break;
				else
					;
			}
			++nit;
		}
		if (nit == nei.end())  // no existing bottomTriangle
			bottomTriangle = n;
		else
			bottomTriangle = nit->triangle;
	}
	// check for an already existing vertex
	topVertex = -1;
	int edge = -1;
	float param;
	// don't allow an existing vertex as it invalidates sequential triangle edge movements
	if (uv[1] < 0.0001f) {
		uv[1] = 0.0001f;
		if (uv[0] < 0.0001f)
			uv[0] = 0.0001f;
		else if (uv[0] > 0.9999f)
			uv[0] = 0.9999f;
		else
			;
		edge = 0;
		param = uv[0];
	}
	else if (uv[0] < 0.0001f) {
		uv[0] = 0.0001f;
		if (uv[1] > 0.9999f)
			uv[1] = 0.9999f;
		edge = 2;
		param = 1.0f - uv[1];
	}
	else if (uv[0] + uv[1] > 0.9999f) {
		edge = 1;
		param = uv[1];
	}
	else
		;
	if (topVertex > -1){
		auto dbit = _deepBed.find(topVertex);
		assert(dbit != _deepBed.end());
		if (dbit->second.deepMtVertex > -1)
			bottomVertex = dbit->second.deepMtVertex;
		else
			bottomVertex = createDeepBedVertex(dbit);
		return;
	}
	// now either a splitEdge() or an addNewPointInMidTriangle() creation for top and bottom
	deepPoint dp;
	if (edge > -1){  // splitEdge() op
		tr = _mt->triangleVertices(topTriangle);
		auto dbit0 = _deepBed.find(tr[edge]);
		auto dbit1 = _deepBed.find(tr[(edge + 1) % 3]);
		assert(dbit0 != _deepBed.end() && dbit1 != _deepBed.end());
		Vec3f gridLocus;
		int* tr = _mt->triangleVertices(topTriangle);
		int tet = _vbt->parametricEdgeTet(tr[edge], tr[(edge + 1) % 3], param, gridLocus);
		topVertex = _mt->splitTriangleEdge(topTriangle, edge, param);
//		assert(topVertex == _vbt->_vertexTets.size());  No longer necessarily true with skin-border deepCut
		_vbt->_vertexTets.push_back(tet);
		_vbt->_barycentricWeights.push_back(Vec3f());
		_vbt->gridLocusToBarycentricWeight(gridLocus, _vbt->tetCentroid(tet), _vbt->_barycentricWeights.back());

		dp.gridLocus = dbit0->second.gridLocus*(1.0f - param) + dbit1->second.gridLocus*param;
		tet = flapBottomTet(tet, dp.gridLocus);
		assert(tet > -1);
		if (bottomTriangle < n){
			tr = _mt->triangleVertices(bottomTriangle);
			for (j = 0; j < 3; ++j){
				if (tr[j] == dbit1->second.deepMtVertex){
					assert(tr[(j + 1) % 3] == dbit0->second.deepMtVertex);
					dp.deepMtVertex = _mt->splitTriangleEdge(bottomTriangle, j, 1.0f - param);
					bottomVertex = dp.deepMtVertex;
					_vbt->_vertexTets.push_back(tet);
					_vbt->_barycentricWeights.push_back(Vec3f());
					_vbt->gridLocusToBarycentricWeight(dp.gridLocus, _vbt->tetCentroid(tet), _vbt->_barycentricWeights.back());
					break;
				}
			}
			assert(j < 3);
		}
		else if (dbit0->second.deepMtVertex > -1 && dbit1->second.deepMtVertex > -1) {
			std::vector<materialTriangles::neighborNode> nei;
			_mt->getNeighbors(dbit0->second.deepMtVertex, nei);
			auto nit = nei.begin();
			while (nit != nei.end()) {
				if (nit->vertex == dbit1->second.deepMtVertex)
					break;
				++nit;
			}
			assert(nit != nei.end());
			tr = _mt->triangleVertices(nit->triangle);
			for (j = 0; j < 3; ++j) {
				if (tr[j] == dbit1->second.deepMtVertex) {
					assert(tr[(j + 1) % 3] == dbit0->second.deepMtVertex);
					dp.deepMtVertex = _mt->splitTriangleEdge(nit->triangle, j, 1.0f - param);
					bottomVertex = dp.deepMtVertex;
					_vbt->_vertexTets.push_back(tet);
					_vbt->_barycentricWeights.push_back(Vec3f());
					_vbt->gridLocusToBarycentricWeight(dp.gridLocus, _vbt->tetCentroid(tet), _vbt->_barycentricWeights.back());
					break;
				}
			}
			assert(j < 3);
		}
		else  // no existing flap bottom
			dp.deepMtVertex = bottomVertex = addSurfaceVertex(tet, dp.gridLocus);
		_deepBed.insert(std::make_pair(topVertex, dp));
		return;
	}
	// New vertex to be created in mid triangle
	int tet = _vbt->parametricTriangleTet(_mt->triangleVertices(topTriangle), uv, dp.gridLocus);
	topVertex = _mt->addNewVertexInMidTriangle(topTriangle, uv);
	int topTexture = _mt->numberOfTextures() - 1;  // texture of this new vertex
	assert(topVertex == _vbt->_vertexTets.size());
	_vbt->_vertexTets.push_back(tet);
	Vec3f bw;
	_vbt->gridLocusToBarycentricWeight(dp.gridLocus, _vbt->_tetCentroids[tet], bw);
	_vbt->_barycentricWeights.push_back(bw);
	float uvDeep[2] = {0.0f, 1.0f};;
	if (bottomTriangle == n){
		dp.gridLocus = deepLocus;
		dp.deepMtVertex = -1;
		auto pr = _deepBed.emplace(topVertex, dp);
		createDeepBedVertex(pr.first);
		// vertices no longer have a unique texture
		bottomVertex = pr.first->second.deepMtVertex;
		return;
	}
	int *trBot = _mt->triangleVertices(bottomTriangle);
	for (j = 0; j < 3; ++j){
		if (trBot[j] == deepVerts[0])
			break;
	}
	assert(j < 3);
	assert(trBot[(j+1)%3] == deepVerts[2]);
	if (j < 1){
		uvDeep[0] = uv[1];
		uvDeep[1] = uv[0];
	}
	else if (j < 2){
		uvDeep[0] = 1.0f - uv[0] - uv[1];
		uvDeep[1] = uv[1];
	}
	else if (j < 3){
		uvDeep[0] = uv[0];
		uvDeep[1] = 1.0f - uv[0] - uv[1];
	}
	else
		assert(false);
	tet = parametricMTtriangleTet(bottomTriangle, uvDeep, dp.gridLocus);
	bottomVertex = _mt->addNewVertexInMidTriangle(bottomTriangle, uvDeep);
	dp.deepMtVertex = bottomVertex;
	assert(_vbt->_vertexTets.size() == bottomVertex);
	_vbt->_vertexTets.push_back(tet);
	_vbt->gridLocusToBarycentricWeight(dp.gridLocus, _vbt->_tetCentroids[tet], bw);
	_vbt->_barycentricWeights.push_back(bw);
	// not necessary to set position or texture.  Done in addNewVertexInMidTriangle()
	_deepBed.emplace(topVertex, dp);
}

bool skinCutUndermineTets::topDeepSplit(std::vector<int> &topV, std::vector<int> &deepV, bool frontSplit, bool backSplit)
{
	std::list<int> newTopVerts, newDeepVerts, tV, dV;
	for (int n = (int)topV.size(), i = 1; i < n; ++i){
		tV.push_back(topV[i - 1]);
		dV.push_back(deepV[i - 1]);
		if (!planeCutSurfaceLine(topV[i - 1], topV[i], deepV[i - 1], deepV[i], newTopVerts, newDeepVerts)){
			tV.splice(tV.end(), newTopVerts);
			dV.splice(dV.end(), newDeepVerts);
			return false;
		}
		tV.splice(tV.end(), newTopVerts);
		dV.splice(dV.end(), newDeepVerts);
	}
	tV.push_back(topV.back());
	dV.push_back(deepV.back());
	if (backSplit && tV.front() == tV.back()) {  // is a closed loop
		tV.pop_back();
		dV.pop_back();
		if(!topDeepSplit_Sub(tV, dV, false, false))
			return false;
		std::list<int> deepV, topV;
		topV.push_back(tV.back());
		deepV.push_back(dV.back());
		topV.push_back(tV.front());
		deepV.push_back(dV.front());
		return topDeepSplit_Sub(topV, deepV, true, true);
	}
	return topDeepSplit_Sub(tV, dV, frontSplit, backSplit);
}

bool skinCutUndermineTets::topDeepSplit_Sub(std::list<int> &topVerts, std::list<int> &deepVerts, bool frontSplit, bool backSplit){
	// broken into subroutine for later deep cut use
	// split top
	std::vector<int> oppVerts, oppTxs, oppBotVerts;
	int startVertex = -1, endVertex = -1;
	if (frontSplit)
		startVertex = deepVerts.front();
	if (backSplit)
		endVertex = deepVerts.back();
	flapSurfaceSplitter(startVertex, endVertex, topVerts, oppVerts);
	// unlike top surface, irregular undermining can cause alternating or multiple deep flap surface cuts
	startVertex = -1;
	endVertex = -1;
	oppBotVerts.assign(deepVerts.size(), -1);
	if (frontSplit)
		startVertex = topVerts.front();
	if (backSplit)
		endVertex = topVerts.back();
	auto dvit = deepVerts.begin();
	int i = 0;
	while (dvit != deepVerts.end()) {
		if (frontSplit && dvit == deepVerts.begin()) {
			assert(*_mt->vertexFaceTriangle(*dvit) != 0x80000000);
			std::vector<materialTriangles::neighborNode> nei;
			_mt->getNeighbors(*dvit, nei);
			auto nit = nei.begin();
			assert(nit->triangle > -1);
			while (nit != nei.end()) {
				if (_mt->triangleMaterial(nit->triangle) == 4)
					break;
				++nit;
			}
			if (nit == nei.end()) {  // bottom vert of a not undermined area where a split will occur
				// annoying case where a different texture must be put on either side of the split in a Tin
				// // COURT now done in Tin routine
				oppBotVerts[i] = *dvit;;
				++dvit; ++i;
				assert(dvit != deepVerts.end());
			}
			else {  // another possibility for a non-split bottom front end
				++dvit;
				if (*_mt->vertexFaceTriangle(*dvit) == 0x80000000) {  // non-split bottom front
					--dvit;
					oppBotVerts[i] = *dvit;;
					++dvit; ++i;
				}
				else
					--dvit;
			}
		}
		while (*_mt->vertexFaceTriangle(*dvit) == 0x80000000) {
			oppBotVerts[i] = *dvit;;
			++dvit; ++i;
			if (dvit == deepVerts.end())
				break;
		}
		if (!backSplit && dvit != deepVerts.end() && *dvit == deepVerts.back()) { // possible unconnected end condition
			--dvit;
			if (*_mt->vertexFaceTriangle(*dvit) == 0x80000000) {
				++dvit;
				oppBotVerts[i] = *dvit;;
				++i;
			}
			++dvit;
		}
		if (backSplit && *dvit == deepVerts.back()) {
			// code to handle undermine crossover
			--dvit;
			if (*_mt->vertexFaceTriangle(*dvit) == 0x80000000) {  // no deep split to last point possible
				++dvit;
				oppBotVerts[i] = *dvit;;
				++dvit; ++i;
			}
			else
				++dvit;
			if (dvit != deepVerts.end()) {  // possible deep split run
				assert(*_mt->vertexFaceTriangle(*dvit) != 0x80000000);
				std::vector<materialTriangles::neighborNode> nei;
				_mt->getNeighbors(*dvit, nei);
				auto nit = nei.begin();
				assert(nit->triangle > -1);
				while (nit != nei.end()) {
					if (_mt->triangleMaterial(nit->triangle) == 4)
						break;
					++nit;
				}
				if (nit == nei.end()) {
					oppBotVerts[i] = *dvit;;
					++dvit; ++i;
				}
			}
		}
		if (dvit == deepVerts.end())
			break;
		// just finished a string of unconnected bottom vertices. If not done look for a flap bottom to cut.
		if (dvit != deepVerts.begin())
			startVertex = -1;
		std::list<int> dVsub;
		while (true) {
			dVsub.push_back(*dvit);
			++dvit;
			if (dvit == deepVerts.end() || *_mt->vertexFaceTriangle(*dvit) == 0x80000000)
				break;
		}
		if (dVsub.size() > 1)
			_solidRecutRequired = true;
		int endV;  // COURT - not yet using endV
		if (dvit != deepVerts.end())
			endV = -1;
		else
			endV = endVertex;
		std::vector<int> oppBotSub;
		flapSurfaceSplitter(startVertex, endV, dVsub, oppBotSub);  // endVertex
		for (auto& obs : oppBotSub)
			oppBotVerts[i++] = obs;
	}
	// create all the new texture coords for incision edges here. Texture seams OK
	float pathLen = 0.0f;
	std::vector<int> topTx, botTx, oppTopTx, oppBotTx;
	topTx.reserve(oppVerts.size()); botTx.reserve(oppVerts.size()); oppTopTx.reserve(oppVerts.size()); oppBotTx.reserve(oppVerts.size());
	int txId = _mt->numberOfTextures();
	float tx[2] = { 0.0f, pathLen };
	Vec3f lV, V;
	_mt->getVertexCoordinate(oppVerts[0], lV.xyz);
	for (int n = oppVerts.size(), i = 0; i < n; ++i) {
		if (i > 0) {
			_mt->getVertexCoordinate(oppVerts[i], V.xyz);
			pathLen += (V - lV).length() * 0.5f;  // COURT note model specific fudge factor
			lV = V;
		}
		tx[1] = pathLen;
		tx[0] = 0.0f;
		topTx.push_back(_mt->addTexture());
		_mt->setTexture(txId++, tx);
		oppTopTx.push_back(_mt->addTexture());
		_mt->setTexture(txId++, tx);
		tx[0] = 1.0f;
		botTx.push_back(_mt->addTexture());
		_mt->setTexture(txId++, tx);
		oppBotTx.push_back(_mt->addTexture());
		_mt->setTexture(txId++, tx);
	}
	auto tvit = topVerts.begin();
	dvit = deepVerts.begin();
	// for new flapCutter()
	int triV[3], triTx[3], lastTv = *tvit, lastDv = *dvit;
	++tvit;  ++dvit;
	i = 1;
	// see incision convention in the header file.
	// Two opposing quads will be triangulated in that format.
	auto triangulateQuad = [&](int lt, int t) {
		Vec3f lv, vn;
		_mt->getVertexCoordinate(lt, lv.xyz);
		_mt->getVertexCoordinate(t, vn.xyz);
		pathLen += (vn - lv).length();
		auto ltit =_deepBed.find(lt);
		assert(ltit != _deepBed.end());
		auto tit = _deepBed.find(t);
		assert(tit != _deepBed.end());
		Vec3f tet[3];
		_vbt->barycentricWeightToGridLocus(_vbt->tetCentroid(_vbt->_vertexTets[lt]), _vbt->_barycentricWeights[lt], tet[0]);
		_vbt->barycentricWeightToGridLocus(_vbt->tetCentroid(_vbt->_vertexTets[t]), _vbt->_barycentricWeights[t], tet[1]);
		tet[2] = tit->second.gridLocus;
		for (int k = 0; k < 3; ++k)
			tet[k] -= ltit->second.gridLocus;
		if ((tet[0] ^ tet[1])*tet[2] < 0.0f) {  // tesselate both sides concave.  Bottom side triangle number always immediately follows top.
			triV[0] = oppVerts[i];
			triV[1] = oppVerts[i - 1];
			triV[2] = ltit->second.deepMtVertex;
			triTx[0] = oppTopTx[i];
			triTx[1] = oppTopTx[i - 1];
			triTx[2] = botTx[i - 1];
			_mt->addTriangle(triV, 3, triTx);
			triV[0] = ltit->second.deepMtVertex;
			triV[1] = tit->second.deepMtVertex;
			triV[2] = oppVerts[i];
			triTx[0] = botTx[i - 1];
			triTx[1] = botTx[i];
			triTx[2] = oppTopTx[i];
			_mt->addTriangle(triV, 3, triTx);
			// do opposite side
			triV[0] = lt;
			triV[1] = t;
			triV[2] = oppBotVerts[i];
			triTx[0] = topTx[i - 1];
			triTx[1] = topTx[i];
			triTx[2] = oppBotTx[i];
			_mt->addTriangle(triV, 3, triTx);
			triV[0] = oppBotVerts[i];
			triV[1] = oppBotVerts[i - 1];
			triV[2] = lt;
			triTx[0] = oppBotTx[i];
			triTx[1] = oppBotTx[i - 1];
			triTx[2] = topTx[i - 1];
			_mt->addTriangle(triV, 3, triTx);
		}
		else {
			triV[0] = oppVerts[i];
			triV[1] = oppVerts[i - 1];
			triV[2] = tit->second.deepMtVertex;
			triTx[0] = oppTopTx[i];
			triTx[1] = oppTopTx[i - 1];
			triTx[2] = botTx[i];
			_mt->addTriangle(triV, 3, triTx);
			triV[0] = ltit->second.deepMtVertex;
			triV[1] = tit->second.deepMtVertex;
			triV[2] = oppVerts[i - 1];
			triTx[0] = botTx[i - 1];
			triTx[1] = botTx[i];
			triTx[2] = oppTopTx[i - 1];
			_mt->addTriangle(triV, 3, triTx);
			// do opposite side
			triV[0] = lt;
			triV[1] = t;
			triV[2] = oppBotVerts[i - 1];
			triTx[0] = topTx[i - 1];
			triTx[1] = topTx[i];
			triTx[2] = oppBotTx[i - 1];
			_mt->addTriangle(triV, 3, triTx);
			triV[0] = oppBotVerts[i];
			triV[1] = oppBotVerts[i - 1];
			triV[2] = t;
			triTx[0] = oppBotTx[i];
			triTx[1] = oppBotTx[i - 1];
			triTx[2] = topTx[i];
			_mt->addTriangle(triV, 3, triTx);
		}
	};
	while (tvit != topVerts.end()){
		triangulateQuad(lastTv, *tvit);
		lastTv = *tvit;
		++tvit;
		++i;
	}
	assert(i == topVerts.size());
	_mt->findAdjacentTriangles(true);
	// get all triangles on the edge of a skin cut.
	_inExCisionTriangles.clear();
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
		if (_mt->triangleMaterial(i) != 2)
			continue;
		int at[3], ae[3];
		_mt->triangleAdjacencies(i, at, ae);
		for (int j = 0; j < 3; ++j) {
			if (_mt->triangleMaterial(at[j]) == 3) {
				_inExCisionTriangles.push_back(i);
				break;
			}
		}
	}
	return true;
}

void skinCutUndermineTets::flapSurfaceSplitter(const int startVertex, const int endVertex, std::list<int> &vertexCutLine, std::vector<int> &oppositeVertices)
{  // divides an _mt surface along a continuous vertexCutLine into two free edge sides in places where the cutting surface is already connected to the _mt surface.
	// Produces an oppositeVertices array of new or repeat vertices corresponding to vertexCutLine.  If cut successful, oppositeVertex will contain the corresponding
	// newly created vertices.  If not cut the vertex from vertexCutLine is simply repeated in oppositeVertex. startVertex and endVertex are present to allow Tin and Tout.
	// If one of these is -1 the first or last vertex in vertexCutLine is assumed to be the start or end vertex and these vertices will be endpoints in the split.
	// In contrast if start/endVertex is > -1 that vertex is the end of the split allowing the first or last point to be doubled for Tin or Tout.
	oppositeVertices.clear();
	oppositeVertices.assign(vertexCutLine.size(), -1);
	auto tvit = vertexCutLine.begin();
	std::vector<materialTriangles::neighborNode> nei;
	int tri = -1, vertIdx, lastV, ovn = 0;
	if (startVertex > -1)
		lastV = startVertex;
	else{
		lastV = *tvit;
		oppositeVertices[ovn++] = lastV;
		++tvit;
	}
	if (endVertex > -1)
		vertexCutLine.push_back(endVertex);
	_mt->findAdjacentTriangles(true);
	assert(*_mt->vertexFaceTriangle(*tvit) != 0x80000000);  // would signal an unconnected startVertex
	_mt->getNeighbors(*tvit, nei);
	auto nit = nei.begin();
	while ( nit != nei.end() ){
		if (nit->vertex == lastV)
			break;
		++nit;
	}
	assert(nit != nei.end());
	int *tr = _mt->triangleVertices(nit->triangle);
	for (int j = 0; j < 3; ++j){
		if (tr[j] == lastV && tr[(j + 1) % 3] == *tvit){
			tri = nit->triangle;
			vertIdx = (j+1) % 3;
			break;
		}
	}
	assert(tri > -1);
	do{
		oppositeVertices[ovn] = _mt->addVertices(1);
		Vec3f v;
		_mt->getVertexCoordinate(*tvit, v.xyz);
		_mt->setVertexCoordinate(oppositeVertices[ovn], v.xyz);
		int lastTx, txId, materialNow = _mt->triangleMaterial(tri);
		int* ttx = _mt->triangleTextures(tri);
		lastTx = txId = ttx[vertIdx];;
		const float* fp = _mt->getTexture(txId);
		float tx[2] = { *fp, fp[1] };
		txId = _mt->addTexture();
		_mt->setTexture(txId, tx);
		auto dbit = _deepBed.find(*tvit);
		if (dbit != _deepBed.end()) {  // is a flap top vertex
			assert(dbit->second.deepMtVertex > -1);
			_deepBed.insert(std::make_pair(oppositeVertices[ovn], dbit->second));
		}
		assert(oppositeVertices[ovn] == _vbt->_vertexTets.size());
		_vbt->_vertexTets.push_back(_vbt->_vertexTets[*tvit]);
		_vbt->_barycentricWeights.push_back(_vbt->_barycentricWeights[*tvit]);
		tr = _mt->triangleVertices(tri);
		tr[vertIdx] = oppositeVertices[ovn];
		_mt->triangleTextures(tri)[vertIdx] = txId;
		unsigned int adj;
		lastV = *tvit;
		++tvit;
		while (tr[(vertIdx + 1) % 3] != *tvit) {
			int at[3], ae[3];
			_mt->triangleAdjacencies(tri, at, ae);
//			adj = _mt->triangleAdjacencies(tri)[vertIdx];
//			assert(adj != 3);
//			tri = adj >> 2;
			tri = at[vertIdx];
			vertIdx = (ae[vertIdx] + 1) % 3;
			tr = _mt->triangleVertices(tri);
			int* ttx = _mt->triangleTextures(tri);
			assert(tr[vertIdx] == lastV);
			tr[vertIdx] = oppositeVertices[ovn];
			if (ttx[vertIdx] != lastTx) {  // crossed a texure seam
				int twoTx[2];
				twoTx[0] = txId;
				fp = _mt->getTexture(ttx[vertIdx]);
				tx[0] = fp[0]; tx[1] = fp[1];
				lastTx = txId = _mt->addTexture();
				_mt->setTexture(txId, tx);
				twoTx[1] = txId;
				if (_mt->triangleMaterial(tri) != materialNow)
					materialNow = _mt->triangleMaterial(tri);
//				else  // this is a single material texture seam vertex  // COURT fix me
//					_mt->addOneMaterialTextureSeamVertex(lastV, twoTx);
			}
			_mt->triangleTextures(tri)[vertIdx] = txId;
		}
		++ovn;
		vertIdx = (vertIdx + 1) % 3;
	} while (*tvit != vertexCutLine.back());
	if (endVertex > -1)
		vertexCutLine.pop_back();
	else
		oppositeVertices[ovn] = vertexCutLine.back();
}

bool skinCutUndermineTets::planeCutSurfaceLine(const int startTopV, const int endTopV, const int startDeepV, const int endDeepV, std::list<int> &newTopVerts, std::list<int> &newDeepVerts)
{
	newTopVerts.clear();
	newDeepVerts.clear();
	Vec3f startP, endP, tempP, planeN;
	float planeD;
	startP.set((float(&)[3])*_mt->vertexCoordinate(startTopV));
	endP.set((float(&)[3])*_mt->vertexCoordinate(endTopV));
	tempP.set((float(&)[3])*_mt->vertexCoordinate(startDeepV));
	planeN = startP - tempP;
	tempP.set((float(&)[3])*_mt->vertexCoordinate(endDeepV));
	planeN += endP - tempP;
	planeN = planeN ^ (endP - startP);
	if (planeN.normalize() < 1e-24)
		return false;
	planeD = planeN*endP;
	auto intersectTriEdge = [&](const int triNum, const int edgeNum, float &edgeParam, Vec3f &intersect) ->bool{
		int *tr = _mt->triangleVertices(triNum);
		Vec3f V, W;
		V.set((const float(&)[3])*_mt->vertexCoordinate(tr[edgeNum]));
		W.set((const float(&)[3])*_mt->vertexCoordinate(tr[(edgeNum + 1) % 3]));
		float vd = planeN*V - planeD, wd = planeN*W - planeD;
		if (!(signbit(vd) ^ signbit(wd)))
			return false;
		edgeParam = vd / (vd - wd);
		intersect = W*edgeParam + V*(1.0f - edgeParam);
		return true;
	};
	unsigned int i, j, n = _mt->numberOfTriangles();
	float edgeParam;
	for (i = 0; i < n; ++i){
		if (_mt->triangleMaterial(i) != 2)
			continue;
		int *tr = _mt->triangleVertices(i);
		for (j = 0; j < 3; ++j)
			if (tr[j] == startTopV)
				break;
		if (j > 2)
			continue;
		if (tr[(j + 1) % 3] == endTopV || tr[(j + 2) % 3] == endTopV)
			return true;
		++j;
		j %= 3;
		if (!intersectTriEdge(i, j, edgeParam, tempP))
			continue;
		if ((startP - tempP).length2() < (startP - endP).length2() && (endP - tempP).length2() < (startP - endP).length2())
			break;
	}
	if (i == n)  // no path start
		return false;
	std::vector<unsigned int> triEdges;
	std::vector<float> edgeParams;
	triEdges.push_back(_mt->triAdjs(i)[j]);
	assert(triEdges.back() != 3);
	edgeParams.push_back(1.0f - edgeParam);
	do{
		i = triEdges.back() >> 2;
		j = triEdges.back() & 3;
		int k;
		if (_mt->triangleVertices(i)[(j + 2) % 3] == endTopV)
			break;
		for (k = 1; k < 3; ++k) {
			if (intersectTriEdge(i, (j + k) % 3, edgeParam, tempP))
				break;
		}
		if (k > 2)
			return false;
		j = (j + k) % 3;
		i = _mt->triAdjs(i)[j];
		if (i == 3)
			return false;
		if (i == triEdges.front()){
			return false;
		}
		else{
			triEdges.push_back(i);
			edgeParams.push_back(1.0f - edgeParam);
		}
	} while (true);
	n = triEdges.size();
	for (i = 0; i < n; ++i){
		j = triEdges[i];
		float uv[2] = {0.0f, 0.0f};
		if ((j & 3) < 1)
			uv[0] = (float)edgeParams[i];
		else if ((j & 3) < 2){
			uv[1] = (float)edgeParams[i];
			uv[0] = 1.0f - uv[1];
		}
		else if ((j & 3) < 3)
			uv[1] = 1.0f - (float)edgeParams[i];
		else
			assert(false);
		int topV, botV;
		createFlapTopBottomVertices(j >> 2, uv, topV, botV);
		newTopVerts.push_back(topV);
		newDeepVerts.push_back(botV);
	}
	return true;
}

bool skinCutUndermineTets::trianglePath(const int triStart, const int endTriangle, const int searchMaterial,  std::vector<int> &triPath)
{  // used in undermining
	if (triStart == endTriangle){
		triPath.assign(1,triStart);
		return true;
	}
	Vec3f vf0, vf1;
	_mt->getTriangleNormal(triStart, vf0, true);
	_mt->getTriangleNormal(endTriangle, vf1, true);
	Vec3f planeNrm(vf0+vf1);
	if (planeNrm.length2() < 1e-16)
		return false;
	float uv[2] = {0.33333f, 0.33333f};
	_mt->getBarycentricPosition(triStart, uv, vf0.xyz);
	_mt->getBarycentricPosition(endTriangle, uv, vf1.xyz);
	planeNrm = planeNrm^(vf0 - vf1);
	float planeD = planeNrm*vf1;
	std::vector<unsigned int> triEdges, try2;
	std::vector<float> params, try2P;
	int mat;
	auto planeDist = [&](int v) ->float{
		Vec3f v3((float(&)[3])*_mt->vertexCoordinate(v));
		return v3*planeNrm - planeD;
	};
	int count = 0, *tr = _mt->triangleVertices(triStart);
	float dNow, dNext, d;
	// two possible triangle paths in a closed manifold surface.
	for (int i = 0; i < 3; ++i) {
		dNow = planeDist(tr[i]);
		dNext = planeDist(tr[(i + 1) % 3]);
		if (std::signbit(dNow) == std::signbit(dNext))
			continue;
		++count;
		unsigned int te = _mt->triAdjs(triStart)[i];
		while ((te >> 2) != triStart && te != 3 && (te >> 2) != endTriangle) {
			mat = _mt->triangleMaterial(te >> 2);
			if (mat != searchMaterial && mat != 10) {
				te = 3;
				break;
			}
			if (count < 2) {
				triEdges.push_back(te);
				assert(fabs(-dNow + dNext) > 1e-16f);
				params.push_back(dNext / (-dNow + dNext));
			}
			else {
				try2.push_back(te);
				assert(fabs(-dNow + dNext) > 1e-16f);
				try2P.push_back(dNext / (-dNow + dNext));
			}
			tr = _mt->triangleVertices(te >> 2);
			assert(std::signbit(planeDist(tr[te & 3])) != std::signbit(dNow));
			d = planeDist(tr[((te & 3) + 2) % 3]);
			if (std::signbit(d) != std::signbit(dNow)) {
				te = _mt->triAdjs(te >> 2)[((te & 3) + 1) % 3];
				dNext = d;
			}
			else {
				te = _mt->triAdjs(te >> 2)[((te & 3) + 2) % 3];
				dNow = d;
			}
		}
		if (te == 3 || (te >> 2) == triStart) { // bad path
			if (count < 2) {
				if (triEdges.empty())
					triEdges.push_back(0xffffffff);
				else
					triEdges.front() = 0xffffffff;
			}
			else {
				if (try2.empty())
					try2.push_back(0xffffffff);
				else
					try2.front() = 0xffffffff;
			}
		}
		else {
			if (count < 2) {
				triEdges.push_back(te);
				assert(fabs(-dNow + dNext) > 1e-16f);
				params.push_back(dNext / (-dNow + dNext));
			}
			else {
				try2.push_back(te);
				assert(fabs(-dNow + dNext) > 1e-16f);
				try2P.push_back(dNext / (-dNow + dNext));
			}
		}
		tr = _mt->triangleVertices(triStart);
	}
	if (count < 2 || (try2.front() == 0xffffffff && triEdges.front() == 0xffffffff))
		return false;
	triPath.clear();
	if (triEdges.front() == 0xffffffff) {
		if (try2.front() == 0xffffffff)
			return false;
		else {
			triPath.reserve(try2.size() + 1);
			triPath.push_back(triStart);
			for (auto tp : try2)
				triPath.push_back(tp >> 2);
			return true;
		}
	}
	if (try2.front() == 0xffffffff){
		triPath.reserve(triEdges.size() + 1);
		triPath.push_back(triStart);
		for (auto tp : triEdges)
			triPath.push_back(tp >> 2);
		return true;
	}
	auto pathLength = [&](std::vector<unsigned int> &tre, std::vector<float> & prm) ->float {
		float dist = 0.0f;
		Vec3f vtx0, vtx1, vtxL;
		int *trV;
		for (int n = (int)tre.size(), i = 0; i < n; ++i) {
			trV = _mt->triangleVertices(tre[i] >> 2);
			_mt->getVertexCoordinate(trV[tre[i] & 3], vtx0.xyz);
			_mt->getVertexCoordinate(trV[((tre[i] & 3) + 1) % 3], vtx1.xyz);
			vtx1 -= vtx0;
			vtx0 += vtx1*(float)prm[i];
			if (i >0)
				dist += (vtx0 - vtxL).length();
			vtxL = vtx0;
		}
		return dist;
	};
	if (pathLength(try2, try2P) < pathLength(triEdges, params)) {
		triPath.reserve(try2.size() + 1);
		triPath.push_back(triStart);
		for (auto tp : try2)
			triPath.push_back(tp >> 2);
	}
	else{
		triPath.reserve(triEdges.size() + 1);
		triPath.push_back(triStart);
		for (auto tp : triEdges)
			triPath.push_back(tp >> 2);
	}
	return true;
}

bool skinCutUndermineTets::setDeepBed(materialTriangles *mt, const std::string &deepBedPath, vnBccTetrahedra *activeVnt)
{
	_mt = mt;
	_vbt = activeVnt;
	std::ifstream istr(deepBedPath.c_str());
	if (!istr.is_open()) {
		istr.close();
		return false;
	}
	_deepBed.clear();
	_deepBed.reserve((size_t)(activeVnt->vertexNumber()*1.2f));  // max number of deep verts after very complex lip repair
	_deepBed.max_load_factor(1.2f);
	deepPoint dp;
	dp.deepMtVertex = -1;
	char s[400];
	while (!istr.eof())
	{
		int topVert;
		istr.getline(s, 399);
		sscanf(s, "%ld %f %f %f", &topVert, &dp.gridLocus.X, &dp.gridLocus.Y, &dp.gridLocus.Z);
		// deep point guaranteed to be inside tet grid
		dp.gridLocus -= activeVnt->getMinimumCorner();
		dp.gridLocus *= (float)activeVnt->_unitSpacingInv;
		_deepBed.emplace(topVert, dp);
	}
	istr.close();
	return true;
}

int skinCutUndermineTets::deepPointTetWeight(const std::unordered_map<int, deepPoint>::iterator &dit, Vec3f &baryWeight)
{  // return deepPoint tet number and baryweight from its grid locus
	bccTetCentroid tc;
	_vbt->gridLocusToLowestTetCentroid(dit->second.gridLocus, tc);
	if (_vbt->tetCentroid(_vbt->getVertexTetrahedron(dit->first)) == tc) {
		_vbt->gridLocusToBarycentricWeight(dit->second.gridLocus, tc, baryWeight);
		return _vbt->getVertexTetrahedron(dit->first);
	}
	std::list<int> lt, tp;
	_vbt->centroidTets(tc, lt);
	int tetOut = -1, count = 0;
	while (lt.empty()) {
		tc = _vbt->centroidUpOneLevel(tc);
		_vbt->centroidTets(tc, lt);
		++count;
	}
	if (count > 15)
		throw(std::logic_error("Modelling error.  Deep bed point not inside solid.\n"));
	if (lt.size() < 2)
		tetOut = lt.front();
	else {
		for (auto tet : lt) {
			assert(false);  // COURT debug me for multires
			if (_vbt->decreasingCentroidPath(_vbt->getVertexTetrahedron(dit->first), tetOut, tp))
				break;
		}
	}
	_vbt->gridLocusToBarycentricWeight(dit->second.gridLocus, tc, baryWeight);
	return tetOut;
}

int skinCutUndermineTets::flapBottomTet(const int topTet, const Vec3f &bottomGridLocus)
{  // material coord flap bottom tet finder. Return -1 signals error in deep bed data input
	bccTetCentroid tc;
	_vbt->gridLocusToLowestTetCentroid(bottomGridLocus, tc);
	if (_vbt->tetCentroid(topTet) == tc) {
		return topTet;
	}
	std::list<int> lt, tp;
	_vbt->centroidTets(tc, lt);
	int tetOut = -1, count = 0;
	while (lt.empty()) {
		tc = _vbt->centroidUpOneLevel(tc);
		_vbt->centroidTets(tc, lt);
		++count;
	}
	if (count > 15)
		throw(std::logic_error("Modelling error.  Deep bed point not inside solid.\n"));
	if (lt.size() < 2)
		tetOut = lt.front();
	else {
		for (auto tet : lt) {
			assert(false);  // COURT debug me for multires
//			if (_vbt->decreasingCentroidPath(_vbt->getVertexTetrahedron(dit->first), tetOut, tp))
//				break;
		}
	}
	return tetOut;
}

void skinCutUndermineTets::collectOldUndermineData()
{
	_prevUnd2.clear();
	_prevBot4.clear();
	_prevEdge3.clear();
	_prevBedSingles.clear();
	// flap bottom has its reverse mirror image topology on the top
	std::vector<int> bot4;
	bot4.reserve(300);
	for (int n = _mt->numberOfTriangles(), j, i = 0; i < n; ++i) {
		if (_mt->triangleMaterial(i) != 2)
			continue;
		int* tr = _mt->triangleVertices(i);
		int deepV[3];
		for (j = 0; j < 3; ++j){
			auto dvit = _deepBed.find(tr[j]);
			if (dvit == _deepBed.end() || dvit->second.deepMtVertex < 0)
				break;
			deepV[j] = dvit->second.deepMtVertex;
		}
		if (j < 3)
			continue;
		std::vector<materialTriangles::neighborNode> nei;
		_mt->getNeighbors(deepV[0], nei);
		auto nit = nei.begin();
		while (nit != nei.end()) {
			if (nit->vertex == deepV[2])
				break;
			++nit;
		}
		if(nit == nei.end())
			continue;
		++nit;
		if (nit == nei.end())
			nit = nei.begin();
		if (nit->vertex == deepV[1]) {
			bot4.push_back(nit->triangle);
			_prevUnd2.push_back(i);
		}
	}
	if (bot4.empty())
		return;
	std::vector<int> f3;
	std::set<int> f4;
	std::vector<int> emptyTex;
	emptyTex.reserve(4);
	auto addSingle = [&](int v, int tx) {
		auto pr = _prevBedSingles.insert(std::make_pair(v, emptyTex));
		auto tex = pr.first->second.begin();
		while (tex != pr.first->second.end()) {
			if (*tex == tx)
				break;
			++tex;
		}
		if (tex == pr.first->second.end())
			pr.first->second.push_back(tx);
	};
	for (auto& bt : bot4) {
		unsigned int* adjs = _mt->triAdjs(bt);
		for (int j = 0; j < 3; ++j) {
			if (_mt->triangleMaterial(adjs[j] >> 2) == 5) {  // potentially alterable flap bottom edge as contains unduplicated vertices
				f4.insert(bt);
				int* tr = _mt->triangleVertices(bt);
				int* ttx = _mt->triangleTextures(bt);
				addSingle(tr[j], ttx[j]);
				addSingle(tr[(j+1)%3], ttx[(j + 1) % 3]);
			}
			if (_mt->triangleMaterial(adjs[j] >> 2) == 3) {  // active incision edges
				f3.push_back(adjs[j] >> 2);
				f3.push_back((adjs[j] >> 2) - 1);  // incision convention
			}
		}
	}
	_prevEdge3.reserve(f3.size());
	for (auto& ef : f3) {
		int* tr = _mt->triangleVertices(ef);
		for (int j = 0; j < 3; ++j) {
			if (_prevBedSingles.find(tr[j]) != _prevBedSingles.end())
				_prevEdge3.push_back(ef);
		}
	}
	for (auto& bt : bot4) {
		int* tr = _mt->triangleVertices(bt);
		for (int j = 0; j < 3; ++j) {
			if (_prevBedSingles.find(tr[j]) != _prevBedSingles.end()) {
				f4.insert(bt);
				break;
			}
		}
	}
	_prevBot4.assign(f4.begin(), f4.end());
	_prevEdge3.shrink_to_fit();
}

bool skinCutUndermineTets::addUndermineTriangle(const int triangle, const int undermineMaterial, bool incisionConnect)
{
	std::vector<int> *edgeTriangles;
	if (undermineMaterial == 2)
		edgeTriangles = &_inExCisionTriangles;
	else{
		assert(undermineMaterial == 7);
		if (_periostealCutEdgeTriangles.empty()) {
			for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
				int mat;
				if ((mat = _mt->triangleMaterial(i)) != 7 && mat != 8)
					continue;
				unsigned int* adjs = _mt->triAdjs(i);
				for (int j = 0; j < 3; ++j) {
					mat = _mt->triangleMaterial(adjs[j] >> 2);
					if (mat != 7 && mat != 8 && mat != 1) {
						_periostealCutEdgeTriangles.push_back(i);
						break;
					}
				}
			}
		}
		edgeTriangles = &_periostealCutEdgeTriangles;
	}
	float uv[2] = {0.33f, 0.33f}, d, minD = FLT_MAX;
	std::vector<int> triPath;
	if (incisionConnect) {
		Vec3f newP, now;
		int closeT = -1;
		_mt->getBarycentricPosition(triangle, uv, newP.xyz);
		for (auto &tIdx : *edgeTriangles) {
			_mt->getBarycentricPosition(tIdx, uv, now.xyz);
			d = (now - newP).length2();
			if (d < minD) {
				minD = d;
				closeT = tIdx;
			}
		}
		if (closeT < 0)
			return false;
		if (!trianglePath(closeT, triangle, undermineMaterial, triPath))
			return false;
	}
	if (_prevUndermineTriangle > -1) {
		std::vector<int> tep2;
		if (!trianglePath(triangle, _prevUndermineTriangle, undermineMaterial, tep2))
			return false;
		if(!triPath.empty())
			triPath.pop_back();
		triPath.insert(triPath.end(), tep2.begin(), tep2.end());
		for (auto &tp : triPath) {
			showPriorUndermine(tp);
			_mt->setTriangleMaterial(tp, 10);
		}
		// now have closed path from cut edge to cut edge
		closeUndermineHoles(triPath, undermineMaterial);
	}
	else {
		if(undermineMaterial == 2)  // COURT - should probably do this for periosteal undermines too. Later.
			collectOldUndermineData();
		_trisUnderminedNow.clear();
		_trisUnderminedNow.assign(_mt->numberOfTriangles(), false);
		for (auto &tp : triPath) {
			showPriorUndermine(tp);
			_mt->setTriangleMaterial(tp, 10);
		}
	}
	_prevUndermineTriangle = triangle;
	return true;
}

void skinCutUndermineTets::undermineSkin() {
	std::vector<int> undTris, newTris;
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
		if (_mt->triangleMaterial(i) != 10)
			continue;
		_mt->setTriangleMaterial(i, 2);
		undTris.push_back(i);
	}
	struct deepVtx {
		int deepV;
		std::vector<int> deepTx;
	}empDv;
	empDv.deepV = -1;
	empDv.deepTx.clear();
	empDv.deepTx.reserve(4);
	std::unordered_map<int, deepVtx> undV;
	undV.reserve(undTris.size() << 1);
	auto findPrevTexture = [&](int tx, std::vector<int> &prevTex) ->int {
		int ret = -1;
		float *pt, *fp = _mt->getTexture(tx);
		for (auto t : prevTex) {
			pt = _mt->getTexture(t);
			if (pt[0] == fp[0] && pt[1] == fp[1]) {
				ret = t;
				break;
			}
		}
		return ret;
	};
	for(auto &t : undTris){
		if (!std::binary_search(_prevUnd2.begin(), _prevUnd2.end(), t)) {
			newTris.push_back(t);
			int *top = _mt->triangleVertices(t);
			int* topTx = _mt->triangleTextures(t);
			for (int j = 0; j < 3; ++j) {
				auto pr = undV.insert(std::make_pair(top[j], empDv));
				if (pr.second) {
					auto dit = _deepBed.find(top[j]);
					if (dit->second.deepMtVertex > -1)
						pr.first->second.deepV = dit->second.deepMtVertex;
					else {
						assert(dit != _deepBed.end());
						pr.first->second.deepV = createDeepBedVertex(dit); // puts result in dit
					}
				}
				int pt = -1;  // previous texture index
				auto pb = _prevBedSingles.find(pr.first->second.deepV);
				if (pb != _prevBedSingles.end()) {
					pt = findPrevTexture(topTx[j], pb->second);
					if (pt > -1) {
						int pt2 = findPrevTexture(pt, pr.first->second.deepTx);
						if (pt2 < 0)
							pr.first->second.deepTx.push_back(pt);
						else
							assert(pt2 == pt);
					}
				}
				if(pt < 0) {
					pt = findPrevTexture(topTx[j], pr.first->second.deepTx);
					if (pt < 0)
						pr.first->second.deepTx.push_back(cloneTexture(topTx[j]));
				}
			}
		}
	}
	// all top vertices are doubled (for bed and flap bottom) except vertices on an edge with both sides material 2 and one triangle undermined and the other not.
	// because all of undTris searched, previous undermines taken into account
	std::set<int> oneDeepV;
	std::vector<int> underminedEdgeTris;
	for (auto &tri : newTris) {
		for (int k = 0; k < 3; ++k) {
			int t = (_mt->triAdjs(tri)[k] >> 2);
			if (_mt->triangleMaterial(t) == 2 && !std::binary_search(undTris.begin(), undTris.end(), t)) {
				oneDeepV.insert(_mt->triangleVertices(tri)[k]);
				oneDeepV.insert(_mt->triangleVertices(tri)[(k + 1) % 3]);
			}
			if (_mt->triangleMaterial(t) == 3)  // there are undermined edge tris
				underminedEdgeTris.push_back(t);  // only need top one
		}
	}
	for (auto nt : newTris) {  // create flap bed
		int vd[3], dTx[3], * top = _mt->triangleVertices(nt), *tx = _mt->triangleTextures(nt);
		for (int i = 0; i < 3; ++i) {
			auto dv = undV.find(top[i]);
			assert(dv != undV.end());
			vd[i] = dv->second.deepV;
			int pt = findPrevTexture(tx[i], dv->second.deepTx);
			assert(pt > -1);
			dTx[i] = pt;
			auto sit = _collisionSpokes.find(top[i]);
			if (sit != _collisionSpokes.end())
				_deepSpokesNow.emplace(vd[i], sit->second);
		}
		_mt->addTriangle(vd, 5, dTx);
	}
	std::map<int, deepVtx*> bottomDoubles;
	for (auto& uv : undV) {
		if (oneDeepV.find(uv.first) == oneDeepV.end()) {  // should be doubled
			float vtx[3];
			int oldV = uv.second.deepV;
			_mt->getVertexCoordinate(oldV, vtx);
			uv.second.deepV = _mt->addVertices(1);
			_mt->setVertexCoordinate(uv.second.deepV, vtx);
			assert(_vbt->_vertexTets.size() == uv.second.deepV);
			_vbt->_vertexTets.push_back(_vbt->getVertexTetrahedron(oldV));
			_vbt->_barycentricWeights.push_back(*_vbt->getVertexWeight(oldV));
			auto dbit = _deepBed.find(uv.first);
			assert(dbit != _deepBed.end());
			dbit->second.deepMtVertex = uv.second.deepV;
			for (auto& dt : uv.second.deepTx) {
				dt = cloneTexture(dt);
			}
			bottomDoubles.insert(std::make_pair(oldV, &uv.second));  // need these to fix old border tris from a previous undermine
		}
	}
	for (auto nt : newTris) {
		int vb[3], bTx[3], * top = _mt->triangleVertices(nt), *ttx = _mt->triangleTextures(nt);
		for (int i = 0; i < 3; ++i) {
			auto dv = undV.find(top[i])->second;  // no need to test find. already done in bed creation
			int pt = findPrevTexture(ttx[i], dv.deepTx);
			assert(pt > -1);
			if (i < 2) {
				vb[1 - i] = dv.deepV;
				bTx[1 - i] = pt;
			}
			else {
				vb[i] = dv.deepV;
				bTx[i] = pt;
			}
		}
		_mt->addTriangle(vb, 4, bTx);
	}
	// now fix old undermined tris with a single deep vertex needing to be doubled
	auto fixEdgeTris = [&](int topT) {  // will also fix its lower successor per incision convention. Textures already assigned by incision process.
		int* tr = _mt->triangleVertices(topT);
		int* trBot = _mt->triangleVertices(topT + 1);
		int v0 = undV[tr[0]].deepV, v1 = undV[tr[1]].deepV;
		trBot[0] = v1;
		trBot[1] = v0;
		if (_mt->triAdjs(topT)[1] >> 2 == topT + 1)
			tr[2] = v0;
		else {
			assert(_mt->triAdjs(topT)[2] >> 2 == topT + 1);
			tr[2] = v1;
		}
	};
	auto fixOldTris = [&](int tri) {
		int *tr = _mt->triangleVertices(tri);
		int* ttx = _mt->triangleTextures(tri);
		for (int i = 0; i < 3; ++i) {
			auto dd = bottomDoubles.find(tr[i]);
			if (dd != bottomDoubles.end()) {
				tr[i] = dd->second->deepV;
				if (_mt->triangleMaterial(tri) == 4) {  // leave fat textures as they are
					// texture seams are a terrible pain
					int pt = findPrevTexture(ttx[i], dd->second->deepTx);
					if (pt < 0) {
						pt = cloneTexture(ttx[i]);
						dd->second->deepTx.push_back(pt);
					}
					ttx[i] = pt;
				}
			}
		}
	};
	for (auto &t : underminedEdgeTris)
		fixEdgeTris(t);
	for (auto t : _prevBot4)
		fixOldTris(t);
	for (auto t : _prevEdge3)
		fixOldTris(t);

//	testIncisionsDeepBed();  // COURT - nuke after debug

	_mt->findAdjacentTriangles(true);
	_prevUndermineTriangle = -1;
}

void skinCutUndermineTets::clearCurrentUndermine(const int underminedTissue){
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i){
		if (_mt->triangleMaterial(i) == 10)
			_mt->setTriangleMaterial(i, underminedTissue);
	}
	_prevUndermineTriangle = -1;
}

void skinCutUndermineTets::showPriorUndermine(int priorTriangle)
{
	int undermineMaterial = _mt->triangleMaterial(priorTriangle);
	std::deque<int> prevUnd;
	prevUnd.push_front(priorTriangle);
	// non-recursively flood fill this hole with previously undermined top triangles
	while (!prevUnd.empty()) {
		for (int j = 0; j < 3; ++j) {
			int tri = _mt->triAdjs(prevUnd.front())[j] >> 2;
			if (_trisUnderminedNow[tri] || _mt->triangleMaterial(tri) != undermineMaterial || !std::binary_search(_prevUnd2.begin(), _prevUnd2.end(), tri))
				continue;
			_trisUnderminedNow[tri] = true;
			prevUnd.push_back(tri);
		}
		_mt->setTriangleMaterial(prevUnd.front(), 10);
		prevUnd.pop_front();
	}
}

bool skinCutUndermineTets::closeUndermineHoles(std::vector<int> &trianglePath, const int undermineMaterial)
{
	std::list<std::vector<int> > holes;
	for (auto &t : trianglePath)
		_trisUnderminedNow[t] = true;
	// find hole seeds on either side of path
	auto tunCopy = _trisUnderminedNow;
	holes.push_back(std::vector<int>());
	for (auto &t : trianglePath){
		for (int i = 0; i < 3; ++i) {
			int tri = _mt->triAdjs(t)[i] >> 2;
			if (tunCopy[tri] || _mt->triangleMaterial(tri) != undermineMaterial)
				continue;
			tunCopy[tri] = true;
			holes.back().push_back(tri);
			std::deque<int> holeQ;
			holeQ.push_front(tri);
			// non-recursively flood fill this hole
			while (!holeQ.empty()) {
				for (int j = 0; j < 3; ++j) {
					tri = _mt->triAdjs(holeQ.front())[j] >> 2;
					if (tunCopy[tri] || _mt->triangleMaterial(tri) != undermineMaterial)  //  || _underminedTriangles.find(t) != _underminedTriangles.end()
						continue;
					tunCopy[tri] = true;
					holes.back().push_back(tri);
					holeQ.push_back(tri);
				}
				holeQ.pop_front();
			}
			if (!holes.back().empty())  // seed grew. start another hole
				holes.push_back(std::vector<int>());
		}
	}
	if (holes.back().empty())
		holes.pop_back();
	if(holes.empty())
		return true;
	auto hitMax = holes.begin();
	auto hit = hitMax;
	++hit;
	while (hit != holes.end()){
		if (hit->size() > hitMax->size())
			hitMax = hit;
		++hit;
	}
	holes.erase(hitMax);  // remove region outside undermine
	for (auto &h : holes){
		auto hit = h.begin();
		while (hit != h.end()) {
			_trisUnderminedNow[*hit] = true;
			_mt->setTriangleMaterial(*hit, 10);
			++hit;
		}
	}
	return true;
}

float skinCutUndermineTets::closestSkinIncisionPoint(const Vec3f xyz, int& triangle, int& edge, float& param)
{	// Input xyz, then overwrites all 4 with the point on the nearest incision edge.  Incision edge facing in opposite direction not selected.
	_mt->findAdjacentTriangles();
	triangle = -1;
	edge = -1;
	param = FLT_MAX;
	float minDsq = FLT_MAX;
	float ret = FLT_MAX;
	for (int n=_mt->numberOfTriangles(), i = 0; i < n; ++i) {
		if (_mt->triangleMaterial(i) != 3)
			continue;
		unsigned int adj = _mt->triAdjs(i)[0];
		if (_mt->triangleMaterial(adj >> 2) != 2)  // incision convention
			continue;
		Vec3f W, P;
		int* tr = _mt->triangleVertices(i);
		_mt->getVertexCoordinate(tr[0], P.xyz);
		_mt->getTriangleNormal(i, W, false);
		if (W * (P - xyz) < 0.0f)  // edge facing wrong direction
			continue;
		_mt->getVertexCoordinate(tr[1], W.xyz);
		W -= P;
		float lenSq, p = (W * (xyz - P)) / (W * W);
		if (p < 0.0f)
			p = 0.0f;
		if (p > 1.0f)
			p = 1.0f;
		W = W * p + P;
		lenSq = (xyz - W).length2();
		if (lenSq < minDsq) {
			minDsq = lenSq;
			param = 1.0f - p;
			triangle = adj >> 2;
			edge = adj & 3;
			ret = true;
		}
	}
	if (triangle > -1)
		ret = sqrt(minDsq);
	return ret;
}

int skinCutUndermineTets::addTinEdgeVertex(const Vec3f &closePoint, const Vec3f &nextConnectedPoint)
{  // find the correct side of incision from usually 2 choices; 1 if part of a moved flap
	// on 7_21_2021 added T out possibility of finding first point in this incision under construction
	int minTri = -1;
	float t, dsq, minD = 1e32f, minParam = 1e32f;
	Vec3f P, dir, closeDir = closePoint - nextConnectedPoint;
	closeDir.normalize();
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i){
		if (_mt->triangleMaterial(i) != 3)
			continue;
		// by new incision convention top cut edge is 0
		if (_mt->triangleMaterial(_mt->triAdjs(i)[0] >> 2) != 2)
			continue;
		// this is a top incision edge
		_mt->getTriangleNormal(i, dir, true);
		if (closeDir*dir < 0.0f) // opposite side
			continue;
		_mt->getVertexCoordinate(_mt->triangleVertices(i)[0], P.xyz);
		_mt->getVertexCoordinate(_mt->triangleVertices(i)[1], dir.xyz);
		dir -= P;
		t = dir*dir;
		if (t < 1e-16f){
			dsq = (closePoint - P).length2();
			t = 0.0f;
		}
		else{
			t = ((closePoint - P)*dir) / t;
			if (t < 0.0f)
				t = 0.0f;
			if (t > 1.0f)
				t = 1.0f;
			dsq = (dir*t + P - closePoint).length2();
		}
		if (dsq < minD){
			minD = dsq;
			minTri = i;
			minParam = t;
		}
	}
	if (_firstTopVertex > -1) {  // possible closed polygon as this is a T out
		_mt->getVertexCoordinate(_firstTopVertex, P.xyz);
		dsq = (closePoint - P).length2();
		if (dsq < minD)
			return _firstTopVertex;
	}
	if (minParam < 0.001f)
		return _mt->triangleVertices(minTri)[0];
	if (minParam > 0.999f)
		return _mt->triangleVertices(minTri)[1];
	if (minTri < 0)
		return -1;
	return TinSub(minTri, minParam);
}

int skinCutUndermineTets::TinSub(const int edgeTriangle, const float edgeParam)
{  // from incision convention we know
	assert(_mt->triangleMaterial(edgeTriangle) == 3 && _mt->triangleMaterial(_mt->triAdjs(edgeTriangle)[0]>>2) == 2);
	// get this incision box. Use new incision convention.
	bool tess12;
	if ((_mt->triAdjs(edgeTriangle)[1] >> 2) == edgeTriangle + 1)  // 1-2 quad tesselation
		tess12 = true;
	else {
		assert((_mt->triAdjs(edgeTriangle)[2] >> 2) == edgeTriangle + 1);
		tess12 = false;
	}
	int layers = 1;
	if (_mt->triangleMaterial((_mt->triAdjs(edgeTriangle + 1)[0] >> 2)) == 3)  // surface groove, not flap bottom, so must cut other side
		layers = 2;
	// if minParam is 0 or 1 no need to split edge, but bottom Tx on this side must be doubled.
	auto oneVertexT = [&](int vId) ->int {
		int oldTx, newTx, topV = _mt->triangleVertices(edgeTriangle)[vId], deepV = _mt->triangleVertices(edgeTriangle + 1)[1- vId];
		std::vector<materialTriangles::neighborNode> nei;
		_mt->getNeighbors(deepV, nei);
		auto nit = nei.begin();
		while (nit != nei.end()) {
			if (nit->vertex == topV)
				break;
			++nit;
		}
		assert(nit != nei.end());
		int i, * tr = _mt->triangleVertices(nit->triangle);
		for (i = 0; i < 3; ++i) {
			if (tr[i] == deepV)
				break;
		}
		assert(i < 3);
		int* ttx = _mt->triangleTextures(nit->triangle);
		oldTx = ttx[i];
		newTx = cloneTexture(oldTx);
		while (ttx[i] == oldTx) {
			ttx[i] = newTx;
			if (nit == nei.begin())
				nit = nei.end();
			--nit;
			tr = _mt->triangleVertices(nit->triangle);
			ttx = _mt->triangleTextures(nit->triangle);
			for (i = 0; i < 3; ++i)
				if (tr[i] == deepV)
					break;
		}
		return topV;
	};
	if (edgeParam < 0.005f) 
		return oneVertexT(0);
	if (edgeParam > 0.995f)
		return oneVertexT(1);
	// must split this incision edge quad
	int tris[4][2], verts[4][3], tex[4][3], *tr = _mt->triangleVertices(edgeTriangle), * ttx = _mt->triangleTextures(edgeTriangle);
	tris[0][0] = edgeTriangle;
	tris[1][0] = edgeTriangle + 1;
	verts[0][0] = tr[1];
	verts[0][2] = tr[0];
	tex[0][0] = ttx[1];
	tex[0][2] = ttx[0];
	tr = _mt->triangleVertices(edgeTriangle + 1);
	ttx = _mt->triangleTextures(edgeTriangle + 1);
	verts[1][0] = tr[0];
	verts[1][2] = tr[1];
	tex[1][0] = ttx[0];
	tex[1][2] = ttx[1];
	if (layers > 1) {
		verts[3][0] = verts[1][2];
		verts[3][2] = verts[1][0];
		unsigned int layer1TE = _mt->triAdjs(edgeTriangle + 1)[0];
		tris[1][1] = layer1TE >> 2;
		tris[0][1] = tris[1][1] - 1;
		ttx = _mt->triangleTextures(layer1TE >> 2);
		tex[3][0] = ttx[0];
		tex[3][2] = ttx[1];
		tr = _mt->triangleVertices(tris[0][1]);
		ttx = _mt->triangleTextures(tris[0][1]);
		verts[2][0] = tr[1];
		verts[2][2] = tr[0];
		tex[2][0] = ttx[1];
		tex[2][2] = ttx[0];
	}
	Vec3f gridLocus, baryWeight;
	tris[2][0] = _mt->numberOfTriangles() + 1;
	int tet = parametricMTedgeTet(edgeTriangle, 0, edgeParam, gridLocus);
	_vbt->gridLocusToBarycentricWeight(gridLocus, _vbt->tetCentroid(tet), baryWeight);
	// must have tris01 + 1 = tris11 for incision convention
	unsigned int topTriEdge = _mt->triAdjs(edgeTriangle)[0];
	verts[0][1] = _mt->splitTriangleEdge(topTriEdge >> 2, topTriEdge & 3, 1.0f - edgeParam);
	tex[0][1] = _mt->triangleTextures(edgeTriangle)[1];
	_vbt->_vertexTets.push_back(tet);
	_vbt->_barycentricWeights.push_back(baryWeight);
	tris[3][0] = _mt->numberOfTriangles();
	tet = parametricMTedgeTet(edgeTriangle + 1, 0, 1.0f - edgeParam, gridLocus);
	_vbt->gridLocusToBarycentricWeight(gridLocus, _vbt->tetCentroid(tet), baryWeight);
	verts[3][1] = verts[1][1] = _mt->splitTriangleEdge(edgeTriangle + 1, 0, 1.0f - edgeParam);
	tex[1][1] = _mt->triangleTextures(edgeTriangle + 1)[1];
	tex[3][1] = _mt->numberOfTextures() - 1;  // even if two layers a new texture will have been created since bottom of a material 3 incision groove will have a texture seam
	assert(verts[1][1] == _vbt->_vertexTets.size());
	_vbt->_vertexTets.push_back(tet);
	_vbt->_barycentricWeights.push_back(baryWeight);
	deepPoint dp;
	dp.deepMtVertex = verts[1][1];
	dp.gridLocus = gridLocus;
	_deepBed.insert(std::make_pair(verts[0][1], dp));
	if (layers > 1) {
		tris[2][1] = tris[3][0] + 1;
		tris[3][1] = _mt->numberOfTriangles();
		tet = parametricMTedgeTet(tris[0][1], 0, 1.0f - edgeParam, gridLocus);
		_vbt->gridLocusToBarycentricWeight(gridLocus, _vbt->tetCentroid(tet), baryWeight);
		verts[2][1] = _mt->splitTriangleEdge(tris[0][1], 0, 1.0f - edgeParam);
		tex[2][1] = _mt->triangleTextures(tris[0][1])[1];
		assert(tex[2][1] == _mt->numberOfTextures() - 2);
		assert(verts[2][1] == _vbt->_vertexTets.size());
		_vbt->_vertexTets.push_back(tet);
		_vbt->_barycentricWeights.push_back(baryWeight);
		_deepBed.insert(std::make_pair(verts[2][1], dp));
	}
	// splitTriangleEdge() calls do all needed for materials 2 & 4. Only material 3 triangles in cut edge need to be reconfigured.
	for (int j, i = 0; i < 2; ++i) {
		for (j = 0; j < layers; ++j) {
			tr = _mt->triangleVertices(tris[i << 1][j]);
			tr[0] = verts[j<<1][i + 1];
			tr[1] = verts[j<<1][i];
			tr[2] = verts[(j<<1) + 1][i +(tess12 ? 1 : 0)];
			ttx = _mt->triangleTextures(tris[i << 1][j]);
			ttx[0] = tex[j<<1][i + 1];
			ttx[1] = tex[j<<1][i];
			ttx[2] = tex[(j<<1) + 1][i + (tess12 ? 1 : 0)];
			tr = _mt->triangleVertices(tris[(i<<1) + 1][j]);
			tr[0] = verts[(j<<1) + 1][i];
			tr[1] = verts[(j << 1) + 1][i + 1];
			tr[2] = verts[j<<1][i + (tess12 ? 0 : 1)];
			ttx = _mt->triangleTextures(tris[(i<<1) + 1][j]);
			ttx[0] = tex[(j << 1) + 1][i];
			ttx[1] = tex[(j << 1) + 1][i + 1];
			ttx[2] = tex[j<<1][i + (tess12 ? 0 : 1)];
		}
		// flapSurfaceSplitter() doesn't dup bottom corner texture, so do it here
		if (i<1) 
			tex[1][1] =cloneTexture(tex[1][1]);
	}
	_mt->findAdjacentTriangles(true);
	return verts[0][1];
}

int skinCutUndermineTets::cloneTexture(int textureIndex) {
	float *fp = _mt->getTexture(textureIndex);
	float tx[2] = { fp[0], fp[1] };
	int ret = _mt->addTexture();
	_mt->setTexture(ret, tx);
	return ret;
}

void skinCutUndermineTets::excise(const int triangle)
{
	if (_mt->triangleMaterial(triangle) == 3 || _mt->triangleMaterial(triangle) == 6)
		return;
	bool notUndermined = false;
	std::vector<bool> xTris;
	xTris.assign(_mt->numberOfTriangles(), false);
	std::list<int> rList;
	rList.push_back(triangle);
	int nConnected = 1, halfTris = _mt->numberOfTriangles()>>1;
	xTris[triangle] = true;
	while (!rList.empty()) {
		unsigned int *adjs = _mt->triAdjs(rList.front());
		rList.pop_front();
		for (int j = 0; j < 3; ++j) {
			if (xTris[adjs[j] >> 2] || adjs[j] == 3)
				continue;
			xTris[adjs[j] >> 2] = true;
			rList.push_back(adjs[j] >> 2);  // recurse
			++nConnected;
			if (nConnected > halfTris) {  // huge excision.  COURT a bit of a hack
				notUndermined = true;
				break;
			}
		}
		if (notUndermined)
			break;
	}
	if (notUndermined) {  // top skin subpatch not yet undermined. Undermine this subpatch then excise
		rList.clear();
		rList.push_back(triangle);
		for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i)
			xTris[i] = false;
		xTris[triangle] = true;
		while (!rList.empty()) {
			unsigned int *adjs = _mt->triAdjs(rList.front());
			rList.pop_front();
			for (int j = 0; j < 3; ++j) {
				int adjTri = adjs[j] >> 2;
				if (adjs[j] == 3 || xTris[adjTri] || _mt->triangleMaterial(adjTri) != 2)
					continue;
				xTris[adjTri] = true;
				rList.push_back(adjTri);  // recurse
			}
		}
		for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
			if (xTris[i])
				_mt->setTriangleMaterial(i, 10);
		}
		undermineSkin();
		excise(triangle);
	}
	else {
		for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
			if (xTris[i]) {
				int *tr = _mt->triangleVertices(i);
				for (int j = 0; j < 3; ++j) {
					_vbt->setVertexTetrahedron(tr[j], -1);  // physics marked deleted
					*_mt->vertexFaceTriangle(tr[j]) = 0x80000000;  // surface vertex marked deleted
					_deepBed.erase(tr[j]);  // remove from deepBed
					_mt->setTriangleMaterial(i, -1);  // surface triangle marked deleted
				}
			}
		}
		// get all triangles on the edge of a skin cut.
		_inExCisionTriangles.clear();
		for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
			if (_mt->triangleMaterial(i) != 2)
				continue;
			unsigned int *adjs = _mt->triAdjs(i);
			for (int j = 0; j < 3; ++j) {
				assert(adjs[j] != 3);
				if (_mt->triangleMaterial(adjs[j] >> 2) == 3) {
					_inExCisionTriangles.push_back(i);
					break;
				}
			}
		}
	}
}

bool skinCutUndermineTets::testIncisionsDeepBed() {  // Looks for intersections of the deep bed with the deep surface of the object.  For debugging.
	_mt->findAdjacentTriangles(true);
	std::map<int, int> fatTop, fatBot, fatBed, cutFatBot, extraFat;
	std::multimap<int, int> topEdge, bedEdge, bottomEdge;
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
		int t, mat = _mt->triangleMaterial(i);
		unsigned int* adjs = _mt->triAdjs(i);
		if (mat == 2) {
			for (int j = 0; j < 3; ++j) {
				t = adjs[j] >> 2;
				if (_mt->triangleMaterial(t) == 3)
					topEdge.insert(std::make_pair(i, t));
			}
		}
		else if (mat == 4) {
			for (int j = 0; j < 3; ++j) {
				t = adjs[j] >> 2;
				if (_mt->triangleMaterial(t) == 3)
					bottomEdge.insert(std::make_pair(i, t));
			}
		}
		else if (mat == 5) {
			for (int j = 0; j < 3; ++j) {
				t = adjs[j] >> 2;
				if (_mt->triangleMaterial(t) == 3)
					bedEdge.insert(std::make_pair(i, t));
			}
		}
		else if (mat == 3) {
			t = adjs[0] >> 2;
			if (_mt->triangleMaterial(t) == 2)
				fatTop.insert(std::make_pair(i, t));
			else if (_mt->triangleMaterial(t) == 4)
				fatBot.insert(std::make_pair(i, t));
			else if (_mt->triangleMaterial(t) == 5)
				fatBed.insert(std::make_pair(i, t));
			else if (_mt->triangleMaterial(t) == 3)
				cutFatBot.insert(std::make_pair(i, t));
			else
				extraFat.insert(std::make_pair(i, t));
		}
		else
			;
	}
	bool success = true;
	if (!extraFat.empty()) {
		std::cout << "Not all fat triangles accounted for. They are (with adjacent 0 edge triangle):\n";
		for (auto& ef : extraFat)
			std::cout << ef.first << " adj[0] tri-" << ef.second << "\n";
	}
	for (auto& fb : cutFatBot) {
		if (cutFatBot.find(fb.second) == cutFatBot.end()) {
			std::cout << "fat bottom edge not undermined triangle " << fb.first << " doesn't have a matching adjacent one-\n";
			success = false;
		}
	}
	for (auto& te : topEdge) {
		if (fatTop.find(te.second) == fatTop.end()) {
			std::cout << "topEdge material 2 triangle " << te.first << " doesn't have a matching fat edge triangle\n";
			success = false;
		}
	}
	for (auto& ft : fatTop) {
		if (topEdge.find(ft.second) == topEdge.end()){
			std::cout << "fatTop triangle " << ft.first << " doesn't have a matching material 2 triangle\n";
			success = false;
		}
	}
	for (auto& be : bedEdge) {
		if (fatBed.find(be.second) == fatBed.end()) {
			std::cout << "bedEdge material 5 triangle " << be.first << " doesn't have a matching fat bed triangle\n";
			success = false;
		}
	}
	for (auto& fb : fatBed) {
		if (bedEdge.find(fb.second) == bedEdge.end()) {
			std::cout << "fatBed triangle " << fb.first << " doesn't have a matching material 5 triangle\n";
			success = false;
		}
	}

	for (auto& be : bottomEdge) {
		if (fatBot.find(be.second) == fatBot.end()) {
			std::cout << "bottomEdge material 4 triangle " << be.first << " doesn't have a matching fat bottom triangle\n";
			success = false;
		}
	}
	for (auto& fb : fatBot) {
		if (bottomEdge.find(fb.second) == bottomEdge.end()) {
			std::cout << "fatBottom triangle " << fb.first << " doesn't have a matching material 4 triangle\n";
			success = false;
		}
	}

	std::set<int> bedEdgeV, botEdgeV;
	for (auto &fb : fatBed) {
		int* tr = _mt->triangleVertices(fb.first);
		bedEdgeV.insert(tr[0]);
		bedEdgeV.insert(tr[1]);
	}
	for (auto& fb : fatBot) {
		int* tr = _mt->triangleVertices(fb.first);
		botEdgeV.insert(tr[0]);
		botEdgeV.insert(tr[1]);
	}
	std::vector<int> iEdgeVerts;
	iEdgeVerts.reserve(4);
	for (auto v : botEdgeV) {
		if (bedEdgeV.find(v) != bedEdgeV.end())
			iEdgeVerts.push_back(v);
	}
	if (iEdgeVerts.size() > 2) {
		for(auto iv : iEdgeVerts)
			std::cout << "Vertex " << iv << " in both bed edge and bottom edge vertex lists\n";  // should be two points at each flap end
	}

	std::vector<int> ut;
	ut.assign(_mt->numberOfTextures(), -1);
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
		int* tr = _mt->triangleVertices(i);
		int* ttx = _mt->triangleTextures(i);
		for (int j = 0; j < 3; ++j) {
			if (ut[ttx[j]] < 0)
				ut[ttx[j]] = tr[j];
			else {
				if (ut[ttx[j]] != tr[j])
					std::cout << "After undermineSkin() for triangle " << i << " material " << _mt->triangleMaterial(i) << " tex coord " << ttx[j] << " has 2 positions " << tr[j] << " and " << ut[ttx[j]] << "\n";
			}
		}
	}

	// now test for any possible deep bed intersections with underside of object.  In material coords.
	auto edgeObjectIntersection = [&](Vec3f& v0, Vec3f& v1, int excludeVertex, int triangle) {
		boundingBox<float> be, bt;
		be.Empty_Box();
		be.Enlarge_To_Include_Point(v0.xyz);
		be.Enlarge_To_Include_Point(v1.xyz);
		float bestI = FLT_MAX;
		int bestTri = -1;
		for (int m, n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
			int mat = _mt->triangleMaterial(i);
			if (mat < 0 || mat == 3 || mat == 4 || mat == 5)
				continue;
			if (triangle > -1 && (mat == 2 || mat == 3))
				continue;
			int* t2 = _mt->triangleVertices(i);
			bt.Empty_Box();
			Vec3f tv[3];
			for (m = 0; m < 3; ++m) {
				if (t2[m] == excludeVertex)
					break;
				_vbt->vertexGridLocus(t2[m], tv[m]);
				bt.Enlarge_To_Include_Point(tv[m].xyz);
			}
			if (m < 3 || !be.Intersection(bt))
				continue;
			tv[1] -= tv[0];
			tv[2] -= tv[0];
			Mat3x3f M;
			M.Initialize_With_Column_Vectors(tv[1], tv[2], v0 - v1);
			Vec3f R = M.Robust_Solve_Linear_System(v0 - tv[0]);
			if (R.X <= 0.0f || R.Y <= 0.0f || R.Z <= 0.0f || R.X > 1.0f || R.Y >= 1.0f || R.Z >= 1.0f || R.X + R.Y >= 1.0f)
				continue;
			if (R.Z < bestI) {
				bestI = R.Z;
				bestTri = i;
			}
			//				std::pair<int, int> ip;
			//				ip.first = k;
			//				ip.second = i;
			//				triangleIntersectPairs.push_back(ip);

		}
		if (bestI < 10000.0f) {
			if (triangle < 0)
				std::cout << "deep bed top vertex " << excludeVertex << " edge intersect with triangle " << bestTri << " at param " << bestI << "\n";
			else
				std::cout << "deep bed triangle " << triangle << " edge intersect with triangle " << bestTri << " at param " << bestI << "\n";
		}
	};
	if (!success)
		std::cout << "Failed testIncisionsUndermines() veracity test.\n";
	else
		std::cout << "Passed testIncisionsUndermines() veracity test.\n";
	return success;
}

bool skinCutUndermineTets::triangleUndermined(int triangle) {
	int dv[3],  *tr = _mt->triangleVertices(triangle);
	for (int i = 0; i < 3; ++i) {
		auto dbit = _deepBed.find(tr[i]);
		if (dbit == _deepBed.end() || dbit->second.deepMtVertex < 0)  // not material 2 would be found here
			return false;
		dv[i] = dbit->second.deepMtVertex;
	}
	std::vector<materialTriangles::neighborNode> nei;
	_mt->getNeighbors(dv[0], nei);
	auto nit = nei.begin();
	while (nit != nei.end()) {
		if (nit->vertex == dv[2])
			break;
		++nit;
	}
	if ((nit == nei.end()))
		return false;
	++nit;
	if ((nit == nei.end()))
		nit = nei.begin();
	if (nit->vertex == dv[1])
		return true;
	return false;
}

skinCutUndermineTets::skinCutUndermineTets()
{
	_inExCisionTriangles.clear();
	_periostealCutEdgeTriangles.clear();
	_prevUndermineTriangle = -1;
	_collisionSpokes.clear();
	_deepSpokesNow.clear();
	_mt = nullptr;
}

skinCutUndermineTets::~skinCutUndermineTets()
{
}
