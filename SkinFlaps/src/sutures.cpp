// File: sutures.cpp
// Author: Court Cutting
// Date: 2/20/12
// Update: 6/10/2019 for bcc tetrahedra
// Purpose: Class for handling sutures and associated graphics

#include <vector>
#include <math.h>
#include <Vec3f.h>
#include <Vec3d.h>
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"
#include "deepCut.h"
#include "GLmatrices.h"
#include <assert.h>
#ifdef linux
#include <stdio.h>
#endif
#include "sutures.h"

float sutures::_sutureSpanGap = 0.03f;
float sutures::_sutureSize=1.0f;
GLfloat sutures::_selectedColor[] = {1.0f, 1.0f, 0.0f, 1.0f};
GLfloat sutures::_unselectedColor[]={0.4f,0.537f,0.984f,1.0f};
GLfloat sutures::_userColor[] = { 0.03f,0.03f,0.99f,1.0f };

int sutures::previousUserSuture(int sutureNumber)
{
	auto s = _sutures.find(sutureNumber);
	if (s == _sutures.end() || s == _sutures.begin())
		return -1;
	do {
		--s;
	}
	while (s->second._type > 1);
	return s->first;
}

int sutures::firstVertexMaterial(int sutureNumber)
{
	auto sut = _sutures.find(sutureNumber);
	if (sut == _sutures.end())
		return -1;
	else
		return sut->second._tri->triangleMaterial(sut->second._tris[0]);
}

void sutures::selectSuture(int sutureNumber)
{
	SUTUREMAP::iterator sit;
	for(sit=_sutures.begin(); sit!= _sutures.end(); ++sit) {
		if(sit->first==sutureNumber)	{
			sit->second.getCylinderShape()->setColor(_selectedColor);
			sit->second.getSphereShape()->setColor(_selectedColor);
			sit->second._selected = true;	}
		else	{
			if (sit->second._type < 2) {
				sit->second.getCylinderShape()->setColor(_userColor);
				sit->second.getSphereShape()->setColor(_userColor);
			}
			else {
				sit->second.getCylinderShape()->setColor(_unselectedColor);
				sit->second.getSphereShape()->setColor(_unselectedColor);
			}
			sit->second._selected = false;
		}
	}
}

void sutures::updateSutureGraphics()
{
	SUTUREMAP::iterator sit=_sutures.begin();
	while(!_sutures.empty() && sit!=_sutures.end())	{
		sutureTets *sut = &(sit->second);
		if (sut->_tris[1] < 0){
			++sit;
			continue;
		}
		GLfloat *mm=sut->getSphereShape()->getModelViewMatrix();
		loadIdentity4x4(mm);
		sut->_selected = true;
		GLfloat *v,*v2;
		float len;
		Vec3f v0,v1,p,dv;
		long *tri;
		tri = sut->_tri->triangleVertices(sut->_tris[0]);
		v = sut->_tri->vertexCoordinate(tri[sut->_edges[0]]);
		v2 = sut->_tri->vertexCoordinate(tri[(sut->_edges[0]+1)%3]);
		for(int i=0; i<3; ++i)
			v0._v[i]=v2[i]*sut->_params[0] + v[i]*(1.0f-sut->_params[0]);
		tri = sut->_tri->triangleVertices(sut->_tris[1]);
		v = sut->_tri->vertexCoordinate(tri[sut->_edges[1]]);
		v2 = sut->_tri->vertexCoordinate(tri[(sut->_edges[1]+1)%3]);
		for(int i=0; i<3; ++i)
			v1._v[i]=v2[i]*sut->_params[1] + v[i]*(1.0f-sut->_params[1]);
		GLfloat sSize = _sutureSize;
		if (sut->_type < 2)
			sSize *= 2.0f;
		scaleMatrix4x4(mm, sSize, sSize, sSize);
		sut->_v1[0]=v1._v[0]; sut->_v1[1]=v1._v[1]; sut->_v1[2]=v1._v[2];
		p = (v0+v1)*0.5f;
		translateMatrix4x4(mm,p._v[0],p._v[1],p._v[2]);
		mm = sut->getCylinderShape()->getModelViewMatrix();
		loadIdentity4x4(mm);
		dv = v0-v1;
		len = dv.length();
		scaleMatrix4x4(mm,sSize*0.5f,sSize*0.5f,len*0.5f);
		if(len<1e-16f)
			;
		else	{
			float qAxis[3]={-dv._v[1]/len,dv._v[0]/len,0.0f};	// normalized dvXj
			len = acos(dv._v[2]/len);
			axisAngleRotateMatrix4x4(mm,qAxis,len);
		}
		translateMatrix4x4(mm,p._v[0],p._v[1],p._v[2]);
		++sit;
	}
}

float* sutures::getSecondVertexPosition(int sutureNumber)
{	// updates entire suture graphics
	SUTUREMAP::iterator sit = _sutures.find(sutureNumber);
	if (sit == _sutures.end())
		return nullptr;
	return sit->second._v1;
}

void sutures::setSecondVertexPosition(int sutureNumber, float *position)
{	// updates entire suture graphics
	SUTUREMAP::iterator sit=_sutures.find(sutureNumber);
	if(sit==_sutures.end())
		return;
	sutureTets *sut = &(sit->second);
	sut->_selected = true;
	GLfloat *v,*v2;
	float p[3],dv[3],len;
	Vec3f v0;
	long *tri = sut->_tri->triangleVertices(sut->_tris[0]);
	v = sut->_tri->vertexCoordinate(tri[sut->_edges[0]]);
	v2 = sut->_tri->vertexCoordinate(tri[(sut->_edges[0]+1)%3]);
	for(int i=0; i<3; ++i)
		v0._v[i]=v2[i]*sut->_params[0] + v[i]*(1.0f-sut->_params[0]);
	GLfloat sSize = _sutureSize;
	if (sut->_type < 2)
		sSize *= 2.0f;
	GLfloat* mm = sut->getSphereShape()->getModelViewMatrix();
	loadIdentity4x4(mm);
	scaleMatrix4x4(mm,sSize, sSize, sSize);
	sut->_v1[0]=position[0]; sut->_v1[1]=position[1]; sut->_v1[2]=position[2];
	p[0]=(v0._v[0]+position[0])*0.5f; p[1]=(v0._v[1]+position[1])*0.5f; p[2]=(v0._v[2]+position[2])*0.5f;
	translateMatrix4x4(mm,p[0],p[1],p[2]);
	mm = sut->getCylinderShape()->getModelViewMatrix();
	loadIdentity4x4(mm);
	dv[0]=v0._v[0]-position[0]; dv[1]=v0._v[1]-position[1]; dv[2]=v0._v[2]-position[2];
	len = dv[0]*dv[0] + dv[1]*dv[1] + dv[2]*dv[2];
	len = sqrt(len);
	scaleMatrix4x4(mm, sSize*0.5f, sSize*0.5f, len*0.5f);
	if(len<1e-16f)
		;
	else	{
		float qAxis[3]={-dv[1]/len,dv[0]/len,0.0f};	// normalized dvXj
		len = acos(dv[2]/len);
		axisAngleRotateMatrix4x4(mm,qAxis,len);
	}
	translateMatrix4x4(mm,p[0],p[1],p[2]);
}

int sutures::setSecondEdge(int sutureNumber, materialTriangles *tri, int triangle, int edge, float param, bool initPhysics)
{
	SUTUREMAP::iterator sit=_sutures.find(sutureNumber);
	if(sit==_sutures.end())
		return 1;
	sit->second._tris[1] = triangle;
	if (sit->second._tris[0] == sit->second._tris[1]){
		--_sutureNow;
		assert(_sutureNow == sutureNumber);
		deleteSuture(sutureNumber);
		return 1;
	}
	sit->second._edges[1] = edge;
	sit->second._params[1] = param;
	if (tri == sit->second._tri)	{	// true for now. May change later.
		Vec3f bw, gridLocus;
		for (int i = 0; i < 2; ++i){
			sit->second._tetIdx[i] = _dc->parametricMTedgeTet(sit->second._tris[i], sit->second._edges[i], sit->second._params[i], gridLocus);
			if (sit->second._tetIdx[i] < 0){
				--_sutureNow;
				assert(_sutureNow == sutureNumber);
				deleteSuture(sutureNumber);
				return 2;
			}
			_vbt->gridLocusToBarycentricWeight(gridLocus, *_vbt->tetCentroid(sit->second._tetIdx[i]), sit->second._baryWeights[i]);
		}
		if (_ptp->solverInitialized()) {
			sit->second._constraintId = _ptp->addSuture(sit->second._tetIdx, reinterpret_cast<const std::array<float, 3>(&)[2]>(sit->second._baryWeights));
			if (initPhysics)
				_ptp->initializePhysics();
		}
		else
			sit->second._constraintId = -1;  // stub until forces applied
		return 0;
	}
	return 3;
}

void sutures::updateSuturePhysics(){  // call after a topology change
	auto sit = _sutures.begin();
	while (sit != _sutures.end()){
		unsigned int sutNow = sit->first;
		Vec3f bw, gridLocus;
		int i;
		for (i = 0; i < 2; ++i){
			sit->second._tetIdx[i] = _dc->parametricMTedgeTet(sit->second._tris[i], sit->second._edges[i], sit->second._params[i], gridLocus);
			_vbt->gridLocusToBarycentricWeight(gridLocus, *_vbt->tetCentroid(sit->second._tetIdx[i]), bw);
			if (sit->second._tetIdx[i] < 0){
				++sit;
				deleteSuture(sutNow);
				break;
			}
			sit->second._baryWeights[i].set(bw._v);  // COURT - old bw was likely fine. Only tetIdx can change
		}
		if (i < 2)
			sit = _sutures.erase(sit);
		else{
			sit->second._constraintId = _ptp->addSuture(sit->second._tetIdx, reinterpret_cast<const std::array<float, 3>(&)[2]>(sit->second._baryWeights));
			++sit;
		}
	}
}

bool sutures::getEdgeAttachment(const int sutureNumber, const bool firstPoint, materialTriangles *(&tri), int &triangle, int &edge, float &param)
{
	SUTUREMAP::iterator sit=_sutures.find(sutureNumber);
	if(sit==_sutures.end())	{
		tri=NULL;
		triangle=-1;
		edge=-1;
		param=-1.0f;
		return false;
	}
	tri=sit->second._tri;
	if (firstPoint)	{
		triangle = sit->second._tris[0];
		edge = sit->second._edges[0];
		param = sit->second._params[0];	}
	else	{
		triangle = sit->second._tris[1];
		edge = sit->second._edges[1];
		param = sit->second._params[1];	}
	return true;
}

int sutures::addUserSuture(materialTriangles *tri, int triangle0, int edge0, float param0)
{
	int sutNum = addSuture(tri, triangle0, edge0, param0);
	_userSutures.insert(std::make_pair(_userSutureNext, sutNum));
	++_userSutureNext;
	return sutNum;
}

int sutures::addSuture(materialTriangles *tri, int triangle0, int edge0, float param0)
{
	std::pair<SUTUREMAP::iterator,bool> hpr;
	hpr = _sutures.insert(std::make_pair(_sutureNow, sutureTets()));
	char name[6];
	sprintf(name,"S_%d",_sutureNow);
	 std::shared_ptr<sceneNode> sh=_shapes->addShape(sceneNode::nodeType::SPHERE,name);
	hpr.first->second.setSphereShape(sh);
	sh->setColor(_selectedColor);
	++_sutureNow;
	hpr.first->second._tris[0] = triangle0;
	hpr.first->second._edges[0] = edge0;
	hpr.first->second._params[0] = param0;
	hpr.first->second._tris[1] = -1;
	hpr.first->second._tri = tri;
	GLfloat *mm = sh->getModelViewMatrix();
	loadIdentity4x4(mm);
	hpr.first->second._selected = true;
	long *t = tri->triangleVertices(triangle0);
	GLfloat *v0,*v1,v[3];
	v0 = tri->vertexCoordinate(t[edge0]);
	v1 = tri->vertexCoordinate(t[(edge0+1)%3]);
	for(int i=0; i<3; ++i)
		v[i] = v1[i]*param0 + v0[i]*(1.0f-param0);
	GLfloat sSize = _sutureSize;
	if (hpr.first->second._type < 2)
		sSize *= 2.0f;
	scaleMatrix4x4(mm, sSize, sSize, sSize);
	translateMatrix4x4(mm,v[0],v[1],v[2]);
	sh =_shapes->addShape(sceneNode::nodeType::CYLINDER,name);
	hpr.first->second.setCylinderShape(sh);
	mm = sh->getModelViewMatrix();
	loadIdentity4x4(mm);
	sh->setColor(_selectedColor);
	scaleMatrix4x4(mm,sSize*0.5f,sSize*0.5f,sSize*0.5f);
	translateMatrix4x4(mm,v[0],v[1],v[2]);
	return _sutureNow-1;
}

int sutures::deleteSuture(int sutureNumber)
{
	int ret = 0;
	SUTUREMAP::iterator sit2, sit=_sutures.find(sutureNumber);
	if(sit==_sutures.end())
		return -1;  // error return
	bool autoSuture = baseToUserSutureNumber(sutureNumber) < 0;
	auto delSut = [&]() {
		if (sit2->second._constraintId > -1)
			_ptp->deleteSuture(sit2->second._constraintId);
		_shapes->deleteShape(sit2->second.getSphereShape());
		_shapes->deleteShape(sit2->second.getCylinderShape());
		sit2 = _sutures.erase(sit2);
	};
	sit2 = sit;
	++sit2;
	while (sit2 != _sutures.end() && sit2->second._type > 1)
		delSut();
	if (!autoSuture && sit2 != _sutures.end() && sit2->second._type == 1) {
		++sit2;
		while (sit2 != _sutures.end() && sit2->second._type > 1)
			delSut();
	}
	sit2 = sit;
	if (autoSuture) {
		while (sit2->second._type > 1) {
			delSut();
			--sit2;
		}
		assert(sit2->second._type == 1);
		sit2->second._type = 0;
		ret = sit2->first;
	}
	else {
		auto usit = _userSutures.find(baseToUserSutureNumber(sit->first));
		assert(usit != _userSutures.end());
		_userSutures.erase(usit);
		delSut();
	}
	_ptp->initializePhysics();
	return ret;
}

void sutures::laySutureLine(int suture2)
{  // creates a line of sutures between 2 sequential user input sutures
	if (suture2 < 1)
		return;
	auto sut1 = _sutures.find(suture2);
	auto sut0 = sut1;
	do{
		--sut0;
	} while (sut0 != _sutures.begin() && sut0->second._type > 1);
	if (sut1 == _sutures.end() || (sut0 == _sutures.begin() && sut0->second._type > 1) || sut0->second._tri != sut1->second._tri)  // second condition could happen after deleting a suture
		return;
	auto s0 = sut0->second;
	auto s1 = sut1->second;
	// at present can only lay suture line between same materials
	materialTriangles *mt = s0._tri;
	Vec3f e0, e1;
	std::vector<float> edgeLengths0, edgeLengths1;
	std::vector<unsigned long> triEdges0, triEdges1;
	unsigned long te, mat;
	float len0, len1, lenMax = 15.0f;
	// due to vagaries in physics, do in material coords
	auto getIncisionSpan = [&](bool forward, int side, std::vector<unsigned long> &tev, std::vector<float> &lens, float &totalLen) {
		totalLen = 0;
		tev.clear();
		lens.clear();
		_vbt->vertexGridLocus(mt->triangleVertices(te >> 2)[((te & 3) + (forward ? 0 : 1)) % 3], e0);
		e0 *= (float)_vbt->getTetUnitSize();
		while ((te >> 2) != s1._tris[0] && (te >> 2) != s1._tris[1] && totalLen < lenMax){
			tev.push_back(te);
			_vbt->vertexGridLocus(mt->triangleVertices(te >> 2)[((te & 3) + (forward ? 1 : 0)) % 3], e1);
			e1 *= (float)_vbt->getTetUnitSize();
			if (lens.empty())
				lens.push_back((e1 - e0).length() * (forward ? 1.0f - s0._params[side] : s0._params[side]));
			else
				lens.push_back((e1 - e0).length());
			totalLen += lens.back();  // consider a bailout length
			e0 = e1;
			while (mt->triangleMaterial(te >> 2) == 2)
				te = mt->triAdjs(te >> 2)[((te & 3) + (forward ? 1 : 2)) % 3];
			te = mt->triAdjs(te >> 2)[te & 3];
		}
		tev.push_back(te);
		_vbt->vertexGridLocus(mt->triangleVertices(te >> 2)[((te & 3) + (forward ? 1 : 0)) % 3], e1);
		e1 *= (float)_vbt->getTetUnitSize();
		float par;
		if ((te >> 2) == s1._tris[0])
			par = s1._params[0];
		else
			par = s1._params[1];
		lens.push_back((e1 - e0).length() * (forward ? par : 1.0f - par));
		totalLen += lens.back();
	};
	mat = mt->triangleMaterial(s0._tris[0]);
	assert(mat != mt->triangleMaterial(mt->triAdjs(s0._tris[0])[s0._edges[0]] >> 2));
	te = (s0._tris[0] << 2) + s0._edges[0];
	getIncisionSpan(true, 0, triEdges0, edgeLengths0, len0);
	te = (s0._tris[0] << 2) + s0._edges[0];
	getIncisionSpan(false, 0, triEdges1, edgeLengths1, len1);
	bool forwardSide0;
	// always span shortest path between sutures
	if (len0 < len1){
		forwardSide0 = true;
		// other side
		te = (s0._tris[1] << 2) + s0._edges[1];
		mat = mt->triangleMaterial(s0._tris[1]);
		assert(mat != mt->triangleMaterial(mt->triAdjs(s0._tris[1])[s0._edges[1]] >> 2));
		getIncisionSpan(false, 1, triEdges1, edgeLengths1, len1);
	}
	else{
		forwardSide0 = false;
		len0 = len1;
		triEdges0 = std::move(triEdges1);
		edgeLengths0 = std::move(edgeLengths1);
		te = (s0._tris[1] << 2) + s0._edges[1];
		mat = mt->triangleMaterial(s0._tris[1]);
		assert(mat != mt->triangleMaterial(mt->triAdjs(s0._tris[1])[s0._edges[1]] >> 2));
		getIncisionSpan(true, 1, triEdges1, edgeLengths1, len1);
	}
	int subdivs = (int)((len0 + len1)* 0.5f / _sutureSpanGap), j0 = 0, j1 = 0, n0 = edgeLengths0.size(), n1 = edgeLengths1.size();
	float lp0=0.0f, lp1=0.0f, d0 = len0 / subdivs, d1 = len1 / subdivs;
	for (int i = 1; i < subdivs; ++i){
		int sutNum;
		while (j0 < n0){
			if (lp0 + edgeLengths0[j0] > d0 * i){
				float param = (d0*i - lp0) / edgeLengths0[j0];
				if (!forwardSide0)
					param = 1.0f - param;
				if (j0 < 1){
					float op = s0._params[0];
					if (forwardSide0)
						param += op - param*op;
					else
						param *= op;
				}
				if(j0 == n0 - 1){
					float op = (triEdges0[j0] >> 2 == s1._tris[0]) ? s1._params[0] : s1._params[1];
					if (forwardSide0)
						param *= op;
					else
						param += op - param*op;
				}
				sutNum = addSuture(mt, triEdges0[j0] >> 2, triEdges0[j0] & 3, param);
				_sutures[sutNum]._type = 2;  // mark as programatically created, not user created
				break;
			}
			else{
				lp0 += edgeLengths0[j0];
				++j0;
			}
		}
		while (j1 < n1){
			if (lp1 + edgeLengths1[j1] > d1 * i){
				float param = (d1*i - lp1) / edgeLengths1[j1];
				if (forwardSide0)
					param = 1.0f - param;
				if (j1 < 1){
					float op = s0._params[1];
					if (forwardSide0)
						param *= op;  // OK start but end is bad
					else
						param += op - param*op;  // OK start but end is bad
				}
				if (j1 == n1 - 1){
					float op = (triEdges1[j1] >> 2 == s1._tris[0]) ? s1._params[0] : s1._params[1];
					if (forwardSide0)
						param += op - param*op;
					else
						param *= op;
				}
				setSecondEdge(sutNum, mt, triEdges1[j1] >> 2, triEdges1[j1] & 3, param, false);  // don't init physics. Do once as group
				break;
			}
			else{
				lp1 += edgeLengths1[j1];
				++j1;
			}
		}
	}
	_ptp->initializePhysics();
}

void sutures::nearestSkinIncisionEdge(const float triUv[2], long &triangle, int &edge, float &param)
{  // problem with this routine is that, in the absence of collisions, overlap can occur so simple nearest spatial incision point won't work with suture always finding same side of skin edge.
	// This overlap impossible in material coords, but since opposing incision edges are identical must disambiguate with incision edge normal.
	// Once closest material coords edge found, do edge neighbor search in spatial coords for correction.
	materialTriangles *mt = _dc->getMaterialTriangles();
	Vec3f gridLocus, startLocus;
	long *tr = mt->triangleVertices(triangle);
	_vbt->vertexGridLocus(tr[0], gridLocus);
	startLocus = gridLocus * (1.0f - triUv[0] - triUv[1]);
	for(int i = 0; i < 2; ++i) {
		_vbt->vertexGridLocus(tr[i+1], gridLocus);
		startLocus += gridLocus * triUv[i];
	}
	unsigned long bestTe=0;
	float bestParam=0.0f, minDist = FLT_MAX;
	for (int n = mt->numberOfTriangles(), i = 0; i < n; ++i) {
		if (mt->triangleMaterial(i) != 3)
			continue;
		unsigned long adj = mt->triAdjs(i)[0];
		if (mt->triangleMaterial(adj >> 2) != 2)
			continue;
		long *tr = mt->triangleVertices(i);
		Vec3f v0, v1;
		_vbt->vertexGridLocus(tr[0], v0);
		_vbt->vertexGridLocus(tr[1], v1);
		v1 -= v0;
		float lSq, t = v1 * (startLocus - v0) / (v1*v1);
		if (t < 0.0f)
			t = 0.0f;
		else if (t > 1.0f)
			t = 1.0f;
		else
			;
		lSq = (v1 * t + v0 - startLocus).length2();
		if (lSq < minDist) {
			Vec3f v2, N;
			_vbt->vertexGridLocus(tr[2], v2);
			v2 -= v0;
			N = v1 ^ v2;
			if (N*(startLocus - v0) > 0)
				continue;
			minDist = lSq;
			bestParam = 1.0f - t;
			bestTe = adj;
		}
	}
	// COURT - started to program correction to spatial coords, but it doesn't appear to be worth the effort
	triangle = bestTe >> 2;
	edge = bestTe & 3;
	param = bestParam;
	return;
}


sutures::sutures()
{
	_sutureNow=0;
	_userSutureNext = 0;
	_sutures.clear();
}

sutures::~sutures()
{
}
