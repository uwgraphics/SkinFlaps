// File: hooks.cpp
// Author: Court Cutting
// Date: 7/17/19
// Purpose: Class for handling tissue hooks and associated graphics using pd constraints

#include <vector>
#include <math.h>
#include "materialTriangles.h"
#include "GLmatrices.h"
#include "Vec3f.h"
#include "vnBccTetrahedra.h"
#include <assert.h>
#ifdef linux
#include <stdio.h>
#endif
#include "hooks.h"

float hooks::_hookSize = 2.5f;
float hooks::_springConstant = 40.0f;
GLfloat hooks::_selectedColor[] = {1.0f, 1.0f, 0.0f, 1.0f};
GLfloat hooks::_unselectedColor[] = {0.043f, 0.898f, 0.102f, 1.0f};

void hooks::deleteHook(int hookNumber)
{
	HOOKMAP::iterator hit = _hooks.find(hookNumber);
	if(hit==_hooks.end())
		return;
	if (hit->second._tri->triangleMaterial(hit->second.triangle) > -1  && hit->second._constraintId > -1){
#ifndef NO_PHYSICS
		_ptp->deleteHook(hit->second._constraintId);
		_ptp->initializePhysics();
#endif
	}
	_shapes->deleteShape(hit->second.getShape());
	_hooks.erase(hit);
}

void hooks::selectHook(int hookNumber)
{
	HOOKMAP::iterator hit;
	for(hit=_hooks.begin(); hit!= _hooks.end(); ++hit) {
		if(hit->first==hookNumber)	{
			hit->second.getShape()->setColor(_selectedColor);
			hit->second._selected = true;
		}
		else	{
			hit->second.getShape()->setColor(_unselectedColor);
			hit->second._selected = false;	}
	}
}

bool hooks::getSelectPosition(unsigned int hookNumber, float(&selectPos)[3])
{
	HOOKMAP::iterator hit = _hooks.find(hookNumber);
	if (hit == _hooks.end())
		return false;
	GLfloat *mvm = hit->second._shape->getModelViewMatrix();
	selectPos[0] = hit->second._selectPosition[0];
	selectPos[1] = hit->second._selectPosition[1];
	selectPos[2] = hit->second._selectPosition[2];
	return true;
}

bool hooks::getHookPosition(unsigned int hookNumber, float(&hookPos)[3])
{
	HOOKMAP::iterator hit = _hooks.find(hookNumber);
	if(hit==_hooks.end())
		return false;
	GLfloat *mvm = hit->second._shape->getModelViewMatrix();
	hookPos[0] = mvm[12];
	hookPos[1] = mvm[13];
	hookPos[2] = mvm[14];
	return true;
}

bool hooks::setHookPosition(unsigned int hookNumber, float(&hookPos)[3])
{
	HOOKMAP::iterator hit = _hooks.find(hookNumber);
	if(hit==_hooks.end())
		return false;
	hit->second.xyz = (Vec3f)hookPos;
#ifndef NO_PHYSICS
	if(hit->second._constraintId > -1)
		_ptp->moveHook(hit->second._constraintId, reinterpret_cast< const std::array<float, 3>(&) >(hit->second.xyz));
	else  // physics not activated yet
		throw(std::logic_error("Attempting to move a hook without physics activation.\n"));
#endif
	GLfloat *mvm = hit->second._shape->getModelViewMatrix();
	mvm[12] = hookPos[0];
	mvm[13] = hookPos[1];
	mvm[14] = hookPos[2];
	return true;
}

int hooks::addHook(materialTriangles *tri, int triangle, float(&uv)[2], bool tiny)
{
	std::pair<HOOKMAP::iterator, bool> hpr;
	hpr = _hooks.insert(std::make_pair(_hookNow, hookConstraint()));
	char name[6];
	sprintf(name,"H_%d",_hookNow);
	std::shared_ptr<sceneNode> sh;
	if(tiny)
		sh = _shapes->addShape(sceneNode::nodeType::SPHERE, name);
	else
		sh = _shapes->addShape(sceneNode::nodeType::CONE, name);
	hpr.first->second.setShape(sh);
	hpr.first->second._strong = tiny;
	++_hookNow;
	hpr.first->second.triangle = triangle;
	hpr.first->second.uv[0] = uv[0];
	hpr.first->second.uv[1] = uv[1];
	hpr.first->second._tri = tri;
	GLfloat *om = sh->getModelViewMatrix();
	loadIdentity4x4(om);
	if(tiny)
		scaleMatrix4x4(om,_hookSize*0.1f,_hookSize*0.1f,_hookSize*0.1f);
	else
		scaleMatrix4x4(om, _hookSize, _hookSize, _hookSize);
	hpr.first->second._selected = true;
	float angle, xyz[3];	// ,*p = hpr.first->second._initialPosition;
	Vec3f n;
	tri->getBarycentricPosition(triangle,uv,xyz);
	hpr.first->second.xyz = (Vec3f)xyz;
	tri->getTriangleNormal(triangle, n, true);
	Vec3f vz(0.0f, 0.0f, 1.0f);
	angle = acos(n*vz);
	Vec3f axis = vz^n;
	if (fabs(axis[0]) < 1e-5 && fabs(axis[1]) < 1e-5 && fabs(axis[2]) < 1e-5) {
		assert(fabs(angle) < 1e-5f || fabs(angle) > 3.1414f);
		axis[1] = 1.0f;  // perpendicular to start of cone
	}
	axisAngleRotateMatrix4x4(om, axis.xyz, angle);
	translateMatrix4x4(om,xyz[0],xyz[1],xyz[2]);
	Vec3f gridLocus, bw;
	if (_vnt->getMaterialTriangles() != nullptr && _ptp->solverInitialized()) {  // COURT - won't need second condition
		int tetIdx = _vnt->parametricTriangleTet(triangle, uv, gridLocus);
		if (tetIdx < 0){
			--_hookNow;
			deleteHook(_hookNow);
			return -1;
		}
		_vnt->gridLocusToBarycentricWeight(gridLocus, _vnt->tetCentroid(tetIdx), bw);
#ifndef NO_PHYSICS
		hpr.first->second._constraintId = _ptp->addHook(tetIdx, reinterpret_cast<const std::array<float, 3>&>(bw), reinterpret_cast<const std::array<float, 3>&>(xyz), tiny);  // COURT - slow. ?put in thread?
#else
		hpr.first->second._constraintId = -1;  // signal that this is a dummy hook that needs a constraint later
#endif
	}
	else
		hpr.first->second._constraintId = -1;  // signal that this is a dummy hook that needs a constraint later
	return _hookNow - 1;
}

bool hooks::updateHookPhysics(){
	auto hit = _hooks.begin();
	while (hit != _hooks.end()){
		if (hit->second._tri->triangleMaterial(hit->second.triangle) < 0){
			int delHook = hit->first;
			++hit;
			deleteHook(delHook);
			continue;
		}
		Vec3f gridLocus, bw;
		int tetIdx = _vnt->parametricTriangleTet(hit->second.triangle, hit->second.uv, gridLocus);
		if (tetIdx < 0){
			--_hookNow;
			++hit;
			deleteHook(_hookNow);
			continue;
		}
		_vnt->gridLocusToBarycentricWeight(gridLocus, _vnt->tetCentroid(tetIdx), bw);
#ifndef NO_PHYSICS
		hit->second._constraintId = _ptp->addHook(tetIdx, reinterpret_cast<const std::array<float, 3>&>(bw), reinterpret_cast<const std::array<float, 3>&>(hit->second.xyz), hit->second._strong);
#endif
		++hit;
	}
	// don't call _ptp->reInitializePhysics() here as done in bccTetScene::updateOldPhysicsLattice() where this routine is called
	return true;
}

hooks::hooks()
{
	_hookNow=0;
	_selectedHook=-1;
}

hooks::~hooks()
{
}
