// File: GLmatrices.cpp
// Author: Court Cutting
// Date: February 3, 2012
// Purpose: Rudimentary matrix handling for openGL purposes.
//	Much of the code is modified from Richard Wright's excellent
//	OpenGL Superbible fifth edition who should get most of the credit for this.

#include <assert.h>
#include "GLmatrices.h"

void GLmatrices::getDragVector(float dScreenX, float dScreenY, float (&position)[3], float (&dragVector)[3])
{
	float fh,fw,height = _mFR[2]*position[0] + _mFR[6]*position[1] + _mFR[10]*position[2] + _mFR[14];	
	height *= -2.0f*tanf(_angleRadians*0.5f);
	fw = dScreenX*height*_screenAspect;
	fh = dScreenY*height;
	dragVector[0] = fw*_mFR[0] + fh*_mFR[1];
	dragVector[1] = fw*_mFR[4] + fh*_mFR[5];
	dragVector[2] = fw*_mFR[8] + fh*_mFR[9];
}

void GLmatrices::setPan(float frX, float frY)
{
	float fh,fw,height = -_zCenter*2.0f*tanf(_angleRadians*0.5f);
	fw = frX*height*_screenAspect;
	fh = frY*height;
	float oldCenter[3];
	oldCenter[0]=_center[0]; oldCenter[1]=_center[1]; oldCenter[2]=_center[2];
	_center[0] += fw*_mFR[0] + fh*_mFR[1];
	_center[1] += fw*_mFR[4] + fh*_mFR[5];
	_center[2] += fw*_mFR[8] + fh*_mFR[9];
	_mFR[12] = -_center[0]*_mFR[0]-_center[1]*_mFR[4]-_center[2]*_mFR[8];
	_mFR[13] = -_center[0]*_mFR[1]-_center[1]*_mFR[5]-_center[2]*_mFR[9];
	_mFR[14] = _zCenter-_center[0]*_mFR[2]-_center[1]*_mFR[6]-_center[2]*_mFR[10];
	// also radius just changed
	float d[3],l;
	d[0]=_center[0]-oldCenter[0]; d[1]=_center[1]-oldCenter[1]; d[2]=_center[2]-oldCenter[2];
	l= (float)sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
	_radius += l*0.5f;
	float zfar,zmin = -_zCenter-_radius*1.5f;
	zfar = zmin + _radius*3.0f;
	if(zmin<_radius*0.01f)
		zmin=_radius*0.01f;
	MakePerspectiveMatrix(_angleRadians,_screenAspect,zmin,zfar);
}

void GLmatrices::changeZoom(float dZoom)
{
	_zCenter += dZoom*_radius;
	if (_zCenter > -0.01f)
		_zCenter = -0.01f;
	for(int i=0; i<16; ++i)
		_mFR[i] = _rotNow[i];
	_mFR[12] = -_center[0]*_mFR[0]-_center[1]*_mFR[4]-_center[2]*_mFR[8];
	_mFR[13] = -_center[0]*_mFR[1]-_center[1]*_mFR[5]-_center[2]*_mFR[9];
	_mFR[14] = _zCenter-_center[0]*_mFR[2]-_center[1]*_mFR[6]-_center[2]*_mFR[10];
	float zfar,zmin = -_zCenter-_radius*2.0f;
	zfar = zmin + _radius*4.0f;
	if(zmin<_radius*0.01f)
		zmin=_radius*0.01f;
	MakePerspectiveMatrix(_angleRadians,_screenAspect,zmin,zfar);
}

void GLmatrices::setFrameAndRotation(GLfloat *rotMatrix)
{	// remember setProjectionMatrix() must be called first or will get undefined results
	for(int i=0; i<16; ++i)
		_rotNow[i]=_mFR[i] = rotMatrix[i];
	_mFR[12] = -_center[0]*_mFR[0]-_center[1]*_mFR[4]-_center[2]*_mFR[8];
	_mFR[13] = -_center[0]*_mFR[1]-_center[1]*_mFR[5]-_center[2]*_mFR[9];
	_mFR[14] = _zCenter-_center[0]*_mFR[2]-_center[1]*_mFR[6]-_center[2]*_mFR[10];
}

void GLmatrices::setView(float angleRadians, float screenAspect)
{	// assumes center and radius have been set
	_screenAspect=screenAspect;
	if(angleRadians<0.1f)
		angleRadians=0.1f;
	_angleRadians=angleRadians;
}

void GLmatrices::resetPerspective()
{	// assumes center and radius have been set
	_zmin = _radius/tanf(_angleRadians * 0.5f);
	_zCenter = -(_zmin +_radius*2.0f);
	MakePerspectiveMatrix(_angleRadians,_screenAspect,_zmin,_zmin+_radius*4.0f);
}

////////////////////////////////////////////////////////////////////////////////////////////
// Create a projection matrix
// Similiar to the old gluPerspective... fov is in radians btw...
void GLmatrices::MakePerspectiveMatrix(float fFov, float fAspect, float zMin, float zMax)
{
	loadIdentity4x4(&_mProj[0]);
    float yMax = zMin * tanf(fFov * 0.5f);
    float yMin = -yMax;
	float xMin = yMin * fAspect;
    float xMax = -xMin; 
	_mProj[0] = (2.0f * zMin) / (xMax - xMin);
	_mProj[5] = (2.0f * zMin) / (yMax - yMin);
	_mProj[8] = (xMax + xMin) / (xMax - xMin);
	_mProj[9] = (yMax + yMin) / (yMax - yMin);
	_mProj[10] = -((zMax + zMin) / (zMax - zMin));
	_mProj[11] = -1.0f;
	_mProj[14] = -((2.0f * (zMax*zMin))/(zMax - zMin));
	_mProj[15] = 0.0f;
}

GLmatrices::GLmatrices()
{
	_zCenter=-4.0f;
	loadIdentity4x4(_mProj);
	loadIdentity4x4(_mFR);
	loadIdentity4x4(_rotNow);
}

GLmatrices::~GLmatrices()
{
}

