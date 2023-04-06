// File: gl3wGraphics.cpp
// Author: Court Cutting, MD
// Date: October 28, 2020
// Purpose: Interface to gl3wGraphics library. This library uses imgui as its
//	GUI library and GLFW as an openGL window.
//	Copyright 2020 - Court Cutting - All rights reserved.

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include "math3d.h"
#include "surgGraphics.h"
#include "gl3wGraphics.h"

bool gl3wGraphics::mouseWheelZoom = true;
float gl3wGraphics::mouseWheelLevel = 500.0f;

bool gl3wGraphics::addCustomSceneNode(std::shared_ptr<sceneNode>& sn, std::vector<int> &txIds, const GLchar *vertexShader, const GLchar *fragmentShader, std::vector<std::string> &attributes)
{  // assumes textures loaded already and input in txIds
	std::sort(txIds.begin(), txIds.end());
	for (int n = (int)txIds.size(), i = 0; i <n; ++i) {
		if (!_texReader.textureExists(txIds[i]))
			return false;
		sn->add2DtextureBufferNumber(_texReader.getOGLtextureNumber(txIds[i]));
	}
	GLuint progNum = 0;  // = sn->getGlslProgramNumber();
	if(!_ls.createCustomProgram(progNum,vertexShader,fragmentShader,attributes))
		return false;
	sn->setGlslProgramNumber(progNum);
	addSceneNode(sn);
	char tName[12];
	sprintf(tName,"Custom_%d",(int)progNum);
	sn->setName(tName);
	return true;
}

std::shared_ptr<sceneNode> gl3wGraphics::loadStaticObjFile(const char *filePath, std::vector<int> &textureIds, bool texturedNotColored)
{
//	_sTris.push_back(staticTriangle(texturedNotColored));
//	addSceneNode(&(_sTris.back()));
//	auto tr = std::make_shared<sceneNode>();
//	addSceneNode(tr);
//	staticTriangle *tr=&_sTris.back();
	_staticTris.setGl3wGraphics(this);
	materialTriangles mt;
	int err = mt.readObjFile(filePath);
	if (err > 0) {
		std::cout << "Couldn't read static .obj file. Error code: " << err << "\n";
		return NULL;
	}
	auto sn = _staticTris.createStaticSceneNode(&mt, textureIds);

	char tName[12];
	sprintf(tName,"Tri_%d", _staticCount++);
	sn->setName(tName);
	addSceneNode(sn);
	return sn;
}

void gl3wGraphics::initializeGraphics()
{
	// initializing GL3 now done by GLFW
	_shapes.setGl3wGraphics(this);
	// Background
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f );
	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	_tBall.computeQuat(_rotQuat,0.0f,0.0f,0.0f,0.0f);
}

void gl3wGraphics::setViewport(unsigned short x, unsigned short y, unsigned short xSize, unsigned short ySize)
{
	assert(x==0 && y==0);
	_xSize=xSize;
	_ySize=ySize;
	    // It's up to the application code to update the OpenGL viewport settings.
    // In order to avoid extensive context switching, consider doing this in
    // OnPaint() rather than here, though.
    glViewport(0, 0, (GLint) xSize, (GLint) ySize);
	_glM.setView(0.7f,(float)xSize/ySize);
}

void gl3wGraphics::mouseButtonEvent(unsigned short screenX, unsigned short screenY, int button, bool dragging)
{
	if(dragging)	{
		_dragging=true;
		if( button<1)	{	// leftMouse
	        /* drag in progress, simulate trackball */
			float spin_quat[4];
			_tBall.computeQuat(spin_quat,
            (2.0f*_lastX - _xSize) / _xSize,
            (_ySize - 2.0f*_lastY) / _ySize,
            (2.0f*screenX - _xSize)    / _xSize,
            (_ySize - 2.0f*screenY)    / _ySize);
			_tBall.add_quats(spin_quat, _rotQuat, _rotQuat);
		}
		else if (button > 1) {	// rightMouse
			if (!mouseWheelZoom)
				_glM.changeZoom((float)(_lastY - screenY) / _ySize);
		}
		else	{	// middleMouse
			_glM.setPan((float)(_lastX-screenX)/_xSize,(float)(screenY-_lastY)/_ySize);
		}
	}
	_lastX = screenX;
	_lastY = screenY;
}

gl3wGraphics::gl3wGraphics() : _staticCount(0)
{
	_dragging=false;
	_lastX=0;
	_lastY=0;
	_ls.setGLmatrices(&_glM);
	_lines.setGl3wGraphics(this);
	_nodes.clear();
}

gl3wGraphics::~gl3wGraphics()
{
}

void gl3wGraphics::drawAll()
{
	// This next line is normally only necessary if there is more than one wxGLCanvas
    // or more than one wxGLContext in the application.
    //SetCurrent(*m_glRC);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    GLfloat m[4][4];
    _tBall.build_rotmatrix( m,_rotQuat);
	// assumes perspective-frame matrix has already been set
	_glM.setFrameAndRotation(&m[0][0]);
//	std::list<sceneNode*>::iterator nit;
	GLuint currentProgram=0;
	for(auto nit = _nodes.begin(); nit != _nodes.end(); ++nit)	{ // textured TRIANGLES will always happen first
		if (!(*nit)->visible)  continue;
		if((*nit)->getGlslProgramNumber()!=currentProgram) {
			currentProgram=(*nit)->getGlslProgramNumber();
			_ls.useGlslProgram(currentProgram);
		}
		_ls.setModelMatrix((*nit)->getModelViewMatrix());	// must reset with new program
		(*nit)->draw();
	}
    glFlush(); // Not really necessary: buffer swapping below implies glFlush()
}

void gl3wGraphics::computeAndSetSceneRadius()
{ // does not change scene center
	float center[3],radius;
	_glM.getSceneCenter(center);
	radius = _glM.getSceneRadius();
	for(auto nit=_nodes.begin(); nit!=_nodes.end(); ++nit) {
		float c[3],r,d[3];
		(*nit)->getBounds(c,r,false);
		d[0]=center[0]-c[0]; d[1]=center[1]-c[1]; d[2]=center[2]-c[2];
		r += sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
		if(r>radius)
			radius = r;
	}
	_glM.setSceneRadius(radius);
}

void gl3wGraphics::getSceneSphere(GLfloat (&center)[3], GLfloat &radius, bool recomputeAll)	{
	auto nit=_nodes.begin();
	(*nit)->getBounds(center,radius,recomputeAll);
	++nit;
	while(nit!=_nodes.end()) {
		float c[3],r,d[3],l;
		(*nit)->getBounds(c,r,recomputeAll);
		d[0]=center[0]-c[0]; d[1]=center[1]-c[1]; d[2]=center[2]-c[2];
		l=sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
		if(l<1e-16)	{
			if(r>radius)
				radius=r;
		}
		else if((l+r)>radius) {	// expand scene sphere
			if(r-l>radius) {
				center[0]=c[0]; center[1]=c[1]; center[2]=c[2];
				radius = r;
			}
			else {
				float p = (l+r-radius)*0.5f/l;
				center[0]-=d[0]*p; center[1]-=d[1]*p; center[2]-=d[2]*p;
				radius += l+r;
				radius *= 0.5f;
			}
		}
		else
			;
		++nit;
	}
}

void gl3wGraphics::frameScene(bool recomputeAll)
{
	float r,c[3];
	getSceneSphere(c,r,recomputeAll);
	_glM.setSceneCenter(c);
	_glM.setSceneRadius(r*0.6f);
	_glM.resetPerspective();
}

void gl3wGraphics::getTrianglePickLine(float(&lineStartPosition)[3], float(&lineDirection)[3])
{
	float sc[3];
	float zCenter, height, verticalAngle, screenAspect;
	_glM.getCameraData(zCenter, height, verticalAngle, screenAspect);
	height = -zCenter*2.0f*tanf(verticalAngle*0.5f);
	sc[0] = (-0.5f + (float)_lastX / _xSize)*height*screenAspect;
	sc[1] = (0.5f - (float)_lastY / _ySize)*height;
	sc[2] = zCenter;
	GLfloat m[16], invM[16], *om;
	const GLfloat *fr = _glM.getFrameAndRotationMatrix();
	for (auto nit = _nodes.begin(); nit != _nodes.end(); ++nit)	{
		if ((*nit)->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)  // was TRISTRIP. Why?
			continue;
		om = (*nit)->getModelViewMatrix();
		for (int i = 0; i < 4; ++i)	{
			for (int j = 0; j < 4; ++j)	// model happens first, then frame-rotation
				m[(i << 2) + j] = fr[j] * om[i << 2] + fr[4 + j] * om[(i << 2) + 1] + fr[8 + j] * om[(i << 2) + 2] + fr[12 + j] * om[(i << 2) + 3];
		}
		invertMatrix4x4(m, invM);
		transformVector3(sc, invM, lineDirection);
		// start point is at origin
		lineStartPosition[0] = invM[12]; lineStartPosition[1] = invM[13]; lineStartPosition[2] = invM[14];
		lineDirection[0] -= lineStartPosition[0]; lineDirection[1] -= lineStartPosition[1]; lineDirection[2] -= lineStartPosition[2];
	}
}

bool gl3wGraphics::pick(unsigned short x, unsigned short y, std::string &name, float (&position)[3], int &triangle, bool excludeShapes, bool excludeStatic)
{
	bool picked=false;
	triangle = -1;
	name = "";
	float sc[3],f[3],b[3],minParam=1e30f;
	float zCenter,height,verticalAngle,screenAspect;
	_glM.getCameraData(zCenter,height,verticalAngle,screenAspect);
	height = -zCenter*2.0f*tanf(verticalAngle*0.5f);
	sc[0] = (-0.5f+(float)x/_xSize)*height*screenAspect;
	sc[1] = (0.5f-(float)y/_ySize)*height;
	sc[2] = zCenter;
	GLfloat m[16],invM[16],*om;	// mfp[16],
	const GLfloat *fr=_glM.getFrameAndRotationMatrix();
	for(auto nit=_nodes.begin(); nit!=_nodes.end(); ++nit)	{
		if(excludeShapes && (*nit)->getType() != sceneNode::nodeType::MATERIAL_TRIANGLES)
			continue;
		if (excludeStatic && (*nit)->getType() == sceneNode::nodeType::STATIC_TRIANGLES)
			continue;
		om = (*nit)->getModelViewMatrix();
		for(int i=0; i<4; ++i)	{
			for(int j=0; j<4; ++j)	// model happens first, then frame-rotation
				m[(i<<2)+j] = fr[j]*om[i<<2] + fr[4+j]*om[(i<<2)+1] + fr[8+j]*om[(i<<2)+2] + fr[12+j]*om[(i<<2)+3];
		}
		invertMatrix4x4(m,invM);
		transformVector3(sc,invM,f);
		// start point is at origin
		b[0]=invM[12]; b[1]=invM[13]; b[2]=invM[14];
		f[0]-=b[0]; f[1]-=b[1]; f[2]-=b[2];
		if((*nit)->getType()==sceneNode::nodeType::CONE) {
			if(_shapes.pickCone(b,f,position,minParam)) {
				picked=true;
				triangle = -1;
				name = (*nit)->getName();
				float vtx[]={position[0],position[1],position[2]};
				transformVector3(vtx,om,position);
			}
		}
		else if((*nit)->getType()==sceneNode::nodeType::SPHERE) {
			if(_shapes.pickSphere(b,f,position,minParam)) {
				picked=true;
				triangle = -1;
				name = (*nit)->getName();
				float vtx[]={position[0],position[1],position[2]};
				transformVector3(vtx,om,position);
			}
		}
		else if((*nit)->getType()==sceneNode::nodeType::CYLINDER) {
			if(_shapes.pickCylinder(b,f,position,minParam)) {
				picked=true;
				triangle = -1;
				name = (*nit)->getName();
				float vtx[]={position[0],position[1],position[2]};
				transformVector3(vtx,om,position);
			}
		}
		else if((*nit)->getType()==sceneNode::nodeType::LINES)
			continue;
		else if ((*nit)->getType() == sceneNode::nodeType::MATERIAL_TRIANGLES) {
			surgGraphics *sg = static_cast<surgGraphics*>((*nit)->getSurgGraphics());
			materialTriangles *tr = sg->getMaterialTriangles();
			float tp[2];
			if(tr->localPick(b,f,position,triangle,tp)) {
				picked=true;
				name = (*nit)->getName();
				float vtx[]={position[0],position[1],position[2]};
				transformVector3(vtx,om,position);
			}
		}
		else {	// (*nit)->getType()==sceneNode::TRIANGLES
			// as of now - do nothing. Static triangles just scenery
/*			staticTriangle *tr = static_cast<staticTriangle*>(*nit);
			if(tr->localPick(b,f,position,triangle,minParam)) {
				picked=true;
				name = (*nit)->getName();
				float vtx[]={position[0],position[1],position[2]};
				transformVector3(vtx,om,position);
			} */
		}
	}
	return picked;
}

void gl3wGraphics::deleteSceneNode(std::shared_ptr<sceneNode> &sn)
{
//	auto typ = sn->getType();
//	if (typ == sceneNode::nodeType::CONE || typ == sceneNode::nodeType::SPHERE || typ == sceneNode::nodeType::CYLINDER) {
//		_shapes.deleteShape(dynamic_cast<shape*>(sn));	// all done in deleteShape() call
//		return;
//	}
	if (!sn)
		return;
	for (auto nit = _nodes.begin(); nit != _nodes.end(); ++nit) {
		if (sn == *nit) {
			_nodes.erase(nit);
			return;
		}
	}
}

sceneNode* gl3wGraphics::getNodePtr(std::string &name)
{
	for(auto nit=_nodes.begin(); nit!=_nodes.end(); ++nit)	{
		if((*nit)->getName() == name)
			return (nit->get());
	}
	return nullptr;
}

void gl3wGraphics::clear() {	// empties all graphics

	assert(false);  // no longer used?

    _texReader.clear();
    _lines.clear();
    _ls.clear();
    _nodes.clear();
}

