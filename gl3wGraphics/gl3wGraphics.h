// File: gl3wGraphics.h
// Author: Court Cutting, MD
// Date: January 21, 2012 - renamed 11/6/2014
// Purpose: Interface to graphics library. This library uses dear_imgui and GLFW as its
//	GUI library and the source of an openGL window.  This library does all the basic graphics
//  operations for the next level up.
//	Copyright 2012 - Court Cutting - All rights reserved.

#ifndef __WX_OPENGL_H__
#define __WX_OPENGL_H__

#include "sceneNode.h"
#include "shapes.h"
#include "lines.h"
#include "textures.h"
#include "lightsShaders.h"
#include "trackball.h"
#include "GLmatrices.h"
#include <list>
#include <memory>
#include <string>
#include <GLFW/glfw3.h>
#include "staticTriangle.h"

class gl3wGraphics
{
public:
	textures _texReader;

	bool addCustomSceneNode(std::shared_ptr<sceneNode> &sn, std::vector<int> &txIds, const GLchar *vertexShader, const GLchar *fragmentShader, std::vector<std::string> &attributes);
	void mouseButtonEvent(unsigned short screenX, unsigned short screenY, int button, bool dragging);
	bool pick(unsigned short x, unsigned short y, std::string &name, float (&position)[3], int &triangle, bool excludeShapes=false, bool excludeStatic = true);
	void getTrianglePickLine(float(&lineStartPosition)[3], float(&lineDirection)[3]);
	void initializeGraphics();
	std::shared_ptr<sceneNode> loadStaticObjFile(const char *filePath, std::vector<int> &textureIds, bool texturedNotColored = true);
	void drawAll();	// draws everything in the scene
	void getSceneSphere(GLfloat (&center)[3], GLfloat &radius, bool recomputeAll=true);
	void frameScene(bool recomputeAll=true);	// sets view data for currently loaded scene
	void computeAndSetSceneRadius();	// alters scene radius without changing center or zoom
	void setViewport(unsigned short x, unsigned short y, unsigned short xSize, unsigned short ySize);
	inline shapes* getShapes() {return &_shapes;}
	inline lines* getLines() {return &_lines;}
	inline staticTriangle* getstaticTriangle() { return &_staticTris; }
	inline lightsShaders* getLightsShaders() {return &_ls;}
	inline textures* getTextures() { return &_texReader; }
	GLmatrices* getGLmatrices() {return &_glM;}
	void addSceneNode(std::shared_ptr<sceneNode> &sn) { sn->visible = true; if (sn->coloredNotTextured()) _nodes.push_back(sn); else _nodes.push_front(sn); }
	void deleteSceneNode(std::shared_ptr<sceneNode> &sn);
	sceneNode* getNodePtr(std::string &name);
	void clear();	// empties all graphics
	void zeroViewRotations() { _rotQuat[0] = 0.0f; _rotQuat[1] = 0.0f; _rotQuat[2] = 0.0f; _rotQuat[3] = 1.0f;}
	gl3wGraphics();
	~gl3wGraphics();

private:
	trackball _tBall;
	float _rotQuat[4];
	unsigned short _xSize,_ySize,_lastX,_lastY;
	int _staticCount;
	std::list<std::shared_ptr<sceneNode> > _nodes;
	lightsShaders _ls;
	shapes _shapes;
	lines _lines;
	staticTriangle _staticTris;
    GLfloat _xrot;
    GLfloat _yrot;
	bool _dragging;
	GLmatrices _glM;
	GLuint  texture;
};

#endif	// __WX_OPENGL_H__
