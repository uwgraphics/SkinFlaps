// staticTriangle.h
// Author: Court Cutting
// Date: 5/20/2014
// Purpose: Triangle object management class for display of static triangle objects.
//       Copyright 2014 - All rights reserved.

#ifndef __staticTriangle__
#define __staticTriangle__

#include <vector>
#include <string>
#include <list>
#include <map>
#include <memory>
#include <set>
#include "sceneNode.h"
#include "materialTriangles.h"

// forward declarations
class gl3wGraphics;

class staticTriangle
{
public:
	std::shared_ptr<sceneNode> createStaticSceneNode(materialTriangles *mt, std::vector<int> &textureIds);  // must be set first before next 2 routines can be called
	void setGl3wGraphics(gl3wGraphics *gl3w) { _gl3w = gl3w; }  // must be set before using this class for graphics output
	staticTriangle() { _snNow = nullptr; }
	~staticTriangle(){}

private:
	materialTriangles *_mt;
	gl3wGraphics *_gl3w;
	std::shared_ptr<sceneNode> _snNow;  // sceneNode currently being worked on
//	GLsizei _triangleArraySize;
//	GLuint _bufferObjects[5];
//	GLuint _vertexArrayBufferObject;

	static GLuint _staticProgram;

	bool createStaticProgram();
	void computeLocalBounds();
	void computeNormalsTangents();

};

#endif    // __staticTriangle__
