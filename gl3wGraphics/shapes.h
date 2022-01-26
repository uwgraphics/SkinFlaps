// File: shapes.h
// Author: Court Cutting
// Date: 3/3/12
// Purpose: Manages untextured cones, cylinders and spheres for graphics.

#ifndef __SHAPES_H__
#define __SHAPES_H__

#include "sceneNode.h"
#include <memory>
#include <list>

// forward declaration
class gl3wGraphics;

class shapes
{
public:
	std::shared_ptr<sceneNode> addShape(const sceneNode::nodeType &type, char *name);
	void deleteShape(std::shared_ptr<sceneNode> shapePtr);
	void drawCylinder();
	bool pickCylinder(const float *lineStart, float *lineDirection, float (&position)[3], float &param);
	void drawSphere();
	bool pickSphere(const float *lineStart, float *lineDirection, float (&position)[3], float &param);
	void drawCone();
	bool pickCone(const float *lineStart, float *lineDirection, float (&position)[3], float &param);
	void getLocalBounds(std::shared_ptr<sceneNode>& sn);
	void setGl3wGraphics(gl3wGraphics *gl3w) { _gl3w = gl3w; }
	shapes();
	~shapes();

private:
	gl3wGraphics *_gl3w;
	void getOrCreateConeGraphic(); // creates cone point at origin with base at z=1.0 with 1.0 diameter
	void getOrCreateSphereGraphic(); // creates sphere at origin with diameter 1.0
	void getOrCreateCylinderGraphic(); // creates cylinder at origin with diameter 1.0, from z=-1 to z=1
	static GLuint _coneBufferObjects[2];
	static GLuint _coneVertexArrayBufferObject;
	static GLuint _sphereBufferObjects[3];
	static GLuint _sphereVertexArrayBufferObject;
	static GLuint _cylinderBufferObjects[2];
	static GLuint _cylinderVertexArrayBufferObject;

};

#endif	// __SHAPES_H__
