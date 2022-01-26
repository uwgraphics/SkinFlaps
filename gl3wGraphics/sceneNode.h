// File: sceneNode.h
// Author: Court Cutting, MD
// Date: February 7, 2012
// Purpose: Basic data and methods every sceneNode must have.

#ifndef __SCENENODE_H__
#define __SCENENODE_H__

#include "GLmatrices.h"
#include <string>
#include <memory>
#include <vector>

// forward declaration
class gl3wGraphics;
class surgGraphics;

class sceneNode
{
public:
	enum class nodeType{ CONE, SPHERE, CYLINDER, TRISTRIP, LINES, STATIC_TRIANGLES, MATERIAL_TRIANGLES };
	bool visible;
	void setType(nodeType type);
	inline nodeType getType() {return _type;}
	void getLocalBounds(GLfloat(&localCenter)[3], GLfloat& Radius);
	void setLocalBounds(GLfloat(&localCenter)[3], GLfloat& Radius);
	void draw(void);
	void getBounds(GLfloat(&center)[3], GLfloat& radius, bool recomputeAll);
	float getRadius() { return _radius; }  // COURT - fix me if radius < 0 compute it
	void setRadius(float& radius) { _radius = radius; }
	void setName(const char *name) {_name=name;}
	const std::string& getName() {return _name;}
	void add2DtextureBufferNumber(GLuint texBufNum) { textureBuffers.push_back(texBufNum); }
	void setGlslProgramNumber(GLuint progNum) {_glslProgram=progNum;}
	inline GLuint getGlslProgramNumber() {return _glslProgram;}
	inline GLfloat* getModelViewMatrix() {return _pat;}
	inline bool coloredNotTextured() {return _coloredNotTextured; }
	GLfloat* getColor() {return _color;}
	inline void setColorLocation(GLint	locObjColor) { _locObjColor = locObjColor; }
	inline void setColor(float (&color)[4]) {_color[0]=color[0]; _color[1]=color[1]; _color[2]=color[2]; _color[3]=color[3];}
	static void setGl3wGraphics(gl3wGraphics *gl3w) { _gl3w = gl3w; }
	static void setSurgGraphics(surgGraphics* sg) { _sg = sg; }
	static surgGraphics* getSurgGraphics() { return _sg; }
	sceneNode();
	~sceneNode();

	GLuint vertexArrayBufferObject;
	std::vector<GLuint> bufferObjects;
	GLsizei elementArraySize;
	std::vector<GLuint> textureBuffers;

protected:
	static surgGraphics* _sg;
	static gl3wGraphics* _gl3w;
	bool _coloredNotTextured;
	nodeType _type;
	std::string _name;
	GLfloat _pat[16];	// GL matrix for local position attitude transform
	bool _boundsComputed;
	GLfloat _localCenter[3],_radius,_color[4];
	GLuint _glslProgram;
	GLint	_locObjColor;		// For colored, non-textured objects

};

#endif	// __SCENENODE_H__
