//////////////////////////////////////////////////////////
// File: surgGraphics.h
// Author: Court Cutting, MD
// Date: Originally 3/3/2015 but a major overhaul to allow texture seams on 3/24/21
// Purpose: Takes as input a materialTriangles class and makes it visible
//    on an openGL canvas.  It creates hard normal edges on
//    the model given an input texture file for the top skin. Lip side
//    texturing, bottom fat and deep muscle are textured procedurally in
//    the fragment shader.
//////////////////////////////////////////////////////////

#ifndef __SURG_GRAPHICS__
#define __SURG_GRAPHICS__

#include <vector>
#include <list>
#include <set>
#include <memory>
#include "sceneNode.h"
#include "materialTriangles.h"

// forward declarations
class gl3wGraphics;
class Vec3f;

class incisionLines
{
public:
	bool isInitialized(){ return _sn && !_sn->bufferObjects.empty(); }
	void sendVertexCoordBuffer(GLuint vertCoordBuf){ _incisionBufferObjects[1] = vertCoordBuf; }  // must be called before addIncisions()
	void addUpdateIncisions(const std::vector<GLuint> &lines);  // 0xffffffff is primitive restart index
	void setGl3wGraphics(gl3wGraphics *gl3w) { _gl3w = gl3w; }
	void setColor(GLfloat(&color)[4]) { for (int i = 0; i < 4; ++i) _color[i] = color[i]; }
	incisionLines() : _sn(nullptr) {}
	~incisionLines(){}

private:
	std::shared_ptr<sceneNode> _sn;
	GLfloat _color[4];
	gl3wGraphics *_gl3w;
	static GLuint _incisionBufferObjects[2];
	static GLuint _incisionVertexArrayBufferObject;
};


class surgGraphics
{
public:
	void draw(void);
	void computeLocalBounds();
	bool setTextureFilesCreateProgram(std::vector<int> &textureIds, const char *vertexShaderFile, const char *fragmentShaderFile);  // must be set first before next 2 routines can be called
	void setNewTopology();
	void updatePositionsNormalsTangents();
	inline 	incisionLines* getIncisionLines() { return &_incis; }
	inline materialTriangles* getMaterialTriangles() {return & _mt;}  // gets the material triangles data class
	inline sceneNode* getSceneNode() { return _sn.get(); }
	void setGl3wGraphics(gl3wGraphics *gl3w) { _gl3w = gl3w; }  // must be set before using this class for graphics output
	surgGraphics(void);
	~surgGraphics(void);

private:
	static const GLchar *skinVertexShader;
	static const GLchar *skinFragmentShader;
	materialTriangles _mt;
	gl3wGraphics *_gl3w;
	std::shared_ptr<sceneNode> _sn;
	std::vector<GLuint> _tris;  // 0xffffffff signals a deleted triangle
	std::vector<GLfloat> _xyz1;
	std::vector<GLfloat> _uv;
	std::vector<int> _uvPos;
	const std::vector<long> *_undermineTriangles;  // if not NULL shade these with material 10
	std::vector<GLuint> _incisionLines;  // indexes into incision lines. 0xffffffff is primitive restart index.
	incisionLines _incis;
	void getSkinIncisionLines();

};

#endif  // __SURG_GRAPHICS__
