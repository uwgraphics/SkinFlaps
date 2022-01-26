// staticTriangle.cpp
// Author: Court Cutting
// Date: 5/20/2014
// Purpose: Triangle object management class for static triangulated objects.
//    Meant to complement triangleUVW class for management of dynamic objects.
//    Topology and textures should be sent anly once. Vertex positions positions
//    however, may be changed and normals and tangents can be recomputed.
//       Copyright 2014 - All rights reserved.

#pragma warning(disable : 4996)
#include <fstream>
#include <set>
#include <assert.h>
#ifdef linux
#include <stdlib.h>
#include <string.h>
#endif
#include "math3d.h"
#include "Vec3f.h"
#include "boundingBox.h"
#include "gl3wGraphics.h"
#include "staticTriangle.h"

//////////////////////// TEMPORARY TEMPORARY TEMPORARY - On SnowLeopard this is suppored, but GLEW doens't hook up properly
//////////////////////// Fixed probably in 10.6.3
#ifdef __APPLE__
#define glGenVertexArrays glGenVertexArraysAPPLE
#define glDeleteVertexArrays  glDeleteVertexArraysAPPLE
#define glBindVertexArray	glBindVertexArrayAPPLE
#endif

GLuint staticTriangle::_staticProgram = 0;

static const GLchar* staticVertexShader = "#version 130 \n"
"in vec4 vVertex;"
"in vec3 vTangent;"
"in vec3 vNormal;"
"in vec2 vTexture;"
"uniform mat4   mvpMatrix;"
"uniform mat4   mvMatrix;"
"uniform mat3   normalMatrix;"
"uniform vec3   vLightPosition;"
"smooth out vec3 vLightDir;"
"smooth out vec3 vEyeDir;"
"smooth out vec2 vTexCoords;"
"void main(void)"
"{"
"	vec3 n = normalize(normalMatrix * vNormal);"
"	vec3 t = normalize(normalMatrix * vTangent);"
"	vec3 b = cross(n, t);"
"	vec3 v;"
"	v.x = dot(vLightPosition, t);"
"	v.y = dot(vLightPosition, b);"
"	v.z = dot(vLightPosition, n);"
"	vLightDir = normalize(v);"
"	vEyeDir = vec3(mvMatrix * vVertex);"
"	v.x = dot(vEyeDir, t);"
"	v.y = dot(vEyeDir, b);"
"	v.z = dot(vEyeDir, n);"
"	vEyeDir = normalize(v);"
"	vTexCoords = vTexture;"
"	gl_Position = mvpMatrix * vVertex;"
"}";

static const GLchar* staticFragmentShader =
"#version 150 \n"
"out vec4 vFragColor; "
"uniform vec4 ambientColor; "
"uniform vec4 diffuseColor; "
"uniform vec4 objectColor; "
"uniform sampler2D colorMap; "
"uniform sampler2D normalMap; "
"uniform sampler2D texture2; "
"uniform sampler2D texture3; "
"uniform int material; "
"uniform mat3   normalMatrix; "
"smooth in vec3 vLightDir; "
"smooth in vec3 vEyeDir; "
"smooth in vec2 vTexCoords; "
"void main(void) "
"{"
"	const float ambientVal = 0.1;"
"	float lightVal = ambientVal;"
"	float specMult = 0.11;"
"	const float diffuseVal = 0.9;"
"	vec3 normDelta = vec3(0.0, 0.0, 1.0);"
"	vec3 litColor = vec3(1.0, 1.0, 1.0);"
"	vec4 tx1 = texture(normalMap, vTexCoords.st);"
"	tx1.rgb -= vec3(0.5);"
"	normDelta = tx1.rgb * 2.0;"
"	vFragColor = texture(colorMap, vTexCoords.st);"
"	lightVal += diffuseVal * max(dot(normDelta, vLightDir), 0.0);"
"	vFragColor *= lightVal;"
"	vec3 reflectDir = reflect(vLightDir, normDelta);"
"	float spec = max(dot(vEyeDir, reflectDir), 0.0);"
"	spec = pow(spec, 40.0);"
"	spec *= specMult;"
"	litColor = min(vFragColor.rgb + spec, vec3(1.0));"
"	vFragColor = vec4(litColor, 1.0); "
"}";

/* void staticTriangle::draw(void)
{
	glBindVertexArray(_vertexArrayBufferObject);
	if(_textured) {
		int n = (int)_2DtextureBuffers.size();
		for (int i = 0; i < n; ++i) {
			glActiveTexture(GL_TEXTURE0 + i);
			glBindTexture(GL_TEXTURE_2D, _2DtextureBuffers[i]);
		}
	}
	glDrawElements(GL_TRIANGLES, _triangleArraySize, GL_UNSIGNED_INT, 0);
    // Never unbind a GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER inside an active vertex array buffer object
	glBindVertexArray(0);

// GLenum errCode;
//if((errCode=glGetError())!=GL_NO_ERROR)
//	errCode=errCode;

}  */  

void staticTriangle::computeLocalBounds()
{
	boundingBox<float> bb;
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i)
		bb.Enlarge_To_Include_Point((const float(&)[3])(*_mt->vertexCoordinate(i)));
	Vec3f vmin, vmax;
	bb.Minimum_Corner(vmin._v);
	bb.Maximum_Corner(vmax._v);
	vmin += vmax;
	vmin *= 0.5f;
	GLfloat lc[3], radius;
	lc[0] = vmin.X;
	lc[1] = vmin.Y;
	lc[2] = vmin.Z;
	radius = (vmax - vmin).length();
	assert(_snNow.use_count() > 0);
	_snNow->setLocalBounds(lc, radius);
}

bool staticTriangle::createStaticProgram() {
	if (_staticProgram > 0)
		return true;
	std::vector<std::string> att;
	att.assign(4, std::string());
	att[0] = "vVertex";
	att[1] = "vNormal";
	att[2] = "vTangent";
	att[3] = "vTexture";
	if (!_gl3w->getLightsShaders()->createCustomProgram(_staticProgram, staticVertexShader, staticFragmentShader, att))
		return false;
	return true;
}

std::shared_ptr<sceneNode> staticTriangle::createStaticSceneNode(materialTriangles* mt, std::vector<int> &textureIds)
{  // must be set first
	_mt = mt;
	if (_staticProgram < 1) {
		if (!createStaticProgram())
			return false;
	}
	_snNow = std::make_shared<sceneNode>();
	_snNow->setType(sceneNode::nodeType::STATIC_TRIANGLES);
	// assumes textures loaded already and input in txIds
	std::sort(textureIds.begin(), textureIds.end());
	for (int n = (int)textureIds.size(), i = 0; i < n; ++i) {
		if (!_gl3w->getTextures()->textureExists(textureIds[i]))
			return false;
		_snNow->add2DtextureBufferNumber(_gl3w->getTextures()->getOGLtextureNumber(textureIds[i]));
	}
	_snNow->setGlslProgramNumber(_staticProgram);
//	_gl3w->getLightsShaders()->useGlslProgram(_glslProgram);  // must be current program. This routine sets other uniforms.
	if (_snNow->bufferObjects.empty()) {
		_snNow->bufferObjects.assign(5, 0);
		glGenBuffers(5, &_snNow->bufferObjects[0]);
	}
	if(!_snNow->vertexArrayBufferObject)
		glGenVertexArrays(1,&_snNow->vertexArrayBufferObject);
	// now make vertex array
	glBindVertexArray(_snNow->vertexArrayBufferObject);
    // Position data
    glBindBuffer(GL_ARRAY_BUFFER, _snNow->bufferObjects[0]);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    // Normal data
    glBindBuffer(GL_ARRAY_BUFFER, _snNow->bufferObjects[1]);
	glEnableVertexAttribArray(1);
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
	// Tangent data
    glBindBuffer(GL_ARRAY_BUFFER, _snNow->bufferObjects[2]);
	glEnableVertexAttribArray(2);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, 0);
    // Texture coordinates
	glBindBuffer(GL_ARRAY_BUFFER, _snNow->bufferObjects[3]);	// TEXTURE_DATA
	glEnableVertexAttribArray(3);
	glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 0, 0);
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _snNow->bufferObjects[4]);	// INDEX_DATA
	// never unbind a GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER inside a vertexArrayBuffer
	glBindVertexArray(0);

	computeNormalsTangents();
	computeLocalBounds();

	return _snNow;
}

void staticTriangle::computeNormalsTangents()
{
	_mt->findAdjacentTriangles(true, false);  // when teeth.obj fixed do next line instead
//	_mt.collectCreateTextureSeams();  // calls findAdjacentTriangles()
	std::vector<GLfloat> xyz1, uv, normals, tangents;
	std::vector<GLuint> tris;
	std::vector<int> uvPos;
	tris.clear();
	xyz1.clear();
	uv.clear();
	std::vector<float>* mtta = _mt->getTextureArray();
	uv.assign(mtta->begin(), mtta->end());
	int n = (int)uv.size() >> 1;
	xyz1.clear();
	xyz1.assign(n << 2, 1.0f);
	uvPos.clear();
	uvPos.assign(n, -1);
	const materialTriangles::matTriangle* trArr = _mt->getTriangleArray(n);
	tris.reserve(n * 3);
	for (int i = 0; i < n; ++i) {
		// include possible deleted triangles so numbering matches up.
		bool valid = true;
		if (trArr[i].material < 0) {
			tris.push_back(0xffffffff);
			valid = false;
		}
		else
			tris.push_back(trArr[i].tex[0]);
		tris.push_back(trArr[i].tex[1]);
		tris.push_back(trArr[i].tex[2]);
		if (valid) {
			for (int j = 0; j < 3; ++j) {
#ifdef DEBUG
				if (uvPos[trArr[i].tex[j]] > -1) {

				//	if (uvPos[trArr[i].tex[j]] != trArr[i].v[j])
				//		int junk = 0;

//					assert(uvPos[trArr[i].tex[j]] == trArr[i].v[j]);
				}
#endif
				uvPos[trArr[i].tex[j]] = trArr[i].v[j];
			}
		}
	}
	// Texture coordinates
	glBindBuffer(GL_ARRAY_BUFFER, _snNow->bufferObjects[3]);	// TEXTURE_DATA
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * uv.size(), &(uv[0]), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _snNow->bufferObjects[4]);	// INDEX_DATA
	// Eliminate deleted triangles from viewing, but to keep the numbering send to graphics card
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * tris.size(), &(tris[0]), GL_STATIC_DRAW);
	_snNow->elementArraySize = (GLsizei)(sizeof(GLuint) * tris.size());
	for (int m = (int)uvPos.size(), i = 0; i < m; ++i) {
		if (uvPos[i] < 0)
			continue;
		float* fp = _mt->vertexCoordinate(uvPos[i]);
		for (int j = 0; j < 3; ++j)
			xyz1[(i << 2) + j] = fp[j];
	}
	normals.assign((uv.size() >> 1) * 3, 0.0f);
	tangents.assign(normals.size(), 0.0f);
	int i = 0, j, k;
	GLfloat* gv[3], * tv[3];
	n = (unsigned int)tris.size();
	Vec3f nrmV, tanV, dXyz[2];
	float d2, dTx[2][2];
	for (i = 0; i < n; i += 3) {
		if (tris[i] > 0xfffffffe)
			continue;
		for (j = 0; j < 3; ++j) {
			nrmV[j] = 0.0f;
			gv[j] = &xyz1[tris[i + j] << 2];
			tv[j] = &uv[tris[i + j] << 1];
		}
		for (j = 0; j < 3; ++j) {
			dXyz[0][j] = gv[1][j] - gv[0][j];
			dXyz[1][j] = gv[2][j] - gv[0][j];
		}
		for (j = 0; j < 2; ++j) {
			dTx[0][j] = tv[1][j] - tv[0][j];
			dTx[1][j] = tv[2][j] - tv[0][j];
		}
		d2 = dTx[0][0] * dTx[1][1] - dTx[1][0] * dTx[0][1];
		if (fabs(d2) < 1e-16f)
			tanV.set(0.0f, 0.0f, 0.0f);
		else
			tanV = (dXyz[0] * dTx[1][1] - dXyz[1] * dTx[0][1]) / d2;
		nrmV = dXyz[0] ^ dXyz[1];
		for (j = 0; j < 3; ++j) {
			k = tris[i + j] * 3;
			normals[k] += nrmV[0];
			normals[k + 1] += nrmV[1];
			normals[k + 2] += nrmV[2];
			tangents[k] += tanV[0];
			tangents[k + 1] += tanV[1];
			tangents[k + 2] += tanV[2];
		}
	}
	auto oit = _mt->_oneMaterialSeams.begin();
	while (oit != _mt->_oneMaterialSeams.end()) {
		auto start = oit;
		do {
			++oit;
		} while (oit != _mt->_oneMaterialSeams.end() && oit->first == start->first);
		GLfloat ns[3] = { 0.0f, 0.0f, 0.0f }, ts[3] = { 0.0f, 0.0f, 0.0f };
		for (auto it = start; it != oit; ++it) {
			for (int j = 0; j < 3; ++j) {
				ns[j] += normals[it->second * 3 + j];
				ts[j] += tangents[it->second * 3 + j];
			}
		}
		while (start != oit) {
			for (int j = 0; j < 3; ++j) {
				normals[start->second * 3 + j] = ns[j];
				tangents[start->second * 3 + j] = ts[j];
			}
			++start;
		}
	}
	auto invSqrt = [](float x) ->float { // Steve Pizer's version of the Quake algorithm
		GLuint i = 0x5F1F1412 - (*(GLuint*)&x >> 1);
		float tmp = *(float*)&i;
		return tmp * (1.69000231f - 0.714158168f * x * tmp * tmp);
	};
	n = (int)normals.size();
	for (i = 0; i < n; i += 3) {
		d2 = normals[i] * normals[i] + normals[i + 1] * normals[i + 1] + normals[i + 2] * normals[i + 2];
		if (d2 < 1e-16f) {
			normals[i] = 0.0f; normals[i + 1] = 0.0f; normals[i + 2] = 1.0f;
		}
		else {
			d2 = invSqrt(d2);
			normals[i] *= d2; normals[i + 1] *= d2; normals[i + 2] *= d2;
		}
		d2 = tangents[i] * tangents[i] + tangents[i + 1] * tangents[i + 1] + tangents[i + 2] * tangents[i + 2];
		if (d2 < 1e-16f) {
			tangents[i] = 1.0f; tangents[i + 1] = 0.0f; tangents[i + 2] = 0.0f;
		}
		else {
			d2 = invSqrt(d2);
			tangents[i] *= d2; tangents[i + 1] *= d2; tangents[i + 2] *= d2;
		}
	}
	// Vertex data
	glBindBuffer(GL_ARRAY_BUFFER, _snNow->bufferObjects[0]);	// VERTEX_DATA */
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * xyz1.size(), &(xyz1[0]), GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, _snNow->bufferObjects[1]);	// NORMAL_DATA
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)* normals.size(), &(normals[0]), GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, _snNow->bufferObjects[2]);	// TANGENT_DATA
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * tangents.size(), &(tangents[0]), GL_STATIC_DRAW);
}


