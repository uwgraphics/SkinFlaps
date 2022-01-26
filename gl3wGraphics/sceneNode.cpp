// File: sceneNode.cpp
// Author: Court Cutting, MD
// Date: February 7, 2012
// Updated: June 12, 2021
// Purpose: Basic data and methods every sceneNode must have.

#include "gl3wGraphics.h"
#include "surgGraphics.h"
#include "sceneNode.h"

gl3wGraphics* sceneNode::_gl3w = nullptr;
surgGraphics* sceneNode::_sg = nullptr;

//	virtual void computeLocalBounds() {}
//	virtual void getLocalBounds(GLfloat (&localCenter)[3], GLfloat &Radius) {Radius=0; localCenter[0]=0;}
void sceneNode::getLocalBounds(GLfloat (&localCenter)[3], GLfloat &Radius) {
	localCenter[0]=_localCenter[0]; localCenter[1]=_localCenter[1]; localCenter[2]=_localCenter[2];
	Radius = _radius;
}

void sceneNode::setLocalBounds(GLfloat(&localCenter)[3], GLfloat& Radius) {
	_localCenter[0] = localCenter[0]; _localCenter[1] = localCenter[1]; _localCenter[2] = localCenter[2];
	_radius = Radius;
	_boundsComputed = true;
}

void sceneNode::draw(void) {
	//	assumes glUseProgram(_program) has already been called
	if (_type == nodeType::MATERIAL_TRIANGLES)  // complex dynamic object. Draw externally
		_sg->draw();
	else if (_type == nodeType::STATIC_TRIANGLES) {
		glBindVertexArray(vertexArrayBufferObject);
		int n = (int)textureBuffers.size();
		for (int i = 0; i < n; ++i) {
			glActiveTexture(GL_TEXTURE0 + i);
			glBindTexture(GL_TEXTURE_2D, textureBuffers[i]);
		}
		glDrawElements(GL_TRIANGLES, elementArraySize, GL_UNSIGNED_INT, 0);
		// Never unbind a GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER inside an active vertex array buffer object
		glBindVertexArray(0);
	}
	else if (_type == nodeType::TRISTRIP) {
		glUniform4fv(_locObjColor, 1, (GLfloat*)_color);
		glBindVertexArray(vertexArrayBufferObject);
		glDrawArrays(GL_TRIANGLE_STRIP, 0, elementArraySize);
		glBindVertexArray(0);
	}
	else if (_type == nodeType::LINES) {
		glPrimitiveRestartIndex(0xffffffff);
		glEnable(GL_PRIMITIVE_RESTART);
	//	glLineWidth(2.0);
		glUniform4fv(_locObjColor, 1, (GLfloat*)_color);
		glBindVertexArray(vertexArrayBufferObject);
		glDrawElements(GL_LINE_STRIP, elementArraySize, GL_UNSIGNED_INT, 0);
		glBindVertexArray(0);
	//	glLineWidth(1.0);
		glDisable(GL_PRIMITIVE_RESTART);
	}
	else if (_type == nodeType::CONE) {
		glUniform4fv(_locObjColor, 1, (GLfloat*)_color);
		_gl3w->getShapes()->drawCone();
	}
	else if (_type == nodeType::CYLINDER) {
		glUniform4fv(_locObjColor, 1, (GLfloat*)_color);
		_gl3w->getShapes()->drawCylinder();
	}
	else if (_type == nodeType::SPHERE) {
		glUniform4fv(_locObjColor, 1, (GLfloat*)_color);
		_gl3w->getShapes()->drawSphere();
	}
	else {
		std::cout << "Trying to draw a sceneNode of unknown type.";
	}

	GLenum errCode;
	if((errCode=glGetError())!=GL_NO_ERROR)
		std::cout << "Graphics draw error code number " << errCode << "\n";
}

void sceneNode::setType(nodeType type) {
	_type = type;
	if (type == sceneNode::nodeType::STATIC_TRIANGLES || type == sceneNode::nodeType::MATERIAL_TRIANGLES)
		_coloredNotTextured = false;
	else
		_coloredNotTextured = true;
}

void sceneNode::getBounds(GLfloat(&center)[3], GLfloat& radius, bool recomputeAll)
{
	if (recomputeAll || !_boundsComputed)
		getLocalBounds(_localCenter, _radius);
	transformVector3(_localCenter, _pat, center);
	GLfloat r = _pat[0] * _radius; r *= r;	// assume no unequal scaling. In some apps this will be a bug. Done for speed.
	r += _pat[4] * _radius * _pat[4] * _radius;
	r += _pat[8] * _radius * _pat[8] * _radius;
	radius = (float)sqrt(r);
}

sceneNode::sceneNode() : _radius(-1.0f)
{
	_boundsComputed = false;
	loadIdentity4x4(_pat);
	_coloredNotTextured = true;
	_glslProgram = 0;
	textureBuffers.clear();
	_color[0] = 1.0f; _color[1] = 1.0f; _color[2] = 1.0f; _color[3] = 1.0f;
	visible = true;
	vertexArrayBufferObject = 0;
	bufferObjects.clear();
}

sceneNode::~sceneNode() {
	if (_type == nodeType::STATIC_TRIANGLES) {  // since multiple instances, do here.  For others do in parent.
		if (!bufferObjects.empty()) {
			glDeleteBuffers((GLsizei)bufferObjects.size(), &bufferObjects[0]);
			bufferObjects.clear();
		}
		if (vertexArrayBufferObject > 0) {
			glDeleteVertexArrays(1, &vertexArrayBufferObject);
			vertexArrayBufferObject = 0;
		}
		if(!textureBuffers.empty())
			glDeleteTextures((GLsizei)textureBuffers.size(), &(textureBuffers[0]));
	}
}

