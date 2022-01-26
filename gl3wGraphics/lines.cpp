#include <assert.h>
#include "gl3wGraphics.h"
#include "boundingBox.h"
#include "lines.h"

void lines::drawLines()
{
	glPrimitiveRestartIndex(0xffffffff);
	glEnable(GL_PRIMITIVE_RESTART);
	glBindVertexArray(_sn->vertexArrayBufferObject);
	//	assumes glUseProgram(_program) has already been called
	if (!_colorOffsets.empty()){
		GLsizei start = 0;
		for(int n = (int)_colors.size(), i = 0; i < n; ++i){
			_gl3w->getLightsShaders()->setColor(_colors[i].data());
			glDrawElements(GL_LINE_STRIP, _colorOffsets[i]-start, GL_UNSIGNED_INT, (char *)NULL + start*sizeof(GL_UNSIGNED_INT));
			start = _colorOffsets[i];
		}
	}
	else{
		_gl3w->getLightsShaders()->setColor(_sn->getColor());
		glDrawElements(GL_LINE_STRIP, _linesSize, GL_UNSIGNED_INT, 0);
	}
	// Unbind to anybody
	glBindVertexArray(0);
	glDisable(GL_PRIMITIVE_RESTART);
}

void lines::computeLocalBounds()
{
	glBindBuffer(GL_ARRAY_BUFFER, _sn->bufferObjects[0]);	// VERTEX_DATA
	GLint size;
	glGetBufferParameteriv(GL_ARRAY_BUFFER, GL_BUFFER_SIZE, &size);
	std::vector<GLfloat> vtx;
	vtx.assign(size / sizeof(GLfloat), 0.0f);
	glGetBufferSubData(GL_ARRAY_BUFFER, 0, size, &vtx[0]);
	boundingBox<GLfloat> bb;
	bb.Empty_Box();
	size /= sizeof(GLfloat);
	for (int i = 0; i < size; i += 4)
		bb.Enlarge_To_Include_Point((GLfloat (&)[3])vtx[i]);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
//	bb.Center(_localCenter);
	GLfloat len[3];
	bb.Edge_Lengths(len);
	float radius;
	radius = 0.0f;
	for (int i = 0; i < 3; ++i)
		radius += len[i] * len[i] * 0.25f;
	radius = (float)sqrt(radius);
	_sn->setRadius(radius);
}

void lines::clear()
{
	if (!_sn)
		return;
	if (_sn->bufferObjects[0]>0) {
		glBindBuffer(GL_ARRAY_BUFFER, _sn->bufferObjects[0]);	// VERTEX_DATA
		glBufferData(GL_ARRAY_BUFFER, 0, NULL, GL_DYNAMIC_DRAW);
		// Indexes
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _sn->bufferObjects[1]);	// INDEX_DATA
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, NULL, GL_STATIC_DRAW);
	}
	// release for next use
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void lines::remove()
{
	if (!_sn)
		return;
	_gl3w->deleteSceneNode(_sn);
	if(_sn->vertexArrayBufferObject)
		glDeleteVertexArrays(1,&_sn->vertexArrayBufferObject);
	_sn->vertexArrayBufferObject = 0;
	if (_sn->bufferObjects[0]) {
		glDeleteBuffers(2, &_sn->bufferObjects[0]);
		_sn->bufferObjects.clear();
	}
	_sn.reset();
}

bool lines::updatePoints(const std::vector<GLfloat> &points)  // updates point positions related to initial addLines() call
{
	if(_pointsSize!=(GLsizei)points.size())
		return false;
	glBindBuffer(GL_ARRAY_BUFFER, _sn->bufferObjects[0]);	// VERTEX_DATA
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GLfloat)*_pointsSize, &(points[0]));
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	return true;
}

void lines::addLines(const std::vector<GLfloat> &points, const std::vector<GLuint> &lines, const std::vector<std::array<float, 4> > &colors, const std::vector<int> &colorOffsets)
{
	assert(colors.size() == colorOffsets.size());
	_colors.assign(colors.size(), std::array<GLfloat, 4>());
	for (int n = (int)colors.size(), i = 0; i < n; ++i){
		_colors[i][0] = (GLfloat)colors[i][0];
		_colors[i][1] = (GLfloat)colors[i][1];
		_colors[i][2] = (GLfloat)colors[i][2];
		_colors[i][3] = (GLfloat)colors[i][3];
	}
	_colorOffsets.assign(colorOffsets.begin(), colorOffsets.end());
	addLines(points, lines);
}

void lines::addLines(const std::vector<GLfloat> &points, const std::vector<GLuint> &lines)
{	// each element of points array has 4 floats as xyz1. lines uses 0xffffffff for
	// line strip primitive restart index.
//	if(lines.empty())
//		return;
/*	if (!_sn) {
		_sn = std::make_shared<sceneNode>();
		_sn->setGl3wGraphics(_gl3w);
	}
	if (!_sn->vertexArrayBufferObject) {
		glGenVertexArrays(1, &_sn->vertexArrayBufferObject);
		_gl3w->addSceneNode(_sn);
		GLuint program = _gl3w->getLightsShaders()->getOrCreateLineProgram();
		_sn->setGlslProgramNumber(program);
		_sn->setColorLocation(glGetUniformLocation(program, "objectColor"));
		GLfloat color[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
		_sn->setColor(color);
		_sn->setType(sceneNode::nodeType::LINES);
		GLfloat* mm = _sn->getModelViewMatrix();
		loadIdentity4x4(mm);
		_sn->setName("tetLines");
		assert(_sn->bufferObjects.empty()); // _linesBufferObjects[0] < 1);
		_sn->bufferObjects.assign(2, 0);
		glGenBuffers(2, &_sn->bufferObjects[0]);
	} */
	if(!_visible)
		setLinesVisible(true);
	_linesSize = (GLsizei)lines.size();
	_pointsSize = (GLsizei)points.size();
	// Vertex and normal data
    glBindBuffer(GL_ARRAY_BUFFER, _sn->bufferObjects[0]);	// VERTEX_DATA
    glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*points.size(), &(points[0]), GL_DYNAMIC_DRAW);
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _sn->bufferObjects[1]);	// INDEX_DATA
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*lines.size(), &(lines[0]), GL_STATIC_DRAW);
	_sn->elementArraySize = (GLsizei)lines.size();
	// Create the master vertex array object
	glBindVertexArray(_sn->vertexArrayBufferObject);
    // Vertex data
    glBindBuffer(GL_ARRAY_BUFFER, _sn->bufferObjects[0]);	// VERTEX DATA
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
    // Indexes
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _sn->bufferObjects[1]);	// INDEX_DATA
    // Unbind to anybody
	glBindVertexArray(0);
	// release for next use
    glBindBuffer( GL_ARRAY_BUFFER, 0);
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0);
}

void lines::setLinesVisible(bool visible){
	_visible = visible;
	if (!visible) {
		if (_sn)
			_sn->visible = false;
		return;
	}
	if (!_sn) {
		_sn = std::make_shared<sceneNode>();
		_sn->setGl3wGraphics(_gl3w);
	}
	_sn->visible = true;
	if (!_sn->vertexArrayBufferObject) {
		glGenVertexArrays(1, &_sn->vertexArrayBufferObject);
		GLuint program = _gl3w->getLightsShaders()->getOrCreateLineProgram();
		_sn->setGlslProgramNumber(program);
		_sn->setColorLocation(glGetUniformLocation(program, "objectColor"));
		GLfloat color[4] = { 1.0f, 1.0f, 1.0f, 1.0f };
		_sn->setColor(color);
		_sn->setType(sceneNode::nodeType::LINES);
		GLfloat* mm = _sn->getModelViewMatrix();
		loadIdentity4x4(mm);
		_sn->setName("tetLines");
		assert(_sn->bufferObjects.empty()); // _linesBufferObjects[0] < 1);
		_sn->bufferObjects.assign(2, 0);
		glGenBuffers(2, &_sn->bufferObjects[0]);
		_gl3w->addSceneNode(_sn);
	}
}


lines::lines() : _visible(false)
{
	_colors.clear();
	_colorOffsets.clear();
	_sn = nullptr;
}


lines::~lines()
{
	_colors.clear();
	_colorOffsets.clear();
	remove();  // delete any openGL buffers
}
