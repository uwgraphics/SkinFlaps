#ifndef __LINES_H__
#define __LINES_H__

#include <vector>
#include <array>
#include <memory>
#include "sceneNode.h"

// forward declaration
class gl3wGraphics;

class lines
{
public:
	void drawLines();
	// addLines() follows. Each element of points array has 4 floats as xyz1. lines uses 0xffffffff for
	// line strip primitive restart index.
	void addLines(const std::vector<GLfloat> &points, const std::vector<GLuint> &lines);  // 0xffffffff is primitive restart index
	void addLines(const std::vector<GLfloat> &points, const std::vector<GLuint> &lines, const std::vector<std::array<float, 4> > &colors, const std::vector<int> &colorOffsets);
	bool updatePoints(const std::vector<GLfloat> &points);  // updates point positions related to initial addLines() call
	std::shared_ptr<sceneNode>& getSceneNode() { return _sn; }
	void clear();
	void remove();
	inline bool linesVisible() { return _visible; }
	inline void setLinesVisible(bool visible);
	//	bool empty() {return _linesVertexArrayBufferObject<1;	}
	void setGl3wGraphics(gl3wGraphics *gl3w) { _gl3w = gl3w; }
	void computeLocalBounds();
	lines();
	~lines();

private:
	gl3wGraphics *_gl3w;
	bool _visible;
	std::shared_ptr<sceneNode> _sn;
//	GLuint _linesBufferObjects[2];
//	GLuint _linesVertexArrayBufferObject;
	GLsizei _linesSize;
	GLsizei _pointsSize;
	std::vector<std::array<GLfloat, 4> > _colors;
	std::vector<GLuint> _colorOffsets;

};

#endif	// __LINES_H__