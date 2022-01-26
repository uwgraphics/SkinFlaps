// File: fence.h
// Author: Court Cutting, MD
// Date: June 7, 2012
// Purpose: User interface for creating a fence on a glslTriangle object.
//     This will be used to specify a desired incision line.

#ifndef __FENCE_H__
#define __FENCE_H__

#include <list>
#include <set>
#include <memory>
#include "Vec3f.h"
#include "gl3wGraphics.h"

// forward declarations
class materialTriangles;

class fence
{
public:
	void addPost(materialTriangles *tri, int triangle, float(&xyz)[3], float(&normal)[3], bool connectToNearestEdge, bool openEnd = false);
	int getPostData(std::vector<Vec3f> &positions,std::vector<Vec3f> &normals,std::vector<int> &triangles,std::vector<float> &uv,bool &edgeStart,bool &edgeEnd, bool& startOpen, bool& endOpen);
	void updatePosts(const std::vector<Vec3f> &positions, const std::vector<Vec3f> &normals);
	void deleteLastPost();
	void selectPost(int postNumber);
	inline void getSpherePos(int postNumber, Vec3f &xyz) { xyz = _posts[postNumber].spherePos; }
	inline void getPostPos(int postNumber, Vec3f &xyz) { xyz = _posts[postNumber].xyz; }
	inline void getPostNormal(int postNumber, Vec3f &nrm) { nrm = _posts[postNumber].nrm; }
	inline int numberOfPosts() { return (int)_posts.size(); }
	void setSpherePos(int postNumber, Vec3f& xyz);
	void setGl3wGraphics(gl3wGraphics *gl3w) {
		_gl3w = gl3w; _shapes = gl3w->getShapes(); _glm = gl3w->getGLmatrices();
		if (fence::_fenceSize < 10000.0f) _initialized = true;
	}
	void setFenceSize(float size) {fence::_fenceSize=size; if(_gl3w !=NULL) _initialized=true;}
	void clear();	// deletes this fence COURT - ?nuke as no longer used
	bool isInitialized()	{return _initialized;}
	fence();
	~fence();

private:
	struct fencePost{
		bool connectToEdge;
		bool openEnd;
		Vec3f xyz;
		Vec3f nrm;
		int triangle;
		float uv[2];
		std::shared_ptr<sceneNode>  cylinderShape;
		std::shared_ptr<sceneNode> sphereShape;
		Vec3f spherePos;
	};
	std::shared_ptr<sceneNode> _wall;
	std::vector<fencePost> _posts;
	gl3wGraphics *_gl3w;
	shapes *_shapes;
	GLmatrices *_glm;
	std::vector<Vec3f> _xyz, _norms;
	bool _initialized;
	static float _fenceSize;
	static GLfloat _selectedColor[4],_unselectedColor[4];

	void displayRemoveWall();
	void computeLocalBounds();

};

#endif // #ifndef __FENCE_H__
