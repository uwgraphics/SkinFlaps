// File: fences.cpp
// Author: Court Cutting, MD
// Date: June 7, 2012
// Purpose: User interface for creating a fence on a glslTriangle object.
//     This will be used to specify a desired incision line.

#include "shapes.h"
#include "GLmatrices.h"
#include "Vec3f.h"
#include <stdio.h>
#include <assert.h>
#include "boundingBox.h"
#include "materialTriangles.h"
#include "fence.h"

float fence::_fenceSize=10000.0f;
GLfloat fence::_selectedColor[]={1.0f,1.0f,0.0f,1.0f};
GLfloat fence::_unselectedColor[]={0.0f,1.0f,0.0f,1.0f};

void fence::clear()	// deletes current fence
{
	std::vector<fencePost>::iterator pit;
	for(pit=_posts.begin(); pit!=_posts.end(); ++pit)	{
		_gl3w->deleteSceneNode(pit->cylinderShape);
		_gl3w->deleteSceneNode(pit->sphereShape);
	}
	_posts.clear();
	_xyz.clear();
	_norms.clear();
	if(_wall!= nullptr)
		_wall->visible = false;  // don't delete if already made.  Will reuse.
}

void fence::updatePosts(const std::vector<Vec3f> &positions, const std::vector<Vec3f> &normals)
{  // this routine only called by the deep cutter.  Make posts longer to ease user alteration of normals
	if (positions.size() != _posts.size())
		return;
	_xyz.assign(positions.size() * 2, Vec3f());
	for (int n = (int)_posts.size(), i = 0; i < n; ++i) {
		_posts[i].xyz = positions[i];
		_posts[i].nrm = normals[i];
		_posts[i].nrm.normalize();
		std::shared_ptr<sceneNode> sh = _posts[i].cylinderShape;
		GLfloat *mm = sh->getModelViewMatrix();
		loadIdentity4x4(mm);
		scaleMatrix4x4(mm, _fenceSize*0.1f, _fenceSize*0.1f, _fenceSize * 3.0f);
		Vec3f vz(0.0f, 0.0f, 1.0f), vn;
		float angle = acos(_posts[i].nrm*vz);
		vn = vz^_posts[i].nrm;
		axisAngleRotateMatrix4x4(mm, vn._v, angle);
		translateMatrix4x4(mm, positions[i].x(), positions[i].y(), positions[i].z());
		sh = _posts[i].sphereShape;
		mm = sh->getModelViewMatrix();
		loadIdentity4x4(mm);
		scaleMatrix4x4(mm, _fenceSize*0.5f, _fenceSize*0.5f, _fenceSize*0.5f);
		vn.set(_posts[i].nrm);
		_xyz[i*2] = positions[i] - vn * 2.5f * _fenceSize;
		_xyz[i*2 + 1] = positions[i] + vn * 2.5f * _fenceSize;
		vn *= _fenceSize * 3.0f;
		vn += positions[i];
		translateMatrix4x4(mm, vn._v[0], vn._v[1], vn._v[2]);
	}
	displayRemoveWall();  // now makeWall always true
}

void fence::deleteLastPost()
{
	if (_posts.empty())
		return;
	std::vector<fencePost>::reverse_iterator pit = _posts.rbegin();
	_gl3w->deleteSceneNode(pit->cylinderShape);
	pit->cylinderShape.reset();
	_gl3w->deleteSceneNode(pit->sphereShape);
	pit->sphereShape = nullptr;
	_posts.pop_back();
	_xyz.pop_back();
	_xyz.pop_back();
	displayRemoveWall();
}

void fence::addPost(materialTriangles *tri, int triangle, float(&xyz)[3], float(&normal)[3], bool connectToNearestEdge, bool openEnd)
{
	char name[8];
	sprintf(name,"P_%d",(int)_posts.size());
	_posts.push_back(fencePost());
	if(connectToNearestEdge && _posts.empty())	{
		int edge;
		float param;
		tri->getNearestHardEdge(xyz,triangle,edge,param);
		_posts.back().triangle = triangle;
		if(edge<1)	{
			_posts.back().uv[0] = param;
			_posts.back().uv[1] = 0.0f;
		}
		else if(edge<2)	{
			_posts.back().uv[0] = 1.0f - param;
			_posts.back().uv[1] = param;
		}
		else	{
			_posts.back().uv[1] = 1.0f - param;
			_posts.back().uv[0] = 0.0f;
		}
	}
	else	{
		_posts.back().triangle = triangle;
		tri->getBarycentricProjection(triangle,xyz,_posts.back().uv);
	}
	_posts.back().openEnd = openEnd;
	_posts.back().xyz[0]=xyz[0]; _posts.back().xyz[1]=xyz[1]; _posts.back().xyz[2]=xyz[2];
	_posts.back().nrm[0]=normal[0]; _posts.back().nrm[1]=normal[1]; _posts.back().nrm[2]=normal[2];
	_posts.back().connectToEdge = connectToNearestEdge;
	 std::shared_ptr<sceneNode> sh = _posts.back().cylinderShape = _shapes->addShape(sceneNode::nodeType::CYLINDER,name);
	GLfloat *mm = sh->getModelViewMatrix();
	loadIdentity4x4(mm);
	sh->setColor(_selectedColor);
	scaleMatrix4x4(mm,_fenceSize*0.1f,_fenceSize*0.1f,_fenceSize*2.0f);
	Vec3f vz(0.0f,0.0f,1.0f),vn(normal[0], normal[1], normal[2]);
	vn.normalize();
	float angle = acos(vn*vz);
	vz = vz^vn;
	axisAngleRotateMatrix4x4(mm, vz._v, angle);
	translateMatrix4x4(mm,xyz[0],xyz[1],xyz[2]);

	sprintf(name, "NP_%d", (int)_posts.size()-1);
	sh = _posts.back().sphereShape = _shapes->addShape(sceneNode::nodeType::SPHERE, name);
	mm = sh->getModelViewMatrix();
	loadIdentity4x4(mm);
	sh->setColor(_selectedColor);
	scaleMatrix4x4(mm, _fenceSize*0.5f, _fenceSize*0.5f, _fenceSize*0.5f);
	//	vn already done
	vn *= _fenceSize * 2.0f;
	translateMatrix4x4(mm, xyz[0] + vn.X, xyz[1] + vn.Y, xyz[2] + vn.Z);
	_posts.back().spherePos = Vec3f(xyz) + vn;

	vn.set(normal);
	vn *= _fenceSize * 2;
	_xyz.push_back(Vec3f(xyz) - vn * 0.75f);
	_xyz.push_back(Vec3f(xyz) + vn * 0.75f);
	displayRemoveWall();  // now makeWall always true
}

void fence::setSpherePos(int postNumber, Vec3f& xyz){
	fencePost* fp = &_posts[postNumber];
	GLfloat* mm = fp->sphereShape->getModelViewMatrix();
	loadIdentity4x4(mm);
	fp->sphereShape->setColor(_selectedColor);
	scaleMatrix4x4(mm, _fenceSize * 0.5f, _fenceSize * 0.5f, _fenceSize * 0.5f);
	Vec3f vn = xyz - fp->xyz;
	vn.normalize();
	fp->nrm = vn;
	vn *= _fenceSize * 2.0f;
	vn += fp->xyz;
	translateMatrix4x4(mm, vn.X, vn.Y, vn.Z);
	fp->spherePos = vn;

	mm = fp->cylinderShape->getModelViewMatrix();
	loadIdentity4x4(mm);
	fp->cylinderShape->setColor(_selectedColor);
	scaleMatrix4x4(mm, _fenceSize * 0.1f, _fenceSize * 0.1f, _fenceSize * 2.0f);
	Vec3f vz(0.0f, 0.0f, 1.0f);
	vn -= fp->xyz;
	vn.normalize();
	float angle = acos(vn * vz);
	vz = vz ^ vn;
	axisAngleRotateMatrix4x4(mm, vz._v, angle);
	translateMatrix4x4(mm, fp->xyz[0], fp->xyz[1], fp->xyz[2]);

	_xyz[postNumber*2] = fp->xyz - vn * 0.75f;
	_xyz[postNumber * 2 + 1] = fp->xyz + vn * 0.75f;
	displayRemoveWall();
}


void fence::displayRemoveWall() {
	if (_xyz.size() > 3) {  // valid wall to display
		if (!_wall) {
			_wall = std::make_shared<sceneNode>();
			_wall->setType(sceneNode::nodeType::TRISTRIP);
			GLfloat color[4] = { 1.0f, 0.4f, 0.0f, 1.0f };
			_wall->setColor(color);
			GLfloat* mm = _wall->getModelViewMatrix();
			loadIdentity4x4(mm);
			GLuint program = _gl3w->getLightsShaders()->getOrCreateColorProgram();
			_wall->setGlslProgramNumber(program);
			_wall->setColorLocation(glGetUniformLocation(program, "objectColor"));
			if (!_wall->vertexArrayBufferObject)
				glGenVertexArrays(1, &_wall->vertexArrayBufferObject);
			if (_wall->bufferObjects.size() != 2) {
				_wall->bufferObjects.assign(2, 0);
				glGenBuffers(2, &_wall->bufferObjects[0]);
			}
			// Create the master vertex array object
			glBindVertexArray(_wall->vertexArrayBufferObject);
			// Vertex data
			glBindBuffer(GL_ARRAY_BUFFER, _wall->bufferObjects[0]);	// VERTEX DATA
			glEnableVertexAttribArray(0);
			glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 0, 0);
			// Normal data
			glBindBuffer(GL_ARRAY_BUFFER, _wall->bufferObjects[1]);	// NORMAL_DATA
			glEnableVertexAttribArray(1);
			glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
			// Unbind to anybody
			glBindVertexArray(0);
			_gl3w->addSceneNode(_wall);
		}
		_wall->visible = true;
		int n = (int)_xyz.size();
		_norms.clear();
		_norms.reserve(n/2);
		Vec3f v0, v1, v2, N0, N1, lastN(0.0f, 0.0f, 0.0f);
		for (int i = 2; i < n; i += 2) {
			v0 = _xyz[i - 2];
			v1 = _xyz[i - 1];
			v2 = _xyz[i];
			N0 = (v1 - v0) ^ (v2 - v0);
			N0.normalize();
			v1 = _xyz[i + 1];
			N1 = (v0 - v2) ^ (v1 - v2);
			N1.normalize();
			if (i < 3) {
				_norms.push_back(N0);
			}
			else {
				lastN += N0;
				lastN.normalize();
				_norms.push_back(lastN);
			}
			if (i == n - 2) {
				_norms.push_back(N1);
			}
			lastN = N1;
		}
		std::vector<GLfloat> grPos;
		grPos.reserve((n/2-1) * 400 + 16);
		auto pushBackV = [&](Vec3f &v) {
			grPos.push_back(v.X);
			grPos.push_back(v.Y);
			grPos.push_back(v.Z);
			grPos.push_back(1.0f);
		};
		std::vector<Vec3f> grNorm;
		grNorm.reserve(grPos.capacity()>>2);
		if (_posts.front().openEnd) {
			v0 = _xyz[0];
			v1 = _xyz[1];
			v2 = v0 + v1 - _xyz[2] -  _xyz[3];
			v2.normalize();
			v2 *= _fenceSize * 2.2f;
			v0 += v2;
			v1 += v2;
			pushBackV(v0);
			pushBackV(v1);
			grNorm.push_back(_norms[0]);
			grNorm.push_back(_norms[0]);
		}
		for (int i = 2; i < n; i+=2) {
			for (int j = 0; j < 50; ++j) {
				v0 = _xyz[i-2]*((49.0f - j) / 49.0f);
				v0 += _xyz[i] * (j / 49.0f);
				pushBackV(v0);
				v0 = _xyz[i - 1] * ((49.0f - j) / 49.0f);
				v0 += _xyz[i + 1] * (j / 49.0f);
				pushBackV(v0);
				N0 = _norms[i/2 - 1] * ((49.0f - j) / 49.0f);
				N0 += _norms[i/2] * (j / 49.0f);
				grNorm.push_back(N0);
				grNorm.push_back(N0);
			}
		}
		if (_posts.back().openEnd) {
			v0 = _xyz[n-2];
			v1 = _xyz[n-1];
			v2 = v0 + v1 - _xyz[n-3] - _xyz[n-4];
			v2.normalize();
			v2 *= _fenceSize * 2.2f;
			v0 += v2;
			v1 += v2;
			pushBackV(v0);
			pushBackV(v1);
			grNorm.push_back(_norms.back());
			grNorm.push_back(_norms.back());
		}
		// Vertex and normal data
		glBindBuffer(GL_ARRAY_BUFFER, _wall->bufferObjects[0]);	// VERTEX_DATA
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * grPos.size(), &(grPos[0]), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, _wall->bufferObjects[1]);	// NORMAL_DATA
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * grNorm.size() * 3, &(grNorm[0]), GL_STATIC_DRAW);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		_wall->elementArraySize = (GLsizei)(grPos.size() >> 2);
	}
	else {
		if (_wall)
			_wall->visible = false;
	}
}

void fence::computeLocalBounds() {
	boundingBox<float> bb;
	bb.Empty_Box();
	for (int n = (int)_xyz.size(), i = 0; i < n; ++i) 
		bb.Enlarge_To_Include_Point((const float(&)[3])_xyz[i]);
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
	_wall->setLocalBounds(lc, radius);
}

void fence::selectPost(int postNumber)
{
	for (int n = (int)_posts.size(), i = 0; i < n; ++i) {
		if (i == postNumber) {
			_posts[i].cylinderShape->setColor(_selectedColor);
			_posts[i].sphereShape->setColor(_selectedColor);
		}
		else {
			_posts[i].cylinderShape->setColor(_unselectedColor);
			_posts[i].sphereShape->setColor(_unselectedColor);
		}
	}
}

int fence::getPostData(std::vector<Vec3f> &positions,std::vector<Vec3f> &normals,std::vector<int> &triangles,std::vector<float> &uv,bool &edgeStart,bool &edgeEnd, bool& startOpen, bool& endOpen)
{
	positions.clear();
	normals.clear();
	triangles.clear();
	uv.clear();
	int n=(int)_posts.size();
	positions.reserve(n*3);
	normals.reserve(n*3);
	triangles.reserve(n);
	uv.reserve(n<<1);
	edgeStart = _posts.front().connectToEdge;
	startOpen = _posts.front().openEnd;
	edgeEnd = _posts.back().connectToEdge;
	endOpen = _posts.back().openEnd;
	std::vector<fencePost>::iterator fit;
	for(fit=_posts.begin(); fit!=_posts.end(); ++fit)	{
		positions.push_back(fit->xyz);
		normals.push_back(fit->nrm);
		triangles.push_back(fit->triangle);
		uv.push_back(fit->uv[0]);
		uv.push_back(fit->uv[1]);
	}
	return n;
}

fence::fence() :_initialized(false)
{
	_gl3w =NULL;
	_wall = nullptr;
	_xyz.clear();
	_norms.clear();
}

fence::~fence()
{
	if (!_wall)
		return;
	if (!_wall->bufferObjects.empty()) {
		glDeleteBuffers((GLsizei)_wall->bufferObjects.size(), &_wall->bufferObjects[0]);
		_wall->bufferObjects.clear();
	}
	if (_wall->vertexArrayBufferObject > 0) {
		glDeleteVertexArrays(1, &_wall->vertexArrayBufferObject);
		_wall->vertexArrayBufferObject = 0;
	}
}

