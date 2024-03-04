// File: hooksConstraint.h
// Author: Court Cutting
// Date: 1/17/19
// Purpose: Class for handling tissue hooks and associated graphics using projective dynamics constraints rather than forces

#ifndef __HOOKS_H__
#define __HOOKS_H__

#include <map>
#include <memory>
#include "shapes.h"
#include "Vec3f.h"

// #include "deepCut.h"
#include "skinCutUndermineTets.h"

#include "pdTetPhysics.h"

// forward declarations
class materialTriangles;
class GLmatrices;
class shapes;
class vnBccTetrahedra;

class hookConstraint
{
public:
	inline void setShape(std::shared_ptr<sceneNode> &shape) {_shape=shape;}
	inline std::shared_ptr<sceneNode>  getShape() {return _shape;}
	hookConstraint() : _shape(nullptr), _constraintId(-1) {}
	~hookConstraint() {}
protected:
	materialTriangles *_tri;
	int triangle;
	float uv[2];
	Vec3f xyz, _selectPosition;  // xyz is current hook position
	bool _selected, _strong;
	std::shared_ptr<sceneNode> _shape;
	int _constraintId;
	friend class hooks;
};

class hooks
{
public:
	int addHook(materialTriangles *tri, int triangle, float(&uv)[2], bool tiny = false);
	bool getHookPosition(unsigned int hookNumber, float (&hookPos)[3]);
	bool setHookPosition(unsigned int hookNumber, float (&hookPos)[3]);
	bool getSelectPosition(unsigned int hookNumber, float(&selectPos)[3]);
	bool updateHookPhysics();  // must do after a physics lattice change
	int getNumberOfHooks() {return (int)_hooks.size();}
	float getHookSize() {return _hookSize;}
	void setHookSize(float size) {_hookSize=size;}
	void selectHook(int hookNumber);
	void deleteHook(int hookNumber);
	static inline void setSpringConstant(float k) { _springConstant = k; }
	static inline float getSpringConstant() { return _springConstant; }
	void setShapes(shapes *shps) {_shapes=shps;}
	void setGLmatrices(GLmatrices *GLm) {_glm=GLm;}
	void setPhysicsLattice(pdTetPhysics *pdtp) { _ptp = pdtp; }
	void setVnBccTetrahedra(vnBccTetrahedra *vnt) { _vnt = vnt; }
//	void setIncisions(deepCut *dc) { _deepCut = dc; }
	bool empty() { return _hooks.empty(); }
	hooks();
	~hooks();

private:
	GLmatrices *_glm;
	pdTetPhysics *_ptp;
	vnBccTetrahedra *_vnt;

//	deepCut *_deepCut;
	skinCutUndermineTets* _deepCut;  // replace later with above

	shapes *_shapes;
	typedef std::map<unsigned int, hookConstraint> HOOKMAP;
	HOOKMAP _hooks;
	unsigned int _hookNow;
	int _selectedHook;
	static float _springConstant;
	static float _hookSize;
	static GLfloat _selectedColor[4], _unselectedColor[4];  // , _insideSkullColor[4];

//	int parametricMTtriangleTet(const int mtTriangle, const float(&uv)[2], Vec3f& gridLocus);
};

#endif	// __HOOKS_H__
