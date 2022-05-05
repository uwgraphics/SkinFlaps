// File: sutures.h
// Author: Court Cutting
// Date: 2/20/12
// Update: 6/10/2019 for bcc tetrahedra
// Purpose: Class for handling sutures and associated graphics

#ifndef __SUTURES_H__
#define __SUTURES_H__

#include <map>
#include <memory>
#include "Vec3f.h"
#include "pdTetPhysics.h"
#include "shapes.h"

// forward declarations
class materialTriangles;
class GLmatrices;
class lightsShaders;
class vnBccTetrahedra;
class surgicalActions;
class deepCut;

class sutureTets
{
public:
	inline void setSphereShape(std::shared_ptr<sceneNode> &sphere) {_sphereShape=sphere;}
	inline std::shared_ptr<sceneNode> getSphereShape() {return _sphereShape;}
	inline void setCylinderShape(std::shared_ptr<sceneNode> &cylinder) {_cylinderShape= cylinder;}
	inline std::shared_ptr<sceneNode> getCylinderShape() {return _cylinderShape;}
	sutureTets() : _tri(NULL), _type(0), _sphereShape(nullptr), _cylinderShape(nullptr) {
		_edges[0] = -1; _tris[0] = -1; _params[0] = -1.0f;
		_edges[1] = -1; _tris[1] = -1; _params[1] = -1.0f; }
	~sutureTets() {}
protected:
	materialTriangles *_tri;
	int _edges[2], _tris[2], _tetIdx[2];
	float _params[2];
	GLfloat _v1[3];
	bool _selected;
	int _type;  // 0=user entered without link to previous, 1=user entered with link to previous, 2=programatically created in suture strip between a type 0-1 pair
	std::shared_ptr<sceneNode> _sphereShape;
	std::shared_ptr<sceneNode> _cylinderShape;
	int _constraintId;  // index of the suture constraint
	Vec3f _baryWeights[2];
	friend class sutures;
};

class sutures
{
public:
	int addUserSuture(materialTriangles *tri, int triangle0, int edge0, float param0);
	int setSecondEdge(int sutureNumber, materialTriangles *tri, int triangle, int edge, float param);  // return 0= normal, 1=one sided suture, 2=one tet suture, 3=different soft bodies
	void setSecondVertexPosition(int sutureNumber, float *position);	// updates this sutures second graphics. position[3]
	float* getSecondVertexPosition(int sutureNumber);
	void updateSuturePhysics();  // must call after a physics lattice topo change
	int deleteSuture(int sutureNumber);  // return 0 is user suture deletion, else autoSuture deletion returning number of the linked suture that generated it.
	void updateSutureGraphics();
	int getNumberOfSutures() {return (int)_sutures.size();}
	void nearestSkinIncisionEdge(const float triUv[2], int &triangle, int &edge, float &param);
	float getSutureSize() {return _sutureSize;}
	void setSutureSize(float size) {_sutureSize=size;}
	int firstVertexMaterial(int sutureNumber);
	bool getEdgeAttachment(const int sutureNumber, const bool firstPoint, materialTriangles *(&tri), int &triangle, int &edge, float &param);
	void selectSuture(int sutureNumber);
	inline bool isLinked(int sutureNumber) { return _sutures[sutureNumber]._type == 1; }
	inline void setLinked(int sutureNumber, bool link) { _sutures[sutureNumber]._type = (link ? 1 : 0); }
	int previousUserSuture(int sutureNumber);
	void laySutureLine(int suture2);  // creates a line of sutures between 2 input sutures
	inline void setShapes(shapes *shps) { _shapes = shps; }
	inline void setGLmatrices(GLmatrices *GLm) {_glm=GLm;}
	inline void setPhysicsLattice(pdTetPhysics *ptp) { _ptp = ptp; }
	inline void setVnBccTetrahedra(vnBccTetrahedra *vbt) { _vbt = vbt; }
	inline void setDeepCut(deepCut *dc){ _dc = dc; }
	inline void setSurgicalActions(surgicalActions* sa) { _surgAct = sa; }
	inline static void setAutoSutureSpacing(float spacing) { _sutureSpanGap = spacing; }
	inline bool empty() { return _sutures.empty(); }
	void clear()	{_sutures.clear(); }

	inline int userToBaseSutureNumber(const int userSutureNumber) {
		auto us = _userSutures.find(userSutureNumber);
		if (us == _userSutures.end())
			return -1;
		return us->second;
	}

	inline int baseToUserSutureNumber(const int baseSutureNumber) {
		auto us = _userSutures.begin();
		while (us != _userSutures.end()) {
			if (us->second == baseSutureNumber)
				return us->first;
			++us;
		}
		return -1;
	}

	sutures();
	~sutures();

private:
	deepCut *_dc;
	pdTetPhysics *_ptp;
	vnBccTetrahedra *_vbt;
	surgicalActions* _surgAct;
	GLmatrices *_glm;
	shapes *_shapes;
	typedef std::map<unsigned int, sutureTets> SUTUREMAP;
	SUTUREMAP _sutures;
	std::map<int, int> _userSutures;
	unsigned int _sutureNow;  // must keep unique number and not decrement with deletions.
	unsigned int _userSutureNext;  // must keep unique number and not decrement with deletions.
	static float _sutureSpanGap;
	static float _sutureSize;
	static GLfloat _selectedColor[4], _unselectedColor[4], _userColor[4];
	int addSuture(materialTriangles *tri, int triangle0, int edge0, float param0);
};

#endif	// __SUTURES_H__
