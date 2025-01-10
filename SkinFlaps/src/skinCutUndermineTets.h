////////////////////////////////////////////////////////////////////////////
// File: skinCutUndermineTets.h
// Author: Court Cutting
// Date: 9/10/2017
// Purpose: Interactive skin-mucosal incision class which operates in spatial coordinates. Requires setup with a triangulated subdermal surface corresponding to skin-mucosal surface of
//    closed manifold surface surrounding elastic solid modelled as virtual noded cubes.  Both superficial and deep surfaces in material coordinates on setup.
//    Allows user to make depth and normal agnostic incisions on skin-mucosa in spatial coordinates.  Also allows undermining of flaps along an incision line
//    away from the deep bed.  Updates both the virtual noded cubes and the surface model.
////////////////////////////////////////////////////////////////////////////

#ifndef __SKIN_CUT_UNDERMINE_TETS__
#define __SKIN_CUT_UNDERMINE_TETS__

#include <vector>
#include <array>
#include <list>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include "Vec3f.h"
#include "Vec3d.h"
#include "Vec2d.h"
#include "Vec2f.h"
#include "Mat3x3d.h"

#pragma warning (disable : 4267)

// forward declarations
class materialTriangles;
class vnBccTetrahedra;
class gl3wGraphics;

class skinCutUndermineTets
{
public:

	void setGl3wGraphics(gl3wGraphics *gl3w) { _gl3w = gl3w; }  // for debug - nuke later
	bool skinCut(std::vector<Vec3f> &topCutPoints, std::vector<Vec3f> &topNormals, bool startOpen, bool endOpen);  // history version
	float closestSkinIncisionPoint(const Vec3f xyz, int& triangle, int& edge, float& param);  // Input xyz, returns all 4
	bool addUndermineTriangle(const int triangle, const int undermineMaterial, bool incisionConnect);
	void undermineSkin();  // completes a user specified undermine
	void clearCurrentUndermine(const int underminedTissue);
	bool triangleUndermined(int triangle);
	void excise(const int triangle);
	bool physicsRecutRequired(){ return _solidRecutRequired; }
	bool setDeepBed(materialTriangles *mt, const std::string &deepBedPath, vnBccTetrahedra *activeVnt);
	inline static void setVnBccTetrahedra(vnBccTetrahedra *activeVnt) { _vbt = activeVnt;  }
	inline void setMaterialTriangles(materialTriangles *mt) { _mt = mt; }
	inline materialTriangles* getMaterialTriangles(){ return _mt; }
	skinCutUndermineTets();
	skinCutUndermineTets(const skinCutUndermineTets&) = delete;
	skinCutUndermineTets& operator=(const skinCutUndermineTets&) = delete;
	~skinCutUndermineTets();

protected:
	static gl3wGraphics *_gl3w;
	static materialTriangles *_mt;  // embedded surface
	static vnBccTetrahedra *_vbt;  // above surface embedded in these current cut tets.
	struct deepPoint{
		Vec3f gridLocus;
		int deepMtVertex;  // get tet & barycentrics from here when > -1
	};
	static std::unordered_map<int, deepPoint> _deepBed;
	// next is data of previously undermined triangles. _prevUnd2 are all previouslu undermined top triangles. Rest are previous undermines containing a non-duplicated deep vertex.  All are sorted vectors except _prevBot5.
	// filled before each undermine by collectOldUndermineData()
	std::vector<int> _prevUnd2, _prevEdge3;
	std::vector<std::pair<int, int> > _prevBot4;  // first is the mat 4 bottom triangle, second is the mat 2 top triangle that created it either by an incision or an undermine. Remember vert order reversed.
	std::map<int, std::vector<int> > _prevBedSingles;
	std::vector<bool> _trisUnderminedNow;  // triangles turned on in current undermine
	std::vector<int> _inExCisionTriangles;  // material 2 triangles at the edge of an incision or excision
	std::vector<int> _periostealCutEdgeTriangles;  // periosteal triangles at the edge of a deepCut
	std::unordered_map<int, Vec3f> _collisionSpokes, _deepSpokesNow;
	int _prevUndermineTriangle, _firstTopVertex;
	bool _startOpen, _endOpen, _solidRecutRequired;

	int deepPointTetWeight(const std::unordered_map<int, deepPoint>::iterator &dit, Vec3f &baryWeight);  // return deepPoint tet number and baryweight from its grid locus
	int addSurfaceVertex(const int tet, const Vec3f &gridLocus);  // ? nuke
	int createDeepBedVertex(std::unordered_map<int, deepPoint>::iterator &dit);
	int addTinEdgeVertex(const Vec3f &closePoint, const Vec3f &nextConnectedPoint);
	int TinSub(const int edgeTriangle, const float edgeParam);
	int flapBottomTet(const int topVertex, const Vec3f &bottomGridLocus);
	bool topDeepSplit(std::vector<int> &topV, std::vector<int> &deepV, bool frontSplit, bool backSplit);
	bool topDeepSplit_Sub(std::list<int> &topVerts, std::list<int> &deepVerts, bool frontSplit, bool backSplit);
	bool planeCutSurfaceLine(const int startTopV, const int endTopV, const int startDeepV, const int endDeepV, std::list<int> &newTopVerts, std::list<int> &newDeepVerts);
	void createFlapTopBottomVertices(const int topTriangle, float(&uv)[2], int &topVertex, int &bottomVertex);
	void flapSurfaceSplitter(const int startVertex, const int endVertex, std::list<int> &vertexCutLine, std::vector<int> &oppositeVertices);
	// next set are for undermining
	bool trianglePath(const int triStart, const int endTriangle, const int searchMaterial, std::vector<int> &triPath);
	bool closeUndermineHoles(std::vector<int> &trianglePath, const int undermineMaterial);
	void showPriorUndermine(int priorTriangle);
	void collectOldUndermineData();
	int cloneTexture(int textureIndex);
	bool testIncisionsDeepBed();  // Looks for intersections of the deep bed with the deep surface of the object.  For debugging.

	// What follows describes an incision convention which must be adhered to throughout these routines to facilitate not only this code
	// but other interrogations of cut skin edges as well:
	// Two immediate neighbor surface points and their corresponding deep points form a quad which will be converted into 2 triangles on either side.
	// This quad will be duplicated with one side bounding one side of the incision solid and the other side bounding the opposite solid side.
	// They must be triangulated such that both sides are concave, or at worst perfectly planar in material coordinates.  This is done so there is no collision
	// in the tetrahedron cutter.  As a result the zero edge of a material 3 triangle will always be on top (adjacent to material 2) or on the bottom (adjacent to material 4, 5 or 6).
	// Determining which tesselation was done requires checking number of triangle adjacent to edge 1 or 2. Bottom tri number is always top number + 1.
};

#endif // __SKIN_CUT_UNDERMINE_TETS__