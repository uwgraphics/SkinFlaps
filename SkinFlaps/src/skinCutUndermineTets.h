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
	bool skinCut(std::vector<Vec3f> &topCutPoints, std::vector<Vec3f> &topNormals, bool startOpen, bool endOpen);
	bool addUndermineTriangle(const int triangle, const int undermineMaterial, bool incisionConnect);
	void undermineSkin();  // completes a user specified undermine
	void clearCurrentUndermine(const int underminedTissue);
	void excise(const long triangle);
	long parametricMTedgeTet(const int triangle, const int edge, const float param, Vec3f &gridLocus);
	long parametricMTtriangleTet(const int mtTriangle, const float(&uv)[2], Vec3f &gridLocus);
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
		long deepMtVertex;  // get tet & barycentrics from here when > -1
	};
	static std::unordered_map<long, deepPoint> _deepBed;
	// next is data of previously undermined triangles. _prevUnd2 are all previouslu undermined top triangles. Rest are previous undermines containing a non-duplicated deep vertex.  All are sorted vectors except _prevBot5.
	// filled before each undermine by collectOldUndermineData()
	std::vector<int> _prevUnd2, _prevBot4, _prevEdge3;
	std::map<long, std::vector<long> > _prevBedSingles;
	std::vector<bool> _trisUnderminedNow;  // triangles turned on in current undermine
	std::vector<long> _inExCisionTriangles;  // material 2 triangles at the edge of an incision or excision
	std::vector<long> _periostealCutEdgeTriangles;  // periosteal triangles at the edge of a deepCut
	std::unordered_map<int, Vec3f> _collisionSpokes, _deepSpokesNow;
	long _prevUndermineTriangle, _firstTopVertex;
	bool _startOpen, _endOpen, _solidRecutRequired;

	long deepPointTetWeight(const std::unordered_map<long, deepPoint>::iterator &dit, Vec3f &baryWeight);  // return deepPoint tet number and baryweight from its grid locus
	long addSurfaceVertex(const long tet, const Vec3f &gridLocus);  // ? nuke
	long createDeepBedVertex(std::unordered_map<long, deepPoint>::iterator &dit);
	long addTinEdgeVertex(const Vec3f &closePoint, const Vec3f &nextConnectedPoint);
	long TinSub(const long edgeTriangle, const float edgeParam);
	long flapBottomTet(const long topTet, const Vec3f &bottomGridLocus);
	bool topDeepSplit(std::vector<long> &topV, std::vector<long> &deepV, bool frontSplit, bool backSplit);
	bool topDeepSplit_Sub(std::list<long> &topVerts, std::list<long> &deepVerts, bool frontSplit, bool backSplit);
	bool planeCutSurfaceLine(const long startTopV, const long endTopV, const long startDeepV, const long endDeepV, std::list<long> &newTopVerts, std::list<long> &newDeepVerts);
	void createFlapTopBottomVertices(const long topTriangle, float(&uv)[2], long &topVertex, long &bottomVertex);
	void flapSurfaceSplitter(const long startVertex, const long endVertex, std::list<long> &vertexCutLine, std::vector<long> &oppositeVertices);
	// next set are for undermining
	bool trianglePath(const long triStart, const long endTriangle, const int searchMaterial, std::vector<long> &triPath);
	bool closeUndermineHoles(std::vector<long> &trianglePath, const int undermineMaterial);
	void showPriorUndermine(long priorTriangle);
	void collectOldUndermineData();
	long cloneTexture(long textureIndex);
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