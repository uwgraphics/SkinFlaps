///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// File: deepCut.h
// Author: Court Cutting MD
// Date: 1/27/2020
// Purpose: Yet another version of deepCut tool.  This version is designed to deal with transiently inverted tets which arise in the physics.
//    Previous versions found the cut tets in spatial coords.  Unfortunately when inverted tets were righted this pulled parts of the triangulated
//    surface inside out.  The resulting occasional self intersections crashed the tet cutter.  For this reason the interior of this deep cutter
//    will have its vertices found in material coordinates so no inversions are possible.  This will result in the following loss of realism:
//    In previous versions pulling up the central part of a surface while holding the perifery down would result in one cut side being concave
//    and the other convex as a result of a single planar cut.  This version will lack that reality and both sides will be relatively planar.
//    This version will also dispense with the deep side of incisions being planar when cutting material 2.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __DEEP_CUT__
#define __DEEP_CUT__

#include <vector>
#include <array>
#include "Vec3f.h"
#include "materialTriangles.h"
#include "skinCutUndermineTets.h"

#pragma warning (disable : 4267)

// forward declarations
class vnBccTetrahedra;
class fence;
struct rayTriangleIntersect;

class deepCut : public skinCutUndermineTets
{
public:
	void setGl3wGraphics(gl3wGraphics *gl3w) { _gl3w = gl3w; }  // for debug - nuke later

	bool inputCorrectFence(fence* fp);
	int addDeepPost(const int triangle, const float(&uv)[2], const Vec3d& rayDirection, bool closedEnd);
	inline void popLastDeepPost() { if(!_deepPosts.empty()) _deepPosts.pop_back(); }
	inline int numberOfDeepPosts() { return (int)_deepPosts.size(); }
	bool preventPreviousCrossover(const int postNum);
	void getDeepPosts(std::vector<Vec3f>& xyz, std::vector<Vec3f>& nrm);
	bool cutDeep();  // data already loaded in _deepPosts in this updated version
	void clearDeepCutter(){_deepPosts.clear();}
	long addPeriostealUndermineTriangle(const int topTriangle, const Vec3f &linePickDirection, const bool incisionConnect);  // can only follow a deepCut through periosteum.
	deepCut() { _deepXyz.clear(); _deepPosts.clear(); }
	deepCut(const deepCut&) = delete;
	deepCut& operator=(const deepCut&) = delete;
	~deepCut(){}

protected:
	bool _startOpen, _endOpen;
	std::vector<std::array<long, 3> > _firstSideTriangles;  // new triangles produced on first side of the cut. in sort order so can do binary search, but not usually consecutive.
	struct surfaceCutLine {  // a cut line intersecting a cutting triangle with the embedded surface
		int rtiPostTo;
		int rtiIndexTo;
		double lowestV;  // smallest v value of the bilinear surface in an interpost scl
		std::list<long> deepVerts;  // sequential deep mtTriangle vertices for this intersect line
	};
	struct rayTriangleIntersect {
		int triangle;
		double uv[2];
		Vec3d intersect;
		double rayParam;
		int postNum;
		int rayIndex;
		bool solidDown;
		int mat2Vert;
		int deepVert;
		surfaceCutLine scl;
		std::vector<int> dcl;  // deep cut line vertices through interior solid. Always top to bottom.
	};
	struct bilinearPatch {  // patch data for nVidia intersection routine
		Vec3d e10;
		Vec3d e01;
		Vec3d e11;
		Vec3d e00;
		Vec3d qn;
		Vec3d P00;
		Vec3d P10;
		Vec3d P01;
		Vec3d P11;
	};
	struct deepPost {
		bool closedEnd;
		std::vector<rayTriangleIntersect> triIntersects;
		Vec3d rayDirection;
		bilinearPatch bl;
		std::vector<materialTriangles::matTriangle> quadTriangles;
	};
	std::vector< deepPost> _deepPosts;

	std::vector<Vec3d> _deepXyz;  // deep spatial coords for each mt vertex. material 2 vertices use deepBed coords.  rayIntersectSolids() repeatedly use these
	float _maxSceneSize;
	static float _cutSpacingInv;  // spacing between interior cut points inverted
	int _preDeepCutVerts;
	int _previousSkinTopEnd, _loopSkinTopBegin;
	std::list<std::list<long> > _holePolyLines;  // pair first is deepVert, second topVert

	void bilinearNormal(const double& u, const double& v, const bilinearPatch& bl, Vec3d& normal) {
		normal = bl.e11 * u + bl.e00 * (1.0 - u);
		Vec3d V = bl.e01 * v + bl.e10 * (1.0 - v);
		normal = normal ^ V;
	}

	bool getDeepSpatialCoordinates();  // used in new version.  Must have physics paused until deepCut complete or will be invalid.
	bool updateDeepSpatialCoordinates();
	bool rayIntersectMaterialTriangles(const Vec3d& rayStart, const Vec3d& rayDirection, std::vector<rayTriangleIntersect>& intersects);
	bool connectToPreviousPost(int postNum);
	double surfacePath(rayTriangleIntersect& from, const rayTriangleIntersect& to, const bilinearPatch* holeBl, const bool cutPath, double& minimumBilinearV);
	void cutSkinLine(int startV, int endV, std::vector<unsigned long>& te, std::vector<float>& params, bool Tin, bool Tout, surfaceCutLine& scl);
	void cutDeepSurface(int startV, int endV, std::vector<unsigned long>& te, std::vector<float>& params, surfaceCutLine& scl);
	bool connectOpenEnd(int postNum, int &interpostEnd);
	bool deepCutQuad(int postNum);
	void getDeepCutLine(rayTriangleIntersect& top, rayTriangleIntersect& bot);  // in material loci, not yet projected onto spatial plane
	bool interiorSpatialTet(const Vec3f pos, int& tet, Vec3f& baryWeight);
	void rayBilinearPatchIntersection(Vec3d& rayStart, Vec3d& rayN, const Vec3d& P00, const Vec3d& P10, const Vec3d& P01, const Vec3d& P11, double(&rayParam)[2], double(&faceParam)[2][2]);
	void makeBilinearPatch(const Vec3d& P00, const Vec3d& P10, const Vec3d& P01, const Vec3d& P11, bilinearPatch& bl);
	int bilinearRayIntersection(const Vec3d& rayStart, const Vec3d& rayDir, const bilinearPatch& bl, double (&rayParam)[2], Vec2d (&faceParams)[2]);
	void findCutInteriorHoles(const bilinearPatch& bl, const std::vector<Vec2d> &bUv, const std::vector<std::pair<long, Vec2d> >& deepOuterPolygon, std::list< std::vector<std::pair<long, Vec2d> > >& holes);

};

#endif  // __DEEP_CUT__





