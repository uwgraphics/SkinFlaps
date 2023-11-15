//////////////////////////////////////////////////////////
// File: materialTriangles.h
// Author: Court Cutting, MD
// Date: 2/26/2015
// Purpose: Triangle storage class using only uniquely linked xyz position
//    and uv texture data for vertices.  Normals are not used.
//    Triangle listings include 3 vertex indices, 3 vertex textures, and an integer material identifier.
//    An auxilliary graphics class, surgGraphics,
//    with its included vertex shaders provides vertex normal
//    doubling and tripling for graphics and user purposes.
//    The aux class surgGraphics provides
//    a tissue specific fragment shader to procedurally texture
//    the model for graphics purposes.
//////////////////////////////////////////////////////////

#ifndef __MATERIAL_TRIANGLES__
#define __MATERIAL_TRIANGLES__

#include <vector>
#include <array>
#include <map>
#include <unordered_map>
#include <set>
#include <string>
#include <iostream>
#include "Vec2f.h"
#include "Vec3f.h"

// forward declarations

class materialTriangles
{
public:
	void clear();
	int readObjFile(const char *fileName);  // uses shading group separators to separate materials
	bool writeObjFile(const char *fileName, const char* materialFileName=nullptr);
	void getVertexCoordinate(unsigned int vertex, float (&xyz)[3]) const;
	bool getBarycentricProjection(const int triangle, const float (&xyz)[3], float(&uv)[2]);
	void getBarycentricPosition(const int triangle, const float (&uv)[2], float (&xyz)[3]);
	void getBarycentricNormal(const int triangle, const float(&uv)[2], float(&nrm)[3]);
	void getBarycentricTexture(const int triangle, const float(&uv)[2], float(&texture)[2]);  // Court fix me
	void reserveVertices(int n) { _xyz.reserve(n); }
	void reserveTextures(int n) { _uv.reserve(n); }
	int addVertices(int numberToAdd = 1);  // Warning these three routines invalidate all pointers and iterators    // Court fix me
	inline int addTexture() { _uv.push_back(Vec2f(0.0f, 0.0f)); return (int)_uv.size() - 1; }
	inline void getTexture(const int txIndx, float(&tx)[2]) { tx[0] = _uv[txIndx].X; tx[1] = _uv[txIndx].Y; }
	inline void setTexture(const int txIndex, const float(&tx)[2]) { _uv[txIndex].X = tx[0]; _uv[txIndex].Y = tx[1]; }
	void reserveTriangles(int n) { _triPos.reserve(n); _triTex.reserve(n); _triMat.reserve(n);}
	int addTriangle(const int(&vertices)[3], const int material, const int(&textures)[3]);    // newer version
	void deleteTriangle(const int triangle) { _triMat[triangle] = -1; _triPos[triangle][0] = -1; }  // invalidate, but leave data in place.
	// ray inputs below are 3 element array pointers. Outputs triangles intersected and parameters along line.
	int findAdjacentTriangles(bool forceCompute=false);    // builds adjacency array for rapid neighbor searches
	void triangleVertexNeighbors(const int triangle, const int vertexNumber, std::vector<int>& neighborTriangles, std::vector<int>& neighborVertices);
	struct neighborNode{  // COURT - nuke this local variant
		int	vertex;
		int triangle;
	};
	void getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors);
	// next routine given a triangle and an edge(0-2), returns adjacent triangle and edge #.  If none, returns -1.
	inline float* vertexCoordinate(int vertex) {return (float*)&_xyz[vertex].xyz;}	// next 4 calls have no error checking for speed. Careful.
	inline float* getTexture(int txIndex) {return (float*)&_uv[txIndex].xy;}
	inline const float* getTexture(int txIndex) const { return (float*)&_uv[txIndex].xy; }
	inline int* triangleVertices(int triangle) {return _triPos[triangle].data(); }
	inline const int* triangleVertices(int triangle) const { return _triPos[triangle].data(); }
	inline int* triangleTextures(int triangle) { return _triTex[triangle].data(); }
	void triangleAdjacencies(int triangle, int(&adjTris)[3], int (&adjEdges)[3]);
	inline const int* triangleTextures(int triangle) const { return &(_triTex[triangle][0]); }
	inline int triangleMaterial(int triangle) const { return _triMat[triangle]; }
	inline void setTriangleMaterial(int triangle, int material) {_triMat[triangle] = material;}
	void setVertexCoordinate(int vertex, const float(&newCoord)[3]);
//	int getVertexTriangle(int vertexNumber){return _vertexFace[vertexNumber]&0x3fffffff;}	// gets triangle vertex is a member of
	void getTriangleNormal(int triangle, Vec3f& normal, bool normalized=true);
//	void getAreaNormal(const int triangle, const float (&uv)[2], const float radius, float(&normal)[3], bool normalized = true);
//	void getMeanVertexNormal(int vertex, float(&normal)[3], int onlyMaterial = -1);
	void getMeanVertexNormal(const int triangle, const int index, float(&normal)[3], int onlyMaterial = -1, bool normalize = true);  // if onlyMaterial>-1 only use neighbor triangles with material == onlyMaterial
	inline int numberOfTriangles() const { return (int)_triPos.size(); }
	inline int numberOfVertices() { return (int)_xyz.size(); }
	inline int numberOfTextures() { return (int)_uv.size(); }
	inline unsigned int* triAdjs(int triangleNumber) { return _adjs[triangleNumber].data(); }  // this routine is for legacy code and should no longer be used.
	inline unsigned int* vertexFaceTriangle(int vertex) {return &(_vertexFace[vertex]);} // careful if you don't know what you're doing
	inline void setName(const char *name) { _name.assign(name); }
	inline const std::string* getName() { return &_name; }
	inline const Vec3f triangleNormalNotNormalized(int triangle) {
		int* tr = _triPos[triangle].data();
		Vec3f v0, v1, N;
		v0 = _xyz[tr[1]] - _xyz[tr[0]];
		v1 = _xyz[tr[2]] - _xyz[tr[0]];
		return v0 ^ v1; }

//	matTriangle* getTriangleArray(int &numberOfTriangles);
//	float* getPositionArray(int &numberOfVertices);
//	float* getTextureArray(int &numberOfTextures);
//	std::vector<matTriangle>* getTriangleArray() { return &_tris; }

	inline const std::vector<std::array<int, 3> >& getTrianglePositionArray() { return _triPos; }
	inline const std::vector<std::array<int, 3> >& getTriangleTextureArray() { return  _triTex; }
	inline const std::vector<int>& getTriangleMaterialArray() { return _triMat; }
	std::vector<Vec3f>* getPositionArrayPtr() { return &_xyz; }
	std::vector<Vec3f>& getPositionArray() {return _xyz;}
	std::vector<Vec2f>& getTextureArray() { return _uv; }
	bool localPick(const float *lineStart, const float *lineDirection, float(&position)[3], int &triangle, float(&triangleParam)[2], const int onlyMaterial = -1);
	int linePick(const Vec3f& lineStart, const Vec3f& lineDirection, std::vector<Vec3f> &positions, std::vector<int> &triangles, std::vector<float> &params, const int onlyMaterial=-1);
//	int linePick(const float *lineStart, const float *lineDirection, std::vector<float> &rayParams, std::vector<int> &triangles, std::vector<Vec2f> &triangleParams);
	int splitTriangleEdge(int triangle, int edge, const float parameter);
	int addNewVertexInMidTriangle(int triangle, const float (&uvParameters)[2]);
//	int isValidClosedManifold();  // return -1 if inconsistent and # of topological handles if consistent
	bool deleteEdge(int triangle, int edge);  // always leaves triangle vertex[edge] behind and deletes vertex[edge+1] as well as the 2 triangles on either side of edge
	void closestPoint(const float(&xyz)[3], int& triangle, float(&uv)[2], int onlyMaterial = -1);

	// default behavior of deleteEdge() is remaining vertex is an average of the initial two. If you want asomething else (e.g. volume preservation) compute externally.
//	void cleanAndPack();  // warning - invalidates all triangle and vertex indices.
//	void cleanAndPack(std::vector<int> &newVertexMap, std::vector<int> &newTriangleMap); // returns mapping of old indices to new
//	bool hasSelfIntersection(const bool isClosed, std::vector<std::pair<int, int> >& triangleIntersectPairs); // slow. should only be used to debug surface
//	void partitionTriangleMaterials();  // key is material, second is index into triangle array of the next element beyond key material
	float getDiameter();

	materialTriangles(void);
	~materialTriangles(void);
	materialTriangles(const materialTriangles& x);

private:
	std::vector<std::array<int, 3> > _triPos;
	std::vector<std::array<int, 3> > _triTex;
	std::vector<int> _triMat;
	std::vector<Vec3f> _xyz;    // 3 float per vertex position data.
	std::vector<Vec2f> _uv;    // 2 float per vertex texture data.  Now changed such that a single vertex _xyz can be linked to more than one _uv to allow for texture seams.
    bool _adjacenciesComputed;
	std::vector<std::array<unsigned int, 3> > _adjs; 	// for each triangle edge low 2 bits are the edge number of the adjacent triangle.
        // If low 2 bits==3 and high order 30 bits==0, there is no adjacent triangle.
        // high 30 bits are the triangle number which must be bit shifted down 2
	std::string _name;
	struct edge {
		unsigned int reversed : 1;
		unsigned int vtxMin : 31;
		unsigned int matched : 1;  // check for non-manifold surface
		unsigned int vtxMax : 31;
		unsigned int adjCode;	// adjacency code for this edge
	};
	struct edgeTest {		// must be a less than operator
		bool operator()(const edge &e1,const edge &e2) const
		{ 
			if( e1.vtxMin < e2.vtxMin )
				return true;
			else if( e1.vtxMin > e2.vtxMin )
				return false;
			else	// vtxMin values equal
				return (e1.vtxMax < e2.vtxMax);
		}
	};
	std::vector<unsigned int> _vertexFace;

	void makeVertexToTriangleMap();
	bool parseNextInputFileLine(std::ifstream *infile, std::string &unparsedLine, std::vector<std::string> &parsedLine);
	void interpolateEdgeTextures(int triangle, int edge, int newVert, float param);
//	void recurseTriangleNormals(const int triangle, std::set<int> &trisDone, float(&center)[3], float radiusSq, float(&normalSum)[3]);
	bool rayTriangleIntersection(const Vec3f &rayOrigin, const Vec3f &rayDirection, const int triangle, float &rayParam, float(&triParam)[2], Vec3f &intersect);
	// be careful of next routine if you aren't expert. While local correction is faster, findAdjacentTriangles() is much less error prone.
//	void correctLocalNeighborArrays(std::vector<int> &changedTriangles);  // does local patching of adjacency and vertexFace arrays for input changedTriangles and the triangles adjacent to them
	struct lineHit{
		int triangle;
		Vec2f uv;  // parametric position on triangle
		Vec3f v;  // coordinate hit
	};
	int rayHits(const float *rayStart, const float *rayDirection, std::map<float, lineHit> &hits);

};

#endif  // __MATERIAL_TRIANGLES__
