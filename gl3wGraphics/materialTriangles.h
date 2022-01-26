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
	bool writeObjFile(const char *fileName);
	void collectCreateTextureSeams();
	void addOneMaterialTextureSeamVertex(long vertex, long(&textures)[2]);
	void getVertexCoordinate(unsigned int vertex, float (&xyz)[3]) const;
	bool getBarycentricProjection(const int triangle, const float (&xyz)[3], float(&uv)[2]);
	void getBarycentricPosition(const int triangle, const float (&uv)[2], float (&xyz)[3]);
	void getBarycentricNormal(const int triangle, const float(&uv)[2], float(&nrm)[3]);
	void getBarycentricTexture(const int triangle, const float(&uv)[2], float(&texture)[2]);  // Court fix me
	void reserveVertices(int n) { _xyz.reserve(n * 3); }
	void reserveTextures(int n) { _uv.reserve(n << 1); }
	long addVertices(int numberToAdd = 1);  // Warning these three routines invalidate all pointers and iterators    // Court fix me
	long cloneVertex(long sourceVertex);  // makes a new vertex which is a copy of sourceVertex    // Court fix me
	inline int addTexture() { _uv.push_back(0.0f); _uv.push_back(0.0f); return (int)(_uv.size() >> 1) - 1; }
	inline void getTexture(const int txIndx, float(&tx)[2]) { tx[0] = _uv[txIndx << 1]; tx[1] = _uv[(txIndx << 1) + 1]; }
	inline void setTexture(const int txIndex, const float(&tx)[2]) { _uv[txIndex << 1] = tx[0]; _uv[(txIndex << 1) + 1] = tx[1]; }
	void reserveTriangles(int n) { _tris.reserve(n); }
	long addTriangle(long(&vertices)[3], long material);  // for backward compatibility
	long addTriangle(const long(&vertices)[3], const long material, const long(&textures)[3]);    // newer version
	// ray inputs below are 3 element array pointers. Outputs triangles intersected and parameters along line.
	int rayIntersect(const float *rayStart, const float *rayDirection, std::vector<int> &triangles, std::vector<float> &params);
	struct neighborNode{  // COURT - nuke this local variant
		long	vertex;
		long triangle;
	};
	int findAdjacentTriangles(bool forceCompute=false, bool fullManifoldTest = false);    // builds adjacency array for rapid neighbor searches    // Court fix me
	void getNeighbors(unsigned int vertex, std::vector<neighborNode> &neighbors);    // Court fix me
	// next routine given a triangle and an edge(0-2), returns adjacent triangle and edge #.  If none, returns -1.
	inline void edgeAdjacency(int &triangle, int &edge) {unsigned long adj = _adjs[triangle * 3 + edge]; if (adj == 0x03) { triangle = -1; edge = -1; } else { triangle = adj >> 2; edge = adj & 0x03; }}
	inline float* vertexCoordinate(int vertex) {return (float*)&_xyz[(vertex<<1)+vertex];}	// next 4 calls have no error checking for speed. Careful.
	inline float* getTexture(int txIndex) {return (float*)&_uv[txIndex <<1];}
	inline const float* getTexture(int txIndex) const { return (float*)&_uv[txIndex << 1]; }
	inline bool vertexDisconected(int vertex) { return _vertexFace[vertex] == 0x80000000; }
	inline long* triangleVertices(long triangle) {return &_tris[triangle].v[0];}
	inline const long* triangleVertices(long triangle) const { return &_tris[triangle].v[0]; }
	inline long* triangleTextures(long triangle) { return &(_tris[triangle].tex[0]); }
	inline const long* triangleTextures(long triangle) const { return &(_tris[triangle].tex[0]); }
	inline long triangleMaterial(long triangle) const { return _tris[triangle].material; }
	inline void setTriangleMaterial(long triangle, int material) {_tris[triangle].material=material;}
	void setVertexCoordinate(int vertex, const float(&newCoord)[3]);
	int getVertexTriangle(int vertexNumber){return _vertexFace[vertexNumber]&0x3fffffff;}	// gets triangle vertex is a member of
	void getTriangleNormal(int triangle, float (&normal)[3], bool normalized=true);
	void getAreaNormal(const int triangle, const float (&uv)[2], const float radius, float(&normal)[3], bool normalized = true);
	void getMeanVertexNormal(int vertex, float(&normal)[3], long onlyMaterial = -1);  // if onlyMaterial>-1 only use neighbor triangles with material == onlyMaterial
	inline int numberOfTriangles() const { return (int)_tris.size(); }
	inline long numberOfVertices() { return (long)(_xyz.size() / 3); }
	inline long numberOfTextures() { return (long)(_uv.size()>>1); }
	inline unsigned long* triAdjs(long triangleNumber) { return &(_adjs[(triangleNumber << 1) + triangleNumber]); }
	inline unsigned long* vertexFaceTriangle(int vertex) {return &(_vertexFace[vertex]);} // careful if you don't know what you're doing
	inline void setName(const char *name) { _name.assign(name); }
	inline const std::string* getName() { return &_name; }

	struct alignas(32) matTriangle{
		long v[3];
		long material;
		long tex[3];  // could use pad if needed   long pad;
	};
	matTriangle* getTriangleArray(int &numberOfTriangles);
	float* getPositionArray(int &numberOfVertices);
	float* getTextureArray(int &numberOfTextures);
	std::vector<matTriangle>* getTriangleArray() { return &_tris; }
	std::vector<float>* getPositionArray() {return &_xyz;}
	std::vector<float>* getTextureArray() { return &_uv; }
	void getNearestHardEdge(float(&xyz)[3], int &triangle, int &edge, float &param, int materialLimit = -1);
	bool localPick(const float *lineStart, const float *lineDirection, float(&position)[3], int &triangle, float(&triangleParam)[2], const int onlyMaterial = -1);
	int linePick(const float *lineStart, const float *lineDirection, std::vector<float> &positions, std::vector<int> &triangles, std::vector<float> &params, const int onlyMaterial=-1);
	int linePick(const float *lineStart, const float *lineDirection, std::vector<float> &rayParams, std::vector<int> &triangles, std::vector<Vec2f> &triangleParams);
	int splitTriangleEdge(int triangle, int edge, const float parameter);
	int addNewVertexInMidTriangle(int triangle, const float (&uvParameters)[2]);
	bool inside(const float(&xyz)[3]);  // if a closed manifold surface, returns if xyz is inside  COURT - get rid of this
	int isManifoldConsistent();  // return -1 if inconsistent and # of topological handles if consistent
	bool deleteEdge(int triangle, int edge);  // always leaves triangle vertex[edge] behind and deletes vertex[edge+1] as well as the 2 triangles on either side of edge

	// default behavior of deleteEdge() is remaining vertex is an average of the initial two. If you want asomething else (e.g. volume preservation) compute externally.
	void cleanAndPack();  // warning - invalidates all triangle and vertex indices.
	void cleanAndPack(std::vector<int> &newVertexMap, std::vector<int> &newTriangleMap); // returns mapping of old indices to new
	bool topoCheck(); // checks current topology versus recomputed from scratch
	bool hasSelfIntersection(const bool isClosed, std::vector<std::pair<int, int> >& triangleIntersectPairs); // slow. should only be used to debug surface
	void partitionTriangleMaterials();  // key is material, second is index into triangle array of the next element beyond key material
	void closestPoint(const float(&xyz)[3], long &triangle, float(&uv)[2], long onlyMaterial = -1);  // closest barycentric position to point xyz. Can limit search to onlyMaterial if desired(-1 is no limit)
	bool textureFind(const float(&txIn)[2], const long materialIn, int &triangle, float(&uv)[2]);

	materialTriangles(void);
	~materialTriangles(void);
	materialTriangles(const materialTriangles& x);

private:
	std::vector<matTriangle> _tris;	// 3 indices into the vertex arrays and the material that make a triangle. Material -1 indicates deletion.
	std::vector<float> _xyz;    // 3 float per vertex position data.
	std::vector<float> _uv;    // 2 float per vertex texture data.  Now changed such that a single vertex _xyz can be linked to more than one _uv to allow for texture seams.
    bool _adjacenciesComputed;
    std::vector<unsigned long> _adjs;	// low 2 bits are the edge number of the adjacent triangle.
        // If low 2 bits==3 and high order 30 bits==0, there is no adjacent triangle.
        // high 30 bits are the triangle number which must be bit shifted down 2
	std::string _name;
//	std::unordered_multimap<long, long> _oneMaterialSeams;  // 2 textures at same position with same material
	std::multimap<long, long> _oneMaterialSeams;  // 2 textures at same position with same material
	struct edge	{
		unsigned long reversed : 1;
		unsigned long vtxMin : 31;
		unsigned long matched : 1;  // check for non-manifold surface
		unsigned long vtxMax : 31;
		unsigned long adjCode;	// adjacency code for this edge
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
	void makeVertexToTriangleMap();
	std::vector<unsigned long> _vertexFace;
	bool parseNextInputFileLine(std::ifstream *infile, std::string &unparsedLine, std::vector<std::string> &parsedLine);
	void interpolateEdgeTextures(int triangle, int edge, int newVert, float param);
	void recurseTriangleNormals(const long triangle, std::set<long> &trisDone, float(&center)[3], float radiusSq, float(&normalSum)[3]);
	bool rayTriangleIntersection(const Vec3f &rayOrigin, const Vec3f &rayDirection, const int triangle, float &rayParam, float(&triParam)[2], Vec3f &intersect);
	// be careful of next routine if you aren't expert. While local correction is faster, findAdjacentTriangles() is much less error prone.
	void correctLocalNeighborArrays(std::vector<long> &changedTriangles);  // does local patching of adjacency and vertexFace arrays for input changedTriangles and the triangles adjacent to them
	bool checkTopology();
	struct lineHit{
		int triangle;
		Vec2f uv;  // parametric position on triangle
		Vec3f v;  // coordinate hit
	};
	int rayHits(const float *rayStart, const float *rayDirection, std::map<float, lineHit> &hits);

	friend class mtConvert;	// this class benefits tremendously from interior access
	friend class incisionUv;	// interior access vital here
	friend class annealSkinBed;  // COURT - nuke a lot of these friends
	friend class growSkinBed;
	friend class deepCut;
	friend class deepCutTet;
	friend class incisionVnCubes;
	friend class skinCutUndermine;
	friend class skinCutUndermineTets;
	friend class surgGraphics;
	friend class staticTriangle;

};

#endif  // __MATERIAL_TRIANGLES__
