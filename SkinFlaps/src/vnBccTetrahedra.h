////////////////////////////////////////////////////////////////////////////
// File: vnBccTetrahedra.h
// Author: Court Cutting
// Date: 3/1/2019
// Purpose: Basic virtual noded bcc tet class where tets in space are not unique, but may be duplicated by the use of virtual nodes.
//     Full description of this concept is given in original work by Molino N,  Bao Z, and Fedkiw R: http://physbam.stanford.edu/~fedkiw/papers/stanford2004-01.pdf
////////////////////////////////////////////////////////////////////////////

#ifndef __VN_BCC_TETS__
#define __VN_BCC_TETS__

#include <vector>
#include <forward_list>
#include <list>
#include <set>
#include <array>
#include <memory>
#include <unordered_map>
#include <atomic>
#include <stdexcept>
#include "Vec3f.h"
#include "Mat3x3f.h"

#pragma warning (disable : 4267)

// forward declarations
class materialTriangles;

// each bcc tet has a unique centroid where only one of its axes is a non-integer containing 0.5.  This is the axis not used in forming its
// "Delaunay" tetrahedralization of a Cartesian line segment with a neighbor from the other lattice immediately perpendicular to it.  Thus an odd-even Cartesian
// line segment pair do not require one Cartesian axis to define it.  In this axis the centroid of the tet will be at xx.5.  For this reason
// all centroid coordinates are multiplied by two so that an array of 3 shorts will hold it.
typedef std::array<unsigned short, 3> bccTetCentroid;

class vnBccTetrahedra
{
public:
	void clear();
	void materialCoordsToNodeSpatialVector();
	inline bool empty(){ return _tetNodes.empty() || _nodeGridLoci.empty(); }
	inline int nodeNumber() { return (int)_nodeGridLoci.size(); }
	inline int tetNumber() { return (int)_tetNodes.size(); }
	inline int vertexNumber() { return (int)_vertexTets.size(); }
	inline const int* tetNodes(int tetIndex){ return _tetNodes[tetIndex].data(); }
	const std::vector<std::array<int, 4> >& getTetNodeArray() { return _tetNodes; }
	const std::vector<std::array<unsigned short, 3> >& getTetCentroidArray() { return _tetCentroids; }  // remember actual material coord centroids are half of each value to enable integer packing.
	inline void centroidTets(const bccTetCentroid &tc, std::list<int> &tets){ auto pr = _tetHash.equal_range(tc); tets.clear(); while (pr.first != pr.second){ tets.push_back(pr.first->second); ++pr.first; } }
	inline const int getVertexTetrahedron(const int vertex) const {return _vertexTets[vertex];}
	inline void setVertexTetrahedron(const int vertex, const int newTetIndex){ _vertexTets[vertex] = newTetIndex; }
	inline const Vec3f* getVertexWeight(const int vertex) const { return &_barycentricWeights[vertex]; }
	inline const Vec3f &getMinimumCorner() { return _minCorner; }
	inline const Vec3f &getMaximumCorner() { return _maxCorner; }
	inline const Vec3f &nodeSpatialCoordinate(const int nodeIndex) { return _nodeSpatialCoords[nodeIndex]; }
	inline const Vec3f nodeMaterialCoordinate(const int nodeIndex) { return _minCorner + Vec3f(reinterpret_cast<short (&)[3]>(*_nodeGridLoci[nodeIndex].data())) * (float) _unitSpacing; }
	inline const float *nodeSpatialCoordinatePtr(const int nodeIndex) { return _nodeSpatialCoords[nodeIndex].xyz; }
	inline double getTetUnitSize() { return _unitSpacing; }
	inline double getTetUnitSizeInv() { return _unitSpacingInv; }
	// next set of routines do transformations from one coordinate system to another
	inline void spatialToGridCoords(const Vec3f &spatialCoords, Vec3f &gridCoords){ gridCoords = spatialCoords - _minCorner; gridCoords *= (float)_unitSpacingInv; }  // only for MATERIAL spatial coord input
	void gridLocusToTetCentroid(const Vec3f &gridLocus, bccTetCentroid &tetCentroid);
	void gridLocusToBarycentricWeight(const Vec3f &gridLocus, const bccTetCentroid &tetCentroid, Vec3f &barycentricWeight);
	void barycentricWeightToGridLocus(const int tet, const Vec3f& barycentricWeight, Vec3f& gridLocus);
	void barycentricWeightToGridLocus(const bccTetCentroid &tetCentroid, const Vec3f &barycentricWeight, Vec3f &gridLocus);
	void vertexGridLocus(const int vertex, Vec3f &gridLocus);  // always material coords
	void vertexMaterialCoordinate(const int vertex, std::array<float, 3> &matCoord);
	int firstInteriorTet() { return _firstInteriorTet; }  // all tets before this are surface tets possibly virtual noded and with out unique centroid. Here on are unique interior tets.
	inline const bccTetCentroid& tetCentroid(int tet) { return _tetCentroids[tet]; }  // do we want to keep these when they are so easily computed?

	inline const bccTetCentroid computeTetCentroid(int tet) {
		bccTetCentroid tc;
		std::array<short, 3> center;
		auto& tn = _tetNodes[tet];
		center = _nodeGridLoci[tn[0]];
		for (int i = 1; i < 4; ++i) {
			auto gl = _nodeGridLoci[tn[i]];
			for(int j=0; j<3; ++j)
				center[j] += gl[j];
		}
		for (int j = 0; j < 3; ++j)
			tc[j] = center[j] >> 1;
		return tc;
	}

	inline void centroidToXyzHalfAxis(const bccTetCentroid cntd, std::array<short, 3>& xyz, int &halfAxis) {
		for (int i = 0; i < 3; ++i) {  // no error check for speed
			if (cntd[i] & 1)
				halfAxis = i;
			xyz[i] = cntd[i] >> 1;
		}
	}

	inline void getBarycentricTetPosition(const int tet, const Vec3f &barycentricWeight, Vec3f &position)
	{
		const int *n = _tetNodes[tet].data();
		position.set(_nodeSpatialCoords[n[0]] * (1.0f - barycentricWeight.X - barycentricWeight.Y - barycentricWeight.Z));
		for (int i = 1; i < 4; ++i)
			position += _nodeSpatialCoords[n[i]] * barycentricWeight[i - 1];
	}

	inline void vertexBarycentricPosition(const int vertex, Vec3f &position)
	{
		const int *n = _tetNodes[_vertexTets[vertex]].data();
		float *bw = _barycentricWeights[vertex].xyz;
		position.set(_nodeSpatialCoords[n[0]] * (1.0f - *bw - bw[1] - bw[2]));
		for (int i = 1; i < 4; ++i)
			position += _nodeSpatialCoords[n[i]] * bw[i - 1];
	}

	void setNodeSpatialCoordinatePointer(std::array<float, 3> *spatialCoordPtr) { _nodeSpatialCoords = reinterpret_cast<Vec3f*>(spatialCoordPtr); }  // assumes vector of coords created elsewhere
	void setNodeSpatialCoordinatePointer(Vec3f *spatialCoordPtr) { _nodeSpatialCoords = spatialCoordPtr; }  // assumes vector of coords created elsewhere
	const Vec3f* getNodeSpatialCoordPointer() { if (_nodeSpatialCoords == nullptr) throw(std::logic_error("Trying to access nodeSpatialCoordinate vector before it has been allocated andassigned")); return _nodeSpatialCoords; }
	// next set of routines traverse topological paths through the bcc data
	int faceAdjacentTet(const bccTetCentroid tc, const int face, bccTetCentroid& tcAdj);  // fundamental code for all topological path routines.
	// Above returns adjacent face # and adjacent tet centroid. a -1 return signals an illegal centroid.
	int faceAdjacentTetNodeIndices(const bccTetCentroid tc, const int face, bccTetCentroid& tcAdj, int(&adjNodeIndices)[3]);  // adds node indices of adjTet if you need them

	int faceAdjacentTets(const int tet, const int face, std::list<int> &adjTets);  // return adjacent face index 0-3
	
	inline void faceNodes(const int tet, const int face, int(&nodes)[3])
	{
		const int *tn = tetNodes(tet);
		for (int i = 0; i < 3; ++i)
			nodes[i] = (tn[(face + i) & 3]);
	}

	inline bool adjacentCentroids(const bccTetCentroid &tc0, const bccTetCentroid &tc1){
		int d, n = 0;
		for (int i = 0; i < 3; ++i) {
			d = abs(tc0[i] - tc1[i]);
			if (d > 1)
				return false;
			else if (d > 0)
				++n;
			else
				;
		}
		return n == 2;
	}

	void edgeAdjacentTets(const int tet, const int edge, std::list<int> &adjTets);  // input one of six edges in permutation order 0-123, 1-23, and 2-3
	void edgeNodes(const int tet, const int edge, int &n0, int &n1);  // same edge numbering as above
	bool decreasingCentroidPath(int startTet, const int targetTet, std::list<int> &tetPath);  // true if constantly decreasing distance centroid path exists.
	int parametricTriangleTet(const int mtTriangle, const float(&uv)[2], Vec3f& gridLocus);  // returns grid locus and tetrahedron at parametric location uv in input triangle

//	inline void setNodeFixationState(const int node, bool fixed){ if (fixed) _fixedNodes.insert(node);  else  _fixedNodes.erase(node); }
//	inline bool nodeFixed(const int node){ return _fixedNodes.find(node) != _fixedNodes.end(); }
//	inline const short* nodeGridLocation(const int tetNode){ return _nodeGridLoci[tetNode].data(); }
	inline const std::array<short, 3>& nodeGridLocation(const int tetNode) { return _nodeGridLoci[tetNode]; }
	void centroidToNodeLoci(const bccTetCentroid& centroid, short (&gridLoci)[4][3]);
	void centroidToNodeLocus(const bccTetCentroid& centroid, const int nodeIndex, short(&gridLocus)[3]);

	void nodeCentroids(std::array<short, 3>& node, bccTetCentroid cntrd[24]);  // for these routines a coordinate greater than 65533 indicates a centroid out of positive grid bounds
	void CartesianEdgeCentroids(const short(&edgeMidpoint)[3], bccTetCentroid(&centroidLoci)[4]);  // centroids surrounding a Cartesian edge
	void CartesianEdgeCentroids(const short(&edgeMidpoint)[3], Vec3f(&centroidLoci)[4]); 
	void unitCubeCentroids(const short(&minimumCorner)[3], bccTetCentroid(&cntrd)[6]);  // centroids whose tet is part of a unit cube
	void unitCubeCentroids(const short(&minimumCorner)[3], Vec3f(&centroidLoci)[6]);
	inline int firstInteriorTetrahedron() { return _firstInteriorTet;  }
	bool insideTet(const bccTetCentroid& tc, const std::array<short, 3>& nodeLocus);  // material coords, not spatial
	bool insideTet(const bccTetCentroid& tc, const Vec3f& gridLocus);

	vnBccTetrahedra(const vnBccTetrahedra&) = delete;
	vnBccTetrahedra& operator=(const vnBccTetrahedra&) = delete;
	vnBccTetrahedra();
	~vnBccTetrahedra();

protected:
	std::vector<std::array<short, 3> > _nodeGridLoci;  // can be negative at min planes
	// The convention adopted in this class for tetrahedral node ordering is positive as in CGAL 4.14.  This means when viewing the tet from the outside that faces 012 and 230 are listed clockwise
	// and 123 and 301 are counterclockwise.  Each BCC tet can be defined by two line segments in two of the Cartesian axes.  The axis not used is the non-integer halfCoordinateAxis of the tets centroid.
	// The remaining two axis line segments are given from minimum to maximum in vertices 0 to 1 in the axis immediately following the halfCoordAxis (x to y, y to z, and z to x).
	// with vertices 2 and 3 ordered such that the CGAL positive ordering convention is consistent.  With this convention enforced, it isn't necessary to save 3-integer node grid locations.
	// The bcc tet centroid specifies its four 3-short integer node locations implicitly. See centroidToNodeLoci().
	std::vector<std::array<int, 4> > _tetNodes;
	std::vector<bccTetCentroid> _tetCentroids;
	struct bccTetCentroidHasher {
		std::size_t operator()(const bccTetCentroid& k) const
		{  // hash function
			long long ll = k[0];
			ll <<= 16;
			ll += k[1];
			ll <<= 16;
			ll += k[2];
			std::hash<long long> hash_funct;
			return hash_funct(ll);
		}
	};
	std::unordered_multimap<bccTetCentroid, int, bccTetCentroidHasher> _tetHash;  // bccTetCenter and index into _tetNodes
	materialTriangles *_mt;  // embedded surface
	std::vector<int> _vertexTets;
	std::vector<Vec3f> _barycentricWeights;
	Vec3f *_nodeSpatialCoords;
	Vec3f _minCorner, _maxCorner;
	double _unitSpacing, _unitSpacingInv;
	int _gridSize[3];  // Number of subdivisions in x,y, and z.
	int _firstInteriorTet;  // first tet of all unique interior tets. Goes to end of list.
	static Mat3x3f _barycentricInverses[6];

	friend class vnBccTetCutter;
	friend class vnBccTetCutter_omp;
	friend class skinCutUndermineTets;
};

#endif // __VN_BCC_TETS__