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
//#include "Vec3d.h"  // COURT - nuke

#pragma warning (disable : 4267)

// forward declarations
class materialTriangles;

// each bcc tet has a unique centroid where only one of its axes is a non-integer containing 0.5.  This is the axis not used in forming its
// "Delaunay" tetrahedralization of a Cartesian line segment with a neighbor from the other lattice immediately perpendicular to it.  Thus an odd-even Cartesian
// line segment pair do not require one Cartesian axis to define it.  In this axis the centroid of the tet will be at xx.5.  This half coordinate axis
// is identified in the structure below.
class bccTetCentroid{
public:
	union{
		struct {
			short halfCoordAxis; // 0-2 the axis where the centroid is its xyz[] + 0.5
			std::array<short, 3> xyz;
		};
		long long ll;  // use this for hash table indexing
	};
	inline bool operator == (const bccTetCentroid btc) const { return ll == btc.ll; }
	inline bool operator < (const bccTetCentroid btc) const { if (halfCoordAxis < btc.halfCoordAxis) return true; else if (halfCoordAxis > btc.halfCoordAxis) return false; else return xyz < btc.xyz; }
};

class vnBccTetrahedra
{
public:
	void clear();
	void materialCoordsToSpatialVector();
	inline bool empty(){ return _tetNodes.empty() || _nodeGridLoci.empty(); }
	inline int nodeNumber() { return (int)_nodeGridLoci.size(); }
	inline int tetNumber() { return (int)_tetNodes.size(); }
	inline int vertexNumber() { return (int)_vertexTets.size(); }
	inline const long* tetNodes(int tetIndex){ return _tetNodes[tetIndex].data(); }
	const std::vector<std::array<long, 4> >& getTetNodeArray() { return _tetNodes; }
	inline const bccTetCentroid* tetCentroid(int tet){ return &_tetCentroids[tet]; }
	inline void centroidTets(const bccTetCentroid &tc, std::list<long> &tets){ auto pr = _tetHash.equal_range(tc.ll); tets.clear(); while (pr.first != pr.second){ tets.push_back(pr.first->second); ++pr.first; } }
	inline const long getVertexTetrahedron(const int vertex) const {return _vertexTets[vertex];}
	inline void setVertexTetrahedron(const int vertex, const long newTetIndex){ _vertexTets[vertex] = newTetIndex; }
	inline const Vec3f* getVertexWeight(const int vertex) const { return &_barycentricWeights[vertex]; }
	inline const Vec3f &getMinimumCorner() { return _minCorner; }
	inline const Vec3f &getMaximumCorner() { return _maxCorner; }
	inline const Vec3f &nodeSpatialCoordinate(const long nodeIndex) { return _nodeSpatialCoords[nodeIndex]; }
	inline const Vec3f nodeMaterialCoordinate(const long nodeIndex) { return _minCorner + Vec3f(reinterpret_cast<short (&)[3]>(*_nodeGridLoci[nodeIndex].data())) * (float) _unitSpacing; }
	inline float *nodeSpatialCoordinatePtr(const long nodeIndex) { return _nodeSpatialCoords[nodeIndex]._v; }
	inline double getTetUnitSize() { return _unitSpacing; }
	inline double getTetUnitSizeInv() { return _unitSpacingInv; }
	// next set of routines do transformations from one coordinate system to another
	inline void spatialToGridCoords(const Vec3f &spatialCoords, Vec3f &gridCoords){ gridCoords = spatialCoords - _minCorner; gridCoords *= (float)_unitSpacingInv; }  // only for MATERIAL spatial coord input
	void gridLocusToTetCentroid(const Vec3f &gridLocus, bccTetCentroid &tetCentroid);
	void gridLocusToBarycentricWeight(const Vec3f &gridLocus, const bccTetCentroid &tetCentroid, Vec3f &barycentricWeight);
	void barycentricWeightToGridLocus(const bccTetCentroid &tetCentroid, const Vec3f &barycentricWeight, Vec3f &gridLocus);
	void vertexGridLocus(const long vertex, Vec3f &gridLocus);  // always material coords
	void vertexMaterialCoordinate(const long vertex, std::array<float, 3> &matCoord);
	int firstInteriorTet() { return _firstInteriorTet; }  // all tets before this are surface tets possibly virtual noded and with out unique centroid. Here on are unique interior tets.

	inline void getBarycentricTetPosition(const long tet, const Vec3f &barycentricWeight, Vec3f &position)
	{
		const long *n = _tetNodes[tet].data();
		position.set(_nodeSpatialCoords[n[0]] * (1.0f - barycentricWeight.X - barycentricWeight.Y - barycentricWeight.Z));
		for (int i = 1; i < 4; ++i)
			position += _nodeSpatialCoords[n[i]] * barycentricWeight[i - 1];
	}

	inline void vertexBarycentricPosition(const long vertex, Vec3f &position)
	{
		const long *n = _tetNodes[_vertexTets[vertex]].data();
		float *bw = _barycentricWeights[vertex]._v;
		position.set(_nodeSpatialCoords[n[0]] * (1.0f - *bw - bw[1] - bw[2]));
		for (int i = 1; i < 4; ++i)
			position += _nodeSpatialCoords[n[i]] * bw[i - 1];
	}

	void setNodeSpatialCoordinatePointer(std::array<float, 3> *spatialCoordPtr) { _nodeSpatialCoords = reinterpret_cast<Vec3f*>(spatialCoordPtr); }  // assumes vector of coords created elsewhere
	void setNodeSpatialCoordinatePointer(Vec3f *spatialCoordPtr) { _nodeSpatialCoords = spatialCoordPtr; }  // assumes vector of coords created elsewhere
	const Vec3f* getNodeSpatialCoordPointer() { if (_nodeSpatialCoords == nullptr) throw(std::logic_error("Trying to access nodeSpatialCoordinate vector before it has been allocated andassigned")); return _nodeSpatialCoords; }
	// next set of routines traverse topological paths through the bcc data
	int faceAdjacentTet(const bccTetCentroid tc, const int face, bccTetCentroid& tcAdj);  // fundamental code for all topological path routines. Returns adjacent face # and adjacent tet centroid.
	int faceAdjacentTets(const long tet, const int face, std::list<long> &adjTets);  // return adjacent face index 0-3
	
	inline void faceNodes(const long tet, const int face, long(&nodes)[3])
	{
		const long *tn = tetNodes(tet);
		for (int i = 0; i < 3; ++i)
			nodes[i] = (tn[(face + i) & 3]);
	}

	void edgeAdjacentTets(const long tet, const int edge, std::list<long> &adjTets);  // input one of six edges in permutation order 0-123, 1-23, and 2-3
	void edgeNodes(const long tet, const int edge, long &n0, long &n1);  // same edge numbering as above
	bool decreasingCentroidPath(long startTet, const long targetTet, std::list<long> &tetPath);  // true if constantly decreasing distance centroid path exists.

	inline void setNodeFixationState(const long node, bool fixed){ if (fixed) _fixedNodes.insert(node);  else  _fixedNodes.erase(node); }
	inline bool nodeFixed(const long node){ return _fixedNodes.find(node) != _fixedNodes.end(); }
	inline const  short* nodeGridLocation(const long tetNode){ return _nodeGridLoci[tetNode].data(); }
	inline const  std::array<short, 3> nodeGridLocation(const bccTetCentroid &tetCentroid, const int nodeIndex)
	{	// given a bccTetCentroid and its nodeIndex (0-3) return material coord location of node
		bccTetCentroid tc;
		tc = tetCentroid;
		bool below01 = (tc.xyz[tc.halfCoordAxis] + tc.xyz[(tc.halfCoordAxis + 2) % 3]) & 1;
		if (nodeIndex < 2){
			if (below01) ++tc.xyz[tc.halfCoordAxis];
			nodeIndex & 1 ? ++tc.xyz[(tc.halfCoordAxis + 1) % 3] : --tc.xyz[(tc.halfCoordAxis + 1) % 3];
		}
		else{
			if (!below01) ++tc.xyz[tc.halfCoordAxis];
			below01 == (nodeIndex < 3) ? --tc.xyz[(tc.halfCoordAxis + 2) % 3] : ++tc.xyz[(tc.halfCoordAxis + 2) % 3];
		}
		return tc.xyz;
	}
	const std::unordered_multimap<long long, long>* getTetHash() { return _tetHash.empty() ? nullptr : &_tetHash; }  // bccTetCenter and index into _tetNodes

	vnBccTetrahedra(const vnBccTetrahedra&) = delete;
	vnBccTetrahedra& operator=(const vnBccTetrahedra&) = delete;
	vnBccTetrahedra();
	~vnBccTetrahedra();

private:
	std::set<long> _fixedNodes;
	std::vector<std::array<short, 3> > _nodeGridLoci;
	// The convention adopted in this class for tetrahedral node ordering is positive as in CGAL 4.14.  This means when viewing the tet from the outside that faces 012 and 230 are listed clockwise
	// and 123 and 301 are counterclockwise.  Each BCC tet can be defined by two line segments in two of the Cartesian axes.  The axis not used is the non-integer halfCoordinateAxis of the tets centroid.
	// The remaining two axis line segments are given from minimum to maximum in vertices 0 to 1 in the axis immediately following the halfCoordAxis (x to y, y to z, and z to x).
	// with vertices 2 and 3 ordered such that the CGAL positive ordering convention is consistent.  With this convention enforced, it isn't necessary to save 3-integer node grid locations.
	// The bcc tet centroid specifies its four 3-integer node locations implicitly. See nodeGridLocation().
	std::vector<std::array<long, 4> > _tetNodes;
	std::vector<bccTetCentroid> _tetCentroids;
	std::unordered_multimap<long long, long> _tetHash;  // bccTetCenter and index into _tetNodes
	materialTriangles *_mt;  // embedded surface
	std::vector<long> _vertexTets;
	std::vector<Vec3f> _barycentricWeights;
	Vec3f *_nodeSpatialCoords;
	Vec3f _minCorner, _maxCorner;
	double _unitSpacing, _unitSpacingInv;
	int _gridSize[3];  // COURT - perhaps should be moved to the cutter.  Just confusing to the user.
	int _firstInteriorTet;  // first tet of all unique interior tets. Goes to end of list.

	friend class vnBccTetCutter;
	friend class vnBccTetCutterTbb;
	friend class IntersectPlanePolygon;
	friend class skinCutUndermineTets;
	friend class deepCut;
	friend class remapTetPhysics;

};

#endif // __VN_BCC_TETS__