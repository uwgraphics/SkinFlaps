////////////////////////////////////////////////////////////////////////////
// File: vnBccTetrahedra.h
// Author: Court Cutting
// Date: 3/1/2019
// Modified: 4/2/2023
// Purpose: Basic virtual noded bcc tet class where tets in space are not unique, but may be duplicated by the use of virtual nodes.
//     Full description of this concept is given in original work by Molino N,  Bao Z, and Fedkiw R: http://physbam.stanford.edu/~fedkiw/papers/stanford2004-01.pdf
//     Latest major modification handles multi resolution tets with sizes that increase by binary steps.  Lowest bit of centroid reveals tet size. 
//     New centroid code left shifts x, y, and z up one bit so half coordinate stored in a short.
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
	void getTJunctionConstraints(std::vector<int>& subNodes, std::vector<std::vector<int> >& macroNodes, std::vector<std::vector<float> >& macroBarycentrics);  // T junctions created in multires cutter
	const std::vector<bccTetCentroid>& getTetCentroidArray() { return _tetCentroids; }  // remember actual material coord centroids are half of each value to enable integer packing.
	inline void centroidTets(const bccTetCentroid &tc, std::list<int> &tets){ auto pr = _tetHash.equal_range(tc); tets.clear(); while (pr.first != pr.second){ tets.push_back(pr.first->second); ++pr.first; } }
	inline const int getVertexTetrahedron(const int vertex) const {return _vertexTets[vertex];}
	inline void setVertexTetrahedron(const int vertex, const int newTetIndex){ _vertexTets[vertex] = newTetIndex; }
	inline const Vec3f* getVertexWeight(const int vertex) const { return &_barycentricWeights[vertex]; }
	inline const Vec3f &getMinimumCorner() { return _minCorner; }
	inline const Vec3f &getMaximumCorner() { return _maxCorner; }
	inline const Vec3f &nodeSpatialCoordinate(const int nodeIndex) { return _nodeSpatialCoords[nodeIndex]; }
	inline const Vec3f nodeMaterialCoordinate(const int nodeIndex) { return _minCorner + Vec3f(reinterpret_cast<short (&)[3]>(*_nodeGridLoci[nodeIndex].data())) * (float) _unitSpacing; }
	inline const float *nodeSpatialCoordinatePtr(const int nodeIndex) { return _nodeSpatialCoords[nodeIndex].xyz; }
	inline const std::array<short, 3>& nodeGridLocation(const int tetNode) { return _nodeGridLoci[tetNode]; }
	inline double getTetUnitSize() { return _unitSpacing; }
	inline double getTetUnitSizeInv() { return _unitSpacingInv; }

	inline materialTriangles* getMaterialTriangles() { return _mt; }

	// next set of routines do transformations from one coordinate system to another
	inline void spatialToGridCoords(const Vec3f &spatialCoords, Vec3f &gridCoords){ gridCoords = spatialCoords - _minCorner; gridCoords *= (float)_unitSpacingInv; }  // only for MATERIAL spatial coord input
	void gridLocusToLowestTetCentroid(const Vec3f &gridLocus, bccTetCentroid &tetCentroid);
	void gridLocusToBarycentricWeight(const Vec3f &gridLocus, const bccTetCentroid &tetCentroid, Vec3f &barycentricWeight);
	void barycentricWeightToGridLocus(const int tet, const Vec3f& barycentricWeight, Vec3f& gridLocus);
	void barycentricWeightToGridLocus(const bccTetCentroid &tetCentroid, const Vec3f &barycentricWeight, Vec3f &gridLocus);
	void vertexGridLocus(const int vertex, Vec3f &gridLocus);  // always material coords
	void vertexMaterialCoordinate(const int vertex, std::array<float, 3> &matCoord);
	inline int firstInteriorTet() { return _firstInteriorTet; }  // all tets before this are surface tets possibly virtual noded and with out unique centroid. Here on are unique interior tets.
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
	const Vec3f* getNodeSpatialCoordPointer() { if (_nodeSpatialCoords == nullptr) throw(std::logic_error("Trying to access nodeSpatialCoordinate vector before it has been allocated and assigned")); return _nodeSpatialCoords; }
	// next set of routines traverse topological paths through the bcc data
	int edgeCircumCentroids(bccTetCentroid tc, int edge, bccTetCentroid(&circumCentroids)[6]);  // gets centroids of same size as tc surrounding edge. Edges listed as sequential pairs from nodes 0 to 3, then 0 to 2, then 1 to 3.
	int faceAdjacentMultiresTet(const bccTetCentroid tc, const int face, bccTetCentroid& tcAdj);
	void nodeMicroCentroids(const std::array<short, 3>& node, bccTetCentroid (&cntrd)[24]);  // for these routines a coordinate greater than 65533 indicates a centroid out of positive grid bounds

	inline void faceNodes(const int tet, const int face, int(&nodes)[3])
	{
		const int *tn = tetNodes(tet);
		for (int i = 0; i < 3; ++i)
			nodes[i] = (tn[(face + i) & 3]);
	}

	inline bool adjacentMicrotetCentroids(const bccTetCentroid &tc0, const bccTetCentroid &tc1){
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

	void edgeNodes(const int tet, const int edge, int &n0, int &n1);  // same edge numbering as above
	int parametricTriangleTet(const int triangle, const float(&uv)[2], Vec3f& gridLocus);  // returns grid locus and tetrahedron at parametric location uv in input triangle
	int parametricEdgeTet(const int vertex0, const int vertex1, const float param, Vec3f& gridLocus);

	void centroidToNodeLoci(const bccTetCentroid& centroid, short (&gridLoci)[4][3]);

	void unitCubeCentroids(const short(&minimumCorner)[3], bccTetCentroid(&cntrd)[6]);  // centroids whose tet is part of a unit cube
	void unitCubeCentroids(const short(&minimumCorner)[3], Vec3f(&centroidLoci)[6]);
	inline int firstInteriorTetrahedron() { return _firstInteriorTet;  }
	bool insideTet(const bccTetCentroid& tc, const std::array<short, 3>& nodeLocus);  // material coords, not spatial
	bool insideTet(const bccTetCentroid& tc, const Vec3f& gridLocus);
	inline int numberOfMegaTets() { return _nMegatets; }

	// multi resolution tet utility routines
	const bccTetCentroid centroidUpOneLevel(const bccTetCentroid& tcIn);
	bool subtetCentroids(const bccTetCentroid& macroCentroid, bccTetCentroid(&subCentroids)[8]);  // invalid subtet outside positive octant labelled as all USHRT_MAX
	inline int centroidLevel(const bccTetCentroid& tc) {
		int bitNow = 1, ored = tc[0] | tc[1] | tc[2];
		for (int i = 1; i < 6; ++i) {
			if (ored & bitNow)
				return i;
			bitNow <<= 1;
		}
		return -1;  // invalid centroid greater than level 5
	}
	inline void centroidType(const bccTetCentroid& tc, int& level, int& halfCoordinate, bool& up) {
		int tbit = 1;
		for (level = 1; level < 11; ++level) {
			for (halfCoordinate = 0; halfCoordinate < 3; ++halfCoordinate) {
				if (tc[halfCoordinate] & tbit)
					break;
			}
			if (halfCoordinate < 3)
				break;
			tbit <<= 1;
		}
		if (level > 10)
			throw(std::logic_error("centroidType() sent a centroid with a level greater than 10.\n"));
		tbit <<= 1;
		if ((tbit & tc[halfCoordinate]) == (tbit & tc[(halfCoordinate + 1) % 3]))
			up = false;
		else
			up = true;
	}

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

	struct unsigned3 {
		std::array<unsigned short, 3> tc;
		unsigned short pad;
	};
	union btHash {
		uint64_t ll;
		unsigned3 us3;
	};
	struct bccTetCentroidHasher {
		std::size_t operator()(const std::array<unsigned short, 3>& k) const
		{  // hash function
			btHash bh;
			bh.us3.tc = k;
			bh.us3.pad = 0;
			std::hash<long long> hash_funct;
			return hash_funct(bh.ll);
		}
	};
	std::unordered_multimap<bccTetCentroid, int, bccTetCentroidHasher> _tetHash;  // bccTetCenter and index into _tetNodes
	typedef std::unordered_multimap<bccTetCentroid, int, bccTetCentroidHasher>::iterator THIT;
	materialTriangles *_mt;  // embedded surface
	std::vector<int> _vertexTets;
	std::vector<Vec3f> _barycentricWeights;
	int _tetSubdivisionLevels;  // if not doing multiresolution, this will be 1.  If using multiresolution, this will be the highest level of subdivision.
	Vec3f *_nodeSpatialCoords;
	Vec3f _minCorner, _maxCorner;
	double _unitSpacing, _unitSpacingInv;
	int _gridSize[3];  // Number of subdivisions in x,y, and z.
	int _firstInteriorTet;  // first tet of all unique interior tets. Goes to end of list.
	int _nMegatets;  // number of unique largest level megatets.  Guaranteed to go from 0 to this number in the tet vectors _tetNodes and _tetCentroids.
	static Mat3x3f _barycentricInverses[6];
	struct decimatedFaceNode {
		std::vector<int> faceNodes;
		std::vector<float> faceBarys;
	};
	std::unordered_map<int, decimatedFaceNode> _tJunctionConstraints;  // first is decimated node, second its parametric location on a super tet face. Data created in cutter.
	inline bool decimatedNode(const int node, const std::vector<int> *faceNodes, const std::vector<float> *faceBarys) {
		auto dn = _tJunctionConstraints.find(node);
		if (dn == _tJunctionConstraints.end())
			return false;
		faceNodes = &dn->second.faceNodes;
		faceBarys = &dn->second.faceBarys;
		return true;
	}

	int vertexSolidLinePath(const int vertex, const Vec3f materialTarget);  // if solid path found, returns tet id containing materialTarget. Else if no path, return -1;

	friend class vnBccTetCutter;
	friend class vnBccTetCutter_tbb;
	friend class remapTetPhysics;
	friend class skinCutUndermineTets;
	friend class deepCut;
};

#endif // __VN_BCC_TETS__