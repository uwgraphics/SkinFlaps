// File: multiresBccTets.h
// Author: Court Cutting MD
// Date: November 26, 2022
// Purpose: Takes a very high resolution vnBccTetrahedra object and does a binary decimation reducing each 8 bcc tet cell
//    into a single lower resolution tet based on various criteria.  This can be done in binary layers as would be done
//    in the inverse of subdivision surface.  The most common criterion preventing decimation is the presence of virtual
//    noding in the sub-tets.  This preserves the separability of solids which are close together.  Another criterion
//    preventing decimation might be a minimum number of sub-tets.  For example if only 2 subtets are present in the next
//    larger binary tet, no decimation would be done.  Other criteria limiting decimation can be introduced.
//    Later work on this class starts with creating large macrotets first.  When the user operates on the mesh, those parts
//    have their macrotets removed and the mesh within is recut at fine resolution then decimated up to join the rest of the mesh.
//    This allows for incremental incision creation while leaving the less important, periferal parts of the model at low res for speed.
//    This open source code is available with a BSD-3 license detailed here: https://opensource.org/licenses/BSD-3-Clause

#ifndef __MULTIRES_BCC_TETS__
#define __MULTIRES_BCC_TETS__

#include "vnBccTetrahedra.h"
//#include "vnBccTetCutter_tbb.h"

class multiresBccTets : public vnBccTetrahedra
{
public:
//	void createMacroTets(materialTriangles* mt, const int nLevels, const int maximumDimensionMacroSubdivs);
//	void subdivideMacroTets(std::vector<int>& macroTets);
//	void decimateOld(int level, int nSubtetsRequired, bool onlyInteriorTets = false);
	void decimate(int level);
	void getInterNodeConstraints(std::vector<int>& subNodes, std::vector<std::vector<int> >& macroNodes, std::vector<std::vector<float> >& macroBarycentrics);
	void barycentricWeightToGridLocusLevels(const bccTetCentroid& tetCentroid, const Vec3f& barycentricWeight, Vec3f& gridLocus);
	int faceAdjacentMacrotet(const bccTetCentroid &macroTc, const int face, bccTetCentroid& tcAdj);
	multiresBccTets() {}
	~multiresBccTets() {}

	inline int getLevel(const bccTetCentroid &tc) {
		int bitNow = 1, ored = tc[0] | tc[1] | tc[2];
		for (int i = 1; i < 32; ++i) {
			if (ored & bitNow)
				return i;
			bitNow <<= 1;
		}
		return -1;  // never happens
	}

private:
//	vnBccTetCutter_tbb _vbtc;

//	struct unsigned3 {
//		std::array<unsigned short, 3> tc;
//		unsigned short pad;
//	};
	struct signed3 {
		std::array<short, 3> tc;
		unsigned short pad;
	};
	union btHash {
		uint64_t ll;
//		unsigned3 us3;
		signed3 ss3;
	};
	struct arrayShort3Hasher {
		std::size_t operator()(const std::array<short, 3>& k) const
		{  // hash function
			btHash bh;
			bh.ss3.tc = k;
			bh.ss3.pad = 0;
			std::hash<long long> hash_funct;
			return hash_funct(bh.ll);
		}
	};

	std::unordered_map<int, std::pair<int, int> > _decimatedNodes;  // first is id of node decimated, second are the two nodes of the edge first node splits.
	typedef std::unordered_map<int, std::pair<int, int> >::iterator DNIT;
	struct decimatedFaceNode {
		std::vector<int> faceNodes;
		std::vector<float> faceBarys;
	};
	std::unordered_map<int, decimatedFaceNode> _decNodeConstraints;  // first is decimated node, second its parametric location on a super tet face.


	void gridLocusToBarycentricWeightLevel(const Vec3f& gridLocus, const bccTetCentroid& tc, Vec3f& barycentricWeight);
	void centroidPromoteOneLevel(const bccTetCentroid& tcLevelIn, bccTetCentroid& tc);

	friend class vnBccTetCutter;
	friend class vnBccTetCutter_tbb;
	friend class remapTetPhysics;
	friend class skinCutUndermineTets;
	friend class deepCut;
};
#endif  // #define __MULTIRES_BCC_TETS__
