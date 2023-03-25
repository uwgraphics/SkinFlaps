// File: bccTetDecimator.h
// Author: Court Cutting MD
// Date: November 26, 2022
// Purpose: Takes a very high resolution vnBccTetrahedra object and does a binary decimation reducing each 8 bcc tet cell
//    into a single lower resolution tet based on various criteria.  This can be done in binary layers as would be done
//    in the inverse of subdivision surface.  The most common criterion preventing decimation is the presence of virtual
//    noding in the sub-tets.  This preserves the separability of solids which are close together.  Another criterion
//    preventing decimation might be a minimum number of sub-tets.  For example if only 2 subtets are present in the next
//    larger binary tet, no decimation would be done.  Other criteria limiting decimation can be introduced.
//    Copyright 2022 Court Cutting.
//    This open source code is available with a BSD-3 license detailed here: https://opensource.org/licenses/BSD-3-Clause

#ifndef __BCC_TET_DECIMATOR__
#define __BCC_TET_DECIMATOR__

#include "vnBccTetrahedra.h"

// forward declarations
class vnBccTetCutter;

class bccTetDecimator : public vnBccTetrahedra
{
public:
	void createMacroTets(vnBccTetCutter* vbtc, int nLevels);
//	void subdivideMacroTets(std::vector<int>& macroTets);
	void decimate(int level, int nSubtetsRequired, bool onlyInteriorTets = false);
	void centroidToNodeLociLevels(const bccTetCentroid& centroid, short(&gridLoci)[4][3]);
	void barycentricWeightToGridLocusLevels(const bccTetCentroid& tetCentroid, const Vec3f& barycentricWeight, Vec3f& gridLocus);

	inline int getLevel(const bccTetCentroid &tc) {
		int bitNow = 1, ored = tc[0] | tc[1] | tc[2];
		for (int i = 1; i < 32; ++i) {
			if (ored & bitNow)
				return i;
			bitNow <<= 1;
		}
		return -1;  // never happens
	}

	struct decimatedFaceNode {
		std::vector<int> faceNodes;
		std::vector<float> faceParams;
	};
	std::unordered_map<int, decimatedFaceNode> decNodeConstraints;  // first is decimated node, second its parametric location on a super tet face.

private:
	vnBccTetCutter* _vbtc;
	struct arrayShort3Hasher {
		std::size_t operator()(const std::array<short, 3>& k) const
		{  // hash function
			long long lli = k[0];
			lli <<= 16;
			lli += k[1];
			lli <<= 16;
			lli += k[2];
			std::hash<long long> hash_funct;
			return hash_funct(lli);
		}
	};
	std::unordered_map<std::array<short, 3>, int, arrayShort3Hasher> _uniqueCornerNodes;
	std::unordered_map<int, std::pair<int, int> > _decimatedNodes;  // first is id of node decimated, second are the two nodes of the edge first node splits.
	typedef std::unordered_map<int, std::pair<int, int> >::iterator DNIT;

	void gridLocusToBarycentricWeightLevel(const Vec3f& gridLocus, const bccTetCentroid& tc, Vec3f& barycentricWeight);
	void centroidPromoteOneLevel(const bccTetCentroid& tcLevelIn, bccTetCentroid& tc);
	int recurseSubtets(const bccTetCentroid& tc, const int& nNeeded, int firstInteriorTet);

};
#endif  // #define __BCC_TET_DECIMATOR__
