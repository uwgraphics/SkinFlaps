///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// File: tetSubset.h
// Author: Court Cutting MD
// Date: 7/15/2021
// Purpose: Inputs a virtual noded BCC tetrahedral structure. This code creates a subset of those
// tetrahedra which have a unique location (i.e. centroid) that are inside an input closed manifold surface.
// Currently it is used to specify a cluster of tets that have different physical properties than the rest.
// It can also be used to specify muscles within the volume that can be assigned contractile properties.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef __TET_SUBSET__
#define __TET_SUBSET__

#include <vector>
#include <list>

// forward declarations
class vnBccTetrahedra;
class materialTriangles;

class tetSubset
{
public:
	bool createSubset(vnBccTetrahedra* vbt, const std::string name,	float lowTetWeight,	float highTetWeight,
		float strainMin, float strainMax, const std::list<std::string> &objFiles);
	void sendTetSubsets(vnBccTetrahedra* vbt, const materialTriangles* mt, pdTetPhysics* ptp);
	tetSubset() { }
	tetSubset(const tetSubset&) = delete;
	tetSubset& operator=(const tetSubset&) = delete;
	~tetSubset() {}

private:
	struct tetSub{
		std::string name;
		float lowTetWeight;
		float highTetWeight;
		float strainMin;
		float strainMax;
		std::vector<long long> subsetCentroids;
	};
	std::list<tetSub> _tetSubs;
};

#endif  // __TET_SUBSET__
