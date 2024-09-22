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
#include "boundingBox.h"
#include "materialTriangles.h"

// forward declarations
class vnBccTetrahedra;

typedef std::array<unsigned short, 3> bccTetCentroid;

class tetSubset
{
public:
	bool createSubset(vnBccTetrahedra* vbt, const std::string objFile, float lowTetWeight, float highTetWeight, float strainMin, float strainMax);
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
		materialTriangles mt;
		std::vector<bccTetCentroid> subsetCentroids;
	};
	std::list<tetSub> _tetSubs;

	struct centLine {
		uint16_t C0;  // value at hc+1 %3
		uint16_t C1;  // value at hc+2 %3
		uint32_t hc;  // half coordinate axis
	};
	union centLineHash {
		uint64_t ll;
		centLine cl;
	};
	struct centroidLineHasher {
		std::size_t operator()(const centLine& cntrL) const
		{  // hash function
			centLineHash clh;
			clh.cl.C0 = cntrL.C0;
			clh.cl.C1 = cntrL.C1;
			clh.cl.hc = cntrL.hc;
			std::hash<uint64_t> hash_funct;
			return hash_funct(clh.ll);
		}
	};
	struct centroidLineEquals {
		bool operator()(const centLine& clL, const centLine& clR) const
		{
			if (clL.C0 != clR.C0)
				return false;
			else if (clL.C1 != clR.C1)
				return false;
			else if (clL.hc != clR.hc)
				return false;
			else
				return true;
		}
	};
	std::unordered_map<centLine, std::multimap<double, char>, centroidLineHasher, centroidLineEquals> _centroidLines;  // each centroid line keeps its intersects with whether or mot solid begins there


	void centroidLineIntersectTriangle(Vec3f(&tri)[3]);
};

#endif  // __TET_SUBSET__
