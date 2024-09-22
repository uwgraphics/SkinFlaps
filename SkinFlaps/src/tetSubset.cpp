///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// File: tetSubset.cpp
// Author: Court Cutting MD
// Date: 8/19/2024
// Purpose: Inputs a virtual noded BCC tetrahedral structure. This code creates a subset of those
// tetrahedra which have a unique location (i.e. centroid) that are inside an input closed manifold surface.
// Currently it is used to specify a cluster of tets that have different physical properties than the rest.
// It can also be used to specify muscles within the volume that can be assigned contractile properties.
// This class has been completely rewritten to work with multiresolution tets.  It currently works only on
// deep bed tets which are level 1.  The rest of the model not in the surgical zone is unaffected - WHICH IS
// AN ERROR that should be fixed later.  This simplification is for execution speed.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <array>
#include <chrono>
#include "Vec3f.h"
#include "Vec3d.h"
#include "Mat2x2d.h"
#include "vnBccTetrahedra.h"
#include "insidePolygon.h"
#include "Mat3x3f.h"
#include "pdTetPhysics.h"
#include "tetSubset.h"

bool tetSubset::createSubset(vnBccTetrahedra* vbt, const std::string objFile, float lowTetWeight, float highTetWeight, float strainMin, float strainMax) {
	// must be careful that vbt unit spacing, etc. has been set correctly to multires grid settings before calling this routine!
	_centroidLines.clear();
	_tetSubs.push_back(tetSub());
	tetSub* ts = &_tetSubs.back();
	if (ts->mt.readObjFile(objFile.c_str()))
		return false;
	if (ts->mt.findAdjacentTriangles(true) > 0) {
		std::cout << "Topological error in tet subset file: " << objFile << "\n";
		return false;
	}
	ts->name = objFile;
	size_t front = objFile.rfind('\\');
	ts->name.erase(0, front+1);
	front = ts->name.rfind(".obj");
	ts->name.erase(front, ts->name.size());
	ts->lowTetWeight = lowTetWeight;
	ts->highTetWeight = highTetWeight;
	ts->strainMin = strainMin;
	ts->strainMax = strainMax;
	ts->subsetCentroids.clear();

//	auto tStartSteady = std::chrono::steady_clock::now();

	// convert vertex coords to grid material coords
	for (int n = ts->mt.numberOfVertices(), i = 0; i < n; ++i) {
		Vec3f spat, mat;
		ts->mt.getVertexCoordinate(i, spat.xyz);
		vbt->spatialToGridCoords(spat, mat);
		ts->mt.setVertexCoordinate(i, mat.xyz);
	}
	for (int n = ts->mt.numberOfTriangles(), i = 0; i < n; ++i) {
		int* tr = ts->mt.triangleVertices(i);
		Vec3f tri[3];
		for (int j = 0; j < 3; ++j)
			ts->mt.getVertexCoordinate(tr[j], tri[j].xyz);
		centroidLineIntersectTriangle(tri);
	}
	ts->subsetCentroids.reserve(_centroidLines.size()*10);
	for (auto& cl : _centroidLines) {
		int c0 = (cl.first.hc + 1) % 3, c1 = (cl.first.hc + 2) % 3;
		bccTetCentroid tc;
		tc[c0] = cl.first.C0 << 1;
		tc[c1] = cl.first.C1 << 1;
		auto iit = cl.second.begin();
		while (iit != cl.second.end()) {
			float bot = (float)iit->first, top;
			if (!iit->second) {  // self intersection.  Should fix model.
				++iit;
				assert(iit->second);
				++iit;
				continue;
			}
			++iit;
			while (iit->second && iit != cl.second.end()) {
				assert(iit->first - bot < 1e-3);
				++iit;
			}
			top = (float)iit->first;
			++iit;
			while (iit != cl.second.end() && iit->second < 1) {
				assert(iit->first - top < 1e-8);
				++iit;
			}
			// all half coordinate centroid values occur at .5
			float start = floor(bot) + 0.5f;
			if (bot > start)
				start += 1.0;
			while (start < top) {
				tc[cl.first.hc] = start * 2.0f + 0.00001f;
				ts->subsetCentroids.push_back(tc);
				start += 1.0;
			}
		}
	}
	ts->subsetCentroids.shrink_to_fit();
	ts->mt.clear();

	// takes 7 milliseconds to imput both cleft lip subsets
//	auto tEndSteady = std::chrono::steady_clock::now();
//	std::chrono::nanoseconds diff = tEndSteady - tStartSteady;
//	std::cout << "Time taken = " << diff.count() * 0.001 << " microseconds \n";

	return true;
}

void tetSubset::sendTetSubsets(vnBccTetrahedra* vbt, const materialTriangles* mt, pdTetPhysics* ptp) {
	if(!_centroidLines.empty())
		_centroidLines.clear();

//	auto tStartSteady = std::chrono::steady_clock::now();

	for (auto& ts : _tetSubs) {
		std::vector<int> tets;
		tets.reserve(ts.subsetCentroids.size());
		for (auto& tc : ts.subsetCentroids) {
			auto pr = vbt->_tetHash.equal_range(tc);
			if (std::distance(pr.first, pr.second) == 1)
				tets.push_back(pr.first->second);
		}
		ptp->tetSubset(ts.lowTetWeight, ts.highTetWeight, ts.strainMin, ts.strainMax, tets);
	}

	// for the 2 subsets in the unilat cleft model takes 9-12 milliseconds
//	auto tEndSteady = std::chrono::steady_clock::now();
//	std::chrono::nanoseconds diff = tEndSteady - tStartSteady;
//	std::cout << "Time taken = " << diff.count() * 0.001 << " microseconds \n";

}

void tetSubset::centroidLineIntersectTriangle(Vec3f(&tri)[3]) {
	// warning not const. Changes values in tri
	Vec3d triN, dT[2] = { tri[1] - tri[0], tri[2] - tri[0] };
	triN = dT[0] ^ dT[1];
	for (int hc = 0; hc < 3; ++hc) {
		int c0 = (hc + 1) % 3, c1 = (hc + 2) % 3;
		int xy[4] = { INT_MAX, -2, INT_MAX, -2 };
		for (int i = 0; i < 3; ++i) {
			if (tri[i][c0] < xy[0])
				xy[0] = ceil(tri[i][c0]);
			if (tri[i][c0] > xy[1])
				xy[1] = floor(tri[i][c0]);
			if (tri[i][c1] < xy[2])
				xy[2] = ceil(tri[i][c1]);
			if (tri[i][c1] > xy[3])
				xy[3] = floor(tri[i][c1]);
		}
		// now get any Z line intersects
		Mat2x2d M;
		M.x[0] = dT[0][c0];  M.x[1] = dT[0][c1];  M.x[2] = dT[1][c0]; M.x[3] = dT[1][c1];
		char solidBegin = triN[hc] < 0.0;  // negative Z starts a solid
		centLine cl;
		cl.hc = hc;
		for (cl.C0 = xy[0]; cl.C0 <= xy[1]; ++cl.C0) {
			bool odd = cl.C0 & 1;
			for (cl.C1 = xy[2]; cl.C1 <= xy[3]; ++cl.C1) {
				if (odd == (bool)(cl.C1 & 1))
					continue;
				Vec2d R = M.Robust_Solve_Linear_System(Vec2d(cl.C0 - tri[0][c0], cl.C1 - tri[0][c1]));
				if (R[0] < -1e-8 || R[0] > 1.0000001 || R[1] < -1e-8 || R[1] > 1.00000001 || R[0] + R[1] >= 1.0000001)  // don't want vertex hits, only edges
					continue;
				double intr = (tri[0][hc] + dT[0][hc] * R[0] + dT[1][hc] * R[1]);
				auto pr = _centroidLines.insert(std::make_pair(cl, std::multimap<double, char>()));
				pr.first->second.insert(std::make_pair(intr, solidBegin));
			}
		}
	}
}

