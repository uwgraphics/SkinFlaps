///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// File: tetSubset.cpp
// Author: Court Cutting MD
// Date: 7/15/2021
// Purpose: Inputs a virtual noded BCC tetrahedral structure. This code creates a subset of those
// tetrahedra which have a unique location (i.e. centroid) that are inside an input closed manifold surface.
// Currently it is used to specify a cluster of tets that have different physical properties than the rest.
// It can also be used to specify muscles within the volume that can be assigned contractile properties.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <array>
#include "Vec3f.h"
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"
#include "boundingBox.h"
#include "insidePolygon.h"
#include "Mat3x3f.h"
#include "pdTetPhysics.h"
#include "tetSubset.h"

bool tetSubset::createSubset(vnBccTetrahedra* vbt, const std::string name, float lowTetWeight, float highTetWeight, float strainMin, float strainMax,
	const std::list<std::string>& objFiles) {
	// creates a subset of the overall tet volume in vbt.  Various uses including tissue subtypes
	std::vector<materialTriangles> mtVec;
	mtVec.reserve(objFiles.size());
	for (auto& of : objFiles) {
		mtVec.push_back(materialTriangles());
		if (mtVec.back().readObjFile(of.c_str()))
			return false;
	}
	_tetSubs.push_back(tetSub());
	tetSub* ts = &_tetSubs.back();
	ts->name = name;
	ts->lowTetWeight = lowTetWeight;
	ts->highTetWeight = highTetWeight;
	ts->strainMin = strainMin;
	ts->strainMax = strainMax;
	ts->subsetCentroids.clear();
	for (auto& mt : mtVec) {
		// convert vertex coords to grid material coords
		boundingBox<float> bb;
		bb.Empty_Box();
		for (int n = mt.numberOfVertices(), i = 0; i < n; ++i) {
			Vec3f spat, mat;
			mt.getVertexCoordinate(i, spat._v);
			vbt->spatialToGridCoords(spat, mat);
			mt.setVertexCoordinate(i, mat._v);
			bb.Enlarge_To_Include_Point(mat._v);
		}
		auto th = vbt->getTetHash();
		long long last;
		last = LONG_MAX;
		last <<= 32;  // impossible centroid
		for (auto hit = th->begin(); hit != th->end(); ++hit) {
			if (hit->first == last)
				continue;
			bccTetCentroid tc;
			tc.ll = hit->first;
			Vec3f V(tc.xyz[0], tc.xyz[1], tc.xyz[2]);
			if (tc.halfCoordAxis < 1)
				V.X += 0.5f;
			else if (tc.halfCoordAxis < 2)
				V.Y += 0.5f;
			else
				V.Z += 0.5f;
			if (!bb.Inside(V._v)) {
				last = hit->first;
				continue;
			}
			int nIntersects = 0;
			for (int n = mt.numberOfTriangles(), i = 0; i < n; ++i) {
				long* tr = mt.triangleVertices(i);
				Vec3f triV[3];
				for (int j = 0; j < 3; ++j) {
					float* fp = mt.vertexCoordinate(tr[j]);
					triV[j].set(fp[0], fp[1], fp[2]);
				}
				if (triV[0].X < V.X && triV[1].X < V.X && triV[2].X < V.X)
					continue;
				if (triV[0].X > V.X && triV[1].X > V.X && triV[2].X > V.X)
					continue;
				if (triV[0].Y < V.Y && triV[1].Y < V.Y && triV[2].Y < V.Y)
					continue;
				if (triV[0].Y > V.Y && triV[1].Y > V.Y && triV[2].Y > V.Y)
					continue;
				triV[1] -= triV[0];
				triV[2] -= triV[0];
				Mat3x3f M;
				M.Initialize_With_Column_Vectors(triV[1], triV[2], Vec3f(0.0f, 0.0f, -1.0f));
				Vec3f R = M.Robust_Solve_Linear_System(V - triV[0]);
				if (R.Z < 0.0f || R.X < 0.0f || R.X > 1.0f || R.Y < 0.0f || R.Y > 1.0f)
					continue;
				++nIntersects;
			}
			if (nIntersects & 1) {
				ts->subsetCentroids.push_back(tc.ll);
			}
			last = hit->first;
		}
	}
	return true;
}

void tetSubset::sendTetSubsets(vnBccTetrahedra* vbt, const materialTriangles* mt, pdTetPhysics* ptp) {
	std::set<int> tets5;
	for (int n = mt->numberOfTriangles(), i = 0; i < n; ++i) {
		if (mt->triangleMaterial(i) != 5)
			continue;
		const long* tr = mt->triangleVertices(i);
		for (int j = 0; j < 3; ++j)
			tets5.insert(vbt->getVertexTetrahedron(tr[j]));
	}
	for (auto& ts : _tetSubs) {
		std::vector<int> tets;
		tets.reserve(ts.subsetCentroids.size());
		for (auto& tc : ts.subsetCentroids) {
			int num = vbt->getTetHash()->count(tc);
			auto tit = vbt->getTetHash()->find(tc);
			if (num == 1)
				tets.push_back(tit->second);
			else if (num == 2) {
				for (int j = 0; j < 2; ++j) {
					if (tets5.find(tit->second) != tets5.end())
						tets.push_back(tit->second);
					++tit;
				}
			}
			else
				;
		}
		ptp->tetSubset(ts.lowTetWeight, ts.highTetWeight, ts.strainMin, ts.strainMax, tets);
	}
}
