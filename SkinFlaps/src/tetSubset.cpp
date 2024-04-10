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

// COURT - with multires tets this approach is flawed.  Perhaps should do only in lowest level centroids and have only subcut volumes exhibit this behavior. Currently inactive. Fix later.

bool tetSubset::createSubset(vnBccTetrahedra* vbt, const std::string objFile, float lowTetWeight, float highTetWeight, float strainMin, float strainMax) {
	materialTriangles mt;
	if (mt.readObjFile(objFile.c_str()))
		return false;
	if (mt.findAdjacentTriangles(true) > 0) {
		std::cout << "Topological error in tet subset file: " << objFile << "\n";
		return false;
	}
	_tetSubs.push_back(tetSub());
	tetSub* ts = &_tetSubs.back();
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
	// convert vertex coords to grid material coords
	boundingBox<float> bb;
	bb.Empty_Box();
	for (int n = mt.numberOfVertices(), i = 0; i < n; ++i) {
		Vec3f spat, mat;
		mt.getVertexCoordinate(i, spat.xyz);
		vbt->spatialToGridCoords(spat, mat);
		mt.setVertexCoordinate(i, mat.xyz);
		bb.Enlarge_To_Include_Point(mat.xyz);
	}
	std::set<bccTetCentroid> cents;
	auto tcArr = vbt->getTetCentroidArray();
	for (int n = tcArr.size(), i = 0; i < n; ++i) {
		if (!cents.insert(tcArr[i]).second)  // != cents.end())  // only unique
			continue;
		Vec3f V((const unsigned short (&)[3]) * tcArr[i].data());
		if (!bb.Inside(V.xyz))
			continue;
		int nIntersects = 0;
		for (int n = mt.numberOfTriangles(), i = 0; i < n; ++i) {
			int* tr = mt.triangleVertices(i);
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
			ts->subsetCentroids.push_back(tcArr[i]);  // COURT - not unique
		}
	}
	return true;
}

void tetSubset::sendTetSubsets(vnBccTetrahedra* vbt, const materialTriangles* mt, pdTetPhysics* ptp) {
	std::set<int> tets5;
	for (int n = mt->numberOfTriangles(), i = 0; i < n; ++i) {
		if (mt->triangleMaterial(i) != 5)  // at present only deep bed surfaces are effected.  Fix this.
			continue;
		const int* tr = mt->triangleVertices(i);
		for (int j = 0; j < 3; ++j)
			tets5.insert(vbt->getVertexTetrahedron(tr[j]));
	}
	for (auto& ts : _tetSubs) {
		std::vector<int> tets;
		tets.reserve(ts.subsetCentroids.size());
		for (auto& tc : ts.subsetCentroids) {
			std::list<int> tetList;
			vbt->centroidTets(tc, tetList);
			tets.insert(tets.end(), tetList.begin(), tetList.end());
		}
		ptp->tetSubset(ts.lowTetWeight, ts.highTetWeight, ts.strainMin, ts.strainMax, tets);
	}
}
