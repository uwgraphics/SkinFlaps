////////////////////////////////////////////////////////////////////////////
// File: vnBccTetrahedra.cpp
// Author: Court Cutting MD
// Date: 9/24/2022 revision of 7/1/2016 original
// Purpose: Basic virtual noded cubes class where cubes in space are not unique, but may be duplicated by the use of virtual nodes.
//     Full description of this concept is given in original work by Molino N,  Bao Z, and Fedkiw R: http://physbam.stanford.edu/~fedkiw/papers/stanford2004-01.pdf
////////////////////////////////////////////////////////////////////////////

#include <tuple>
#include <assert.h>
#include <algorithm>
#include <iterator>
#include <set>
#include <array>
#include <functional>
#include "Mat3x3f.h"

#include "Mat3x3d.h"

#include "boundingBox.h"
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"

// set up inverses for quick barycentric coord computation for the 6 material coord orientations.  Could also be used to compute deformation gradients if start orientation needed.
// Usually need only one if svd makes the starting orientation irrelevant since all tets have exactly the same shape and size.

Mat3x3f vnBccTetrahedra::_barycentricInverses[6] = {
	{-0.5f, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0f, 0.5f, -0.5f},
	{0.0f, 0.5f, -0.5f, -0.5f, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f},
	{0.5f, 0.0f, 0.0f, 0.0f, 0.5f, -0.5f, -0.5f, 0.5f, 0.5f},
	{0.5f, -0.5f, -0.5f, 0.5f, 0.0f, 0.0f, 0.0f, -0.5f, 0.5f},
	{0.0f, -0.5f, 0.5f, 0.5f, -0.5f, -0.5f, 0.5f, 0.0f, 0.0f},
	{0.5f, 0.0f, 0.0f, 0.0f, -0.5f, 0.5f, 0.5f, -0.5f, -0.5f} };

void vnBccTetrahedra::clear()
{
	_nodeGridLoci.clear();
	_tetNodes.clear();
	_tetCentroids.clear();
	_tetHash.clear();
	_vertexTets.clear();
	_barycentricWeights.clear();
}

void vnBccTetrahedra::centroidToNodeLoci(const bccTetCentroid& centroid, short (&gridLoci)[4][3]) {
	// fixed for multires
	int c1, c2, hc, size, levelUpBit;
	centroidHalfAxisSize(centroid, hc, size);
	levelUpBit = (size << 1);
	c1 = hc < 2 ? hc + 1 : 0;
	c2 = hc > 0 ? hc - 1 : 2;
	for (int j = 0; j < 4; ++j) {
		gridLoci[j][0] = centroid[0];
		gridLoci[j][1] = centroid[1];
		gridLoci[j][2] = centroid[2];
	}
	if ((centroid[hc] & levelUpBit) == (centroid[c2] & levelUpBit)) {  // 0-1 Cartesian axis below 2-3 
		gridLoci[0][hc] -= size;
		gridLoci[1][hc] -= size;
		gridLoci[2][hc] += size;
		gridLoci[3][hc] += size;
		gridLoci[2][c2] += levelUpBit;
		gridLoci[3][c2] -= levelUpBit;
	}
	else {
		gridLoci[0][hc] += size;
		gridLoci[1][hc] += size;
		gridLoci[2][hc] -= size;
		gridLoci[3][hc] -= size;
		gridLoci[2][c2] -= levelUpBit;
		gridLoci[3][c2] += levelUpBit;
	}
	gridLoci[0][c1] -= levelUpBit;
	gridLoci[1][c1] += levelUpBit;
	for (int j = 0; j < 4; ++j) {  // don't bit shift as possibly negative
		gridLoci[j][0] *= 0.5f;
		gridLoci[j][1] *= 0.5f;
		gridLoci[j][2] *= 0.5f;
	}
}

void vnBccTetrahedra::gridLocusToLowestTetCentroid(const Vec3f &gridLocus, bccTetCentroid &tetCentroid)
{
	short tc[3]; //  vMin, vNow;
	float dxyz[3];
	for (int i = 0; i < 3; ++i) {
		tc[i] = (short)std::floor(gridLocus[i]);
		dxyz[i] = gridLocus[i] - tc[i];
	}

#ifdef _DEBUG
	// get closest tet using centroids. Of course this works, but next version is faster.
	Vec3f cLoc[6], found, newTC;
	unitCubeCentroids(tc, cLoc);
	float d, dMin = FLT_MAX;
	for (int i = 0; i < 6; ++i) {
		d = (gridLocus - cLoc[i]).length2();
		if (d < dMin) {
			dMin = d;
			found = cLoc[i];
		}
	}
	newTC = found * 2.0001f;
	bccTetCentroid testCentroid = { (unsigned short)newTC[0], (unsigned short)newTC[1], (unsigned short)newTC[2] };
#endif

	// 4 different diagonal vectors in unit cubes give 4 different centroid patterns
	Vec3f newC, center = Vec3f(tc[0] + 0.5f, tc[1] + 0.5f, tc[2] + 0.5f);
	newC = center;
	bool split[3] = { (bool)(tc[0] & 1), (bool)(tc[1] & 1), (bool)(tc[2] & 1) };
	if (split[0] == split[1] && split[0] == split[2]) {
		if (dxyz[0] > dxyz[1] && dxyz[0] > dxyz[2]) {
			newC[0] += 0.5f;
			if (dxyz[1] > dxyz[2]) 
				newC[2] -= 0.5f;
			else
				newC[1] -= 0.5f;
		}
		else if (dxyz[1] >= dxyz[0] && dxyz[1] > dxyz[2]) {
			newC[1] += 0.5f;
			if (dxyz[0] > dxyz[2])
				newC[2] -= 0.5f;
			else 
				newC[0] -= 0.5f;
		}
		else {
			newC[2] += 0.5f;
			if (dxyz[0] > dxyz[1])
				newC[1] -= 0.5f;
			else
				newC[0] -= 0.5f;
		}
	}
	else if (split[0] == split[1]) {
		dxyz[2] = 1.0f - dxyz[2];
		if (dxyz[0] > dxyz[1] && dxyz[0] > dxyz[2]) {
			newC[0] += 0.5f;
			if (dxyz[1] > dxyz[2])
				newC[2] += 0.5f;
			else
				newC[1] -= 0.5f;
		}
		else if (dxyz[1] > dxyz[0] && dxyz[1] > dxyz[2]) {
			newC[1] += 0.5f;
			if (dxyz[0] > dxyz[2])
				newC[2] += 0.5f;
			else
				newC[0] -= 0.5f;
		}
		else {
			newC[2] -= 0.5f;
			if (dxyz[0] > dxyz[1])
				newC[1] -= 0.5f;
			else
				newC[0] -= 0.5f;
		}
	}
	else if (split[0] == split[2]) {
		dxyz[1] = 1.0f - dxyz[1];
		if (dxyz[0] > dxyz[1] && dxyz[0] > dxyz[2]) {
			newC[0] += 0.5f;
			if (dxyz[1] > dxyz[2])
				newC[2] -= 0.5f;
			else
				newC[1] += 0.5f;
		}
		else if (dxyz[1] > dxyz[0] && dxyz[1] > dxyz[2]) {
			newC[1] -= 0.5f;
			if (dxyz[0] > dxyz[2]) 
				newC[2] -= 0.5f;
			else
				newC[0] -= 0.5f;
		}
		else {
			newC[2] += 0.5f;
			if (dxyz[0] > dxyz[1])
				newC[1] += 0.5f;
			else
				newC[0] -= 0.5f;
		}
	}
	else {
		assert(split[1] == split[2]);
		dxyz[0] = 1.0f - dxyz[0];
		if (dxyz[0] > dxyz[1] && dxyz[0] > dxyz[2]) {
			newC[0] -= 0.5f;
			if (dxyz[1] > dxyz[2])
				newC[2] -= 0.5f;
			else
				newC[1] -= 0.5f;
		}
		else if (dxyz[1] > dxyz[0] && dxyz[1] > dxyz[2]) {
			newC[1] += 0.5f;
			if (dxyz[0] > dxyz[2])
				newC[2] -= 0.5f;
			else
				newC[0] += 0.5f;
		}
		else {
			newC[2] += 0.5f;
			if (dxyz[0] > dxyz[1])
				newC[1] -= 0.5f;
			else
				newC[0] += 0.5f;
		}
	}
	newC *= 2.0001f;
	tetCentroid = { (unsigned short)newC[0], (unsigned short)newC[1], (unsigned short)newC[2] };

#ifdef _DEBUG
	tetCentroid = testCentroid;
//	assert(tetCentroid == testCentroid);
#endif

//  veracity test
//	if (newC != found) {
//		assert(dxyz[0] == dxyz[1] || dxyz[0] == dxyz[2] || dxyz[1] == dxyz[2]);
//		bccTetCentroid at, tn = { (unsigned short)newC[0], (unsigned short)newC[1], (unsigned short)newC[2] };
//		int i;
//		for (i = 0; i < 4; ++i) {
//			faceAdjacentTet(tetCentroid, i, at);
//			if (at == tetCentroid)
//				break;
//		}
//		if (i > 3)
//			int junk = 0;
//	}
}

void vnBccTetrahedra::barycentricWeightToGridLocus(const int tet, const Vec3f& barycentricWeight, Vec3f& gridLocus) {
	int* tn = _tetNodes[tet].data();
	gridLocus = Vec3f((const short(&)[3]) * _nodeGridLoci[tn[0]].data()) * (1.0f - barycentricWeight[0] - barycentricWeight[1] - barycentricWeight[2]);
	for (int i = 1; i < 4; ++i)
		gridLocus += Vec3f((const short(&)[3]) * _nodeGridLoci[tn[i]].data()) * barycentricWeight[i - 1];
}

void vnBccTetrahedra::barycentricWeightToGridLocus(const bccTetCentroid &tetCentroid, const Vec3f &barycentricWeight, Vec3f &gridLocus)
{
	short gridLoci[4][3];
	centroidToNodeLoci(tetCentroid, gridLoci);
	gridLocus = Vec3f((const short(&)[3])gridLoci[0]) * (1.0f - barycentricWeight.X - barycentricWeight.Y - barycentricWeight.Z);
	for (int i = 1; i < 4; ++i)
		gridLocus += Vec3f((const short(&)[3])gridLoci[i]) * barycentricWeight[i-1];
}

void vnBccTetrahedra::vertexGridLocus(const int vertex, Vec3f &gridLocus)  // always material coords
{
	float *bw = _barycentricWeights[vertex].xyz;
	short gl[4][3];
	centroidToNodeLoci(_tetCentroids[_vertexTets[vertex]], gl);
	gridLocus = Vec3f((const short(&)[3]) gl[0]) * (1.0f - *bw - bw[1] - bw[2]);
	for (int i = 1; i < 4; ++i)
		gridLocus += Vec3f((const short(&)[3]) gl[i]) * bw[i - 1];
}

void vnBccTetrahedra::vertexMaterialCoordinate(const int vertex, std::array<float, 3> &matCoord) {
	Vec3f gridLocus;
	vertexGridLocus(vertex, gridLocus);
	gridLocus *= (float)_unitSpacing;
	gridLocus += _minCorner;
	matCoord[0] = gridLocus.X;
	matCoord[1] = gridLocus.Y;
	matCoord[2] = gridLocus.Z;
}

void vnBccTetrahedra::gridLocusToBarycentricWeight(const Vec3f &gridLocus, const bccTetCentroid &tetCentroid, Vec3f &barycentricWeight)
{  // fixed for multires
	short gl[4][3];
	centroidToNodeLoci(tetCentroid, gl);
	Vec3f V[4];
	for (int i = 0; i < 4; ++i)
		V[i] = { (float)gl[i][0], (float)gl[i][1], (float)gl[i][2] };
	for (int i = 1; i < 4; ++i)
		V[i] -= V[0];
	Mat3x3f M(V[1], V[2], V[3]);
	barycentricWeight = M.Robust_Solve_Linear_System(gridLocus - V[0]);
	assert(barycentricWeight[0] >= 0.0f && barycentricWeight[0] <= 1.0f && barycentricWeight[1]>=0.0f && barycentricWeight[1] <= 1.0f && barycentricWeight[2] >= 0.0f && barycentricWeight[2] <= 1.0f && barycentricWeight[0] + barycentricWeight[1] + barycentricWeight[2] <= 1.0f);
}

int vnBccTetrahedra::faceAdjacentMultiresTet(const bccTetCentroid tc, const int face, bccTetCentroid& tcAdj)
{  // fundamental code for all topological path routines
	// triangle faces are listed cyclic from the 4 tet nodes. Face 0 and 2 are CW, 1 & 3 CCW.
	// Returns face # of the adjacent tet.  If adjacent centroid would be outside positive octant (illegal centroid), face return is -1.
	// This routine corrected for tc of any size. tcAdj returned will be of the same size.
	int ha, size;
	centroidHalfAxisSize(tc, ha, size);
	int adjFace;
	tcAdj = tc;
	int aha = (ha + 1) % 3;
	if (((tc[ha] + tc[aha]) >> 1) & size) {  // up tet
		if (face < 1 || face > 2) {
			adjFace = 1;
			if (tcAdj[ha] < 1)
				return -1;  // next line would be below bounds tet
			tcAdj[ha] -= size;
			aha = (ha + 2) % 3;
			if (face > 2 && tcAdj[aha] < 1)
				return -1;  // next line would be below bounds tet
			tcAdj[aha] += face < 1 ? size : -size;
		}
		else {
			adjFace = face < 2 ? 3 : 0;
			tcAdj[ha] += size;
			if (face > 1 && tcAdj[aha] < 1)
				return -1;  // next line would be below bounds tet
			tcAdj[aha] += face < 2 ? size : -size;
		}
	}
	else {  // down tet
		if (face < 1 || face > 2) {
			adjFace = 2;
			tcAdj[ha] += size;
			aha = (ha + 2) % 3;
			if (face < 1 && tcAdj[aha] < 1)
				return -1;  // next line would be below bounds tet
			tcAdj[aha] += face < 1 ? -size : size;
		}
		else {
			adjFace = face < 2 ? 0 : 3;
			if (tcAdj[ha] < 1)
				return -1;  // next line would be below bounds tet
			tcAdj[ha] -= size;
			if (face > 1 && tcAdj[aha] < 1)
				return -1;  // next line would be below bounds tet
			tcAdj[aha] += face < 2 ? size : -size;
		}
	}
	return adjFace;
}

void vnBccTetrahedra::nodeMicroCentroids(std::array<short, 3>& node, bccTetCentroid cntrd[24]) {
	int count = 0;
	for (int dim = 0; dim < 3; ++dim) {
		for (int pos = -1; pos < 2; pos += 2) {
			bccTetCentroid tc2, tc = {(unsigned short)node[0], (unsigned short)node[1], (unsigned short)node[2]};
			tc[dim] += pos;
			for (int i = 0; i < 3; ++i)
				tc[i] <<= 1;
			for (int i = 0; i < 4; ++i) {
				tc2 = tc;
				int hc = ((i >> 1) + 1 + dim) % 3;
				if (i & 1)
					++tc2[hc];
				else
					--tc2[hc];
				cntrd[count++] = tc2;
			}
		}
	}
}

void vnBccTetrahedra::unitCubeCentroids(const short(&minimumCorner)[3], bccTetCentroid(&cntrd)[6]) {
	// 4 different diagonal vectors in unit cubes give 4 differentcentroid patterns
	Vec3f cl[6];
	unitCubeCentroids(minimumCorner, cl);
	for (int i = 0; i < 6; ++i) {
		cl[i] *= 2.00001f;
		cntrd[i] = { (unsigned short)cl[i][0], (unsigned short)cl[i][1], (unsigned short)cl[i][2] };
	}
}

void vnBccTetrahedra::unitCubeCentroids(const short(&minimumCorner)[3], Vec3f(&centroidLoci)[6]) {
	// 4 different diagonal vectors in unit cubes give 4 differentcentroid patterns
	bool split[3] = { (bool)(minimumCorner[0] & 1), (bool)(minimumCorner[1] & 1), (bool)(minimumCorner[2] & 1)};
	if (split[0] == split[1] && split[0] == split[2]) {
		split[0] = true;  split[1] = true; split[2] = true;
	}
	else if(split[0] == split[1]){
		split[0] = false;  split[1] = false; split[2] = true;
	}
	else if (split[0] == split[2]) {
		split[0] = false;  split[1] = true; split[2] = false;
	}
	else{
		assert(split[1] == split[2]);
		split[0] = true;  split[1] = false; split[2] = false;
	}
	Vec3f center = Vec3f(minimumCorner[0] + 0.5f, minimumCorner[1] + 0.5f, minimumCorner[2] + 0.5f);
	for (int i = 0; i < 3; ++i) {
		int c2, c1 = i < 2 ? i + 1 : 0;
		c2 = c1 < 2 ? c1 + 1 : 0;
		for (int j = 0; j < 2; ++j) {
			Vec3f& c = centroidLoci[(i << 1) + j];
			c = center;
			if (j) {
				c[c1] -= 0.5f;
				c[c2] += split[i] ? 0.5f : -0.5f;
			}
			else {
				c[c1] += 0.5f;
				c[c2] += split[i] ? -0.5f : 0.5f;
			}
		}
	}
}

void vnBccTetrahedra::edgeNodes(const int tet, const int edge, int &n0, int &n1)
{ // input one of six edges in permutation order 0-123, 1-23, and 2-3
	int *tn = _tetNodes[tet].data();
	if (edge < 3) {
		n0 = tn[0];
		n1 = tn[edge + 1];
	}
	else if (edge < 5) {
		n0 = tn[1];
		n1 = tn[edge - 1];
	}
	else {
		n0 = tn[2];
		n1 = tn[3];
	}
}

int vnBccTetrahedra::edgeCircumCentroids(bccTetCentroid tc, int edge, bccTetCentroid (&circumCentroids)[6]) {
	// edges sequential pairs from nodes 0 to 3, then 0 to 2, then 1 to 3
	short gl[4][3];
	centroidToNodeLoci(tc, gl);
	int ret;
	int N[2];
	if (edge == 0 || edge == 2) {
		N[0] = edge;  N[1] = edge + 1;
		circumCentroids[5] = { (unsigned short)(gl[edge][0] + gl[edge + 1][0]), (unsigned short)(gl[edge][1] + gl[edge + 1][1]), (unsigned short)(gl[edge][2] + gl[edge + 1][2]) };
	}
	else if (edge < 2) {
		N[0] = 1;  N[1] = 2;
	}
	else if (edge < 4) {
		N[0] = 0;  N[1] = 3;
	}
	else if (edge < 5) {
		N[0] = 0; N[1] = 2;
	}
	else {
		N[0] = 1; N[1] = 3;
	}
	if (edge < 1 || edge == 2) {  // Cartesian edges
		int size, ha;
		centroidHalfAxisSize(tc, ha, size);
		ret = 2;
		int  offset = edge < 1 ? 2 : 1;
		circumCentroids[0] = circumCentroids[5];
		circumCentroids[0][ha] += size;
		circumCentroids[1] = circumCentroids[5];
		circumCentroids[1][(ha + offset) % 3] += size;
		if (circumCentroids[5][ha] >= size) {
			circumCentroids[ret] = circumCentroids[5];
			circumCentroids[ret++][ha] -= size;
		}
		if (circumCentroids[5][(ha + offset) % 3] >= size) {
			circumCentroids[ret] = circumCentroids[5];
			circumCentroids[ret++][(ha + offset) % 3] -= size;
		}
	}
	else {
		// Diagonal edge
		ret = 0;
		int m0[3][3] = { gl[N[1]][0], gl[N[0]][1], gl[N[0]][2],
					gl[N[0]][0], gl[N[1]][1], gl[N[0]][2],
					gl[N[0]][0], gl[N[0]][1], gl[N[1]][2] };
		int m1[3][3] = { gl[N[0]][0], gl[N[1]][1], gl[N[1]][2],
					gl[N[1]][0], gl[N[0]][1], gl[N[1]][2],
					gl[N[1]][0], gl[N[1]][1], gl[N[0]][2] };
		for (int i = 0; i < 3; ++i) {
			for (int k, j = 1; j < 3; ++j) {
				int C[3];
				for (k = 0; k < 3; ++k) {
					C[k] = m0[i][k] + m1[(i + j) % 3][k];
					if (C[k] < 0)
						break;
				}
				if (k > 2)
					circumCentroids[ret++] = { (unsigned short)C[0], (unsigned short)C[1], (unsigned short)C[2] };
			}
		}
	}

#ifdef _DEBUG
	for (int i = 0; i < ret; ++i){
		short g[4][3];
		centroidToNodeLoci(circumCentroids[i], g);
		int count = 0;
		for (int j = 0; j < 4; ++j) {
			if (g[j][0] == gl[N[0]][0] && g[j][1] == gl[N[0]][1] && g[j][2] == gl[N[0]][2])
				++count;
			if (g[j][0] == gl[N[1]][0] && g[j][1] == gl[N[1]][1] && g[j][2] == gl[N[1]][2])
				++count;
		}
		if (count != 2)
			throw(std::logic_error("Program error in edgeCircumCentroids().\n"));
	}
#endif

	return ret;
}

int vnBccTetrahedra::vertexSolidLinePath(const int vertex, const Vec3f materialTarget) {  // if linear material coord path from vertex through solid to target found,
	// returns tet idx containing materialTarget. Else if no path, return -1. VN y junction hit making search impossible returns -2, but (sh/c)ould have checked all paths.
	// remember that only microtets can virtual node and be duplicated.  This facilitates search through tets whose level is above 1.
	// Any existing adjacent macrotet is guaranteed to link. Must only check duplicated level 1 tets.

	if (vertex == 3486)
		int junk = 0;

	Vec3f tmp;
	vertexGridLocus(vertex, tmp);
	Vec3d N = materialTarget - tmp;
	Vec3d vLoc(tmp);
	struct tetLink {
		int tet;
		int level;
		int face;
		double p;
	}prevTet, tetNow;
	prevTet.tet = _vertexTets[vertex];
	prevTet.level = centroidLevel(_tetCentroids[prevTet.tet]);
	prevTet.p = 0.0;
	auto nodesConnect = [&]() ->bool{
		const int *topNodes = (int *)_tetNodes[tetNow.tet].data();
		std::set<int> nodeSet;
		if (prevTet.level > tetNow.level) {
			nodeSet.insert(_tetNodes[tetNow.tet].begin(), _tetNodes[tetNow.tet].end());
			topNodes = (int*)_tetNodes[prevTet.tet].data();
		}
		else
			nodeSet.insert(_tetNodes[prevTet.tet].begin(), _tetNodes[prevTet.tet].end());
		for (int i = 0; i < 4; ++i) {
			if (nodeSet.find(topNodes[i]) != nodeSet.end())
				return true;
		}
		assert(abs(prevTet.level - tetNow.level) < 2);
		return false;
	};
	int adjFace = -1;
	auto tetIntersect = [&](const bccTetCentroid &tetCent, int& tetFace, double& lineParam) ->bool {
		short gl[4][3];
		centroidToNodeLoci(tetCent, gl);
		int maxFace = -1;
		double maxP = -1.0;
		for (tetFace = 0; tetFace < 4; ++tetFace) {
			if (tetFace == adjFace)
				continue;
			Vec3d V0(gl[tetFace]), V1(gl[(tetFace + 1) & 3]), V2(gl[(tetFace + 2) & 3]);
			V1 -= V0;
			V2 -= V0;
			Mat3x3d M(-N, V1, V2);
			Vec3d R = M.Robust_Solve_Linear_System(vLoc - V0);
			if (R[0] >= prevTet.p && R[1] >= 0.0 && R[1] <= 1.0 && R[2] >= 0.0 && R[2] <= 1.0 && R[1] + R[2] <= 1.0) {  // for single point hit would want R[0] >= prevTet.p, but adjFace would have to be always good.
				if (adjFace > -1) {
					lineParam = R[0];
					return true;
				}
				else {  // if adjFace unknown get both face intersects and get the larger p.
					if (R[0] > maxP) {
						maxP = R[0];
						maxFace = tetFace;
					}
				}
			}
		}
		if (maxFace > -1) {
			tetFace = maxFace;
			lineParam = maxP;
			return true;
		}
		return false;
	};
	std::list<tetLink> tree;
	do {
		while (prevTet.p < 1.0) {

			if (prevTet.tet == 10685)
				int junk = 0;

			if (!tetIntersect(_tetCentroids[prevTet.tet], prevTet.face, prevTet.p))
				throw(std::logic_error("Program error in vertexSolidLinePath()\n"));
			if (prevTet.p >= 1.0)
				return prevTet.tet;
			bccTetCentroid tcAdj;
			adjFace = faceAdjacentCentroid(_tetCentroids[prevTet.tet], prevTet.face, tcAdj);
			// get next tet in path and assure node link to previous, otherwise return -1
			auto pr = _tetHash.equal_range(tcAdj);
			if (pr.first == pr.second) {
				bccTetCentroid tc = tcAdj;
				tetNow.level = prevTet.level;
				tetNow.tet = -1;
				adjFace = -1;  // no longer valid if not a corner subtet. Use tetIntersect routine to find face with maximum p.
				while (tetNow.level < _tetSubdivisionLevels) {
					tc = centroidUpOneLevel(tc);
					++tetNow.level;
					auto tit = _tetHash.find(tc);  // looking up will only find unique tets 
					if (tit != _tetHash.end()) {
						tetNow.tet = tit->second;
						break;
					}
				}
				if (tetNow.tet < 0) {  // not found looking up
					// Not found up. Now look down
					if (prevTet.level < 2)  // there is no down path
						break;
					tetNow.level = prevTet.level;
					bccTetCentroid subC8[8];
					double minP;
					while (pr.first == pr.second && tetNow.level > 1) {
						if (!subtetCentroids(tcAdj, subC8))  // can't look down any further
							break;;
						--tetNow.level;
						// can't be sure which is first hit so run all 8
						int firstSubtet = 8;
						minP = DBL_MAX;
						for (int i = 0; i < 8; ++i) {
							if (subC8[i][0] == USHRT_MAX)  // subtet out of positive octant
								continue;
							double p;  // smaller tets will have a smaller p than last run
							int face;
							if (tetIntersect(subC8[i], face, p)) {
								if (p < minP) {
									minP = p;
									firstSubtet = i;
								}
							}
						}
						if (firstSubtet > 7)  // line does not pass through this subtet
							throw(std::logic_error("Program error in vertexSolidLinePath()\n"));
						tcAdj = subC8[firstSubtet];
						pr = _tetHash.equal_range(tcAdj);
					}
					if (pr.first == pr.second)  // no solid path found
						break;
					if (std::distance(pr.first, pr.second) == 1)
						tetNow.tet = pr.first->second;
					else {
						assert(tetNow.level == 1);
						while (pr.first != pr.second) {  // this must be a virtual noded multi tet level 1. No need to reset p.
							tetNow.tet = pr.first->second;
							tetNow.level = 1;
							if (nodesConnect()) {  // valid branch
								tetLink branch;
								branch.level = 1;
								branch.p = prevTet.p;
								branch.tet = pr.first->second;
								tree.push_back(branch);
							}
							++pr.first;
						}
						if (tree.empty())
							break;
						else {
							tetNow = tree.back();
							tree.pop_back();
						}
					}
				}
				// Skipping next check allowing for permissive connection here due to limited tJunction creation algorithm.  This may need fix in the future. 
//				if (!nodesConnect())  // no common node. Must look through tree for another path.
//					break;
			}
			else if (std::distance(pr.first, pr.second) == 1) {
				tetNow.level = prevTet.level;  // must be
				tetNow.tet = pr.first->second;
				// tetNow.p = p;  // to do later
				if (!nodesConnect())  // no common node. Must look through tree for another path.
					break;
			}
			else {
				assert(prevTet.level == 1);  // otherwise faceAdjacentCentroid() couldn't generate this result
				while (pr.first != pr.second) {  // this must be a virtual noded multi tet level 1. No need to reset p.
					tetNow.tet = pr.first->second;
					tetNow.level = 1;
					if (nodesConnect()) {  // valid branch
						tetLink branch;
						branch.level = 1;
						branch.p = prevTet.p;
						branch.tet = pr.first->second;
						tree.push_back(branch);
					}
					++pr.first;
				}
				if (tree.empty())
					break;
				else {
					tetNow = tree.back();
					tree.pop_back();
				}
			}
			// Successful solid path move from prevTet to tetNow
			prevTet.tet = tetNow.tet;
			prevTet.level = tetNow.level;
		}
		if (prevTet.p < 1.0f) {  // no path out on this branch
			if (!tree.empty()) { // try from this new branch point
				while (!tree.empty()) {
					tetNow = tree.back();
					tree.pop_back();
					if (!nodesConnect())  // no common node. Must look through tree for another path.
						break;
				}
			}
			else
				break;
		}
		else
			break;
	} while (true);
	if (prevTet.p > 1.0)
		return prevTet.tet;
	return -1;
}

void vnBccTetrahedra::getTJunctionConstraints(std::vector<int>& subNodes, std::vector<std::vector<int> >& macroNodes, std::vector<std::vector<float> >& macroBarycentrics) {
	size_t snSize = _tJunctionConstraints.size();
	subNodes.clear();
	macroNodes.clear();
	macroBarycentrics.clear();
	subNodes.reserve(snSize);
	macroNodes.reserve(snSize);
	macroBarycentrics.reserve(snSize);
	for (auto& dn : _tJunctionConstraints) {
		subNodes.push_back(dn.first);
		macroNodes.push_back(dn.second.faceNodes);  // don't std::move() as need these later for incisions
		macroBarycentrics.push_back(dn.second.faceBarys);
	}
}

void vnBccTetrahedra::materialCoordsToNodeSpatialVector()
{
	if (_nodeSpatialCoords == nullptr) {
		throw(std::logic_error("Trying to fill spatialCoordinateVector before it has been allocated and assigned."));
		exit(0);
	}
	for (int n = _nodeGridLoci.size(), i = 0; i < n; ++i) {
		const short* np = _nodeGridLoci[i].data();
		auto &vp = _nodeSpatialCoords[i];
		vp.set((float)np[0], (float)np[1], (float)np[2]);
		vp *= (float)_unitSpacing;
		vp += _minCorner;
	}
}

bool vnBccTetrahedra::insideTet(const bccTetCentroid& tc, const Vec3f& gridLocus) {
	int dd = 1, hc = -1;
	while (true) {
		if (tc[0] & dd) {
			hc = 0;
			break;
		}
		if (tc[1] & dd) {
			hc = 1;
			break;
		}
		if (tc[2] & dd) {
			hc = 2;
			break;
		}
		dd <<= 1;
	}
	dd <<= 1;
	int c1 = (hc + 1) % 3, c2 = (hc + 2) % 3;
	bool up = (tc[hc] & dd) == (tc[c2] & dd) ? true : false;
	if (up) {
		int V[2];
		V[0] = tc[c2] >> 1;
		if (dd < 3)
			V[1] = (tc[hc] >> 1);
		else
			V[1] = (tc[hc] >> 1) - (dd >> 2);
		if (gridLocus[c2] - gridLocus[hc] > V[0] - V[1])
			return false;
		if (-gridLocus[c2] - gridLocus[hc] > -V[0] - V[1])
			return false;
		if (dd < 3)
			++V[1];
		else
			V[1] += (dd >> 1);
		V[0] = tc[c1] >> 1;
		if (gridLocus[c1] + gridLocus[hc] > V[0] + V[1])
			return false;
		if (-gridLocus[c1] + gridLocus[hc] > -V[0] + V[1])
			return false;
	}
	else {
		int V[2];
		V[0] = tc[c2] >> 1;
		if (dd < 3)
			V[1] = (tc[hc] >> 1) + 1;
		else
			V[1] = (tc[hc] >> 1) + (dd >> 2);
		if (gridLocus[c2] + gridLocus[hc] > V[0] + V[1])
			return false;
		if (-gridLocus[c2] + gridLocus[hc] > -V[0] + V[1])
			return false;
		if (dd < 3)
			--V[1];
		else
			V[1] -= (dd >> 1);
		V[0] = tc[c1] >> 1;
		if (gridLocus[c1] - gridLocus[hc] > V[0] - V[1])
			return false;
		if (-gridLocus[c1] - gridLocus[hc] > -V[0] - V[1])
			return false;
	}
	return true;
}

bool vnBccTetrahedra::insideTet(const bccTetCentroid& tc, const std::array<short, 3>& nodeLocus) {
	int dd = 1, hc = -1;
	while (true) {
		if (tc[0] & dd) {
			hc = 0;
			break;
		}
		if (tc[1] & dd) {
			hc = 1;
			break;
		}
		if (tc[2] & dd) {
			hc = 2;
			break;
		}
		dd <<= 1;
	}
	dd <<= 1;
	int c1 = (hc + 1) % 3, c2 = (hc + 2) % 3;
	bool up = (tc[hc] & dd) == (tc[c2] & dd) ? true : false;
	if (up) {
		int V[2];
		V[0] = tc[c2] >> 1;
		if (dd < 3)
			V[1] = (tc[hc] >> 1);
		else
			V[1] = (tc[hc] >> 1) + (dd >> 2);
		if (nodeLocus[c2] - nodeLocus[hc] > V[0] - V[1])
			return false;
		if (-nodeLocus[c2] - nodeLocus[hc] > -V[0] - V[1])
			return false;
		if (dd < 3)
			++V[1];
		else
			V[1] += (dd >> 1);
		V[0] = tc[c1] >> 1;
		if (nodeLocus[c1] + nodeLocus[hc] > V[0] + V[1])
			return false;
		if (-nodeLocus[c1] + nodeLocus[hc] > -V[0] + V[1])
			return false;
	}
	else {
		int V[2];
		V[0] = tc[c2] >> 1;
		if (dd < 3)
			V[1] = (tc[hc] >> 1) + 1;
		else
			V[1] = (tc[hc] >> 1) + (dd >> 2);
		if (nodeLocus[c2] + nodeLocus[hc] > V[0] + V[1])
			return false;
		if (-nodeLocus[c2] + nodeLocus[hc] > -V[0] + V[1])
			return false;
		if (dd < 3)
			--V[1];
		else
			V[1] -= (dd >> 1);
		V[0] = tc[c1] >> 1;
		if (nodeLocus[c1] - nodeLocus[hc] > V[0] - V[1])
			return false;
		if (-nodeLocus[c1] - nodeLocus[hc] > -V[0] - V[1])
			return false;
	}
	return true;
}

int vnBccTetrahedra::parametricEdgeTet(const int vertex0, const int vertex1, const float param, Vec3f& gridLocus)
{  // new multi resolution version
	Vec3f tV[2];
	vertexGridLocus(vertex0, tV[0]);
	vertexGridLocus(vertex1, tV[1]);
	gridLocus = tV[0] * (1.0f - param) + tV[1] * param;
	bccTetCentroid tC;
	gridLocusToLowestTetCentroid(gridLocus, tC);  // new version only return a valid tc from current mesh
	// if this lowest level centroid not found, promote til found.
	int level = 1;
	assert(!_tetHash.empty());
	auto pr = _tetHash.equal_range(tC);
	while (pr.first == pr.second) {  // not found.  Move up centroid hierarchy if possible.
		tC = centroidUpOneLevel(tC);
		++level;
		if (level > 16)  // nobody would decimate this much
			throw(std::logic_error("Surface point chosen not embedded in an existing tetrahedron.\n"));
		pr = _tetHash.equal_range(tC);
	}
	if (level > 1) {  // upper level tets guaranteed not to virtual node so are unique
		assert(std::distance(pr.first, pr.second) == 1);
		return pr.first->second;
	}
	int tetOut = getVertexTetrahedron(vertex0);
	if (tC == _tetCentroids[tetOut])
		return tetOut;
	tetOut = getVertexTetrahedron(vertex1);
	if (tC == _tetCentroids[tetOut])
		return tetOut;
	// find candidate cubes
	int nTets = std::distance(pr.first, pr.second);
	if (nTets < 1) {
		assert(false);
		return -1;
	}
	else if (nTets < 2)
		return pr.first->second;
	else {
		int tetOut;
		if (param < 0.5f)
			tetOut = vertexSolidLinePath(vertex0, gridLocus);
		else
			tetOut = vertexSolidLinePath(vertex1, gridLocus);
		assert(tetOut > -1);
		return tetOut;
	}
	assert(false);
	return -1;
}

int vnBccTetrahedra::parametricTriangleTet(const int triangle, const float (&uv)[2], Vec3f& gridLocus)
{  // new multi resolution version
	Vec3f tV[3];
	const int* tr = _mt->triangleVertices(triangle);
	for (int i = 0; i < 3; ++i)
		vertexGridLocus(tr[i], tV[i]);
	gridLocus = tV[0] * (1.0f - uv[0] - uv[1]) + tV[1] * uv[0] + tV[2] * uv[1];
	bccTetCentroid tC;
	gridLocusToLowestTetCentroid(gridLocus, tC);  // new version only return a valid tc from current mesh
	// if this lowest level centroid not found, promote til found.
	int level = 1;
	assert(!_tetHash.empty());
	auto pr = _tetHash.equal_range(tC);
	while (pr.first == pr.second) {  // not found.  Move up centroid hierarchy if possible.
		tC = centroidUpOneLevel(tC);
		++level;
		if (level > 16)  // nobody would decimate this much
			throw(std::logic_error("Surface point chosen not embedded in an existing tetrahedron.\n"));
		pr = _tetHash.equal_range(tC);
	}
	if (level > 1) {  // upper level tets guaranteed not to virtual node so are unique
		assert(std::distance(pr.first, pr.second) == 1);
		return pr.first->second;
	}
	for (int i = 0; i < 3; ++i) {
		int tetOut = _vertexTets[tr[i]];
		if (tC == _tetCentroids[tetOut])
			return tetOut;
	}
	// find candidate cubes
	int nTets = std::distance(pr.first, pr.second);
	if (nTets < 1) {
		assert(false);
		return -1;
	}
	else if (nTets < 2)
		return pr.first->second;
	else {
		// look for shared neighbor node first
		for (auto tcit = pr.first; tcit != pr.second; ++tcit) {
			std::set<int> lnodes(_tetNodes[tcit->second].begin(), _tetNodes[tcit->second].end());
			for (int j, i = 0; i < 3; ++i) {
				auto& vt = _tetNodes[_vertexTets[tr[i]]];
				for (j = 0; j < 4; ++j) {
					if (lnodes.find(vt[j]) != lnodes.end())
						break;
				}
				if (j < 4)
					return tcit->second;
			}
		}
		for (int j, i = 0; i < 3; ++i) {
			int tet = vertexSolidLinePath(tr[i], gridLocus);
			if (tet > -1)
				return tet;
			else
				int junk = 0;
		}
	}
	return -1;
}

bool vnBccTetrahedra::subtetCentroids(const bccTetCentroid& macroCentroid, bccTetCentroid(&subCentroids)[8]) {
	int hc, level;
	centroidHalfAxisSize(macroCentroid, hc, level);
	if (level < 2)
		return false;
	int levelUp = level << 1, levelDown = level >> 1;
	int c1 = (hc + 1) % 3, c2 = (hc + 2) % 3;

	bool up = (macroCentroid[hc] & levelUp) == (macroCentroid[c2] & levelUp) ? true : false;

	for (int i = 0; i < 8; ++i)
		subCentroids[i] = macroCentroid;
	auto markInvalid = [&](bccTetCentroid& tc) {
		for (int j = 0; j < 3; ++j)
			tc[j] = USHRT_MAX;
	};
	// 4 corners have same hc.  List corner tets in same order as nodes to ease later processing.
	if (up) {
		subCentroids[0][hc] -= levelDown;
		subCentroids[1][hc] -= levelDown;
		subCentroids[2][hc] += levelDown;
		subCentroids[3][hc] += levelDown;
	}
	else {
		subCentroids[0][hc] += levelDown;
		subCentroids[1][hc] += levelDown;
		subCentroids[2][hc] -= levelDown;
		subCentroids[3][hc] -= levelDown;
	}
	if (subCentroids[0][c1] < level)
		markInvalid(subCentroids[0]);
	else
		subCentroids[0][c1] -= level;
	subCentroids[1][c1] += level;
	if(subCentroids[up ? 3 : 2][c2] < level)
		markInvalid(subCentroids[up ? 3 : 2]);
	else
		subCentroids[up ? 3 : 2][c2] -= level;
	subCentroids[up ? 2 : 3][c2] += level;
	// 4 core tets ring around hc axis
	if(subCentroids[4][c1] < levelDown)
		markInvalid(subCentroids[4]);
	else
		subCentroids[4][c1] -= levelDown;
	subCentroids[5][c1] += levelDown;
	if (subCentroids[6][c2] < levelDown)
		markInvalid(subCentroids[6]);
	else
		subCentroids[6][c2] -= levelDown;
	subCentroids[7][c2] += levelDown;
#ifdef _DEBUG
	for (int i = 0; i < 8; ++i) {
		if (subCentroids[i][0] == USHRT_MAX)  // || subCentroids[i][1] == USHRT_MAX || subCentroids[i][2] == USHRT_MAX)
			continue;
		assert(centroidUpOneLevel(subCentroids[i]) == macroCentroid);
	}
#endif
	return true;
}

const bccTetCentroid vnBccTetrahedra::centroidUpOneLevel(const bccTetCentroid& tcIn) {
	bccTetCentroid tcUp;
	int levelBit, levelX2, levelX4, hc;  // levelBit = 1, level = getResolutionLevel(tcIn);
	centroidHalfAxisSize(tcIn, hc, levelBit);
	levelX2 = levelBit << 1;
	levelX4 = levelX2 << 1;
	int c1 = (hc + 1) % 3, c2 = (hc + 2) % 3;
	tcUp = tcIn;
	assert((tcUp[c1] & levelX2) != (tcUp[c2] & levelX2));
	// none of the 4 core subtets have the same hc as the supertet
	// if making tc[hc] a multiple of 4 with a one unit move creates a valid level up tet, is a center core subtet
	tcUp[hc] += (tcUp[hc] & levelX2) ? levelBit : -levelBit;
	if (tcUp[c1] & levelX2) {
		if ((tcUp[hc] & levelX4) != (tcUp[c2] & levelX4))  // valid level 2
			return tcUp;
	}
	if (tcUp[c2] & levelX2) {
		if ((tcUp[hc] & levelX4) != (tcUp[c1] & levelX4))  // valid level 2
			return tcUp;
	}
	// is corner subtet and not core so will have same hc axis in level 2
	tcUp[hc] = tcIn[hc];
	tcUp[hc] += tcUp[hc] & levelX2 ? -levelBit : levelBit;
	if (tcUp[c1] & levelX2) {
		if (tcUp[c2] & levelX4) {
			tcUp[c1] += tcUp[c1] & levelX4 ? levelX2 : -levelX2;
		}
		else {
			tcUp[c1] += tcUp[c1] & levelX4 ? -levelX2 : levelX2;
		}
	}
	else {
		assert(tcUp[c2] & levelX2);
		if (tcUp[c1] & levelX4) {
			tcUp[c2] += tcUp[c2] & levelX4 ? levelX2 : -levelX2;
		}
		else {
			tcUp[c2] += tcUp[c2] & levelX4 ? -levelX2 : levelX2;
		}
	}
	return tcUp;
}

int vnBccTetrahedra::faceAdjacentCentroid(const bccTetCentroid& tc, const int face, bccTetCentroid& tcAdj)
{	// triangle faces are listed cyclic from the 4 tet nodes. Face 0 and 2 are CW, 1 & 3 CCW.
	// Returns face # of the adjacent macrotet.  If adjacent centroid would be outside positive octant (illegal centroid), face return is -1.
	int adjFace;
	tcAdj = tc;
	int ha, aha, size;
	centroidHalfAxisSize(tc, ha, size);
	tcAdj[ha] -= size;
	if (face < 1 || face >2) {
		aha = (ha + 2) % 3;
		tcAdj[aha] += size;
		if (((tc[ha] + tc[aha]) >> 1) & size) {  // down tet
			tcAdj[ha] += 2 * size;
			adjFace = 2;
			if (face < 1) {
				if (tcAdj[aha] < 2 * size)
					return -1;
				tcAdj[aha] -= 2 * size;
			}
		}
		else {
			adjFace = 1;
			if (face > 2) {
				if (tcAdj[aha] < 2 * size)
					return -1;
				tcAdj[aha] -= 2 * size;
			}
		}
	}
	else {
		aha = (ha + 1) % 3;
		tcAdj[aha] += size;
		if (face > 1) {
			if (tcAdj[aha] < 2 * size)
				return -1;
			tcAdj[aha] -= 2 * size;
		}
		if (((tc[ha] + tc[aha]) >> 1) & size) {  // up tet
			tcAdj[ha] += 2 * size;
			adjFace = face > 1 ? 0 : 3;
		}
		else
			adjFace = face > 1 ? 3 : 0;
	}
	return adjFace;
}

vnBccTetrahedra::vnBccTetrahedra() : _nodeSpatialCoords(nullptr), _firstInteriorTet(-1)
{
}


vnBccTetrahedra::~vnBccTetrahedra()
{
}
