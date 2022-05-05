////////////////////////////////////////////////////////////////////////////
// File: vnBccTetrahedra.cpp
// Author: Court Cutting
// Date: 7/1/2016
// Purpose: Basic virtual noded cubes class where cubes in space are not unique, but may be duplicated by the use of virtual nodes.
//     Full description of this concept is given in original work by Molino N,  Bao Z, and Fedkiw R: http://physbam.stanford.edu/~fedkiw/papers/stanford2004-01.pdf
////////////////////////////////////////////////////////////////////////////

#include <tuple>
#include <assert.h>
#include <algorithm>
#include <set>
#include <array>
#include <functional>
#include "Mat3x3f.h"
#include "boundingBox.h"
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"

void vnBccTetrahedra::clear()
{
	_fixedNodes.clear();
	_nodeGridLoci.clear();
	_tetNodes.clear();
	_tetCentroids.clear();
	_tetHash.clear();
	_vertexTets.clear();
	_barycentricWeights.clear();
}

void vnBccTetrahedra::gridLocusToTetCentroid(const Vec3f &gridLocus, bccTetCentroid &tetCentroid)
{
	std::array<short, 3> vMin, vNow;
	// get enclosing tet
	int i;
	for (i = 0; i < 3; ++i)
		vMin[i] = (short)std::floor(gridLocus[i]);
	for (i = 0; i < 3; ++i){
		vNow[i] = vMin[i];
		float dy = gridLocus[i] - vMin[i];
		//  2 permutations
		int c2, c1 = i < 2 ? i + 1 : 0;
		c2 = c1 < 2 ? c1 + 1 : 0;
		int j;
		for (j = 0; j < 2; ++j){
			vNow[c1] = vMin[c1];
			if ((vMin[c1] & 1) == j)
				++vNow[c1];
			vNow[c2] = vMin[c2];
			if ((vMin[c2] & 1) != j)
				++vNow[c2];
			if ((vNow[i] & 1) == (vNow[c1] & 1)){
				if (fabs(gridLocus[c1] - vNow[c1]) > dy)
					continue;
				if (1.0f - fabs(gridLocus[c2] - vNow[c2]) < dy)
					continue;
			}
			else{
				if (fabs(gridLocus[c2] - vNow[c2]) > dy)
					continue;
				if (1.0f - fabs(gridLocus[c1] - vNow[c1]) < dy)
					continue;
			}
			break;
		}
		if (j < 2){
			tetCentroid.halfCoordAxis = i;
			tetCentroid.xyz = vNow;
			break;
		}
	}
	assert(i < 3);
}

void vnBccTetrahedra::barycentricWeightToGridLocus(const bccTetCentroid &tetCentroid, const Vec3f &barycentricWeight, Vec3f &gridLocus)
{
	std::array<short, 3> tet[4] = {tetCentroid.xyz, tetCentroid.xyz, tetCentroid.xyz, tetCentroid.xyz};
	--tet[0][(tetCentroid.halfCoordAxis + 1) % 3];
	++tet[1][(tetCentroid.halfCoordAxis + 1) % 3];
	bool below01 = (tetCentroid.xyz[tetCentroid.halfCoordAxis] + tetCentroid.xyz[(tetCentroid.halfCoordAxis + 2) % 3]) & 1;
	if (below01){
		++tet[0][tetCentroid.halfCoordAxis];
		++tet[1][tetCentroid.halfCoordAxis];
		--tet[2][(tetCentroid.halfCoordAxis + 2) % 3];
		++tet[3][(tetCentroid.halfCoordAxis + 2) % 3];;
	}
	else{
		++tet[2][tetCentroid.halfCoordAxis];
		++tet[3][tetCentroid.halfCoordAxis];
		++tet[2][(tetCentroid.halfCoordAxis + 2) % 3];
		--tet[3][(tetCentroid.halfCoordAxis + 2) % 3];;
	}
	gridLocus = Vec3f((const short(&)[3])*tet[0].data()) * (1.0f - barycentricWeight.X - barycentricWeight.Y - barycentricWeight.Z);
	for (int i = 1; i < 4; ++i)
		gridLocus += Vec3f((const short(&)[3])*tet[i].data()) * barycentricWeight[i-1];
}

void vnBccTetrahedra::vertexGridLocus(const int vertex, Vec3f &gridLocus)  // always material coords
{
	int *tn = _tetNodes[_vertexTets[vertex]].data();
	float *bw = _barycentricWeights[vertex]._v;
	gridLocus = Vec3f((const short(&)[3])*_nodeGridLoci[tn[0]].data()) * (1.0f - *bw - bw[1] - bw[2]);
	for (int i = 1; i < 4; ++i)
		gridLocus += Vec3f((const short(&)[3])*_nodeGridLoci[tn[i]].data()) * bw[i - 1];
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
{
	Vec3f B(gridLocus);
	// set barycentric coordinate within that tet
	Mat3x3f Minv;  // column major order.  The 6 solution matrices inverted for speed.
	if ((tetCentroid.xyz[tetCentroid.halfCoordAxis] + tetCentroid.xyz[(tetCentroid.halfCoordAxis + 1) % 3]) & 1){  // main axis below secondary
		if (tetCentroid.halfCoordAxis < 1){
			Minv.x[0] = -0.5f; Minv.x[1] = 0.5f; Minv.x[2] = 0.5f; Minv.x[3] = 0.5f; Minv.x[4] = 0.0f; Minv.x[5] = 0.0f; Minv.x[6] = 0.0f; Minv.x[7] = 0.5f; Minv.x[8] = -0.5f;
			//			M.x[0] = 0.0f; M.x[1] = 2.0f; M.x[2] = 0.0f; M.x[3] = 1.0f; M.x[4] = 1.0f; M.x[5] = 1.0f; M.x[6] = 1.0f; M.x[7] = 1.0f; M.x[8] = -1.0f;
		}
		else if (tetCentroid.halfCoordAxis < 2){
			Minv.x[0] = 0.0f; Minv.x[1] = 0.5f; Minv.x[2] = -0.5f; Minv.x[3] = -0.5f; Minv.x[4] = 0.5f; Minv.x[5] = 0.5f; Minv.x[6] = 0.5f; Minv.x[7] = 0.0f; Minv.x[8] = 0.0f;
			//			M.x[0] = 0.0f; M.x[1] = 0.0f; M.x[2] = 2.0f; M.x[3] = 1.0f; M.x[4] = 1.0f; M.x[5] = 1.0f; M.x[6] = -1.0f; M.x[7] = 1.0f; M.x[8] = 1.0f;
		}
		else{
			Minv.x[0] = 0.5f; Minv.x[1] = 0.0f; Minv.x[2] = 0.0f; Minv.x[3] = 0.0f; Minv.x[4] = 0.5f; Minv.x[5] = -0.5f; Minv.x[6] = -0.5f; Minv.x[7] = 0.5f; Minv.x[8] = 0.5f;
			//			M.x[0] = 2.0f; M.x[1] = 0.0f; M.x[2] = 0.0f; M.x[3] = 1.0f; M.x[4] = 1.0f; M.x[5] = 1.0f; M.x[6] = 1.0f; M.x[7] = -1.0f; M.x[8] = 1.0f;
		}
		B -= Vec3f((short(&)[3])*tetCentroid.xyz.data());
		B[(tetCentroid.halfCoordAxis + 1) % 3] += 1.0f;
	}
	else{
		if (tetCentroid.halfCoordAxis < 1){
			Minv.x[0] = 0.5f; Minv.x[1] = -0.5f; Minv.x[2] = -0.5f; Minv.x[3] = 0.5f; Minv.x[4] = 0.0f; Minv.x[5] = 0.0f; Minv.x[6] = 0.0f; Minv.x[7] = -0.5f; Minv.x[8] = 0.5f;
			//			M.x[0] = 0.0f; M.x[1] = 2.0f; M.x[2] = 0.0f; M.x[3] = -1.0f; M.x[4] = 1.0f; M.x[5] = -1.0f; M.x[6] = -1.0f; M.x[7] = 1.0f; M.x[8] = 1.0f;
		}
		else if (tetCentroid.halfCoordAxis < 2){
			Minv.x[0] = 0.0f; Minv.x[1] = -0.5f; Minv.x[2] = 0.5f; Minv.x[3] = 0.5f; Minv.x[4] = -0.5f; Minv.x[5] = -0.5f; Minv.x[6] = 0.5f; Minv.x[7] = 0.0f; Minv.x[8] = 0.0f;
			//			M.x[0] = 0.0f; M.x[1] = 0.0f; M.x[2] = 2.0f; M.x[3] = -1.0f; M.x[4] = -1.0f; M.x[5] = 1.0f; M.x[6] = 1.0f; M.x[7] = -1.0f; M.x[8] = 1.0f;
		}
		else{
			Minv.x[0] = 0.5f; Minv.x[1] = 0.0f; Minv.x[2] = 0.0f; Minv.x[3] = 0.0f; Minv.x[4] = -0.5f; Minv.x[5] = 0.5f; Minv.x[6] = 0.5f; Minv.x[7] = -0.5f; Minv.x[8] = -0.5f;
			//			M.x[0] = 2.0f; M.x[1] = 0.0f; M.x[2] = 0.0f; M.x[3] = 1.0f; M.x[4] = -1.0f; M.x[5] = -1.0f; M.x[6] = 1.0f; M.x[7] = 1.0f; M.x[8] = -1.0f;
		}
		B -= Vec3f((short(&)[3])*tetCentroid.xyz.data());
		B[tetCentroid.halfCoordAxis] -= 1.0f;
		B[(tetCentroid.halfCoordAxis + 1) % 3] += 1.0f;
	}
	barycentricWeight = Minv * B;
}

int vnBccTetrahedra::faceAdjacentTet(const bccTetCentroid tc, const int face, bccTetCentroid &tcAdj)
{  // fundamental code for all topological path routines
	// triangle faces are listed cyclic from the 4 tet nodes. Face 0 and 2 are CW, 1 & 3 CCW.
	// Returns face # of the adjacent tet.
	int adjFace;
	tcAdj = tc;
	if (face < 1 || face >2){
		tcAdj.halfCoordAxis = (tc.halfCoordAxis + 2) % 3;
		if ((tc.xyz[tc.halfCoordAxis] + tc.xyz[tcAdj.halfCoordAxis]) & 1){  // down tet
			++tcAdj.xyz[tc.halfCoordAxis];
			adjFace = 2;
			if (face < 1)
				--tcAdj.xyz[tcAdj.halfCoordAxis];
		}
		else{
			adjFace = 1;
			if (face > 2)
				--tcAdj.xyz[tcAdj.halfCoordAxis];
		}
	}
	else{
		tcAdj.halfCoordAxis = (tc.halfCoordAxis + 1) % 3;
		if (face > 1)
			--tcAdj.xyz[tcAdj.halfCoordAxis];
		if ((tc.xyz[tc.halfCoordAxis] + tc.xyz[tcAdj.halfCoordAxis]) & 1){  // up tet
			++tcAdj.xyz[tc.halfCoordAxis];
			adjFace = face > 1 ? 0 : 3;
		}
		else
			adjFace = face > 1 ? 3 : 0;
	}
	return adjFace;
}

int vnBccTetrahedra::faceAdjacentTets(const int tet, const int face, std::list<int> &adjTets)
{  // returns adjacent face index 0-3
	adjTets.clear();
	bccTetCentroid tcAdj;
	int adjFace = faceAdjacentTet(_tetCentroids[tet], face, tcAdj);
	const int *tn = tetNodes(tet);
	std::set<int> faceNodes;
	for (int i = 0; i < 3; ++i)
		faceNodes.insert(tn[(face + i) & 3]);
	auto tr = _tetHash.equal_range(tcAdj.ll);
	while (tr.first != tr.second){
		tn = tetNodes(tr.first->second);
		int i;
		for (i = 0; i < 3; ++i){
			if(faceNodes.find(tn[(adjFace + i) & 3]) == faceNodes.end())
				break;
		}
		if (i > 2)
			adjTets.push_back(tr.first->second);
		++tr.first;
	}
	return adjFace;
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

void vnBccTetrahedra::edgeAdjacentTets(const int tet, const int edge, std::list<int> &adjTets)
{ // input one of six edges in permutation order 0-123, 1-23, and 2-3
	adjTets.clear();
	int n0, n1, *tn = _tetNodes[tet].data();
	edgeNodes(tet, edge, n0, n1);
/*	if (edge < 3){
		n0 = tn[0];
		n1 = tn[edge + 1];
	}
	else if (edge < 5){
		n0 = tn[1];
		n1 = tn[edge - 1];
	}
	else{
		n0 = tn[2];
		n1 = tn[3];
	} */
	bccTetCentroid tc = _tetCentroids[tet];
	bool tetUp = (tc.xyz[(tc.halfCoordAxis + 1) % 3] + tc.xyz[tc.halfCoordAxis]) & 1;
	std::list<bccTetCentroid> nTets1;
	if (edge < 1){  // Cartesian edge
		if (!tetUp)
			++tc.xyz[tc.halfCoordAxis];
		for (int i = 0; i < 4; ++i){
			nTets1.push_back(tc);
			if (i & 1)
				nTets1.back().halfCoordAxis = (tc.halfCoordAxis + 2) % 3;
			if (i > 1)
				--nTets1.back().xyz[nTets1.back().halfCoordAxis];
		}
	}
	else if (edge > 4){  // edge==5 Cartesian edge
		// set midpoint of edge
		if (tetUp)
			++tc.xyz[tc.halfCoordAxis];
		for (int i = 0; i < 4; ++i){
			nTets1.push_back(tc);
			if (i > 1)
				nTets1.back().halfCoordAxis = (tc.halfCoordAxis + 1) % 3;
			if (i & 1)
				--nTets1.back().xyz[nTets1.back().halfCoordAxis];
		}
	}
	else{
		//  Get 6 diagonal centroid neighbors
		int c2, c1, c0 = tc.halfCoordAxis;
		c1 = c0 < 2 ? c0 + 1 : 0;
		c2 = c1 < 2 ? c1 + 1 : 0;
		for (int i = 0; i < 6; ++i){
			nTets1.push_back(tc);
			nTets1.back().halfCoordAxis = (c0 + (i>>1)) % 3;
			if (edge < 2){
				if (tetUp){
					if (i < 1){
						++nTets1.back().xyz[c2];
						--nTets1.back().xyz[c1];
					}
					else if (i>1 && i < 4){
						--nTets1.back().xyz[c1];
						i & 1 ? ++nTets1.back().xyz[c2] : ++nTets1.back().xyz[c0];
					}
					else if (i > 4){
						--nTets1.back().xyz[c1];
						++nTets1.back().xyz[c0];
					}
					else;
				}
				else{
					if (i < 1){
						--nTets1.back().xyz[c2];
						--nTets1.back().xyz[c1];
					}
					else if (i>1 && i < 4){
						--nTets1.back().xyz[c1];
						if (i & 1){
							++nTets1.back().xyz[c0];
							--nTets1.back().xyz[c2];
						}
					}
					else if(i>3){
						--nTets1.back().xyz[c2];
						i & 1 ? --nTets1.back().xyz[c1] : ++nTets1.back().xyz[c0];
					}
					else;
				}
			}
			else if (edge < 3){
				if (tetUp){
					if (i < 1){
						--nTets1.back().xyz[c1];
						--nTets1.back().xyz[c2];
					}
					else if (i>1 && i < 4){
						--nTets1.back().xyz[c1];
						i & 1 ? --nTets1.back().xyz[c2] : ++nTets1.back().xyz[c0];
					}
					else if (i > 3){
						--nTets1.back().xyz[c2];
						if (i & 1){
							--nTets1.back().xyz[c1];
							++nTets1.back().xyz[c0];
						}
					}
					else;
				}
				else{
					if (i < 1){
						++nTets1.back().xyz[c2];
						--nTets1.back().xyz[c1];
					}
					else if (i>1 && i < 4){
						--nTets1.back().xyz[c1];
						if (i & 1){
							++nTets1.back().xyz[c0];
							++nTets1.back().xyz[c2];
						}
					}
					else if (i>3){
						i & 1 ? --nTets1.back().xyz[c1] : ++nTets1.back().xyz[c0];
					}
					else;
				}
			}
			else if (edge < 4){
				if (tetUp){
					if (i < 1){
						++nTets1.back().xyz[c1];
						++nTets1.back().xyz[c2];
					}
					else if (i>1 && i < 4){
						i & 1 ? ++nTets1.back().xyz[c2] : ++nTets1.back().xyz[c0];
					}
					else if (i > 4){
						++nTets1.back().xyz[c1];
						++nTets1.back().xyz[c0];
					}
					else;
				}
				else{
					if (i < 1){
						++nTets1.back().xyz[c1];
						--nTets1.back().xyz[c2];
					}
					else if (i == 2){
						--nTets1.back().xyz[c2];
						++nTets1.back().xyz[c0];
					}
					else if (i>3){
						--nTets1.back().xyz[c2];
						i & 1 ? ++nTets1.back().xyz[c0] : ++nTets1.back().xyz[c1];
					}
					else;
				}
			}
			else{
				if (tetUp){
					if (i < 1){
						++nTets1.back().xyz[c1];
						--nTets1.back().xyz[c2];
					}
					else if (i>1 && i < 4){
						i & 1 ? --nTets1.back().xyz[c2] : ++nTets1.back().xyz[c0];
					}
					else if (i > 3){
						--nTets1.back().xyz[c2];
						if (i & 1){
							++nTets1.back().xyz[c0];
							++nTets1.back().xyz[c1];
						}
					}
					else;
				}
				else{
					if (i < 1){
						++nTets1.back().xyz[c1];
						++nTets1.back().xyz[c2];
					}
					else if (i == 2){
						++nTets1.back().xyz[c0];
						++nTets1.back().xyz[c2];
					}
					else if (i > 3){
						i & 1 ? ++nTets1.back().xyz[c1] : ++nTets1.back().xyz[c0];
					}
					else;
				}

			}
		}
	}
	for (auto &nt : nTets1){
		auto er = _tetHash.equal_range(nt.ll);
		while (er.first != er.second){
			tn = _tetNodes[er.first->second].data();
			int count = 0;
			for (int k = 0; k < 4; ++k){
				if (tn[k] == n0 || tn[k] == n1){
					++count;
					if (count > 1){
						adjTets.push_back(er.first->second);
						break;
					}
				}
			}
			++er.first;
		}
	}
}

bool vnBccTetrahedra::decreasingCentroidPath(const int startTet, const int targetTet, std::list<int> &tetPath)
{  // true if constantly decreasing distance centroid path exists.
	Vec3f loc, target;
	bccTetCentroid *tc = &_tetCentroids[targetTet];
	target.set((float)tc->xyz[0], (float)tc->xyz[1], (float)tc->xyz[2]);
	target[tc->halfCoordAxis] += 0.5f;
	int tetNow = startTet;
	tc = &_tetCentroids[startTet];
	loc.set((float)tc->xyz[0], (float)tc->xyz[1], (float)tc->xyz[2]);
	loc[tc->halfCoordAxis] += 0.5f;
	float d2, d2min = (loc - target).length2();
	bool moved;
	tetPath.clear();
	struct branch{
		int lastTetTried;
		std::list<int> alternateTets;
	};
	int bestAdjFace = -1;
	std::list<branch> branches;
	while (tetNow != targetTet){
		std::list<int> adjTets, bestAdjTets;
		moved = false;
		for (int i = 0; i < 4; ++i){
			int adjFace = faceAdjacentTets(tetNow, i, adjTets);
			if (i == bestAdjFace || adjTets.empty())
				continue;
			tc = &_tetCentroids[adjTets.front()];
			loc.set((float)tc->xyz[0], (float)tc->xyz[1], (float)tc->xyz[2]);
			loc[tc->halfCoordAxis] += 0.5f;
			if ((d2 = (loc - target).length2()) < d2min){
				d2min = d2;
				tetNow = adjTets.front();
				bestAdjTets = std::move(adjTets);
				bestAdjFace = adjFace;
				moved = true;
			}
		}
		if (moved){
			if (bestAdjTets.size() > 1){
				branches.push_back(branch());
				branches.back().lastTetTried = tetNow;
				bestAdjTets.pop_front();
				branches.back().alternateTets = std::move(bestAdjTets);
			}
			tetPath.push_back(tetNow);
		}
		else{
			// try moving up branch list
			if (branches.empty())
				return false;
			while (tetPath.back() != branches.back().lastTetTried)  // COURT - not debugged yet
				tetPath.pop_back();
			tetNow = branches.back().alternateTets.front();
			tetPath.back() = tetNow;
			branches.back().alternateTets.pop_front();
			if (branches.back().alternateTets.empty())
				branches.pop_back();
		}
	}
	return true;
}

void vnBccTetrahedra::materialCoordsToSpatialVector()
{
	if (_nodeSpatialCoords == nullptr) {
		throw(std::logic_error("Trying to fill spatialCoordinateVector before it has been allocated and assigned."));
		exit(0);
	}
	for (int n = _nodeGridLoci.size(), i = 0; i < n; ++i) {
		const short *np = _nodeGridLoci[i].data();
		Vec3f *vp = &_nodeSpatialCoords[i];
		vp->set((float)np[0], (float)np[1], (float)np[2]);
		*vp *= (float)_unitSpacing;
		*vp += _minCorner;
	}
}


vnBccTetrahedra::vnBccTetrahedra() : _nodeSpatialCoords(nullptr)
{
}


vnBccTetrahedra::~vnBccTetrahedra()
{
}
