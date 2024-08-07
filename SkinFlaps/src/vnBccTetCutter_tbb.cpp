#include <assert.h>
#include <limits>
#include <algorithm>
#include <cmath>
#include <functional>
#include "Mat2x2d.h"
#include "Vec3d.h"
#include "Mat2x2f.h"
#include "Mat3x3f.h"
#include "Mat3x3d.h"
#include "triTriIntersect_ShenD.h"

#include <chrono>
#include <ctime>  
#include <iostream>
#include <fstream>

#include "remapTetPhysics.h"  // for high speed spatial coord remapping after topological changes
#include "vnBccTetCutter_tbb.h"

void vnBccTetCutter_tbb::addNewMultiresIncision() {
	// Extend _vMatCoords to new vertices. Things may have already been moved so use old tet locations assigned in incision tool for new vertices
	_lastVertexSize = _vMatCoords.size();
	_vMatCoords.insert(_vMatCoords.end(), _mt->numberOfVertices() - _lastVertexSize, Vec3f());
	_vertexTetCentroids.insert(_vertexTetCentroids.end(), _mt->numberOfVertices() - _lastVertexSize, bccTetCentroid());
	for (int n = _mt->numberOfVertices(), i = _lastVertexSize; i < n; ++i) {  // incision process assigned an old tet and a weight to new incision vertices
		int tet = _vbt->getVertexTetrahedron(i);
		if (tet < -1)  // excised vertex
			continue;
		_vbt->barycentricWeightToGridLocus(tet, *_vbt->getVertexWeight(i), _vMatCoords[i]);
		_vbt->gridLocusToLowestTetCentroid(_vMatCoords[i], _vertexTetCentroids[i]);
	}
	auto getTriangleVertexCentroids = [&](int triangle) {
		if (_mt->triangleMaterial(triangle) < 0)	// signals a deleted triangle
			return;
		int* tr = _mt->triangleVertices(triangle);
		for (int j = 0; j < 3; ++j) {
			if (!(_vertexTetCentroids[tr[j]][0] < USHRT_MAX))
				_vbt->gridLocusToLowestTetCentroid(_vMatCoords[tr[j]], _vertexTetCentroids[tr[j]]);  // All tri vertices must be converted to lowest.  This invalidates next step.
		}
	};
	for (int n = _mt->numberOfTriangles(), i = _lastTriangleSize; i < n; ++i) // all new triangles part of an incision
		getTriangleVertexCentroids(i);
	// COURT 4x faster than single thread using tbb hash container requiring no reduction. tbb version ~30% faster than omp before reduction and reduction using critical section. tbb hash container very helpful.
	// get _centTris only of the new incision triangles
	_centTris.clear();
#if defined( _DEBUG )
	for (int i = _lastTriangleSize; i < _mt->numberOfTriangles(); ++i) {
		if (_mt->triangleMaterial(i) < 0)	// signals a deleted triangle
			continue;
		inputTriangleTetsTbb(i, _centTris);
	}
#else
	tbb::parallel_for(
		tbb::blocked_range<size_t>(_lastTriangleSize, _mt->numberOfTriangles()),
		[&](tbb::blocked_range<size_t> r) {
			for (size_t i = r.begin(); i != r.end(); ++i) {
				if (_mt->triangleMaterial(i) < 0)	// signals a deleted triangle
					continue;
				inputTriangleTetsTbb(i, _centTris);
			}
		});
#endif
	// Look for and delete any new megatets which have new incision triangles penetrating them that will need to be recut.
	std::unordered_set <bccTetCentroid, bccTetCentroidHasher> incisMegaCentroids;
	auto possibleMegatetReduction = [&](bccTetCentroid& tc) {
		auto mtit = _megatetTetTris.find(tc);
		if (mtit != _megatetTetTris.end()) { // &&	_vbt->_tetNodes[mtit->second.tetIdx][0] > -1)
			_vnCentroids.push_back(tc);
			assert(_vbt->_tetNodes[mtit->second.tetIdx][0] > -1);
			_vbt->_tetNodes[mtit->second.tetIdx][0] = -1;  // mark for deletion
			_vnTris.insert(mtit->second.tris.begin(), mtit->second.tris.end());
			_megatetTetTris.erase(mtit);
		}
	};
	for (auto& ct : _centTris) {
		auto tc = ct.first;
		for (int j = 1; j < _vbt->_tetSubdivisionLevels; ++j)
			tc = _vbt->centroidUpOneLevel(tc);
		if(incisMegaCentroids.insert(tc).second)
			possibleMegatetReduction(tc);
	}
	// add a megatet border around the incision to be subdivided to soften any stark junction between microtets and megatets.
	for (auto& ic : incisMegaCentroids) {
		for (int i = 0; i < 6; ++i) {
			bccTetCentroid cC[6];
			int nCC = _vbt->edgeCircumCentroids(ic, i, cC);
			for (int j = 0; j < nCC; ++j) {
				possibleMegatetReduction(cC[j]);
			}
		}
	}
	incisMegaCentroids.clear();
	std::vector<int> borderTris(_vnTris.begin(), _vnTris.end());
	for (auto t : borderTris) {
		assert(t < _lastTriangleSize);
		getTriangleVertexCentroids(t);
	}

#if defined( _DEBUG )
	for (int i = 0; i < borderTris.size(); ++i) {
		if (_mt->triangleMaterial(borderTris[i]) < 0)	// signals a deleted triangle
			continue;
		inputTriangleTetsTbb(borderTris[i], _centTris);
	}
#else
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, borderTris.size()),
		[&](tbb::blocked_range<size_t> r) {
			for (size_t i = r.begin(); i != r.end(); ++i) {
			if (_mt->triangleMaterial(borderTris[i]) < 0)	// signals a deleted triangle
				continue;
			inputTriangleTetsTbb(borderTris[i], _centTris);
		}
	});
#endif

	// save these new triangles for reprocessing next time
	for (int n = _mt->numberOfTriangles(), i = _lastTriangleSize; i < n; ++i) {
		if (_mt->triangleMaterial(i) < 0)
			continue;
		_vnTris.insert(i);
	}

	_surfaceCentroids.clear();
	_surfaceTetTris.clear();
	_vbt->_tetNodes.erase(_vbt->_tetNodes.begin() + _vbt->_nMegatets, _vbt->_tetNodes.end());
	_vbt->_tetCentroids.erase(_vbt->_tetCentroids.begin() + _vbt->_nMegatets, _vbt->_tetCentroids.end());
	_vbt->_nodeGridLoci.erase(_vbt->_nodeGridLoci.begin() + _meganodeSize, _vbt->_nodeGridLoci.end());
	// now perform remainder of recut operation
	macrotetRecutCore();
}

void vnBccTetCutter_tbb::macrotetRecutCore() {
	// reused for multiple incisions
	pack();  // removes all tets and nodes marked for deletion leaving only megatets
	_vbt->_nMegatets = _vbt->_tetNodes.size();  // reduced after pack
	_meganodeSize = _vbt->_nodeGridLoci.size();
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_vbt->_nMegatets * 1.5);
	for (int n = _vbt->_tetCentroids.size(), i = 0; i < n; ++i)
		_vbt->_tetHash.insert(std::make_pair(_vbt->_tetCentroids[i], i));  // at this time only hash unique megatets
	// get unique tet faces at the boundary of object and of the virtual noded tets that were removed in contact with tets that remain.
	_megatetBounds.clear();
	_megatetBounds.reserve(_vnCentroids.size() << 2);  // COURT check rough guess later
	std::unordered_map<int, std::set<tetTris*> > bnTris;
	bnTris.reserve(_vnCentroids.size() << 2);  // COURT check rough guess later
	for (auto& vnc : _vnCentroids) {
		for (int j = 0; j < 4; ++j) {
			bccTetCentroid adjTc;
			int adjFace = _vbt->faceAdjacentMultiresTet(vnc, j, adjTc);
			if (adjFace > -1) {  // -1 is an object boundary face, not a borderFace
				auto mttit = _megatetTetTris.find(adjTc);
				if (mttit != _megatetTetTris.end()) {
					auto& tn = _vbt->_tetNodes[mttit->second.tetIdx];
					megatetFace mf;
					mf.nodes[1] = tn[(adjFace + 1) & 3];
					if (adjFace & 1) {
						mf.nodes[0] = tn[adjFace];
						mf.nodes[2] = tn[(adjFace + 2) & 3];
					}
					else {
						mf.nodes[2] = tn[adjFace];
						mf.nodes[0] = tn[(adjFace + 2) & 3];
					}
					mf.tris = &mttit->second.tris;
					_megatetBounds.push_back(mf);
					for (auto& t : mf.nodes) {
						auto bntIt = bnTris.insert(std::make_pair(t, std::set<tetTris*>())).first;
						if (!mttit->second.tris.empty())
							bntIt->second.insert(&mttit->second);
					}
				}
			}
		}
	}
	_boundingNodeData.clear();
	_boundingNodeData.reserve(bnTris.size() * 1.1f);
	for (auto& bnt : bnTris) {
		auto pr = _boundingNodeData.insert(std::make_pair(_vbt->_nodeGridLoci[bnt.first], boundingNodeTris()));
		pr.first->second.node = bnt.first;
		pr.first->second.megaTetTris.assign(bnt.second.begin(), bnt.second.end());
	}
	bnTris.clear();

	// clear Z intersect arrays for finding interior nodes
	// am reusing these structures for recuts
	for (auto& eXy : evenXy) {
		for (auto& xy : eXy)
			xy.clear();
	}
	for (auto& oXy : oddXy) {
		for (auto& xy : oXy)
			xy.clear();
	}

	// macrotets will be guaranteed not to virtual node.  Subcut any found at this stage.
	// build microtets in deleted solid
	// now need all interior nodes inside the recut volume.  Unfortunately bounding tris of that recut volume may be partially empty.
	// So must get all interior nodes from entire volume and keep only those inside the recut volume.
	for (auto& vnc : _vnCentroids)
		addCentroidMicronodesZ(vnc);
	_zIntr.clear();
	// COURT time this to see if multithreading worth it.
	for (int n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
		if (_mt->triangleMaterial(i) < 0)
			continue;
		Vec3d triVec[3];
		int* tr = _mt->triangleVertices(i);
		for (int i = 0; i < 3; ++i)
			triVec[i] = _vMatCoords[tr[i]];
		zIntersectTriangleTbb(triVec, true, _zIntr);
	}
	for (auto ziv : _zIntr) {
		if (ziv.flags.odd)
			oddXy[ziv.x][ziv.y].insert(std::make_pair(ziv.zInt, ziv.flags));
		else
			evenXy[ziv.x][ziv.y].insert(std::make_pair(ziv.zInt, ziv.flags));
	}
	_zIntr.clear();
	std::vector<tetTriangles> tetTriVec;
	tetTriVec.reserve(_centTris.size());
	_surfaceCentroids.clear();
	_surfaceCentroids.reserve(_centTris.size());
	for (auto ctit = _centTris.begin(); ctit != _centTris.end(); ++ctit) {
		_surfaceCentroids.insert(std::make_pair(ctit->first, std::vector<int>()));
		tetTriVec.push_back(tetTriangles());
		tetTriVec.back().tc = ctit->first;
		tetTriVec.back().tris = std::move(ctit->second);
	}
	_centTris.clear();
	_interiorNodes.clear();
	createInteriorMicronodes();
	// Some of the _tetTris may be invalid if they are outside the recut volume.
	for (int n = tetTriVec.size(), i = 0; i < n; ++i) {
		auto tc = tetTriVec[i].tc;
		for (int j = 1; j < _vbt->_tetSubdivisionLevels; ++j)
			tc = _vbt->centroidUpOneLevel(tc);
		if (_megatetTetTris.find(tc) != _megatetTetTris.end())
			tetTriVec[i].tc[0] = USHRT_MAX;
	}
	_nSurfaceTets.store(_vbt->_tetNodes.size());  // this atomic must not step on any megatets that have already been created. Atomic used to multithread next section
	_newTets.clear();
	_ntsHash.clear();
#if defined( _DEBUG )
	for (int i = 0; i<tetTriVec.size(); ++i) {
		if (tetTriVec[i].tc[0] == USHRT_MAX)
			continue;
		getConnectedComponents(tetTriVec[i], _newTets, _ntsHash);  // for this centroid split its triangles into solid connected components
	}
#else
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, tetTriVec.size()),
		[&](tbb::blocked_range<size_t> r) {
			for (size_t i = r.begin(); i != r.end(); ++i) {
				if (tetTriVec[i].tc[0] == USHRT_MAX)
					continue;
				getConnectedComponents(tetTriVec[i], _newTets, _ntsHash);  // for this centroid split its triangles into solid connected components
			}
		});
#endif
	tetTriVec.clear();

	// vbt->_tetCentroids has the megatets remaining
	int incr = _nSurfaceTets - _vbt->_tetCentroids.size();
	_vbt->_tetCentroids.insert(_vbt->_tetCentroids.end(), incr, bccTetCentroid());
	_vbt->_tetNodes.insert(_vbt->_tetNodes.end(), incr, std::array<int, 4>());
	_surfaceTetTris.assign(incr, tetTris());
	for (auto& nt : _newTets) {
		_vbt->_tetCentroids[nt.tetIdx] = nt.tc;  // COURT don't hash them here
		_vbt->_tetNodes[nt.tetIdx] = std::move(nt.tetNodes);
		_vbt->_tetHash.insert(std::make_pair(nt.tc, nt.tetIdx));
		auto scit = _surfaceCentroids.find(nt.tc);
		assert(scit != _surfaceCentroids.end());
		scit->second.push_back(nt.tetIdx);
		_surfaceTetTris[nt.tetIdx - _vbt->_nMegatets].tetIdx = nt.tetIdx;  // if careful indexing don't need theis
		_surfaceTetTris[nt.tetIdx - _vbt->_nMegatets].tris = std::move(nt.tris);
	}
	_newTets.clear();
	// _centroidTriangles will be used later for vertex tetId and multires version, so don't delete
	struct extNodeLoc {
		std::array<short, 3> loc;
		std::list<nodeTetSegment> tetNodes;
	};
	std::vector<extNodeLoc> extNodeLocs;
	extNodeLocs.reserve(_ntsHash.size());
	for (auto& ns : _ntsHash) {
		extNodeLoc enl;
		enl.loc = ns.first;
		enl.tetNodes = std::move(ns.second);
		extNodeLocs.push_back(enl);
	}
	_ntsHash.clear();

	_firstNewExteriorNode = _vbt->_nodeGridLoci.size();
	oneapi::tbb::concurrent_vector<extNode> eNodes;
#if defined( _DEBUG )
	for (int i = 0; i < extNodeLocs.size(); ++i)
		assignExteriorTetNodes(extNodeLocs[i].loc, extNodeLocs[i].tetNodes, eNodes);
#else
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, extNodeLocs.size()),
		[&](tbb::blocked_range<size_t> r) {
			for (size_t i = r.begin(); i != r.end(); ++i)
				assignExteriorTetNodes(extNodeLocs[i].loc, extNodeLocs[i].tetNodes, eNodes);
		});
#endif
	for (auto& en : eNodes) {
		int eNode = _vbt->_nodeGridLoci.size();
		_vbt->_nodeGridLoci.push_back(std::move(en.loc));
		for (auto& ti : en.tiPairs)
			_vbt->_tetNodes[ti.first][ti.second] = eNode;
	}
	eNodes.clear();

	_vbt->_firstInteriorTet = _vbt->_tetNodes.size();
	fillInteriorMicroTets(_vnCentroids);
	// wed seams between macrotets and recut microtet regions with T junctions
	linkMicrotetsToMegatets();
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_vbt->_tetNodes.size());
	for (int n = _vbt->_tetNodes.size(), i = 0; i < n; ++i)  // firstInteriorTet
		_vbt->_tetHash.insert(std::make_pair(_vbt->_tetCentroids[i], i));
	// now reconnect stranded vertices to their new barycentric tet loci
	for (int n = _vbt->_vertexTets.size(), v = 0; v < n; ++v) {
		// _vertexTetCentroids[v] already converted to lowest microtet centroid values
		if (_vbt->_vertexTets[v] < -1)  // excised vertex
			continue;
		if (_vbt->_vertexTets[v] < 0) {  // only those not in a megatet
			if (!(_vertexTetCentroids[v][0] < USHRT_MAX))
				_vbt->gridLocusToLowestTetCentroid(_vMatCoords[v], _vertexTetCentroids[v]);
			auto vtet = _vbt->_tetHash.equal_range(_vertexTetCentroids[v]);
			int limit = 0;
			while (vtet.first == vtet.second) {
				_vertexTetCentroids[v] = _vbt->centroidUpOneLevel(_vertexTetCentroids[v]);
				vtet = _vbt->_tetHash.equal_range(_vertexTetCentroids[v]);
				++limit;
			}
			if (limit > 10)
				throw(std::logic_error("Stranded vertex can't find its enclosing tetrahedron."));
			if (std::distance(vtet.first, vtet.second) < 2)
				_vbt->_vertexTets[v] = vtet.first->second;
			else {
				auto ct = _surfaceCentroids.find(_vertexTetCentroids[v]);
				assert(ct != _surfaceCentroids.end());
				int i;
				for (auto& ti : ct->second) {
					auto& tt = _surfaceTetTris[ti - _vbt->_nMegatets];
					assert(tt.tetIdx == ti);
					for (auto& tri : tt.tris) {
						const int* tr = _mt->triangleVertices(tri);
						for (i = 0; i < 3; ++i) {
							if (tr[i] == v) {
								_vbt->_vertexTets[v] = tt.tetIdx;
								break;
							}
						}
						if (i < 3)
							break;
					}
					if (i < 3)
						break;
				}
				assert(i < 3);
			}
			// set barycentric coordinate within that tet
			auto& bw = _vbt->_barycentricWeights[v];
			_vbt->gridLocusToBarycentricWeight(_vMatCoords[v], _vertexTetCentroids[v], bw);
			// any vertex on a tet face boundary should be nudged inside
			bool nudge = false;
			for (int j = 0; j < 3; ++j) {
				if (bw[j] < 1e-4f || bw[j] > 0.9999f)
					nudge = true;
			}
			if (bw[0] + bw[1] + bw[2] > 0.9999)
				nudge = true;
			if (nudge) {  // vertices on tet boundarys must be pulled inside to avoid dual identities
				Vec3f nV = (Vec3f(0.25f, 0.25f, 0.25f) - bw) * 0.0002f;
				bw += nV;
				_vbt->barycentricWeightToGridLocus(_vertexTetCentroids[v], bw, _vMatCoords[v]);
			}
		}
	}
	// clean up
	// keep _vMatCoords and _vertexTetCentroids for augmentation with each new incision
	_interiorNodes.clear();
	if (_rtp != nullptr) {  // will need this for fast remapTetPhysics class
		_rtp->clearVnTetTris();
		for (auto& ct : _surfaceCentroids) {
			if (ct.second.size() < 2)
				continue;
			for (auto& ti : ct.second) {
				auto &tt = _surfaceTetTris[ti - _vbt->_nMegatets];
				assert(tt.tetIdx == ti);
				_rtp->insertVnTetTris(tt.tetIdx, tt.tris);
			}
		}
	}
	_surfaceCentroids.clear();
	_surfaceTetTris.clear();  // only reused in makeFirst(), but not again.
	_boundingNodeData.clear();
	_megatetBounds.clear();
	_lastTriangleSize = _mt->numberOfTriangles();
	_lastVertexSize = _mt->numberOfVertices();
	_vbt->_tetNodes.shrink_to_fit();
	_vbt->_tetCentroids.shrink_to_fit();
	_vbt->_nodeSpatialCoords = nullptr;  // this is owned by the physics system. Should be assigned after physics library creates it.
}

void vnBccTetCutter_tbb::createFirstMacroTets(materialTriangles* mt, vnBccTetrahedra* vbt, const int nLevels, const int maximumDimensionMacroSubdivs) {
	_mt = mt;
	_vbt = vbt;
	makeFirstVnTets(_mt, vbt, maximumDimensionMacroSubdivs);
	_vbt->_tetSubdivisionLevels = nLevels;  // Creating nLevels of multiresolution tets.
	int mult = (1 << (nLevels - 1)), shiftUp = nLevels - 1;
	// macrotets guaranteed not to virtual node.  Subcut any found at this stage.
	// collect all triangles in virtual noded tets and nonVN megatets
	_vnCentroids.clear();
	_vnTris.clear();
	_vnTris.reserve(_mt->numberOfTriangles() >> 2);  // COURT revisit this guess
	_megatetTetTris.clear();
	_megatetTetTris.reserve(_vbt->_tetNodes.size());
	for (auto& ct : _surfaceCentroids) {  // left in place after makeFirstVnTets()
		auto tc = ct.first;
		for (int i = 0; i < 3; ++i)
			tc[i] <<= shiftUp;
		if (ct.second.size() > 1) {
			_vnCentroids.push_back(tc);
			for (auto& ti : ct.second) {
				auto& tl = _surfaceTetTris[ti];
				assert(tl.tetIdx == ti);
				_vbt->_tetNodes[tl.tetIdx][0] = -1;  // mark for deletion
				_vnTris.insert(tl.tris.begin(), tl.tris.end());
			}
		}
		else {
			assert(ct.second.size() == 1);
			auto& tt = _surfaceTetTris[ct.second.front()];
			_vbt->_tetCentroids[tt.tetIdx] = tc;
			_megatetTetTris.insert(std::make_pair(tc, std::move(tt)));
		}
	}
	_surfaceCentroids.clear();
	_surfaceTetTris.clear();
	// now add interior megatets not penetrated by a triangle
	for (int n = _vbt->tetNumber(), i = _vbt->_firstInteriorTet; i < n; ++i) {
		tetTris tt;
		tt.tetIdx = i;
		auto& tc = vbt->_tetCentroids[i];
		for (int j = 0; j < 3; ++j)
			tc[j] <<= shiftUp;
		_megatetTetTris.insert(std::make_pair(tc, tt));
	}
	for (auto& mc : _vMatCoords)
		mc *= mult;
	vbt->_unitSpacingInv *= mult;
	vbt->_unitSpacing = 1.0 / vbt->_unitSpacingInv;
	for (int n = _vbt->_nodeGridLoci.size(), i = 0; i < n; ++i) {
		auto& nl = vbt->_nodeGridLoci[i];
		for (int j = 0; j < 3; ++j)
			nl[j] *= mult;  // can be negative
	}

	std::vector<int> vnTriVec(_vnTris.begin(), _vnTris.end());
	std::unordered_set<int> vtVerts;
	vtVerts.reserve(vnTriVec.size());
	for (auto t : vnTriVec) {
		int* tr = _mt->triangleVertices(t);
		for (int j = 0; j < 3; ++j)
			vtVerts.insert(tr[j]);
	}
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i)
		_vertexTetCentroids[i][0] = USHRT_MAX;  // reset all to empty
	for (auto v : vtVerts)
		_vbt->gridLocusToLowestTetCentroid(_vMatCoords[v], _vertexTetCentroids[v]);  // All tri vertices must be converted to lowest.  This invalidates next step.
	vtVerts.clear();
	// COURT 4x faster than single thread using tbb hash container requiring no reduction. tbb version ~30% faster than omp before reduction and reduction using critical section. tbb hash container very helpful.
#if defined( _DEBUG )
	for (int n = vnTriVec.size(), i = 0; i<n; ++i)
		inputTriangleTetsTbb(vnTriVec[i], _centTris);
#else
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, vnTriVec.size()),
		[&](tbb::blocked_range<size_t> r) {
			for (size_t i = r.begin(); i != r.end(); ++i) {
				inputTriangleTetsTbb(vnTriVec[i], _centTris);
			}
		});
#endif
	vnTriVec.clear();
	_interiorNodes.clear();  // COURT perhaps keep this and delete vn tet interiors
	// setup Z intersect arrays for finding interior nodes
	Vec3f maxMaterialCorner = (_vbt->_maxCorner - _vbt->_minCorner) * (float)_vbt->_unitSpacingInv;
	for (int i = 0; i < 3; ++i)
		_vbt->_gridSize[i] = 1 + (int)std::floor(maxMaterialCorner.xyz[i]);
	evenXy.clear();
	oddXy.clear();
	// setup lines parallel with Z axis
	evenXy.assign(_vbt->_gridSize[0] >> 1, std::vector<std::multimap<double, zIntersectFlags> >());  // 0th i always empty
	oddXy.assign(_vbt->_gridSize[0] >> 1, std::vector<std::multimap<double, zIntersectFlags> >());
	int gsy = _vbt->_gridSize[1] >> 1;
	for (int n = _vbt->_gridSize[0] >> 1, i = 0; i < n; ++i) {
		evenXy[i].assign(gsy, std::multimap<double, zIntersectFlags>());  // 0th j always empty
		oddXy[i].assign(gsy, std::multimap<double, zIntersectFlags>());
	}
	macrotetRecutCore();
}

bool vnBccTetCutter_tbb::makeFirstVnTets(materialTriangles* mt, vnBccTetrahedra* vbt, int maximumGridDimension)
{  // initial creation of vbt based only on materialTriangles input amd maxGridDim.
	if (maximumGridDimension > 0x8ffe)
		throw(std::logic_error("Maximum grid dimension requested must be less than 32K."));
	// WARNING - no complete tests are done to check for non-self-intersecting closed manifold triangulated surface input!!
	// This is essential. findAdjacentTriangles() is the closest test this routine provides.  Test externally.
	_mt = mt;
	_vbt = vbt;
	_vbt->_nMegatets = 0;
	_vbt->_mt = mt;
	_vbt->_barycentricWeights.clear();
	_vbt->_barycentricWeights.assign(mt->numberOfVertices(), Vec3f());
	_vbt->_tetHash.clear();
	_vbt->_tetNodes.clear();
	evenXy.clear();
	oddXy.clear();
	if (_mt->findAdjacentTriangles(true))	return false;
	if (!setupBccIntersectionStructures(maximumGridDimension))
		return false;
	_vbt->_tetSubdivisionLevels = 1;  // Not creating multiresolution tets.
	// COURT 4x faster than single thread using tbb hash container requiring no reduction. tbb version ~30% faster than omp before reduction and reduction using critical section. tbb hash container very helpful.
	auto procTri = [&](size_t i) {
		if (_mt->triangleMaterial(i) < 0)
			return;
		inputTriangleTetsTbb(i, _centTris);
		Vec3d triVec[3];
		int* tr = _mt->triangleVertices(i);
		for (int j = 0; j < 3; ++j)
			triVec[j] = _vMatCoords[tr[j]];
		zIntersectTriangleTbb(triVec, true, _zIntr);
		};
#ifdef _DEBUG
	for (size_t i = 0; i < _mt->numberOfTriangles(); ++i)
		procTri(i);
#else
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, _mt->numberOfTriangles()),
		[&](tbb::blocked_range<size_t> r) {
			for (size_t i = r.begin(); i != r.end(); ++i)
				procTri(i);
		});
#endif
	std::vector<tetTriangles> tetTriVec;
	tetTriVec.assign(_centTris.size(), tetTriangles());
	_surfaceCentroids.clear();
	_surfaceCentroids.reserve(_centTris.size());
	int count = 0;
	for (auto ctit = _centTris.begin(); ctit != _centTris.end(); ++ctit) {
		_surfaceCentroids.insert(std::make_pair(ctit->first, std::vector<int>()));
		tetTriVec[count].tc = ctit->first;
		tetTriVec[count++].tris = std::move(ctit->second);
	}
	_centTris.clear();
	for (auto ziv : _zIntr) {
		if (ziv.flags.odd)
			oddXy[ziv.x][ziv.y].insert(std::make_pair(ziv.zInt, ziv.flags));
		else
			evenXy[ziv.x][ziv.y].insert(std::make_pair(ziv.zInt, ziv.flags));
	}
	_zIntr.clear();
	// create and hash all interior nodes.  Very fast (< 0.002 sec) so don't bother multithreading
	createInteriorNodes();
	//	evenXy.clear(); oddXy.clear();  // am reusing these structures for remakeVnTets()
	_interiorNodes.clear();  // only nodes created thus far are interior nodes
	_interiorNodes.reserve(_vbt->_nodeGridLoci.size());
	for (int n = _vbt->_nodeGridLoci.size(), j = 0; j < n; ++j)
		_interiorNodes.insert(std::make_pair(_vbt->_nodeGridLoci[j], j));
	_surfaceTetTris.clear();
	_nSurfaceTets.store(0);  // this atomic must not step on any megatets that have already been created. In this routine there are none.  Atomic used to multithread next section

#ifdef _DEBUG
	for (int i = 0; i < tetTriVec.size(); ++i)
		getConnectedComponents(tetTriVec[i], _newTets, _ntsHash);  // for this centroid split its triangles into solid connected components
#else
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, tetTriVec.size()),
		[&](tbb::blocked_range<size_t> r) {
			for (size_t i = r.begin(); i != r.end(); ++i)
				getConnectedComponents(tetTriVec[i], _newTets, _ntsHash);  // for this centroid split its triangles into solid connected components
		});
#endif

	_vbt->_tetCentroids.assign(_nSurfaceTets, bccTetCentroid());
	_vbt->_tetNodes.assign(_nSurfaceTets, std::array<int, 4>());
	_surfaceTetTris.assign(_nSurfaceTets, tetTris());
	for (auto& nt : _newTets) {
		_vbt->_tetCentroids[nt.tetIdx] = nt.tc;  // COURT don't hash them here
		_vbt->_tetNodes[nt.tetIdx] = std::move(nt.tetNodes);
		auto scit = _surfaceCentroids.find(nt.tc);
		assert(scit != _surfaceCentroids.end());
		scit->second.push_back(nt.tetIdx);
		_surfaceTetTris[nt.tetIdx].tetIdx = nt.tetIdx;  // COURT if careful indexing don't need this
		_surfaceTetTris[nt.tetIdx].tris = std::move(nt.tris);
	}
	_newTets.clear();
	// _centroidTriangles will be used later for vertex tetId and multires version, so don't delete
	struct extNodeLoc {
		std::array<short, 3> loc;
		std::list<nodeTetSegment> tetNodes;
	};
	std::vector<extNodeLoc> extNodeLocs;
	extNodeLocs.reserve(_ntsHash.size());
	for (auto& ns : _ntsHash) {
		extNodeLoc enl;
		enl.loc = ns.first;
		enl.tetNodes = std::move(ns.second);
		extNodeLocs.push_back(enl);
	}
	_ntsHash.clear();
	// get tets where vertices reside
	_vbt->_vertexTets.clear();
	_vbt->_vertexTets.assign(_mt->numberOfVertices(), -1);
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i) {
		auto ct = _surfaceCentroids.find(_vertexTetCentroids[i]);
		assert(ct != _surfaceCentroids.end());
		if(ct->second.size() < 2)
			_vbt->_vertexTets[i] = _surfaceTetTris[ct->second.front()].tetIdx;
		else {
			for (auto& ti : ct->second) {
				auto& tt = _surfaceTetTris[ti];
				assert(tt.tetIdx == ti);
				for (auto& tri : tt.tris) {
					int* tr = _mt->triangleVertices(tri);
					int j;
					for (j = 0; j < 3; ++j) {
						if (tr[j] == i) {
							_vbt->_vertexTets[i] = tt.tetIdx;
							break;
						}
					}
					if (j < 3)
						break;
				}
				if (_vbt->_vertexTets[i] > -1)
					break;
			}
			assert(_vbt->_vertexTets[i] > -1);
		}
	}

	oneapi::tbb::concurrent_vector<extNode> eNodes;
#ifdef _DEBUG
	for (int i = 0; i < extNodeLocs.size(); ++i)
		assignExteriorTetNodes(extNodeLocs[i].loc, extNodeLocs[i].tetNodes, eNodes);
#else
	tbb::parallel_for(
		tbb::blocked_range<size_t>(0, extNodeLocs.size()),
		[&](tbb::blocked_range<size_t> r) {
			for (size_t i = r.begin(); i != r.end(); ++i)
				assignExteriorTetNodes(extNodeLocs[i].loc, extNodeLocs[i].tetNodes, eNodes);
		});
#endif

	for (auto& en : eNodes) {
		int eNode = _vbt->_nodeGridLoci.size();
		_vbt->_nodeGridLoci.push_back(std::move(en.loc));
		for (auto& ti : en.tiPairs)
			_vbt->_tetNodes[ti.first][ti.second] = eNode;
	}
	eNodes.clear();

	_vbt->_firstInteriorTet = _vbt->_tetNodes.size();
	fillNonVnTetCenter();  // fast. Don't bother multithreading
	_interiorNodes.clear();  // COURT perhaps keep this and delete vn tet interiors
	_vbt->_tetNodes.shrink_to_fit();
	_vbt->_tetCentroids.shrink_to_fit();
	_vbt->_tetHash.clear();
	_vbt->_tetHash.reserve(_vbt->_tetCentroids.size());
	for (int n = _vbt->_tetCentroids.size(), i = 0; i < n; ++i)
		_vbt->_tetHash.insert(std::make_pair(_vbt->_tetCentroids[i], i));
	return true;
}

void vnBccTetCutter_tbb::pack(){
	int tnNow = 0;
	std::vector<int> tnArr;
	tnArr.assign(_vbt->_tetNodes.size(), -1);
	for (int n = _vbt->_tetNodes.size(), i = 0; i < n; ++i) {
		if (_vbt->_tetNodes[i][0] > -1) {  // undeleted tet
			if (tnNow < i) {
				_vbt->_tetNodes[tnNow] = _vbt->_tetNodes[i];
				_vbt->_tetCentroids[tnNow] = _vbt->_tetCentroids[i];
			}
			tnArr[i] = tnNow++;
		}
	}
	_vbt->_tetNodes.erase(_vbt->_tetNodes.begin() + tnNow, _vbt->_tetNodes.end());
	_vbt->_tetCentroids.erase(_vbt->_tetCentroids.begin() + tnNow, _vbt->_tetCentroids.end());
	_vbt->_tetHash.clear();  // invalidate now, hash later
	// will repeatedly use _megatetTetTris, so remove any deleted ones
	for (auto& mt : _megatetTetTris) {
		mt.second.tetIdx = tnArr[mt.second.tetIdx];
		assert(mt.second.tetIdx > -1);
	}
	_vbt->_nMegatets = _megatetTetTris.size();
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i) {
		auto& vt = _vbt->_vertexTets[i];
		if (vt < -1)  // excised vertex
			continue;
		if (vt >= tnArr.size())
			vt = -1;
		else {
			vt = tnArr[vt];
			assert(vt < tnNow);
		}
	}
	tnArr.clear();
	tnArr.assign(_vbt->_nodeGridLoci.size(), -1);
	for (int n = _vbt->_tetNodes.size(), i = 0; i < n; ++i) {
		auto& tn = _vbt->_tetNodes[i];
		for (int j = 0; j < 4; ++j)
			tnArr[tn[j]] = 1;
	}
	tnNow = 0;
	for (int n = tnArr.size(), i = 0; i < n; ++i) {
		if (tnArr[i] > -1) {
			if (tnNow < i)  // ? _nodeSpatialCoords too
				_vbt->_nodeGridLoci[tnNow] = _vbt->_nodeGridLoci[i];
			tnArr[i] = tnNow++;
		}
	}
	_vbt->_nodeGridLoci.erase(_vbt->_nodeGridLoci.begin() + tnNow, _vbt->_nodeGridLoci.end());
	for (int n = _vbt->_tetNodes.size(), i = 0; i < n; ++i) {
		auto& tn = _vbt->_tetNodes[i];
		for (int j = 0; j < 4; ++j)
			tn[j] = tnArr[tn[j]];
	}
}

void vnBccTetCutter_tbb::linkMicrotetsToMegatets(){
	// pack decimated tets and nodes
	std::vector<int> nodeMap;
	nodeMap.assign(_vbt->_nodeGridLoci.size(), -1);
	int offset = _vbt->_nMegatets;
	for (int n = _vbt->_tetNodes.size(), i = _vbt->_nMegatets; i < n; ++i) {
		auto& tn = _vbt->_tetNodes[i];
		if (tn[0] > -1) {
			for (int j = 0; j < 4; ++j)
				nodeMap[tn[j]] = 1;
			if (offset < i) {
				_vbt->_tetNodes[offset] = _vbt->_tetNodes[i];
				_vbt->_tetCentroids[offset] = _vbt->_tetCentroids[i];
			}
			++offset;
		}
	}
	// Don't assign _vbt->_vertexTets here.  Leave -1 and assign later
	_vbt->_tetNodes.resize(offset);
	_vbt->_tetCentroids.resize(offset);
	for (int i = 0; i < _meganodeSize; ++i)  // don't want any original meganodes overwritten
		nodeMap[i] = i;
	offset = _meganodeSize;
	for (int n = nodeMap.size(), i = _meganodeSize; i < n; ++i) {
		if (nodeMap[i] > -1) {
			nodeMap[i] = offset;
			if (offset < i)
				_vbt->_nodeGridLoci[offset] = _vbt->_nodeGridLoci[i];
			++offset;
		}
	}
	_firstNewExteriorNode = nodeMap[_firstNewExteriorNode];
	if (_firstNewExteriorNode < 0)
		throw(std::logic_error("Program error in decimateInteriorMicroTets()\n"));
	_vbt->_nodeGridLoci.resize(offset);
	for (int n = _vbt->_tetNodes.size(), i = _vbt->_nMegatets; i < n; ++i) {
		auto& tn = _vbt->_tetNodes[i];
		for (int j = 0; j < 4; ++j) {
			assert(nodeMap[tn[j]] > -1);
			tn[j] = nodeMap[tn[j]];
		}
	}
	_vbt->_tJunctionConstraints.clear();
	struct boundingTri {
		Vec3f N, tv[3];
		float d;
		Mat2x2f M;
		std::vector<int>* tris;
	};
	std::vector < boundingTri> bt;
	//	bt.reserve(boundingTris.size());
	bt.reserve(_megatetBounds.size());
	for (int n = _megatetBounds.size(), i = 0; i < n; ++i) {
		boundingTri t;
		t.tris = _megatetBounds[i].tris;
		for (int j = 0; j < 3; ++j)
			t.tv[j] = (short(&)[3]) * _vbt->_nodeGridLoci[_megatetBounds[i].nodes[j]].data();
		t.tv[1] -= t.tv[0];
		t.tv[2] -= t.tv[0];
		t.N = t.tv[1] ^ t.tv[2];
		t.d = t.N * t.tv[0];
		t.M = { t.tv[1] * t.tv[1], t.tv[1] * t.tv[2], 0.0f, t.tv[2] * t.tv[2] };
		t.M.x[2] = t.M.x[1];
		bt.push_back(std::move(t));
	}
	auto addTJunction = [&](const int node, int& tri, Vec2f& bary) {
		auto dnc = _vbt->_tJunctionConstraints.insert(std::make_pair(node, vnBccTetrahedra::decimatedFaceNode()));
		if (dnc.second) {
			if (bary[1] < 1e-5f) {
				dnc.first->second.faceNodes.assign(_megatetBounds[tri].nodes.begin(), _megatetBounds[tri].nodes.begin() + 2);
				dnc.first->second.faceBarys.reserve(2);
				dnc.first->second.faceBarys.push_back(1.0f - bary[0]);
				dnc.first->second.faceBarys.push_back(bary[0]);
			}
			else if (bary[0] < 1e-5f) {
				dnc.first->second.faceNodes.reserve(2);
				dnc.first->second.faceNodes.push_back(_megatetBounds[tri].nodes.front());
				dnc.first->second.faceNodes.push_back(_megatetBounds[tri].nodes.back());
				dnc.first->second.faceBarys.reserve(2);
				dnc.first->second.faceBarys.push_back(1.0f - bary[1]);
				dnc.first->second.faceBarys.push_back(bary[1]);
			}
			else if (bary[0] + bary[1] > 0.99995f) {
				dnc.first->second.faceNodes.assign(_megatetBounds[tri].nodes.begin() + 1, _megatetBounds[tri].nodes.end());
				dnc.first->second.faceBarys.reserve(2);
				dnc.first->second.faceBarys.push_back(bary[0]);
				dnc.first->second.faceBarys.push_back(bary[1]);
			}
			else {
				dnc.first->second.faceNodes.assign(_megatetBounds[tri].nodes.begin(), _megatetBounds[tri].nodes.end());
				dnc.first->second.faceBarys.reserve(3);
				dnc.first->second.faceBarys.push_back(1.0f - bary[0] - bary[1]);
				dnc.first->second.faceBarys.push_back(bary[0]);
				dnc.first->second.faceBarys.push_back(bary[1]);
			}
		}
	};
	auto isTjunct = [&](const int node, int& tri, Vec2f& bary, std::vector<int>* triTest) ->bool {
		Vec3f P = (short(&)[3]) * _vbt->_nodeGridLoci[node].data();
		for (tri = 0; tri < bt.size(); ++tri) {
			auto& b = bt[tri];
			if (fabs(P * b.N - b.d) < 1e-8f) {
				Vec3f V = P - b.tv[0];
				bary = b.M.Robust_Solve_Linear_System(Vec2f(V * b.tv[1], V * b.tv[2]));
				if (bary[0] < 0.0f || bary[1] < 0.0f || bary[0] > 1.0f || bary[1] > 1.0f || bary[0] + bary[1] > 1.0f)
					continue;
				if (triTest == nullptr)
					return true;
				else {
					auto bit = b.tris->begin();
					auto tit = triTest->begin();
					while (bit != b.tris->end()) {
						if (*bit < *tit)
							++bit;
						else if (*tit < *bit) {
							while ((*tit < *bit)) {
								++tit;
								if (tit == triTest->end())
									return false;
							}
						}
						else {
							assert(*tit == *bit);
							return true;
						}
					}
					return false;
				}
			}
		}
		return false;
	};
	for (int i = _meganodeSize; i < _firstNewExteriorNode; ++i) {  //       can't use exterior nodes since a virtual node could link across to wrong side?  COURT check this
		int tri;
		Vec2f bary;
		if (isTjunct(i, tri, bary, nullptr)) {
			addTJunction(i, tri, bary);
		}
	}

	std::vector<char> ctOpen;
	ctOpen.assign(_vbt->_nodeGridLoci.size() - _firstNewExteriorNode, 0xff);
	for (auto& tt : _surfaceTetTris) {
		auto &tn = _vbt->_tetNodes[tt.tetIdx];
		for (int j = 0; j < 4; ++j) {
			if (tn[j] >= _firstNewExteriorNode && ctOpen[tn[j] - _firstNewExteriorNode]) {
				ctOpen[tn[j] - _firstNewExteriorNode] = 0x00;
				int tri;
				Vec2f bary;
				if (isTjunct(tn[j], tri, bary, &tt.tris)) {
					addTJunction(tn[j], tri, bary);
				}
			}
		}
	}
	// don't rehash _vbt->_tetHash here
}

void vnBccTetCutter_tbb::createInteriorNodes() {
	short z;
	assert(_vbt->_tetSubdivisionLevels < 2);  // this routine only for first pass
	std::array<short, 3> s3;
	double zTol = _vbt->_gridSize[2] * 1e-8;
	auto solidFilter = [&](std::multimap<double, zIntersectFlags>& mm) ->bool {  // isolate solid runs from surface edge and vertex hits and correct switch of coincident surfaces.
		// first correct any coincident or multihit surfaces from triangle surface
		bool noSubdivSolid = true;
		auto mit = mm.begin();
		while (mit != mm.end()) {
			if (mit->second.surfaceTri) {
				noSubdivSolid = false;
				++mit;
				continue;
			}
			auto mNext = mit;
			++mNext;
			while (mNext != mm.end()) {
				if (mNext->second.surfaceTri) {
					noSubdivSolid = false;
					++mNext;
					continue;
				}
				while (mNext != mm.end() && mNext->first - mit->first < zTol) {  // multihit
					if (mNext->second.surfaceTri) {
						noSubdivSolid = false;
						++mNext;
						continue;
					}
					mNext = mm.erase(mNext);
				}
				while (mNext != mm.end() && mNext->second.surfaceTri) {
					noSubdivSolid = false;
					++mNext;
				}
				if (mNext == mm.end())
					break;
				if (mit->second.solidBegin) {
					if (mNext->second.solidBegin)
						mNext = mm.erase(mNext);
					else
						break;
				}
				else {
					if (!mNext->second.solidBegin)
						mit = mm.erase(mit);
					else
						break;
				}
			}
			mit = mNext;
		}
		mit = mm.begin();
		while (mit != mm.end()) {
			assert(mit->second.surfaceTri);
			auto mNext = mit;
			++mNext;
			while (mNext != mm.end() && mNext->first - mit->first < zTol) {  // possible multihit or coincident surfaces
				assert(mNext->second.surfaceTri);
				while (mNext != mm.end() && mNext->first - mit->first < zTol) {  // possible multihit or coincident surfaces
					if (!mNext->second.surfaceTri) {
						if (mit->second.solidBegin != mNext->second.solidBegin)
							throw(std::logic_error("Solid ordering error in createInteriorNodes()\n"));
						else if (!mit->second.solidBegin) {  // make sure triangle solids are inside tet solids
							auto tmp = mit->second;
							mit->second = mNext->second;
							mNext->second = tmp;
						}
						else;
						++mNext;
					}
					else if (mNext->second.solidBegin == mit->second.solidBegin) // multihit
						mNext = mm.erase(mNext);
					else {  // coincident surface
						if (mit->second.solidBegin) {  // incorrect ordering for a coincident surface cutting a solid block
							auto tmp = mit->second;
							mit->second = mNext->second;
							mNext->second = tmp;
						}
						++mNext;
						assert(mNext->second.surfaceTri);
					}
				}
			}
			mit = mNext;
		}
		if (mm.size() < 2) {
			mm.clear();
			return true;
		}
		return false;
	};
	auto runInteriorNodes = [&](std::multimap<double, zIntersectFlags>& mm, bool evenLine) {
		auto mit = mm.begin();
		auto mend = mm.end();
		if (mit == mend)
			return;
		if (!mit->second.solidBegin)
			throw(std::runtime_error("Model is not a closed manifold surface.\n"));
		while (mit != mend) {
			// create runs of interior vertices
			z = (short)std::floor(mit->first);
			if (z < mit->first)
				++z;
			if ((z & 1) == evenLine)
				++z;
			++mit;
			if (mit->second.solidBegin) {  // correct any coincident pairs or micro collisions
				auto mit2 = mit;
				while (mit2 != mend && mit2->second.solidBegin && mit2->first - mit->first < 2e-5f)
					++mit2;
				if (mit2 == mend || mit2->second.solidBegin)
					throw(std::runtime_error("Model has a self collision.\n"));
				mit = mit2;
			}
			while ( z < mit->first) {
				s3[2] = z;
				_vbt->_nodeGridLoci.push_back(s3);
				z += 2;
			}
			++mit;
			if (mit != mend && !mit->second.solidBegin) {  // correct any coincident pairs or micro collisions
				auto mit2 = mit;
				while (mit2 != mend && !mit2->second.solidBegin && mit2->first - mit->first < 2e-5f)
					++mit2;
				if (mit2 == mend)
					return;
				if (!mit2->second.solidBegin)
					throw(std::runtime_error("Model has a self collision.\n"));
				mit = mit2;
			}
		}
	};
	for (short xi = 0; xi < evenXy.size(); ++xi) {
		for (short yi = 0; yi < evenXy[xi].size(); ++yi) {
			if (evenXy[xi][yi].empty())
				continue;
			s3[0] = (xi + 1) * 2;
			s3[1] = (yi + 1) * 2;
			if (solidFilter(evenXy[xi][yi]))
				continue;
			runInteriorNodes(evenXy[xi][yi], true);
		}
	}
	for (short xi = 0; xi < oddXy.size(); ++xi) {
		for (short yi = 0; yi < oddXy[xi].size(); ++yi) {
			if (oddXy[xi][yi].empty())
				continue;
			s3[0] = xi * 2 + 1;
			s3[1] = yi * 2 + 1;
			if (solidFilter(oddXy[xi][yi]))
				continue;;
			runInteriorNodes(oddXy[xi][yi], false);
		}
	}
}

void vnBccTetCutter_tbb::createInteriorMicronodes() {
	short z;
	std::array<short, 3> s3;
	double zTol = _vbt->_gridSize[2] * 1.5e-4;
	auto solidFilter = [&](std::multimap<double, zIntersectFlags>& mm) ->bool {  // isolate solid runs from surface hits and correct any switch of coincident surfaces.
		bool inSolid = false, noInteriorNodes = true, swapHere;
		auto mit = mm.begin();
		while (mit != mm.end()) {
			if (mit->second.surfaceTri) {
				if (inSolid != mit->second.solidBegin) {
					inSolid = !inSolid;
					swapHere = false;
				}
				else {
					assert(mit != mm.begin());  // possible coincident surface hit.
					swapHere = true;
				}
			}
			else {
				if (!inSolid)
					mit = mm.erase(mit);
				else {
					noInteriorNodes = false;
					++mit;
				}
				continue;
			}
			auto mNext = mit;
			++mNext;
			if (swapHere && (mNext == mm.end() || mNext->first - mit->first > zTol))
				throw(std::logic_error("Program logic error in vnBccTetCutter_tbb::createInteriorMicronodes()\n"));
			while (swapHere && mNext != mm.end() && mNext->first - mit->first < zTol) {  // remove all interior nodes and any multihits at surface intersect
				// interior nodes should already have their multihits eliminated
				if (!mit->second.surfaceTri) {  // interior nodes can't be on the surface
					mm.erase(mit);
					mit = mNext;
					++mNext;
					continue;
				}
				else if (!mNext->second.surfaceTri) {
					mNext = mm.erase(mNext);
					continue;
				}
				else {  // coincident surface or multihit
					if (mit->second.solidBegin != mNext->second.solidBegin) {
						mit->second.solidBegin = !mit->second.solidBegin;
						mNext->second.solidBegin = !mNext->second.solidBegin;
						inSolid = mNext->second.solidBegin;
						swapHere = false;
						++mNext;
					}
					else  // delete any multihits
						mNext = mm.erase(mNext);
				}
			}
			if (swapHere) {  // test for valid exit
				mNext = mit;
				++mNext;
				if (mNext == mm.end() || mNext->first - mit->first > zTol)
					throw(std::logic_error("Program logic error in vnBccTetCutter_tbb::createInteriorMicronodes()\n"));
				++mNext;
				swapHere = false;
			}
			mit = mNext;
		}
		if(inSolid)
			throw(std::logic_error("Program logic error in vnBccTetCutter_tbb::createInteriorMicronodes()\n"));
		if (noInteriorNodes)
			return true;
		if (mm.size() < 3) {
			mm.clear();
			return true;
		}
		return false;
	};
	auto runInteriorNodes = [&](std::multimap<double, zIntersectFlags>& mm, bool evenLine) {
		bool inSolid = false;
		auto mit = mm.begin();
		while (mit != mm.end()) {
			if (mit->second.surfaceTri) {
				if (mit->second.solidBegin != inSolid)
					inSolid = !inSolid;
				else
					throw(std::logic_error("Solid bounding error in createInteriorMicronodes()\n"));
			}
			else {
				assert(inSolid);
				s3[2] = (short)mit->first;
				if (mit->second.macroNode) {
					auto iit = _boundingNodeData.find(s3);
					if (iit != _boundingNodeData.end()) {
						_interiorNodes.insert(std::make_pair(s3, iit->second.node));
						++mit;
						continue;
					}
				}
				_interiorNodes.insert(std::make_pair(s3, _vbt->_nodeGridLoci.size()));
				_vbt->_nodeGridLoci.push_back(s3);
			}
			++mit;
		}
	};
	for (short xi = 0; xi < evenXy.size(); ++xi) {
		for (short yi = 0; yi < evenXy[xi].size(); ++yi) {
			if (evenXy[xi][yi].empty())
				continue;
			s3[0] = (xi + 1) * 2;
			s3[1] = (yi + 1) * 2;
			if (solidFilter(evenXy[xi][yi]))
				continue;
			runInteriorNodes(evenXy[xi][yi], true);
		}
	}
	for (short xi = 0; xi < oddXy.size(); ++xi) {
		for (short yi = 0; yi < oddXy[xi].size(); ++yi) {
			if (oddXy[xi][yi].empty())
				continue;
			s3[0] = xi * 2 + 1;
			s3[1] = yi * 2 + 1;
			if (solidFilter(oddXy[xi][yi]))
				continue;;
			runInteriorNodes(oddXy[xi][yi], false);
		}
	}
}

void vnBccTetCutter_tbb::assignExteriorTetNodes(std::array<short, 3>& locus, std::list<nodeTetSegment>& tetNodeIds, oneapi::tbb::concurrent_vector<extNode>& eNodes) {
	std::list<std::list<nodeTetSegment* > > tetPools;
	auto tnit = tetNodeIds.begin();
	while (tnit != tetNodeIds.end()) {
		std::list<std::list<nodeTetSegment* > * > poolLinks;
		for (auto &tp : tetPools) {
			for (auto &pe : tp) {
				assert(tnit->tetIdx != pe->tetIdx);  // obvious nuke me
				if (sortedVectorsIntersect(_surfaceTetTris[tnit->tetIdx - _vbt->_nMegatets].tris, _surfaceTetTris[pe->tetIdx - _vbt->_nMegatets].tris)) {
					poolLinks.push_back(&tp);
					break;
				}
			}
		}
		if (poolLinks.size() < 1) {
			tetPools.push_back(std::list<nodeTetSegment* >());
			tetPools.back().push_back(&(*tnit));
		}
		else if (poolLinks.size() < 2) {
			poolLinks.front()->push_back(&(*tnit));
		}
		else {
			auto pit = poolLinks.begin();
			(*pit)->push_back(&(*tnit));
			++pit;
			while (pit != poolLinks.end()) {
				poolLinks.front()->splice(poolLinks.front()->end(), **pit);
				for (auto tpit = tetPools.begin(); tpit != tetPools.end(); ++tpit) {
					if (&(*tpit) == *pit) {
						tetPools.erase(tpit);
						break;
					}
				}
				++pit;
			}
		}
		++tnit;
	}
	// each tetPool should get its own exterior, possibly virtual, node
	for (auto &tp : tetPools) {
		extNode en;
		en.node = -1;
		en.loc = locus;
		if (!_boundingNodeData.empty()) {
			auto bnit = _boundingNodeData.find(locus);
			if (bnit != _boundingNodeData.end()) {
				for (auto& bt : bnit->second.megaTetTris) {
					auto tpit = tp.begin();
					while (tpit != tp.end()) {
						if (sortedVectorsIntersect(bt->tris, (*tpit)->tetNodeTris)) {
							en.node = bnit->second.node;
							break;
						}
						++tpit;
					}
					if (tpit != tp.end())
						break;
				}
			}
		}
		for (auto& ntsp : tp)
			en.tiPairs.push_back(std::make_pair(ntsp->tetIdx, ntsp->tetNodeIndex));
		eNodes.push_back(en);
	}
}

int vnBccTetCutter_tbb::nearestRayPatchHit(const Vec3d &rayBegin, Vec3d rayEnd, const std::vector<int>& tris, Vec3d &hitP, double& distanceSq) {  // Return -1 is inside hit, 1 is outside hit and 0 is no hit.
	distanceSq = DBL_MAX;
	Vec3d rE(rayEnd - rayBegin);
	Vec3d N;
	for (auto& t : tris) {
		int* tr = _mt->triangleVertices(t);
		Vec3d tv[3];
		for (int i = 0; i < 3; ++i)
			tv[i].set(_vMatCoords[tr[i]]);
		tv[1] -= tv[0];
		tv[2] -= tv[0];
		double dSq;
		Vec3d P, R;
		Mat3x3d M;
		M.Initialize_With_Column_Vectors(tv[1], tv[2], -rE);
		R = M.Robust_Solve_Linear_System(Vec3d(rayBegin) - tv[0]);
		if (R[0] < -1e-8 || R[0] > 1.00000001 || R[1] < -1e-8 || R[1] > 1.00000001 || R[2] < 0.0 || R[2] > 1.0001 || R[0] + R[1] > 1.00000001)
			continue;
		P = rE * R[2];
		dSq = P * P;
		if (distanceSq > dSq) {
			distanceSq = dSq;
			hitP = rayBegin + rE * R[2];
			N = tv[1] ^ tv[2];
		}
	}
	if (distanceSq < DBL_MAX) {
		if(distanceSq < 1e-3)  // coincident surfaces don't connect
			return 1;
		if (N * rE > 0.0f)
			return 1;
		else
			return -1;
	}
	else
		return 0;
}

void vnBccTetCutter_tbb::processHoles(std::list<patch>& patches, std::list<hole>& holes, const Vec3d(&tetNodes)[4], int(&inodes)[4]) {
	std::list<patch*> outerP;
	for (auto& p : patches) {
		auto hit = holes.begin();
		for (; hit != holes.end(); ++hit)
			if (hit->pptr == &p)
				break;
		if (hit == holes.end()) {
			for (int i = 0; i < 4; ++i) {
				if (p.tetNodes[i] > -1)
					inodes[i] = -1;
			}
			outerP.push_back(&p);
		}
	}
	auto surrounds = [&](const patch* outer, const hole* h, const Vec3d &lineN, double& lineT) ->bool {
		lineT = DBL_MAX;
		Vec3d minTriN;
		bool pathEnclosed = false;
		if (outer->tetNodes[h->penetratedFace] > -1) {
			lineT = 1.0;
			pathEnclosed = true;
		}
		for (auto& t : outer->tris) {
			int* tn = _mt->triangleVertices(t);
			Vec3d V[3];
			for (int i = 0; i < 3; ++i)
				V[i].set((const float(&)[3])_vMatCoords[tn[i]]);
			V[1] -= V[0];
			V[2] -= V[0];
			Mat3x3d M(V[1], V[2], -lineN);
			Vec3d R = M.Robust_Solve_Linear_System(h->outerPoint - V[0]);
			if (R[2] <= 0.0 || R[2] >= 1.0 || R[0] < 0.0 || R[0] >= 1.0 || R[1] < 0.0 || R[1] >= 1.0 || R[0] + R[1] >= 1.0)  // coincident surfaces managed in isHole()
				continue;
			if (lineT > R[2]) {
				lineT = R[2];
				minTriN = V[1] ^ V[2];
				pathEnclosed = (minTriN * lineN > 0.0);
			}
		}
		return pathEnclosed;
	};
	auto encloseHole = [&](patch *enclosingPatch, hole *hole) {
		std::vector<int> mergeVec(enclosingPatch->tris.size() + hole->pptr->tris.size());
		std::merge(enclosingPatch->tris.begin(), enclosingPatch->tris.end(), hole->pptr->tris.begin(), hole->pptr->tris.end(), mergeVec.begin());
		enclosingPatch->tris = std::move(mergeVec);
		hole->pptr->tris.clear();
		auto pit = patches.begin();
		while (pit != patches.end()) {
			if (&(*pit) == hole->pptr) {
				patches.erase(pit);
				break;
			}
			++pit;
		}
	};
	auto hit = holes.begin();
	while (hit != holes.end()) {
		double minT = DBL_MAX;
		patch* enclosingP = nullptr;
		Vec3d toCorner = tetNodes[hit->penetratedFace] - hit->outerPoint;
		for (auto& op : outerP) {
			double t;
			if (surrounds(op, &(*hit), toCorner, t)) {
				if (t < minT) {
					minT = t;
					enclosingP = op;
				}
			}
		}
		if (enclosingP) {
			encloseHole(enclosingP, &(*hit));
			hit = holes.erase(hit);
		}
		else
			++hit;
	}
	if (!holes.empty()) {
		// any holes within holes would have been enclosed in an outerP, so these are independent.
		// Any outerP must be inside a hole
		if (inodes[0] < 0) {  // outerP must be enclosing a hole with a coincident surface
			assert(inodes[1] < 0 && inodes[2] < 0 && inodes[3] < 0 && outerP.size() == 1);
		}
		else {
			assert(inodes[1] > -1 && inodes[2] > -1 && inodes[3] > -1);
			auto hit = holes.begin();
			for (int i = 0; i < 4; ++i)
				hit->pptr->tetNodes[i] = inodes[i];
			++hit;
			while (hit != holes.end()) {
				encloseHole(holes.front().pptr, &(*hit));
				hit = holes.erase(hit);
			}
		}
	}
}

bool vnBccTetCutter_tbb::isHole(const patch* patch, const bccTetCentroid& tc, const Vec3d(&tetNodes)[4], Vec3d& patchPoint, int& tetFaceHit) {
	Vec3d tetFaceN[4];
	double tetD[4];
	for (tetFaceHit = 0; tetFaceHit < 4; ++tetFaceHit) {  // will also be intersected face #
		Vec3d F[3];
		F[0] = tetNodes[tetFaceHit];
		F[1] = tetNodes[(tetFaceHit + 1) & 3];
		F[2] = tetNodes[(tetFaceHit + 2) & 3];
		F[1] -= F[0];
		F[2] -= F[0];
		if (tetFaceHit & 1)  // outward facing tet normal. Clockwiseness alternates.
			tetFaceN[tetFaceHit] = F[1] ^ F[2];
		else
			tetFaceN[tetFaceHit] = F[2] ^ F[1];
		tetD[tetFaceHit] = tetFaceN[tetFaceHit] * F[0];
	}
	int startTri;
	for (startTri = 0; startTri < patch->tris.size(); ++startTri) {
		int j, *tn = _mt->triangleVertices(patch->tris[startTri]);
		for (tetFaceHit = 0; tetFaceHit < 4; ++tetFaceHit) {
			for (j = 0; j < 3; ++j) {
				if (tetFaceN[tetFaceHit] * Vec3d(_vMatCoords[tn[j]]) - tetD[tetFaceHit] > 0.0)
					break;
			}
			if (j < 3)
				break;
		}
		if (tetFaceHit < 4)
			break;
	}
	if (tetFaceHit > 3)
		throw(std::logic_error("Program error in isHole().\n"));
	// Can have nested closed polygons (e.g. tip of pipette), so must search all tris for outermost.
	double minDsq = DBL_MAX;
	Vec3d minTriN;
	for (int n = patch->tris.size(), i = startTri; i < n; ++i) {
		double td[3];
		Vec3d tV[3];
		int* tn = _mt->triangleVertices(patch->tris[i]);
		bool allInside = true;
		for (int j = 0; j < 3; ++j) {
			tV[j].set(_vMatCoords[tn[j]]);
			td[j] = tetFaceN[tetFaceHit] * tV[j] - tetD[tetFaceHit];
			if (td[j] > 0.0)
				allInside = false;
		}
		if (allInside)
			continue;
		Vec3d I[2];
		int v2 = 0, count = 0;
		for (int j = 2; j > -1; --j) {
			if ((td[v2] > 0.0) != (td[j] > 0.0)) {
				double s = td[v2] / (td[v2] - td[j]);
				I[count++] = tV[j] * s + tV[v2] * (1.0 - s);
				if (count > 1)
					break;
			}
			v2 = j;
		}
		I[1] -= I[0];
		double den = I[1] * I[1];
		if (den > 1e-20) {
			double s = ((tetNodes[tetFaceHit] - I[0]) * I[1]) / den;
			if (s < 1e-5)
				s = 1e-5;
			else if (s > 0.99999)
				s = 0.99999;
			else;
			I[0] += I[1] * s;
			double dsq = (tetNodes[tetFaceHit] - I[0]).length2();
			if (minDsq > dsq) {
				minDsq = dsq;
				patchPoint = I[0];
				minTriN = (tV[1] - tV[0]) ^ (tV[2] - tV[0]);
			}
		}
	}
	Vec3d toCorner = tetNodes[tetFaceHit] - patchPoint;
	if (toCorner * minTriN > 0.0)
		return false;
	// pull facePoint slightly inside solid to eliminate coincident surface problem between patches.
	toCorner.q_normalize();
	patchPoint += toCorner * 1e-5;
	return true;
}


void vnBccTetCutter_tbb::getConnectedComponents(const tetTriangles& tt, oneapi::tbb::concurrent_vector<newTet>& nt_vec, NTS_HASH& local_nts) {
	// for this centroid split its triangles into single solid connected components
	std::set<int> trSet, ts;
	trSet.insert(tt.tris.begin(), tt.tris.end());
	patch p;
	std::list<patch> patches;
	typedef std::list<patch>::iterator PListIt;
	p.tetNodes = { -1, -1, -1, -1 };
	p.tetIndex = -1;
	patches.push_back(p);
	auto trVec = &patches.back().tris;  // on entry only one triangle list for this centroid
	while (!trSet.empty()) { // find list of connected triangles
		ts.clear();
		int t, at[3], ae[3];
		t = *trSet.begin();
		trSet.erase(trSet.begin());
		ts.insert(t);
		std::forward_list<int> tList;  // unprocessed adjacencies
		tList.push_front(t);
		while (!tList.empty()) {
			t = tList.front();
			tList.pop_front();
			_mt->triangleAdjacencies(t, at, ae);  // COURT consider a local implementation of this
			for (int i = 0; i < 3; ++i) {
				auto ti = trSet.find(at[i]);
				if (ti == trSet.end())
					continue;
				trSet.erase(ti);
				ts.insert(at[i]);
				tList.push_front(at[i]);
			}
		}
		trVec->assign(ts.begin(), ts.end());  // this leaves all patch triangles in a sorted vector
		if (!trSet.empty()) {
			patches.push_back(p);
			trVec = &patches.back().tris;
		}
	}
	// find interior nodes for each component, joining patches as appropriate
	short gl[4][3];
	_vbt->centroidToNodeLoci(tt.tc, gl);
	int inodes[4] = { -1, -1, -1, -1 };
	for (int i = 0; i < 4; ++i) {
		std::array<short, 3> loc = { gl[i][0], gl[i][1], gl[i][2] };
		auto iit = _interiorNodes.find(loc);
		if (iit != _interiorNodes.end())
			inodes[i] = iit->second;
	}
	if (patches.size() < 2) {  // no interpatch computation required
		for (int i = 0; i < 4; ++i)
			patches.front().tetNodes[i] = inodes[i];
	}
	else {  // resolve which patches should be combined and which left separate, based on enclosure of solid section
		Vec3d tetNodes[4];
		for (int i = 0; i < 4; ++i)
			tetNodes[i].set(gl[i]);
		struct edgeIntrsct {
			bool nextOutside;
			patch* pptr;
		};
		std::list<hole> tetFaceHoles;
		std::multimap<double, edgeIntrsct> edges[6];
		auto removePatch = [&](patch* p) {
			for (auto pi = patches.begin(); pi != patches.end(); ++pi)
				if (&(*pi) == p) {
					patches.erase(pi);
					break;
				}
		};
		bool edgeCut = false;  // interior face intersects or edge cut patches present
		for (auto pit = patches.begin(); pit != patches.end(); ++pit) {
			pit->noEdgeCut = true;
			for (auto& t : pit->tris) {
				int* tr = _mt->triangleVertices(t);
				Vec3d T[3], P, Q;
				for (int i = 0; i < 3; ++i)
					T[i] = _vMatCoords[tr[i]];
				T[1] -= T[0];
				T[2] -= T[0];
				Mat3x3d M(T[1], T[2], Vec3d());
				auto getEdgeIntersect = [&](int edge) {
					for (int i = 0; i < 3; ++i)
						M.x[i + 6] = P[i] - Q[i];
					Vec3d R = M.Robust_Solve_Linear_System(P - T[0]);
					if (R[0] < -1e-8 || R[0] > 1.00000001 || R[1] < -1e-8 || R[1] > 1.00000001 || R[2] <= 0.0 || R[2] >= 1.0 || R[0] + R[1] >= 1.00000001)  // don't want triangle vertex hits, only edges
						return;
					edgeIntrsct ei;
					ei.pptr = &(*pit);
					ei.nextOutside = (T[1] ^ T[2]) * (Q - P) >= 0;
					edges[edge].insert(std::make_pair(R[2], ei));
					pit->noEdgeCut = false;
				};
				P = tetNodes[3];
				for (int i = 0; i < 4; ++i) {
					Q = tetNodes[i];
					getEdgeIntersect(i);
					P = Q;
				}
				P = tetNodes[0];
				Q = tetNodes[2];
				getEdgeIntersect(4);
				P = tetNodes[1];
				Q = tetNodes[3];
				getEdgeIntersect(5);
			}
			if (pit->noEdgeCut) {  // This patch is interior to all edges
				// if facing exterior to any tet node, must be independent of all patches with an edge intersection processed above and must be exterior to all tet nodes.
				// if interior could be interior to any existing solid
				hole h;
				if (isHole(&(*pit), tt.tc, tetNodes, h.outerPoint, h.penetratedFace)) {  // may be a hole inside a solid patch, otherwise a spike will always bound its own tet solid.
					h.pptr = &(*pit);
					tetFaceHoles.push_back(h);
				}
			}
			else
				edgeCut = true;
		}
		std::list<std::set<patch*> > combineP;
		bool interiorLinks[2][4] = { false,false,false,false,false,false,false,false }, zeroEmpty = true, secondPair = false;
		auto combinePatches = [&](std::list<patch*>& pl) {
			auto cp = combineP.begin();
			while (cp != combineP.end()) {
				for (auto pp : pl) {
					if (cp->find(pp) != cp->end()) {
						for (auto pa : pl)
							cp->insert(pa);
						return;
					}
				}
				++cp;
			}
			if (cp == combineP.end()) {
				combineP.push_back(std::set<patch*>());
				for (auto pp : pl)
					combineP.back().insert(pp);
			}
		};
		auto getEdgeNodeIndices = [](const int edge, int& idx0, int& idx1) {
			if (edge < 1) { idx0 = 3; idx1 = 0; }
			else if (edge < 4) { idx0 = edge - 1; idx1 = edge; }
			else if (edge < 5) { idx0 = 0; idx1 = 2; }
			else { idx0 = 1; idx1 = 3; }
		};
		if (edgeCut) {  // unnecessary if no edge was cut
			int nbegin, nend;
			double sameHit = std::max({ _vbt->_gridSize[0], _vbt->_gridSize[1], _vbt->_gridSize[2] }) * 6e-4;
			for (int i = 0; i < 6; ++i) {
				getEdgeNodeIndices(i, nbegin, nend);
				if (!edges[i].empty()) {
					bool inSolid = inodes[nbegin] > -1;
					auto eit = edges[i].begin();
					auto enext = eit;
					++enext;
					while (enext != edges[i].end()) {
						if (enext->first - eit->first < sameHit) {  // possible coincident surface or minor collision hit.  Both sides guaranteed inside.
							if (eit->second.nextOutside == enext->second.nextOutside) {  // possible surface edge or vertex multi hit
								if (inSolid != eit->second.nextOutside) {
									edges[i].erase(enext);
									enext = eit;
									++enext;
									if (enext == edges[i].end())
										break;
								}
							}
							else if (inSolid != eit->second.nextOutside) {  // flip if necessary
								auto swap = enext->second;
								enext->second = eit->second;
								eit->second = swap;
							}
							else
								;
						}
						if (inSolid != eit->second.nextOutside)
							throw(std::logic_error("Solid ordering error in getConnectedComponents()\n"));
						inSolid = !inSolid;
						eit = enext;
						++enext;
					}
					eit = edges[i].begin();
					assert(eit->second.nextOutside == inodes[nbegin] > -1);
					enext = eit;
					++enext;
					while (enext != edges[i].end()) {
						if (!eit->second.nextOutside && enext->second.nextOutside) {  // solid interval. Combine patches
							if (eit->second.pptr != enext->second.pptr) {  // can be equal
								std::list<patch*> pl;
								pl.push_back(eit->second.pptr);
								pl.push_back(enext->second.pptr);
								combinePatches(pl);
							}
						}
						eit = enext;
						++enext;
					}
					assert(eit->second.nextOutside == inodes[nend] < 0);
					auto et = edges[i].begin();
					if (et->second.nextOutside) {
						assert(inodes[nbegin] > -1);
						et->second.pptr->tetNodes[nbegin] = inodes[nbegin];
					}
					else
						assert(inodes[nbegin] < 0);
					et = edges[i].end();
					--et;
					if (!et->second.nextOutside) {
						assert(inodes[nend] > -1);
						et->second.pptr->tetNodes[nend] = inodes[nend];
					}
					else
						assert(inodes[nend] < 0);
				}
				else {
					if (inodes[nbegin] > -1) {
						assert(inodes[nend] > -1);
						if (zeroEmpty) {
							interiorLinks[0][nbegin] = true;
							interiorLinks[0][nend] = true;
							zeroEmpty = false;
						}
						else if (interiorLinks[0][nbegin] || interiorLinks[0][nend]) { // 3rd or 4th member
							if (secondPair) {
								for (int k = 0; k < 4; ++k)
									interiorLinks[0][k] = true;
								secondPair = false;
							}
							else {
								interiorLinks[0][nbegin] = true;
								interiorLinks[0][nend] = true;
							}
						}
						else {  // can have 2 pairs
							interiorLinks[1][nbegin] = true;
							interiorLinks[1][nend] = true;
							secondPair = true;
						}
					}
				}
			}
			// add linked interior nodes to patches
			for (auto& p : patches) {
				if (p.noEdgeCut)
					continue;
				for (int i = 0; i < 4; ++i) {
					if (p.tetNodes[i] > -1) {
						if (interiorLinks[0][i]) {
							for (int j = 0; j < 4; ++j) {
								if (i == j)
									continue;
								if (interiorLinks[0][j])
									p.tetNodes[j] = inodes[j];
							}
							if (!secondPair)
								break;
						}
						else if (secondPair && interiorLinks[1][i]) {
							for (int j = 0; j < 4; ++j) {
								if (i == j)
									continue;
								if (interiorLinks[1][j]) {
									p.tetNodes[j] = inodes[j];
									break;
								}
							}
						}
						else
							;
					}
				}
			}
			// look for patches that combine due to shared interior node (ie solid span bridges more than one edge)
			for (int i = 0; i < 4; ++i) {
				if (inodes[i] < 0)
					continue;
				std::list<patch*> pl;
				for (auto& p : patches) {
					if (p.noEdgeCut)
						continue;
					if (p.tetNodes[i] > -1) {
						assert(p.tetNodes[i] == inodes[i]);
						pl.push_back(&p);
					}
				}
				if (pl.size() > 1)
					combinePatches(pl);
			}
			// combine patches that enclose a solid
			for (auto& cp : combineP) {
				assert(cp.size() > 1);
				auto cpfirst = cp.begin();
				while (cpfirst != cp.end() && (*cpfirst)->tris.empty()) {
					removePatch(*cpfirst);
					++cpfirst;
				}
				if (cpfirst == cp.end())
					continue;
				auto cpnext = cpfirst;
				++cpnext;
				while (cpnext != cp.end()) {
					if ((*cpnext)->tris.empty()) {
						removePatch(*cpnext);
						++cpnext;
						continue;
					}
					std::vector<int> mergeVec((*cpfirst)->tris.size() + (*cpnext)->tris.size());
					std::merge((*cpfirst)->tris.begin(), (*cpfirst)->tris.end(), (*cpnext)->tris.begin(), (*cpnext)->tris.end(), mergeVec.begin());
					(*cpfirst)->tris = std::move(mergeVec);
					for (int i = 0; i < 4; ++i) {
						if ((*cpnext)->tetNodes[i] > -1)
							(*cpfirst)->tetNodes[i] = (*cpnext)->tetNodes[i];
					}
					removePatch(*cpnext);
					++cpnext;
				}
			}
		}
		if (!tetFaceHoles.empty())
			processHoles(patches, tetFaceHoles, tetNodes, inodes);
	}
	for (auto& cTet : patches) {
		auto ntit = nt_vec.push_back(newTet());
		ntit->tetIdx = _nSurfaceTets.fetch_add(1);
		cTet.tetIndex = ntit->tetIdx;
		ntit->tc = tt.tc;
		ntit->tetNodes = cTet.tetNodes;
		ntit->tris.assign(cTet.tris.begin(), cTet.tris.end());    // COURT change to move after new data structure
	}
	for (int i = 0; i < 4; ++i) {
		bool enNotEntered = true;
		NTS_HASH::accessor acc;
		for (auto& p : patches) {
			if (p.tetNodes[i] < 0) {
				if (enNotEntered) {
					std::array<short, 3> loc = { gl[i][0], gl[i][1], gl[i][2] };
					local_nts.insert(acc, std::make_pair(loc, std::list<nodeTetSegment>()));
					enNotEntered = false;
				}
				nodeTetSegment nts;
				nts.tetNodeIndex = i;
				nts.tetIdx = p.tetIndex;
				nts.tetNodeTris.assign(p.tris.begin(), p.tris.end());  // COURT nuke afternew data structure
				acc->second.push_back(nts);
			}
		}
	}
}

bool vnBccTetCutter_tbb::setupBccIntersectionStructures(int maximumGridDimension)
{
	boundingBox<float> bbf;
	bbf.Empty_Box();
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i)
		bbf.Enlarge_To_Include_Point((const float(&)[3])(*_mt->vertexCoordinate(i)));
	bbf.Minimum_Corner(_vbt->_minCorner.xyz);
	bbf.Maximum_Corner(_vbt->_maxCorner.xyz);
	// In this model all grid distances are 1 and all odd or even Cartesian distances are 2.
	_vbt->_unitSpacing = -1.0f;
	int bigDim = -1;
	for (int i = 0; i < 3; ++i){
		float dimSize = bbf.val[(i << 1) + 1] - bbf.val[i << 1];
		if (dimSize > _vbt->_unitSpacing){
			_vbt->_unitSpacing = dimSize;
			bigDim = i;
		}
	}
	// object must be inside positive octant
	float offset = (_vbt->getMaximumCorner() - _vbt->getMinimumCorner()).length() * 0.0001f;
	_vbt->_minCorner -= Vec3f(offset, offset, offset);
	_vbt->_maxCorner += Vec3f(offset, offset, offset);
	double cs = _vbt->_maxCorner.xyz[bigDim] - _vbt->_minCorner.xyz[bigDim];
	_vbt->_unitSpacing = cs / maximumGridDimension;
	_vbt->_unitSpacingInv = 1.0 / _vbt->_unitSpacing;
	Vec3f maxMaterialCorner = (_vbt->_maxCorner - _vbt->_minCorner)*(float)_vbt->_unitSpacingInv;
	for (int i = 0; i < 3; ++i)
		_vbt->_gridSize[i] = 1 + (int)std::floor(maxMaterialCorner.xyz[i]);
	_vMatCoords.clear();
	_vMatCoords.assign(_mt->numberOfVertices(), Vec3f());
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i){
		Vec3f Vf;
		Vf.set((const float(&)[3])*_mt->vertexCoordinate(i));
		_vMatCoords[i].set((Vf -_vbt->_minCorner)* (float)_vbt->_unitSpacingInv);
	}
	_vertexTetCentroids.clear();
	_vertexTetCentroids.assign(_mt->numberOfVertices(), bccTetCentroid());
	_vbt->_barycentricWeights.clear();
	_vbt->_barycentricWeights.assign(_mt->numberOfVertices(), Vec3f());
	for (int n = _mt->numberOfVertices(), i = 0; i < n; ++i){
		_vbt->gridLocusToLowestTetCentroid(_vMatCoords[i], _vertexTetCentroids[i]);
		// set barycentric coordinate within that tet
		auto& bw = _vbt->_barycentricWeights[i];
		_vbt->gridLocusToBarycentricWeight(_vMatCoords[i], _vertexTetCentroids[i], bw);
		// any vertex on a tet face boundary should be nudged inside
		bool nudge = false;
		for (int j = 0; j < 3; ++j) {
			if (bw[j] < 1e-4f || bw[j] > 0.9999f)
				nudge = true;
		}
		if (bw[0] + bw[1] + bw[2] > 0.9999)
			nudge = true;
		if (nudge) {  // vertices on tet boundarys must be pulled inside to avoid dual identities
			Vec3f nV = (Vec3f(0.25f, 0.25f, 0.25f) - bw) * 0.0002f;
			bw += nV;
			_vbt->barycentricWeightToGridLocus(_vertexTetCentroids[i], bw, _vMatCoords[i]);
		}
	}
	// setup lines parallel with Z axis
		// setup lines parallel with Z axis
	evenXy.assign(_vbt->_gridSize[0] >> 1, std::vector<std::multimap<double, zIntersectFlags> >());  // 0th i always empty
	oddXy.assign(_vbt->_gridSize[0] >> 1, std::vector<std::multimap<double, zIntersectFlags> >());
	int gsy = _vbt->_gridSize[1] >> 1;
	for (int n = _vbt->_gridSize[0] >> 1, i = 0; i < n; ++i) {
		evenXy[i].assign(gsy, std::multimap<double, zIntersectFlags>());  // 0th j always empty
		oddXy[i].assign(gsy, std::multimap<double, zIntersectFlags>());
	}
	return true;
}

void vnBccTetCutter_tbb::zIntersectTriangleTbb(Vec3d(&tri)[3], const bool surfaceTriangle, oneapi::tbb::concurrent_vector<zIntrsct>& zi_loc) {
	// warning not const. Changes values in tri
	int xy[4] = { INT_MAX, -2, INT_MAX, -2 };
	for (int i = 0; i < 3; ++i) {
		if (tri[i][0] < xy[0])
			xy[0] = ceil(tri[i][0]);
		if (tri[i][0] > xy[1])
			xy[1] = floor(tri[i][0]);
		if (tri[i][1] < xy[2])
			xy[2] = ceil(tri[i][1]);
		if (tri[i][1] > xy[3])
			xy[3] = floor(tri[i][1]);
	}
	if (!surfaceTriangle) {  // tet triangles can be partially out of bounds
		if (xy[0] < 1)
			xy[0] = 1;
		if (xy[1] >= _vbt->_gridSize[0])
			xy[1] = _vbt->_gridSize[0] - 1;
		if (xy[2] < 1)
			xy[2] = 1;
		if (xy[3] >= _vbt->_gridSize[1])
			xy[3] = _vbt->_gridSize[1] - 1;
	}
	else {
		assert(xy[0] > 0);
		assert(xy[1] < _vbt->_gridSize[0]);
		assert(xy[2] > 0);
		assert(xy[3] < _vbt->_gridSize[1]);

	}
	// now get any Z line intersects
	tri[1] -= tri[0];
	tri[2] -= tri[0];
	double zCut = tri[1][0] * tri[2][1] - tri[1][1] * tri[2][0];
	if (!surfaceTriangle && abs(zCut) < 1e-8)  // tet triangle parallel to Z axis. No intersect possible other than planar which can be ignored.
		return;
	Mat2x2d M;
	M.x[0] = tri[1][0];  M.x[1] = tri[1][1];  M.x[2] = tri[2][0]; M.x[3] = tri[2][1];
	zIntrsct zi;
	zi.flags.surfaceTri = surfaceTriangle;
	zi.flags.solidBegin = zCut < 0.0;  // negative Z starts a solid
	auto triIntersectZ = [&](const int x, const int y, float& z) ->bool {
		Vec2d R = M.Robust_Solve_Linear_System(Vec2d(x - tri[0][0], y - tri[0][1]));
		if (R[0] < -1e-8 || R[0] > 1.0000001 || R[1] < -1e-8 || R[1] > 1.00000001 || R[0] + R[1] >= 1.0000001)  // don't want vertex hits, only edges
			return false;
		z = (tri[0][2] + tri[1][2] * R[0] + tri[2][2] * R[1]);
		return true;
	};
	for (int i = xy[0]; i <= xy[1]; ++i) {
		bool odd = i & 1;
		for (int j = xy[2]; j <= xy[3]; ++j) {
			if (odd != (bool)(j & 1))
				continue;
			if (odd) {
				if (triIntersectZ(i, j, zi.zInt)) {
					zi.flags.odd = true;
					zi.x = (i - 1) >> 1;
					zi.y = (j - 1) >> 1;
					zi_loc.push_back(zi);
				}
			}
			else {
				if (triIntersectZ(i, j, zi.zInt)) {
					zi.flags.odd = false;
					zi.x = (i - 2) >> 1;
					zi.y = (j - 2) >> 1;
					zi_loc.push_back(zi);
				}
			}
		}
	}
}

void vnBccTetCutter_tbb::addCentroidMicronodesZ(const bccTetCentroid& tc) {
	short gl[4][3];
	_vbt->centroidToNodeLoci(tc, gl);
	int ha;
	for (ha = 0; ha < 3; ++ha)
		if (gl[0][ha] == gl[1][ha] && gl[2][ha] == gl[3][ha])
			break;
	assert(ha < 3);
	int c1 = ha > 1 ? 0 : ha + 1;
	int c2 = c1 > 1 ? 0 : c1 + 1;
	int rect[2][2] = { gl[0][c1], gl[1][c1], gl[0][c2], gl[0][c2]}, dha = abs(gl[0][ha] - gl[2][ha]);
	int xyz[3], hincr = gl[0][ha] > gl[2][ha] ? -1 : 1;
	xyz[ha] = gl[0][ha];
	auto zMatrixAdd = [&](bool isMacroNode) {
		zIntersectFlags zf;
		zf.macroNode = isMacroNode;
		zf.odd = zf.solidBegin = 1;
		zf.surfaceTri = 0;
		if (xyz[0] & 1) {
			int xi = (xyz[0] - 1) >> 1, yi = (xyz[1] - 1) >> 1;
			if (xi >= oddXy.size() || yi >= oddXy[xi].size())  // out of solid range so can't be an interior node
				return;
			auto& zLine = oddXy[xi][yi];
			auto pr = zLine.equal_range(xyz[2]);
			if (pr.first == pr.second)
				zLine.emplace_hint(pr.first, std::make_pair((double)xyz[2], zf));
			else {
				auto it = pr.first;
				while (it != pr.second) {
					if (!it->second.surfaceTri)
						return;
					++it;
				}
				zLine.emplace_hint(pr.first, std::make_pair((double)xyz[2], zf));
			}
		}
		else {
			int xi = (xyz[0] - 2) >> 1, yi = (xyz[1] - 2) >> 1;
			if (xi >= evenXy.size() || yi >= evenXy[xi].size())  // out of solid range so can't be an interior node
				return;
			auto& zLine = evenXy[xi][yi];
			auto pr = zLine.equal_range(xyz[2]);
			if (pr.first == pr.second)
				zLine.emplace_hint(pr.first, std::make_pair((double)xyz[2], zf));
			else {
				auto it = pr.first;
				while (it != pr.second) {
					if (!it->second.surfaceTri)
						return;
					++it;
				}
				zLine.emplace_hint(pr.first, std::make_pair((double)xyz[2], zf));
			}
		}
	};
	for (int i = 0; i <= dha; ++i) {
		xyz[c1] = rect[0][0];
		while (xyz[c1] <= rect[0][1]) {
			xyz[c2] = rect[1][0];
			while (xyz[c2] <= rect[1][1]) {
				if ((i < 1 && (xyz[c1] == rect[0][0] || xyz[c1] == rect[0][1])) || (i == dha && (xyz[c2] == rect[1][0] || xyz[c2] == rect[1][1])))
					zMatrixAdd(true);
				else
					zMatrixAdd(false);
				xyz[c2] += 2;
			}
			xyz[c1] += 2;
		}
		++rect[0][0];
		--rect[0][1];
		--rect[1][0];
		++rect[1][1];
		xyz[ha] += hincr;
	}
}

void vnBccTetCutter_tbb::inputTriangleTetsTbb(const int& surfaceTriangle, CENTtris& centTris) {
	int* tr = _mt->triangleVertices(surfaceTriangle);
	Vec3d T[3];
	bccTetCentroid tc[3];
	for (int i = 0; i < 3; ++i) {
		T[i] = _vMatCoords[tr[i]];
		tc[i] = _vertexTetCentroids[tr[i]];
	}
	// triangle-triangle intersection used to non-recursively grow tets who have a face intersecting triangle
	// Of 3 algorithms; Moller, Devillers-Guigue, and Shen-Heng-Tang this routine uses the last of these.
	// Routine modified to set one triangle fixed to evaluate against multiple others.
	std::set<bccTetCentroid> triTetsS;
	for (int i = 0; i < 3; ++i)
		triTetsS.insert(tc[i]);
	auto recurseTriangleTets = [&](bccTetCentroid& tetC) {
		std::map<bccTetCentroid, short> doTets;  // tet that need processing and its face which can be ignored
		doTets.insert(std::make_pair(tetC, -1));
		short gridLoci[4][3];
		auto dtit = doTets.begin();
		while (dtit != doTets.end()) {
			triTetsS.insert(dtit->first);
			_vbt->centroidToNodeLoci(dtit->first, gridLoci);
			bccTetCentroid tAdj;
			short adjFace;
			for (int i = 0; i < 4; ++i) {
				if (i == dtit->second)  // face already processed
					continue;
				adjFace = _vbt->faceAdjacentMultiresTet(dtit->first, i, tAdj);
				if (triTetsS.find(tAdj) != triTetsS.end()) // already found
					continue;
				Vec3d V0(gridLoci[i]), V1(gridLoci[(i + 1) & 3]), V2(gridLoci[(i + 2) & 3]);

				if (tri_tri_intersectD(V0.xyz, V1.xyz, V2.xyz, T[0].xyz, T[1].xyz, T[2].xyz))
					doTets.insert(std::make_pair(tAdj, adjFace));
			}
			doTets.erase(dtit);
			dtit = doTets.begin();
		}
	};
	if (triTetsS.size() > 2 || (triTetsS.size() > 1 && !_vbt->adjacentMicrotetCentroids(*triTetsS.begin(), *triTetsS.rbegin()))) {
		triTetsS.clear();
		recurseTriangleTets(tc[0]);
		for (int i = 1; i < 3; ++i) {
			if (triTetsS.find(tc[i]) == triTetsS.end()) // should be there if recurse OK
				recurseTriangleTets(tc[i]);
		}
	}
	for (auto& tt : triTetsS) {
		CENTtris::accessor acc;
		bool result = centTris.insert(acc, std::make_pair(tt, std::vector<int>({ surfaceTriangle })));
		if(!result)
			acc->second.push_back(surfaceTriangle);
	}
}

void vnBccTetCutter_tbb::fillInteriorMicroTets(std::vector<bccTetCentroid> &recutMacrotets) {
	// COURT - this first pass may be able to be done more efficiently
	std::vector<bccTetCentroid> microTets;  // only interested in unique interior ones
	microTets.reserve(recutMacrotets.size() * (2 << _vbt->_tetSubdivisionLevels) * 16);  // COURT check  
	for (auto& mt : recutMacrotets) {
		std::vector<bccTetCentroid> tcIn, nextTc;
		tcIn.push_back(mt);
		bccTetCentroid subtets[8];
		int s = _vbt->_tetSubdivisionLevels;
		while (s > 1) {
			for (auto ti : tcIn) {
				_vbt->subtetCentroids(ti, subtets);
				for(int j=0; j<8; ++j)
					if(subtets[j][0] < USHRT_MAX)  // out of positive octant
						nextTc.push_back(subtets[j]);
			}
			tcIn = std::move(nextTc);
			nextTc.clear();
			--s;
		}
		microTets.insert(microTets.end(), tcIn.begin(), tcIn.end());
	}
	// COURt change reserve size of _vbt->_tetNodes and _tetCentroids
	for (auto &mt : microTets) {
		std::array<short, 3> gl[4];
		_vbt->centroidToNodeLoci(mt, (short (&)[4][3])gl);
		std::array<int, 4> tn;
		int i;
		for (i = 0; i < 4; ++i) {
			auto iit = _interiorNodes.find(gl[i]);
			if (iit == _interiorNodes.end())
				break;
			tn[i] = iit->second;
		}
		if (i < 4)
			continue;
		if (_surfaceCentroids.find(mt) != _surfaceCentroids.end())
			continue;
		_vbt->_tetNodes.push_back(tn);
		_vbt->_tetCentroids.push_back(mt);
	}
}

void vnBccTetCutter_tbb::fillNonVnTetCenter()
{
	_vbt->_firstInteriorTet = _vbt->_tetNodes.size();
	//  Use sequential z to get surrounding frame with less hash table searches.
	auto nextZnode = [&](const int prevNode, std::array<short, 3>& prevLocus) ->int {  // returns next node in z sequence from prevNode at locus. If former is -1 look for locus. If not found return -1
		prevLocus[2] += 2;
		if (prevNode > -1) {
			auto pn = _vbt->_nodeGridLoci[prevNode + 1].data();
			if (prevLocus[0] == pn[0] && prevLocus[1] == pn[1] && prevLocus[2] == pn[2])
				return prevNode + 1;
		}
		auto in = _interiorNodes.find(prevLocus);
		if (in != _interiorNodes.end())
			return in->second;
		else return -1;
	};
	auto findNode = [&](std::array<short, 3>& locus) ->int {
		auto in = _interiorNodes.find(locus);
		if (in != _interiorNodes.end())
			return in->second;
		else
			return -1;

	};
	auto zBlock = [&](int start, int stop) {
		std::array<int, 4> prevQuad, nextQuad;
		std::array<short, 3> locs[4], yLoc, xLoc = _vbt->_nodeGridLoci[start];
		yLoc = xLoc;
		yLoc[1] += 2;
		int nextY = findNode(yLoc);
		xLoc[0] += 2;
		int nextX = findNode(xLoc);
		locs[0] = xLoc;
		locs[0][0] -= 3;
		--locs[0][1];
		--locs[0][2];
		locs[1] = locs[0];
		locs[1][0] += 2;
		locs[2] = locs[0];
		locs[2][1] += 2;
		locs[3] = locs[2];
		locs[3][0] += 2;
		for (int i = 0; i < 4; ++i)
			prevQuad[i] = findNode(locs[i]);
		for (int i = start; i < stop; ++i) {
			for (int j = 0; j < 4; ++j)
				nextQuad[j] = nextZnode(prevQuad[j], locs[j]);
			// make tets
			bccTetCentroid center = { (unsigned short)(xLoc[0] - 2), (unsigned short)xLoc[1], (unsigned short)(xLoc[2] + 1) };
			for (int j = 0; j < 3; ++j)
				center[j] <<= 1;
			std::array<int, 4> tet;
			if (i < stop - 1) { // next even z ray
				tet[0] = i;
				tet[1] = i + 1;
				if (nextQuad[0] > -1 && nextQuad[1] > -1) {
					tet[2] = nextQuad[0];
					tet[3] = nextQuad[1];
					auto tc = center;
					--tc[1];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (nextQuad[2] > -1 && nextQuad[3] > -1) {
					tet[2] = nextQuad[3];
					tet[3] = nextQuad[2];
					auto tc = center;
					++tc[1];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (nextQuad[0] > -1 && nextQuad[2] > -1) {
					tet[0] = nextQuad[0];
					tet[1] = nextQuad[2];
					tet[2] = i + 1;
					tet[3] = i;
					auto tc = center;
					--tc[0];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (nextQuad[1] > -1 && nextQuad[3] > -1) {
					tet[0] = nextQuad[1];
					tet[1] = nextQuad[3];
					tet[2] = i;
					tet[3] = i + 1;
					auto tc = center;
					++tc[0];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
			}
			center[2] -= 2;  // next even x ray
			center[0] += 2;
			if (nextX > -1) {
				tet[0] = i;
				tet[1] = nextX;
				if (nextQuad[3] > -1 && nextQuad[1] > -1) {
					tet[2] = nextQuad[3];
					tet[3] = nextQuad[1];
					auto tc = center;
					++tc[2];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[1] > -1 && prevQuad[3] > -1) {
					tet[2] = prevQuad[1];
					tet[3] = prevQuad[3];
					auto tc = center;
					--tc[2];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[1] > -1 && nextQuad[1] > -1) {
					tet[0] = prevQuad[1];
					tet[1] = nextQuad[1];
					tet[2] = nextX;
					tet[3] = i;
					auto tc = center;
					--tc[1];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[3] > -1 && nextQuad[3] > -1) {
					tet[0] = prevQuad[3];
					tet[1] = nextQuad[3];
					tet[2] = i;
					tet[3] = nextX;
					auto tc = center;
					++tc[1];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
			}
			center[0] -= 2;
			center[1] += 2;
			if (nextY > -1) {
				tet[0] = i;
				tet[1] = nextY;
				if (prevQuad[2] > -1 && nextQuad[2] > -1) {
					tet[2] = prevQuad[2];
					tet[3] = nextQuad[2];
					auto tc = center;
					--tc[0];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[3] > -1 && nextQuad[3] > -1) {
					tet[2] = nextQuad[3];
					tet[3] = prevQuad[3];
					auto tc = center;
					++tc[0];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (prevQuad[2] > -1 && prevQuad[3] > -1) {
					tet[0] = prevQuad[2];
					tet[1] = prevQuad[3];
					tet[2] = nextY;
					tet[3] = i;
					auto tc = center;
					--tc[2];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
				if (nextQuad[2] > -1 && nextQuad[3] > -1) {
					tet[0] = nextQuad[2];
					tet[1] = nextQuad[3];
					tet[2] = i;
					tet[3] = nextY;
					auto tc = center;
					++tc[2];
					if (_surfaceCentroids.find(tc) == _surfaceCentroids.end()) {  // may already be there if penetration by surface. Only create unpenetrated interior cubes.
						_vbt->_tetCentroids.push_back(tc);
						_vbt->_tetNodes.push_back(tet);
					}
				}
			}
			// advance to next frame
			for (int j = 0; j < 4; ++j)
				prevQuad[j] = nextQuad[j];
			nextX = nextZnode(nextX, xLoc);
			nextY = nextZnode(nextX, yLoc);
		}
	};
	// interior tet nodes are created sequentially in z.  Do each z block as nucleus around even lattice
	for (int n = _interiorNodes.size(), i = 0; i < n; ) {
		auto lp0 = _vbt->_nodeGridLoci[i].data();
		if (lp0[0] & 1) {
			++i;
			continue;
		}
		int j = i + 1;
		auto lp1 = _vbt->_nodeGridLoci[j].data();
		while (lp0[0] == lp1[0] && lp0[1] == lp1[1] && lp0[2] + 2 == lp1[2]) { // this z segment exists
			++j;
			lp0 = lp1;
			lp1 = _vbt->_nodeGridLoci[j].data();
		}
		zBlock(i, j);
		i = j;
	}
}

bool vnBccTetCutter_tbb::latticeTest() {

	// COURT - finish me.  This could be more complete

	bool ret = true;
	for (int i = _vbt->_nMegatets; i < _vbt->tetNumber(); ++i) {
		int* tn = _vbt->_tetNodes[i].data();
		std::array<unsigned int, 3> ct = { 0, 0, 0 };
		for (int j = 0; j < 4; ++j) {
			for (int k = 0; k < 3; ++k)
				ct[k] += _vbt->_nodeGridLoci[tn[j]][k];
		}
		bccTetCentroid tc = _vbt->_tetCentroids[i];
		auto pr = _vbt->_tetHash.equal_range(tc);
		while (pr.first != pr.second) {
			if (pr.first->second == i)
				break;
			++pr.first;
		}
		if (pr.first == pr.second) {
			std::cout << "Tet centroid not found for tet " << i << "\n";
			ret = false;
		}
		for (int k = 0; k < 3; ++k) {
			ct[k] >>= 1;
			if (ct[k] != tc[k]) {
				std::cout << "Bad tet centroid at tet " << i << "\n";
				ret = false;
			}
		}
	}
	for (int i = _meganodeSize; i < _firstNewExteriorNode; ++i) {
		bccTetCentroid cntrd[24];
		_vbt->nodeMicroCentroids(_vbt->_nodeGridLoci[i], cntrd);  // for these routines a coordinate greater than 65533 indicates a centroid out of positive grid bounds
		// an interior node should all be sourrounded by tets containing it
		for (int k, j = 0; j < 24; ++j) {
			auto pr = _vbt->_tetHash.equal_range(cntrd[j]);
			while (pr.first != pr.second) {
				int* tn = _vbt->_tetNodes[pr.first->second].data();
				for (k = 0; k < 4; ++k) {
					if (tn[k] == i)
						break;
				}
				if (k<4)
					break;
				++pr.first;
			}
			if(pr.first == pr.second){  // OK if T junction
				if (_vbt->_tJunctionConstraints.find(i) != _vbt->_tJunctionConstraints.end())
					continue;
				else {
					std::cout << "Interior node " << i << " is not surrounded by tc " << cntrd[j][0] << ", " << cntrd[j][1] << ", " << cntrd[j][2] << "\n";
					ret = false;
				}
			}
		}
	}
	return ret;
}
