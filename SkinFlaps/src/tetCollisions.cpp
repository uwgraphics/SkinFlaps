#include <unordered_set>
#include <array>
#include <map>
#include "boundingBox.h"
#include "Mat2x2f.h"
#include "insidePolygon.h"
#include "materialTriangles.h"
#include "vnBccTetrahedra.h"
#include "pdTetPhysics.h"
#include "tbb/tbb.h"
#include "tetCollisions.h"

#include "tbb/tick_count.h"  // for debug nuke later

materialTriangles *tetCollisions::_mt;
vnBccTetrahedra *tetCollisions::_vnt;
pdTetPhysics* tetCollisions::_ptp;

void tetCollisions::initSoftCollisions(materialTriangles* mt, vnBccTetrahedra* vnt) {
	_mt = mt;
	_vnt = vnt;
	if (!_initialized) {
		// generate 6 cardinal material coordinate inverses for bcc tet deformation gradients
		for (int i = 0; i < 3; ++i) {
			Mat3x3f M(0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f, 1.0f);
			M((i + 1) % 3, 0) = 2.0f;
			M((i + 2) % 3, 2) *= -1.0f;
			M *= (float)_vnt->getTetUnitSize();
			_rest[i << 1] = M.Inverse();
			M((i + 2) % 3, 1) *= -1.0f;
			M((i + 2) % 3, 2) *= -1.0f;
			M(i, 1) *= -1.0f;
			M(i, 2) *= -1.0f;
			_rest[(i << 1) + 1] = M.Inverse();
		}
		_initialized = true;
	}
	std::unordered_set<int> tets;
	tets.reserve(2048);
	_flapBottomTris.clear();
	for (auto& bedV : _bedVerts)
		bedV.second.materialNormal.set(0.0f, 0.0f, 0.0f);
	auto inputFlapBottomTriangle = [&](int triangle) {
		for (int j = 0; j < 3; ++j) {
			if (_mt->triangleMaterial(_mt->triAdjs(triangle)[j] >> 2) == 5)  // hinge triangle, don't process
				return;
		}
		_flapBottomTris.push_back(triangle);
		long* tr = _mt->triangleVertices(triangle);
		for (int j = 0; j < 3; ++j)
			tets.insert(_vnt->getVertexTetrahedron(tr[j]));
	};
	std::set<int> hingeVerts;
	auto inputBedTriangle = [&](int triangle) {
		long* tr = _mt->triangleVertices(triangle);
		bottomRay* br[3];
		Vec3f matPos[3];
		for (int j = 0; j < 3; ++j) {
			if (_mt->triangleMaterial(_mt->triAdjs(triangle)[j] >> 2) == 4) {  // hinge triangle.  Do normal computation, but remove hinge verts later.
				hingeVerts.insert(tr[j]);
				hingeVerts.insert(tr[(j + 1) % 3]);
			}
			int tet = _vnt->getVertexTetrahedron(tr[j]);
			auto tc = _vnt->tetCentroid(tet);
			vnt->barycentricWeightToGridLocus(*tc, *vnt->getVertexWeight(tr[j]), matPos[j]);
			auto bv = _bedVerts.emplace(tr[j], bottomRay());
			if (bv.second) {
				bv.first->second.vertex = tr[j];
				bv.first->second.materialNormal.set(0.0f, 0.0f, 0.0f);
				bv.first->second.restIdx = tc->halfCoordAxis << 1;
				if ((tc->xyz[tc->halfCoordAxis] + tc->xyz[(tc->halfCoordAxis + 2) % 3]) & 1)
					++(bv.first->second.restIdx);
			}
			br[j] = &bv.first->second;
		}
		Vec3f v0, v1, N;
		v0 = matPos[2] - matPos[0];
		v1 = matPos[1] - matPos[0];
		N = v0 ^ v1;
		for (int j = 0; j < 3; ++j)
			br[j]->materialNormal += N;
	};
	for (size_t n = _mt->numberOfTriangles(), i = 0; i < n; ++i) {
		if (_mt->triangleMaterial(i) == 3) {
			int adjMat = _mt->triangleMaterial(_mt->triAdjs(i)[0] >> 2);
			if (adjMat == 5) {
				if (_mt->triangleMaterial(_mt->triAdjs(i - 1)[0] >> 2) != 2)
					throw(std::logic_error("Incision convention violated in initializing soft colli8sions."));
				inputBedTriangle(i);
				inputBedTriangle(i - 1);
			}
			else if (adjMat == 4) {  // should only the lower of the 2 edge triangles be processed?
				if (_mt->triangleMaterial(_mt->triAdjs(i - 1)[0] >> 2) != 2)
					throw(std::logic_error("Incision convention violated in initializing soft colli8sions."));
				inputFlapBottomTriangle(i);
//				inputFlapBottomTriangle(i-1);
			}
			else;  // top incion edge tris processed with their bottom edge tri
		}
		else if (_mt->triangleMaterial(i) == 4)
			inputFlapBottomTriangle(i);
		else if (_mt->triangleMaterial(i) == 5)
			inputBedTriangle(i);
		else
			;
	}
	_flapBottomTris.shrink_to_fit();

	// vertices on hinge boundary between bed and bottom tris can be removed from bed set as collision processing not necessary
	for (auto h : hingeVerts)
		_bedVerts.erase(h);
	_bedRays.clear();
	_bedRays.reserve(_bedVerts.size());
	for (auto& bedV : _bedVerts) {
//		tets.insert(vnt->getVertexTetrahedron(bedV.second.vertex));
		if (!tets.insert(vnt->getVertexTetrahedron(bedV.second.vertex)).second)  // for now tets must be unique
			continue;
		bedV.second.materialNormal.normalize();
		bedV.second.materialNormal *= 0.1f;  // ?scale?
		_bedRays.push_back(&bedV.second);
	}
	_bedRays.shrink_to_fit();
	// _bedRay crossover ignored since will only use shortest one

	//	// run test - results verified
	//	const long *nodes = _vnt->tetNodes(tet);
	//	Vec3f cols[3];
	//	for (int i = 1; i < 4; ++i)
	//		cols[i-1] = _vnt->nodeSpatialCoordinate(nodes[i]) - _vnt->nodeSpatialCoordinate(nodes[0]);
	//	Mat3x3f Q, N(cols[0], cols[1], cols[2]);
	//	Q = N * _rest[r.restIdx];

	if (!tets.empty()) {
		std::vector<long> tetras;
		tetras.assign(tets.begin(), tets.end());
		_ptp->addSoftCollisionTets(tetras);
	}
}

void tetCollisions::findSoftCollisionPairs() {
	if (_flapBottomTris.empty())
		return;

//	tbb::tick_count t0 = tbb::tick_count::now();

	std::vector<boundingBox<float> > flapBox, bedBox;
	boundingBox<float> bb;
	bb.Empty_Box();
	flapBox.assign(_flapBottomTris.size(), bb);
	bedBox.assign(_bedRays.size(), bb);
	for (size_t n = _flapBottomTris.size(), i = 0; i < n; ++i) {
		long* tr = _mt->triangleVertices(_flapBottomTris[i]);
		for (int j = 0; j < 3; ++j)
			flapBox[i].Enlarge_To_Include_Point(reinterpret_cast<const float(&)[3]>(*_mt->vertexCoordinate(tr[j])));
	}
	for (size_t n = _bedRays.size(), i = 0; i < n; ++i) {
		bottomRay* b = _bedRays[i];
		const long* nodes = _vnt->tetNodes(_vnt->getVertexTetrahedron(b->vertex));
		const Vec3f* bw = _vnt->getVertexWeight(b->vertex);
		const Vec3f* sc[4];
		sc[0] = &_vnt->nodeSpatialCoordinate(nodes[0]);
		b->P = *sc[0] * (1.0f - bw->_v[0] - bw->_v[1] - bw->_v[2]);
		for (int j = 1; j < 4; ++j) {
			sc[j] = &_vnt->nodeSpatialCoordinate(nodes[j]);
			b->P += *sc[j] * bw->_v[j - 1];
		}
		Mat3x3f M(*sc[1] - *sc[0], *sc[2] - *sc[0], *sc[3] - *sc[0]);
		// only direction of b->materialNormal matters since deformation gradient will change length.  Normalize here.
		b->N = _rest[b->restIdx] * M * b->materialNormal;
		float scale = inverse_rsqrt(b->N * b->N);
		b->N *= scale * 0.1f;  // QISI and YUTIAN - play with fudge factor here.
		bedBox[i].Empty_Box();
		bedBox[i].Enlarge_To_Include_Point(b->P._v);
		bedBox[i].Enlarge_To_Include_Point((b->P + b->N)._v);
	}
	std::vector<long> topTets, bottomTets;
	topTets.assign(_bedRays.size(), -1);
	bottomTets.assign(_bedRays.size(), -1);
	std::vector<std::array<float, 3> > topBarys, bottomBarys, collisionNormals;
	topBarys.assign(_bedRays.size(), std::array<float, 3>());
	bottomBarys.assign(_bedRays.size(), std::array<float, 3>());
	collisionNormals.assign(_bedRays.size(), std::array<float, 3>());  // should I do this unique spec on my side or continue having Qisi do it?

	tbb::parallel_for(tbb::blocked_range<size_t>(0, _bedRays.size()),
		[&](const tbb::blocked_range<size_t>& r) {
			for (size_t j = r.begin(); j != r.end(); ++j) {
//	for (int n = _bedRays.size(), j = 0; j < n; ++j) {  // serial version
				float nearT = FLT_MAX;
				int nearV = -1;
				bottomRay* bestBottom = nullptr;

				//			int bottomTri = -1;
				//			float bottomUv[2];

				for (size_t n = _flapBottomTris.size(), i = 0; i < n; ++i) {
					bottomRay* b = _bedRays[j];
					if (bedBox[j].Intersection(flapBox[i])) {
						long* tr = _mt->triangleVertices(_flapBottomTris[i]);
						Vec3f* tv[3];
						for (int k = 0; k < 3; ++k)
							tv[k] = reinterpret_cast<Vec3f*>(_mt->vertexCoordinate(tr[k]));
						Mat3x3f C(*tv[1] - *tv[0], *tv[2] - *tv[0], -b->N);
						Vec3f R = C.Robust_Solve_Linear_System(b->P - *tv[0]);
						if (R[0] < 1e-6f || R[1] < 1e-6f || R[2] < 2e-1f || R[0] + R[1] > 1.0f || R[0] > 1.0f || R[1] > 1.0f || R[2] > 1.0f)  // R[2] determines how deep the collision must go before processing triggered. Bigger makes less sticky.
							continue;
						if (nearT > R[2]) {

							//						bottomTri = i;
							//						bottomUv[0] = R[0];
							//						bottomUv[1] = R[1];

							bestBottom = b;
							nearT = R[2];

							if (R[0] + R[1] < 0.66667f)
								nearV = tr[0];
							else if (R[0] > R[1])
								nearV = tr[1];
							else
								nearV = tr[2];
						}
					}
				}

				//			if (bottomTri < 0)
				//				continue;

				if (nearV < 0)
					continue;
				// found soft-soft collision pair
				topTets[j] = _vnt->getVertexTetrahedron(nearV);
				const Vec3f* W = _vnt->getVertexWeight(nearV);
				std::array<float, 3> tmp = { W->X, W->Y, W->Z };
				//			Vec3f gridLocus, bw;
				//			bccTetCentroid tc;
				//			topTets[j] = parametricMTtriangleTet(bottomTri, bottomUv, gridLocus, tc);  // problem with this version is all tets intersecting bottom of flaps must be initialized making Schur complement large
				//			if (topTets[j] < 0)
				//				continue;
				//			_vnt->gridLocusToBarycentricWeight(gridLocus, tc, bw);
				//			tmp = { bw.X, bw.Y, bw.Z };
				topBarys[j] = tmp;
				bottomTets[j] = _vnt->getVertexTetrahedron(bestBottom->vertex);
				W = _vnt->getVertexWeight(bestBottom->vertex);
				tmp = { W->X, W->Y, W->Z };
				bottomBarys[j] = tmp;
				Vec3f N = bestBottom->N * nearT;
				tmp = { N.X, N.Y, N.Z };
				collisionNormals[j] = tmp;
			}
		}
	);

	int offset = 0;
	for (int n = topTets.size(), i = 0; i < n; ++i) {
		if (topTets[i] > -1) {
			topTets[offset] = topTets[i];
			bottomTets[offset] = bottomTets[i];
			topBarys[offset] = topBarys[i];
			bottomBarys[offset] = bottomBarys[i];
			collisionNormals[offset++] = collisionNormals[i];
		}
	}
	topTets.resize(offset);
	bottomTets.resize(offset);
	topBarys.resize(offset);
	bottomBarys.resize(offset);
	collisionNormals.resize(offset);

/*	tbb::tick_count t1 = tbb::tick_count::now();
	double time = (t1 - t0).seconds();
	if (_maxTime < time)
		_maxTime = time;
	if (_minTime > time)
		_minTime = time;
	std::cout << "Soft collision minimum processing time was " << _minTime << " and maximum processing time was " << _maxTime << "verts colliding " << topTets.size() << "\n"; */
	// tbb speeds up my 20 thread machine by a factor of 4 to 0.4495 milliseconds so overhead significant and on low core machine may not be worth it

/*	topTets.clear();
	bottomTets.clear();
	topBarys.clear();
	bottomBarys.clear();
	collisionNormals.clear(); */

	_ptp->currentSoftCollisionPairs(topTets, topBarys, bottomTets, bottomBarys, collisionNormals);
}

long tetCollisions::parametricMTtriangleTet(const int mtTriangle, const float(&uv)[2], Vec3f &gridLocus, bccTetCentroid &tC)
{  // in material coords
	long* tr = _mt->triangleVertices(mtTriangle);
	Vec3f tV[3];
	for (int i = 0; i < 3; ++i)
		_vnt->vertexGridLocus(tr[i], tV[i]);
	gridLocus = tV[0] * (1.0f - uv[0] - uv[1]) + tV[1] * uv[0] + tV[2] * uv[1];
	_vnt->gridLocusToTetCentroid(gridLocus, tC);
	for (int i = 0; i < 3; ++i) {
		if (tC.ll == _vnt->tetCentroid(_vnt->getVertexTetrahedron(tr[i]))->ll)
			return _vnt->getVertexTetrahedron(tr[i]);
	}

	return -1;

	// find candidate cubes
	std::list<long> cc, tp;
	auto pr = _vnt->getTetHash()->equal_range(tC.ll);
	while (pr.first != pr.second) {
		cc.push_back(pr.first->second);
		++pr.first;
	}
	if (cc.size() < 1) {
//		assert(false);
		return -1;
	}
	if (cc.size() < 2)
		return cc.front();
	for (auto c : cc) {
		for (int i = 0; i < 3; ++i) {
			if (_vnt->decreasingCentroidPath(c, _vnt->getVertexTetrahedron(tr[i]), tp))
				return c;
		}
	}
//	assert(false);
	return -1;
}



float tetCollisions::inverse_rsqrt(float number)
{  // usual Quake cheat
	const float threehalfs = 1.5F;
	float x2 = number * 0.5F;
	float y = number;
	// evil floating point bit level hacking 
	long i = *(long *)&y;
	// value is pre-assumed 
	i = 0x5f3759df - (i >> 1);
	y = *(float *)&i;
	// 1st iteration 
	y = y * (threehalfs - (x2 * y * y));
	return y;
}

void tetCollisions::addFixedCollisionSet(materialTriangles* mt, const std::string& levelSetFile, const std::vector<Vec2f>& txPoly) {  // call once at load

//	return;

	insidePolygon ip;
	std::set<long> collisionVertices;
	for (int n = mt->numberOfTriangles(), i = 0; i < n; ++i) {
		if (mt->triangleMaterial(i) != 2)
			continue;
		const long* tr = mt->triangleTextures(i);
		for (int j = 0; j < 3; ++j) {
			const float* fp = mt->getTexture(tr[j]);
			Vec2f V(fp[0], fp[1]);
			if (ip.insidePolygon2f(V, txPoly))
				collisionVertices.insert(tr[j]).second;
		}
	}
	fixedCollisionSet fc;
	fc.levelSetFilename = levelSetFile;
	fc.vertices.assign(collisionVertices.begin(), collisionVertices.end());
	_fixedCollisionSets.push_back(fc);
}

void tetCollisions::updateFixedCollisions(materialTriangles *mt, vnBccTetrahedra *vnt) {  // call after every topo change
	_mt = mt;
	_vnt = vnt;
	for (auto& fc : _fixedCollisionSets) {
		std::vector<long> tets;
		std::vector<std::array<float, 3> > weights;
		tets.reserve(fc.vertices.size());
		weights.reserve(fc.vertices.size());
		for (auto &v : fc.vertices) {
			int tet = _vnt->getVertexTetrahedron(v);
			if (tet < 0)  // vertex has been excised
				continue;
			tets.push_back(tet);
			const Vec3f* v3 = _vnt->getVertexWeight(v);
			std::array<float, 3> va = { v3->Y, v3->X,v3->Z };
			weights.push_back(va);
		}
		_ptp->addFixedCollisionSet(fc.levelSetFilename, tets, weights);
	}
}
