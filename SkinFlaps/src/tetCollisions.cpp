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

	// now trim lengths of _botRays to prevent crossover. N*logN
/*	for (int n = (int)_botRays.size(), i = 0; i < n; ++i) {
		Vec3f P = _botRays[i].materialXyz, N = _botRays[i].materialNormal;
		boundingBox<float> bb;
		bb.Empty_Box();
		bb.Enlarge_To_Include_Point(P._v);
		bb.Enlarge_To_Include_Point((P + N)._v);
		for (int j = i + i; j < n; ++j) {
			Vec3f M = _botRays[j].materialNormal * _botRays[j].param;
			Vec3f Q = _botRays[j].materialXyz;
			boundingBox<float> b2;
			b2.Empty_Box();
			b2.Enlarge_To_Include_Point(Q._v);
			b2.Enlarge_To_Include_Point((Q + M)._v);
			if (!b2.Intersection(bb))
				continue;
			Mat2x2f Mat;
			Mat.Initialize_With_Column_Vectors(Vec2f(N * N, N * M), Vec2f(M * N, M * M));
			Q -= P;
			Vec2f R = Mat.Robust_Solve_Linear_System(Vec2f(Q * N, Q * M));
			if (R.X <= 0.0f || R.Y <= 0.0f || R.X >= 1.0f || R.Y >= 1.0f)
				continue;
			_botRays[i].param = R.X;
			N *= R.X;
			_botRays[j].param *= R.Y;
		}
	} */


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

	std::vector<boundingBox<float> > triBox;
	boundingBox<float> bb;
	bb.Empty_Box();
	triBox.assign(_flapBottomTris.size(), bb);
	for (size_t n = _flapBottomTris.size(), i = 0; i < n; ++i) {
		long* tr = _mt->triangleVertices(_flapBottomTris[i]);
		for (int j = 0; j < 3; ++j)
			triBox[i].Enlarge_To_Include_Point(reinterpret_cast<const float(&)[3]>(*_mt->vertexCoordinate(tr[j])));
	}
	std::vector<long> topTets, bottomTets;
	topTets.reserve(_bedRays.size());
	bottomTets.reserve(_bedRays.size());
	std::vector<std::array<float, 3> > topBarys, bottomBarys, collisionNormals;
	topBarys.reserve(_bedRays.size());
	bottomBarys.reserve(_bedRays.size());
	collisionNormals.reserve(_bedRays.size());

//	tbb::parallel_for(tbb::blocked_range<size_t>(0, n),
//		[=](const tbb::blocked_range<size_t>& r) {
//			for (size_t i = r.begin(); i != r.end(); ++i)
//				Foo(a[i]);
//		}
//	);

	for (int n = _bedRays.size(), i = 0; i < n; ++i) {
		bottomRay* b = _bedRays[i];
		Vec3f P, N;
		const long* nodes = _vnt->tetNodes(_vnt->getVertexTetrahedron(b->vertex));
		const Vec3f* bw = _vnt->getVertexWeight(b->vertex);
		const Vec3f* sc[4];
		sc[0] = &_vnt->nodeSpatialCoordinate(nodes[0]);
		P = *sc[0] * (1.0f - bw->_v[0] - bw->_v[1] - bw->_v[2]);
		for (int j = 1; j < 4; ++j) {
			sc[j] = &_vnt->nodeSpatialCoordinate(nodes[j]);
			P += *sc[j] * bw->_v[j - 1];
		}
		Mat3x3f M(*sc[1] - *sc[0], *sc[2] - *sc[0], *sc[3] - *sc[0]);
		// only direction of b->materialNormal matters since deformation gradient will change length.  Normalize here.
		N = _rest[b->restIdx] * M * b->materialNormal;
		float scale = inverse_rsqrt(N * N);
		N *= scale * 0.1f;  // QISI and YUTIAN - play with fudge factor here.
		bb.Empty_Box();
		bb.Enlarge_To_Include_Point(P._v);
		bb.Enlarge_To_Include_Point((P + N)._v);
		float nearT = FLT_MAX;
		int nearV = -1;
		for (size_t n = _flapBottomTris.size(), j = 0; j < n; ++j) {
			if (bb.Intersection(triBox[j])) {
				long* tr = _mt->triangleVertices(_flapBottomTris[j]);
				Vec3f* tv[3];
				for (int k = 0; k < 3; ++k)
					tv[k] = reinterpret_cast<Vec3f*>(_mt->vertexCoordinate(tr[k]));
				Mat3x3f C(*tv[1] - *tv[0], *tv[2] - *tv[0], -N);
				Vec3f R = C.Robust_Solve_Linear_System(P - *tv[0]);
				if (R[0] < 1e-6f || R[1] < 1e-6f || R[2] < 2e-1f || R[0] + R[1] > 1.0f || R[0] > 1.0f || R[1] > 1.0f || R[2] > 1.0f)  // R[2] determines how deep the collision must go before processing triggered. Bigger makes less sticky.
					continue;
				if (nearT > R[2]) {
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
		if (nearV < 0)
			continue;
		// found soft-soft collision pair
		topTets.push_back(_vnt->getVertexTetrahedron(nearV));
		const Vec3f* W = _vnt->getVertexWeight(nearV);
		std::array<float, 3> tmp = { W->X, W->Y, W->Z };
		topBarys.push_back(tmp);
		bottomTets.push_back(_vnt->getVertexTetrahedron(b->vertex));
		tmp = { bw->X, bw->Y, bw->Z };
		bottomBarys.push_back(tmp);
		N *= nearT;
		tmp = { N.X, N.Y, N.Z };
		collisionNormals.push_back(tmp);
	}

/*	tbb::tick_count t1 = tbb::tick_count::now();
	double time = (t1 - t0).seconds();
	if (_maxTime < time)
		_maxTime = time;
	if (_minTime > time)
		_minTime = time;
	std::cout << "Soft collision minimum processing time was " << _minTime << " and maximum processing time was " << _maxTime << "verts colliding " << topTets.size() << "\n"; */

	_ptp->currentSoftCollisionPairs(topTets, topBarys, bottomTets, bottomBarys, collisionNormals);
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
