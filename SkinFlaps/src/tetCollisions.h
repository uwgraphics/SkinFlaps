#ifndef __TET_COLLISIONS__
#define __TET_COLLISIONS__

#include <vector>
#include <array>
#include "Vec2f.h"
#include "Vec3f.h"
#include "Mat3x3f.h"

// forward declarations
class materialTriangles;
class vnBccTetrahedra;
class pdTetPhysics;

class tetCollisions
{
public:
	void initSoftCollisions(materialTriangles *mt, vnBccTetrahedra *vnt);  // call after every topo change
	void findSoftCollisionPairs();  // call every physics iteration
//	void addSphereSet(const float radius, const float (&center)[3], const std::vector<Vec2f> &txPoly );
	void addFixedCollisionSet(materialTriangles* mt, const std::string& levelSetFile, const std::vector<Vec2f>& txPoly);
	void updateFixedCollisions(materialTriangles *mt, vnBccTetrahedra *vnt);  // must be done after every topo change
	bool empty() { return _fixedCollisionSets.empty() && _bedRays.empty(); }
	inline void setPdTetPhysics(pdTetPhysics *ptp) { _ptp = ptp; }
	tetCollisions() : _itCount(0), _initialized(false), _minTime((double)FLT_MAX), _maxTime(0.0){
		_fixedCollisionSets.clear(); _flapBottomTris.clear();  // _bedVerts.clear(); _bedVerts.reserve(1024); 
	}
	~tetCollisions() {}

private:
	int _itCount;
	static materialTriangles *_mt;
	static vnBccTetrahedra *_vnt;
	static pdTetPhysics *_ptp;
	bool _initialized;
	Mat3x3f _rest[6];  // material inverses used to compute deformation gradients
	struct bottomRay {
		int vertex;
		Vec3f P;
		Vec3f N;
		Vec3f materialNormal;
		int restIdx;  // only 6 of these in bcc tets
	};
	std::vector<bottomRay> _bedRays;
//	std::unordered_map<int, bottomRay> _bedVerts;
	std::vector<int> _flapBottomTris;
	std::vector<int> _topTets;
	std::vector<Vec3f> _topBarys;
	struct fixedCollisionSet {
		std::string levelSetFilename;
		std::vector<long> vertices;
	};
	std::list< fixedCollisionSet> _fixedCollisionSets;

	long parametricMTtriangleTet(const int mtTriangle, const float(&uv)[2], Vec3f& gridLocus, bccTetCentroid& tC);

	double _minTime, _maxTime;

	float inverse_rsqrt(float number);
};
#endif  // __TET_COLLISIONS__

