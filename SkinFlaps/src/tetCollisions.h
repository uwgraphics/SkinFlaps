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
	bool empty() { return _fixedCollisionSets.empty() && _botRays.empty(); }
	inline void setPdTetPhysics(pdTetPhysics *ptp) { _ptp = ptp; }
	tetCollisions() : _itCount(0), _initialized(false), _minTime((double)FLT_MAX), _maxTime(0.0){
		_fixedCollisionSets.clear(); _bedVerts.clear(); _bedVerts.reserve(1024);
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
		Vec3f materialXyz;
		Vec3f materialNormal;
		int restIdx;  // only 6 of these in bcc tets
		float param;  // parametric length variable. <0 means not unique tetrahedron so not used for now.
	};
	std::vector<bottomRay> _botRays;
	std::unordered_map<int, bottomRay> _bedVerts;
	struct fixedCollisionSet {
		std::string levelSetFilename;
		std::vector<long> vertices;
	};
	std::list< fixedCollisionSet> _fixedCollisionSets;
	std::vector<int> _topTris;

	std::vector<int> _topTets;
	std::vector<Vec3f> _topBarys;


	double _minTime, _maxTime;

	float inverse_rsqrt(float number);
};
#endif  // __TET_COLLISIONS__

