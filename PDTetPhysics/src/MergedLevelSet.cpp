#include "Utilities.h"
#include "MergedLevelSet.h"

#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>

using namespace PhysBAM;


template<class VectorType>
typename MergedLevelSet<VectorType>::T MergedLevelSet<VectorType>::Extended_Phi(const VectorType& pos) const
{
	T minDist = FLT_MAX;
	for (const auto levelSet : m_levelSet) {
		T dist = levelSet->Extended_Phi(pos);
		minDist = dist < minDist ? dist : minDist;
	}
	return minDist;
}

template<class VectorType>
VectorType PhysBAM::MergedLevelSet<VectorType>::Extended_Normal(const VectorType& pos) const
{
	// TODO: cleanup idx, mIdx
	int idx = 0;
	int mIdx = -1;
	T minDist = FLT_MAX;
	VectorType normal = pos;
	for (const auto levelSet : m_levelSet) {
		T dist = levelSet->Extended_Phi(pos);
		idx++;
		if (dist < minDist) {
			minDist = dist;
			normal = levelSet->Extended_Normal(pos);
			mIdx = idx;
		}
	}
#if 0
	if (mIdx == 0)
		std::cout << "leftEye" << std::endl;
	else if (mIdx == 1)
		std::cout << "rightEye" << std::endl;
	else
		std::cout << "teeth" << std::endl;
#endif

	return normal;
}

template<class VectorType>
PhysBAM::MergedLevelSet<VectorType>::~MergedLevelSet()
{
	for (int i = 0; i < m_levelSet.size(); i++)
		if (m_levelSet[i]) delete m_levelSet[i];
}

template<class VectorType>
void MergedLevelSet<VectorType>::addLevelSet(const std::string& collisionObjPath)
{
	std::cout << collisionObjPath << std::endl;
	using IndexType = VECTOR<int, d>;

	std::vector<std::array<int, 3>> triangles;
	std::vector<std::array<T, d>> particles;

	pdUtilities::readObj<T, d>(collisionObjPath, triangles, particles);
	// triangles here need to use 0-based indices

	const size_t nTris = triangles.size();
	const size_t nVerts = particles.size();


	TRIANGULATED_SURFACE<T>* levelset_surface;
	levelset_surface = TRIANGULATED_SURFACE<T>::Create();

	levelset_surface->particles.array_collection->Add_Elements(static_cast<int>(nVerts));
	levelset_surface->mesh.elements.Exact_Resize(static_cast<int>(nTris));
	for (int p = 0; p < nVerts; p++)
		for (int v = 0; v < d; v++)
			levelset_surface->particles.X(p + 1)(v + 1) = particles[p][v];

	for (int e = 0; e < nTris; e++) {
		for (int i = 0; i < d; i++)
			levelset_surface->mesh.elements(e + 1)(i + 1) = triangles[e][i] + 1;
	}

	levelset_surface->mesh.number_nodes = static_cast<int>(nVerts);
	levelset_surface->mesh.Initialize_Neighbor_Nodes();

	VectorType maxCorner = VectorType(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	VectorType minCorner = VectorType(FLT_MAX, FLT_MAX, FLT_MAX);
	for (int i = 1; i <= levelset_surface->particles.X.Size(); i++) {
		const VectorType& X = levelset_surface->particles.X(i);
		for (int j = 1; j <= d; j++) {
			if (X(j) > maxCorner(j))
				maxCorner(j) = X(j);
			if (X(j) < minCorner(j))
				minCorner(j) = X(j);
		}
	}

	LOG::cout << "maxCorner: " << maxCorner << std::endl;
	LOG::cout << "minCorner: " << minCorner << std::endl;

	VectorType gridOrigin = minCorner;
	IndexType gridSize = IndexType(int((maxCorner(1) - minCorner(1)) / gridDX) + 1, int((maxCorner(2) - minCorner(2)) / gridDX) + 1, int((maxCorner(3) - minCorner(3)) / gridDX) + 1);

	LOG::cout << "gridSize: " << gridSize << std::endl;
	// Compute level set

	size_t idx = m_levelSet.size();
	m_levelSet.push_back(nullptr);
	m_levelSet[idx] = LEVELSET_IMPLICIT_OBJECT<VectorType>::Create();
	
	m_levelSet[idx]->levelset.grid.Initialize(gridSize, RANGE<VectorType>(minCorner, maxCorner));
	m_levelSet[idx]->Update_Box();
	m_levelSet[idx]->Update_Minimum_Cell_Size();
	//GRID<VectorType>& levelSetGrid = *new GRID<VectorType>(gridSize, RANGE<VectorType>(minCorner, maxCorner));
	//ARRAY<T, IndexType>& phi = *new ARRAY<T, IndexType>;
	LEVELSET_MAKER_UNIFORM<T> maker;
	std::cout << "before maker" << std::endl;
	maker.Compute_Level_Set(*levelset_surface, m_levelSet[idx]->levelset.grid, m_levelSet[idx]->levelset.phi);
	std::cout << "after maker" << std::endl;
	
	delete levelset_surface;
}

template
class MergedLevelSet<VECTOR<float, 3>>;
