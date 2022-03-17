#pragma once

#include "GridDeformerTet.h"
#ifdef USE_CUDA
#include "CudaSolver.h"
#endif
#include "SchurSolver.h"
#include "MergedLevelSet.h"

namespace PhysBAM {
	template<class VectorType>
	class MergedLevelSet;
}
template <class T, int d>
class PDTetSolver
{
private:
	using DeformerType = PhysBAM::GridDeformerTet<std::vector<PhysBAM::VECTOR<T, d>>>;
	using VectorType = typename DeformerType::VectorType;
	using DiscretizationType = typename DeformerType::DiscretizationType;
	using IntType = int;

	T m_collisionStiffness;
	T m_selfCollisionStiffness;
	int m_nInner;
	T m_weightProportion;
	T m_rangeMin;
	T m_rangeMax;
	T m_uniformMu = 0;

	// int m_nCollisionNodes;

	DeformerType m_gridDeformer;
	PhysBAM::SchurSolver<DiscretizationType, IntType> m_solver_d; // use this when no collision node
#ifdef USE_CUDA
	PhysBAM::CudaSolver<DiscretizationType, IntType> m_solver_c; // use this when there are collision nodes
#else
	PhysBAM::SchurSolver<DiscretizationType, IntType> m_solver_c;
#endif
	// PhysBAM::LEVELSET_IMPLICIT_OBJECT<VectorType>* m_softLevelSet;
	PhysBAM::MergedLevelSet<VectorType>* m_levelSet;
	std::vector<std::string> m_levelSetPaths;

	bool hasCollision = false;

public:

	inline T* getPositionPtr() {
		return &m_gridDeformer.m_X[0](1);
	}

	inline const T* getPositionPtr() const {
		return &m_gridDeformer.m_X[0](1);
	}

	inline void addSubset(const float uniformMu, const float weightProportion, const float strainMin, const float strainMax, const std::vector<int>& tets) {
		for (const int i : tets) {
			m_gridDeformer.m_muLow[i] = uniformMu*weightProportion*weightProportion;
			m_gridDeformer.m_muHigh[i] = uniformMu;
			m_gridDeformer.m_rangeMin[i] = strainMin;
			m_gridDeformer.m_rangeMax[i] = strainMax;
		}
	}

	inline void setParameters(const T mu, const T weightProportion, const T rangeMin, const T rangeMax, const T collisionStiffness, const T selfCollisionStiffness) {
		using namespace PhysBAM;
		// m_gridDeformer.m_uniformMu = mu;
		m_uniformMu = mu;
		m_weightProportion = weightProportion;
		m_rangeMin = rangeMin;
		m_rangeMax = rangeMax;
		m_collisionStiffness = collisionStiffness;
		m_selfCollisionStiffness = selfCollisionStiffness;


		std::cout << "    m_weightProportion = " << m_weightProportion << std::endl;
		std::cout << "    m_rangeMin         = " << m_rangeMin << std::endl;
		std::cout << "    m_rangeMax         = " << m_rangeMax << std::endl;
		std::cout << "    m_collisionStiff   = " << m_collisionStiffness << std::endl;
	}
	void initializeDeformer(const long (*elements)[d+1], const T (*x)[d], const size_t nEls, const size_t nNodes);

	void initializeDeformer(const long(*elements)[d + 1], const size_t nEls, const T tetSize);
	
	int addConstraint(const long tet, const T(&barycentricWeight)[d], const T(&hookPosition)[d], const T stiffness, const T limit = std::numeric_limits<T>::max());  // returns constraint index

	int addConstraint(const int (&index)[d+1], const T(&barycentricWeight)[d], const T(&hookPosition)[d], const T stiffness, const T limit = std::numeric_limits<T>::max());  // returns constraint index

	inline void moveConstraint(const int hookHandle, const T(&newPosition)[d]) {
		for (int v = 0; v < d; v++)
			m_gridDeformer.m_constraints[hookHandle].m_xT(v + 1) = newPosition[v];
	}

	inline void deleteConstraint(const int hookHandle) {
		m_gridDeformer.m_constraints[hookHandle].m_stiffness = 0;
	}

	int addSuture(const long (&tets)[2], const T (&barycentricWeights)[2][d], const T stiffness);  // returns constraint index

	void deleteSuture(const int sutureHandle) {
		if (sutureHandle < m_gridDeformer.m_sutures.size())
		{
			m_gridDeformer.m_sutures[sutureHandle].m_stiffness = 0;
		}
		else {
			const size_t fakeSutureHandle = sutureHandle - m_gridDeformer.m_sutures.size();
			m_gridDeformer.m_fakeSutures[fakeSutureHandle * 2].m_stiffness = 0;
			m_gridDeformer.m_fakeSutures[fakeSutureHandle * 2 + 1].m_stiffness = 0;
		}
	}

	void initializeSolver();  // After constraints have changed computes ATA and does its LDLT()

	void reInitializeSolver();  

	void addCollisionProxies(const long *tets, const T (*weights)[d], size_t length);
	void addSelfCollisionElements(const long* tets, size_t length);

	inline void addLevelSet(const std::string s) { m_levelSetPaths.push_back(s); }

	inline void initializeLevelSet(const T dx) { m_levelSet->initializeLevelSet(m_levelSetPaths, dx); }
	// void initializeLevelSet(const int(*triangles)[d], const T(*vertices)[d], const size_t nTris, const size_t nVerts);

	void solve();  // do least squares solve and process collisions

	PDTetSolver() : m_nInner(1), m_rangeMin(1), m_rangeMax(1), m_weightProportion(0), m_collisionStiffness(0), m_selfCollisionStiffness(0) { m_levelSet = new PhysBAM::MergedLevelSet<VectorType>; }
	~PDTetSolver();

	void premoteSutures();

	inline const std::array<int, 4>& getTetIndices(int tet) { return m_gridDeformer.m_elements[tet]; }  // COURT added

private:
	void updateCollisionConstraints();
public:
	void updateCollisionSutures(const long length, const long* topI, const long* botI, const T* topW, const T* botW, const T* normal); // this should be private and handled by PDSolver it self in future iterations

	inline void releaseSolver() {
		if (m_gridDeformer.m_collisionConstraints.size()) {
#ifdef USE_CUDA
			m_solver_c.releaseCuda();
#endif
			m_solver_c.releasePardiso();
			m_solver_c.deallocate();
		}
		else {
			// TODO: garbage collection for schurSolver
		}
	}
	
	inline void releaseDeformer() {
		m_gridDeformer.deallocateAuxiliaryStructures();
	}
};

