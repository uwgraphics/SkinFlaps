#include "PDTetSolver.h"
#include "Algebra.h"

#include "MergedLevelSet.h"
#include <fstream>
#include <sstream>
#include <set>

namespace {
	// parameter for bcc tet with dual latice of side length 1
const std::array<int, 9> bccDmInv{ 1, 0, 0,
								   0, 1,-1,
								  -1, 1, 1 };
const double vol = 1 / 12.;
}

template<class T, int d>
int PDTetSolver<T, d>::addSuture(const int(&tets)[2], const T(&barycentricWeights)[2][d], const T stiffness)
{
#if 0
	typename DeformerType::Suture suture;

	for (int v = 0; v < d + 1; v++) {
		suture.m_elementIndex1[v] = m_gridDeformer.m_elements[tets[0]][v];
		suture.m_elementIndex2[v] = m_gridDeformer.m_elements[tets[1]][v];
	}
	suture.m_weights1[0] = T(1);
	suture.m_weights2[0] = T(1);
	for (int v = 0; v < d; v++) {
		suture.m_weights1[0] -= barycentricWeights[0][v];
		suture.m_weights1[v + 1] = barycentricWeights[0][v];
		suture.m_weights2[0] -= barycentricWeights[1][v];
		suture.m_weights2[v + 1] = barycentricWeights[1][v];
	}
	suture.m_stiffness = stiffness;
	suture.m_restLength = 0;
	m_gridDeformer.m_sutures.push_back(suture);
	return (int)m_gridDeformer.m_sutures.size();
#endif

	typename DeformerType::Constraint c1{}, c2{};

	for (int v = 0; v < d + 1; v++) {
		c1.m_elementIndex[v] = m_gridDeformer.m_elements[tets[0]][v];
		c2.m_elementIndex[v] = m_gridDeformer.m_elements[tets[1]][v];
	}
	c1.m_weights[0] = T(1);
	c2.m_weights[0] = T(1);
	for (int v = 0; v < d; v++) {
		c1.m_weights[0] -= barycentricWeights[0][v];
		c1.m_weights[v + 1] = barycentricWeights[0][v];
		c2.m_weights[0] -= barycentricWeights[1][v];
		c2.m_weights[v + 1] = barycentricWeights[1][v];
	}
	c1.m_stiffness = 2 * stiffness;
	c2.m_stiffness = 2 * stiffness;
	m_gridDeformer.m_fakeSutures.push_back(c1);
	m_gridDeformer.m_fakeSutures.push_back(c2);
	return (int) (m_gridDeformer.m_sutures.size() - 1 + m_gridDeformer.m_fakeSutures.size()/2);

}

template<class T, int d>
void PDTetSolver<T, d>::initializeSolver()
{
	using IteratorType = typename DeformerType::IteratorType;
	m_gridDeformer.deallocateAuxiliaryStructures();
	m_gridDeformer.initializeCollisionElements();
	m_gridDeformer.initializeAuxiliaryStructures();
	if (m_gridDeformer.m_collisionConstraints.size()||m_gridDeformer.m_collisionSutures.size()) {
		hasCollision = true;
#ifdef USE_CUDA
		m_solver_c.releaseCuda();
#endif
		m_solver_c.releasePardiso();
		m_solver_c.deallocate();

		m_solver_c.initialize(m_gridDeformer.m_nodeType); // initialzie
		m_solver_c.computeTensor(m_gridDeformer.m_elements, m_gridDeformer.m_gradientMatrix, m_gridDeformer.m_elementRestVolume, m_gridDeformer.m_muLow, m_gridDeformer.m_muHigh, m_gridDeformer.m_sutures); // computeTensor
#ifdef USE_CUDA
		m_solver_c.computeE2Tensor(m_gridDeformer.m_elements, m_gridDeformer.m_elementFlags, m_gridDeformer.m_gradientMatrix, m_gridDeformer.m_elementRestVolume, m_gridDeformer.m_muHigh[0] * (1 + m_weightProportion * m_weightProportion)); // computeE2Tensor
#endif
		m_solver_c.initializePardiso(m_gridDeformer.m_constraints, m_gridDeformer.m_sutures, m_gridDeformer.m_fakeSutures); // init pardiso
#ifdef USE_CUDA
		m_solver_c.initializeCuda(m_gridDeformer.m_collisionConstraints, m_gridDeformer.m_collisionSutures); // init Cuda
		std::cout << "using CudaSolver with nInner = " << m_nInner << std::endl;
#else
		std::cout << "using MKLSolver"<< std::endl;
#endif
	}
	else {
		hasCollision = false;
		m_solver_d.releasePardiso();
		m_solver_d.deallocate();

		m_solver_d.initialize(m_gridDeformer.m_nodeType);
		m_solver_d.computeTensor(m_gridDeformer.m_elements, m_gridDeformer.m_gradientMatrix, m_gridDeformer.m_elementRestVolume, m_gridDeformer.m_muHigh[0] * (1 + m_weightProportion * m_weightProportion), m_gridDeformer.m_sutures);
		m_solver_d.initializePardiso(m_gridDeformer.m_constraints, m_gridDeformer.m_sutures, m_gridDeformer.m_fakeSutures);
		std::cout << "using DirectSolver" << std::endl;
	}
}

template<class T, int d>
void PDTetSolver<T, d>::reInitializeSolver()
{
	if (hasCollision) {
		m_solver_c.reInitializePardiso(m_gridDeformer.m_constraints, m_gridDeformer.m_sutures, m_gridDeformer.m_fakeSutures);
#ifdef USE_CUDA
		m_solver_c.reInitializeCuda(m_gridDeformer.m_collisionConstraints, m_gridDeformer.m_collisionSutures);
#endif
	}
	else {
		m_solver_d.reInitializePardiso(m_gridDeformer.m_constraints, m_gridDeformer.m_sutures, m_gridDeformer.m_fakeSutures);
	}
}

template<class T, int d>
void PDTetSolver<T, d>::addCollisionProxies(const int * tets, const T (*weights)[d], size_t length)
{
	//m_gridDeformer.m_collisionConstraints.resize(length);
	// int offset = m_gridDeformer.m_collisionConstraints.size();
	for (int i = 0; i < length; i++) {
		typename DeformerType::Constraint constraint{};
		constraint.m_weights[0] = T(1);
		for (int v = 0; v < d + 1; v++) {
			constraint.m_elementIndex[v] = m_gridDeformer.m_elements[tets[i]][v];
			m_gridDeformer.m_nodeType[constraint.m_elementIndex[v]] = NodeType::Collision;
		}
		for (int v = 0; v < d; v++) {
			constraint.m_weights[0] -= weights[i][v];
			constraint.m_weights[v + 1] = weights[i][v];
		}
		constraint.m_stiffness = 0;
		// set target to origin, this doesn't matter since stiffness is 0
		constraint.m_xT = VectorType();
		m_gridDeformer.m_collisionConstraints.push_back(constraint);
	}
	// IteratorType iterator(m_gridDeformer.m_X);

}

template<class T, int d>
void PDTetSolver<T, d>::addSelfCollisionElements(const int* tets, size_t length)
{
	using VectorType = PhysBAM::VECTOR<T, d>;
	// temporary implementation without any assumption on the collision tets
	m_gridDeformer.m_collisionSutures.resize(length);
	for (int i = 0; i < length; i++) {
		auto& constraint = m_gridDeformer.m_collisionSutures[i];
		constraint.m_weights1 = { T(0) };
		constraint.m_weights2 = { T(0) };

		for (int v = 0; v < d + 1; v++) {
			constraint.m_elementIndex1[v] = m_gridDeformer.m_elements[tets[i]][v];
			constraint.m_elementIndex2[v] = m_gridDeformer.m_elements[tets[i]][v];
			m_gridDeformer.m_nodeType[constraint.m_elementIndex1[v]] = NodeType::Collision;
		}
		
		constraint.m_stiffness = 0;

		constraint.m_normal = VectorType::All_Ones_Vector().Normalized();
	}
}

template<class T, int d>
void PDTetSolver<T, d>::solve()
{
	using StateVariableType = typename DiscretizationType::StateVariableType;
	using IteratorType = typename DeformerType::IteratorType;
	using AlgebraType = PhysBAM::Algebra<StateVariableType>;

	StateVariableType delta_X{};
	StateVariableType f{};
	
	IteratorType iterator(m_gridDeformer.m_X);
	iterator.resize(delta_X);
	iterator.resize(f);
	

	m_gridDeformer.updatePositionBasedState(ElementFlag::unCollisionEl/*, m_rangeMin, m_rangeMax*/ ); // updateR1
	m_gridDeformer.addElasticForce(f, ElementFlag::unCollisionEl /*, m_rangeMin, m_rangeMax, m_weightProportion */); //addR1Force
	m_gridDeformer.addConstraintForce(f); //addConstraintForec

	if (hasCollision) {
#ifdef USE_CUDA
		StateVariableType u{};
		iterator.resize(u);

		for (int v = 0; v < d; v++) {
			m_solver_c.copyIn(f, v); //updateR1
			m_solver_c.forwardSubstitution(); //forwardSubstitution
			m_solver_c.setTemp(v); //setTemp
			m_solver_c.copyOut(f, v); //copyOut

			// exit(1);
		}

		// inner loop
		for (int inner_i = 0; inner_i < m_nInner; inner_i++) {
			//if (frame != 1)
			updateCollisionConstraints();     // updateCollision

		//m_solver_c.updatePardiso(deformer.m_collisionConstraints, deformer.m_collisionSutures);
			m_solver_c.updateCuda(m_gridDeformer.m_collisionConstraints, m_gridDeformer.m_collisionSutures);

			StateVariableType f_temp{};
			iterator.resize(f_temp);

			for (IteratorType iterator(f); !iterator.isEnd(); iterator.next())
				if (iterator.value(m_gridDeformer.m_nodeType) != NodeType::Collision) {
					iterator.value(f_temp) = iterator.value(f);
				}

			m_gridDeformer.updatePositionBasedState(ElementFlag::CollisionEl/*, m_rangeMin, m_rangeMax*/); // updateR2

			m_gridDeformer.addElasticForce(f_temp, ElementFlag::CollisionEl/*, m_rangeMin, m_rangeMax, m_weightProportion*/); // addR2Force

			m_gridDeformer.addCollisionForce(f_temp);     // addCollisionForce

			for (int v = 0; v < d; v++) {
				m_solver_c.copyIn(f_temp, v); // copyIn
				m_solver_c.updateForce(v); // updateForce
				m_solver_c.copyOut(f, v);//copyOut
			}

			//m_boxTest.clearDirichlet(m_boxTest.m_geometry, deformer.m_nodeType, f);
			for (int v = 0; v < d; v++) {
				m_solver_c.copyIn(f, v); //copyIn
				m_solver_c.diagSolve(); //diagSolve
				m_solver_c.updateTemp(v); // updateTemp
				m_solver_c.copyOut(delta_X, v);//copyOutTime
			}

			// update x2
			for (IteratorType iterator(delta_X); !iterator.isEnd(); iterator.next())
				if (iterator.value(m_gridDeformer.m_nodeType) == NodeType::Collision)
					iterator.value(m_gridDeformer.m_X) += iterator.value(delta_X);

			// accum to u
			for (IteratorType iterator(u); !iterator.isEnd(); iterator.next())
				if (iterator.value(m_gridDeformer.m_nodeType) == NodeType::Collision)
					iterator.value(u) += iterator.value(delta_X);

			for (IteratorType iterator(delta_X); !iterator.isEnd(); iterator.next())
				if (iterator.value(m_gridDeformer.m_nodeType) == NodeType::Inactive)
					iterator.value(delta_X) = VectorType();

		}
		// copy in x1 part
		for (IteratorType iterator(u); !iterator.isEnd(); iterator.next())
			if (iterator.value(m_gridDeformer.m_nodeType) != NodeType::Collision)
				iterator.value(u) = iterator.value(delta_X);

		for (int v = 0; v < d; v++) {
			m_solver_c.copyIn(u, v); // copyIn
			m_solver_c.backwardSubstitution(); //backwardSubstitution
			m_solver_c.copyOut(delta_X, v); // copyOut
		}

		// update x1
		for (IteratorType iterator(delta_X); !iterator.isEnd(); iterator.next())
			if (iterator.value(m_gridDeformer.m_nodeType) != NodeType::Collision)
				iterator.value(m_gridDeformer.m_X) += iterator.value(delta_X);
#else
		StateVariableType u{};
		iterator.resize(u);

		updateCollisionConstraints();     // updateCollision
		m_solver_c.updatePardiso(m_gridDeformer.m_collisionConstraints, m_gridDeformer.m_collisionSutures);
			m_gridDeformer.updatePositionBasedState(ElementFlag::CollisionEl /*, m_rangeMin, m_rangeMax*/); // updateR2
			m_gridDeformer.addElasticForce(f, ElementFlag::CollisionEl /*, m_rangeMin, m_rangeMax, m_weightProportion */ ); // addR2Force
			m_gridDeformer.addCollisionForce(f);     // addCollisionForce

			for (int v = 0; v < d; v++) {
				m_solver_c.copyIn(f, v); //copyIn
				m_solver_c.solve(); //diagSolve
				m_solver_c.copyOut(delta_X, v);//copyOutTime
			}

			for (IteratorType i(delta_X); !i.isEnd(); i.next())
				if (i.value(m_gridDeformer.m_nodeType) == NodeType::Inactive)
					i.value(delta_X) = VectorType();

		// update x1
		for (IteratorType i(delta_X); !i.isEnd(); i.next())
				i.value(m_gridDeformer.m_X) += i.value(delta_X);
#endif
	}
	else {
		//m_boxTest.clearDirichlet(m_boxTest.m_geometry, deformer.m_nodeType, f);

		for (int v = 0; v < d; v++) {
			m_solver_d.copyIn(f, v);
			m_solver_d.solve();
			m_solver_d.copyOut(delta_X, v);
		}
		AlgebraType::addTo(m_gridDeformer.m_X, delta_X);
	}

	for (IteratorType i(delta_X); !i.isEnd(); i.next())
		if (i.value(m_gridDeformer.m_nodeType) == NodeType::Inactive)
			i.value(delta_X) = VectorType();
}

template<class T, int d>
PDTetSolver<T, d>::~PDTetSolver()
{
	// TODO: distruct everthing here
}


template<class T, int d>
void PDTetSolver<T, d>::premoteSutures()
{
	for (int i = 0; i < m_gridDeformer.m_fakeSutures.size(); i += 2) {
		typename DeformerType::Suture suture{};
		const DeformerType::Constraint& c1 = m_gridDeformer.m_fakeSutures[i];
		const DeformerType::Constraint& c2 = m_gridDeformer.m_fakeSutures[i+1];
		suture.m_elementIndex1 = c1.m_elementIndex;
		suture.m_elementIndex2 = c2.m_elementIndex;
		suture.m_weights1 = c1.m_weights;
		suture.m_weights2 = c2.m_weights;
		suture.m_stiffness = c1.m_stiffness / 2;
		m_gridDeformer.m_sutures.push_back(suture);
	}
	m_gridDeformer.m_fakeSutures.clear();
#if 0
	dumper::writeElements(m_gridDeformer.m_elements);
	dumper::writePositions(m_gridDeformer.m_X);
	dumper::writeCollisionConstraints(m_gridDeformer.m_collisionConstraints);
	dumper::writeConstraints(m_gridDeformer.m_constraints);
	dumper::writeSutures(m_gridDeformer.m_sutures);
	dumper::writeNodeTypes<d>(m_gridDeformer.m_nodeType);

	{
		std::ofstream myout("DmInverse.txt");
		for (int v = 0; v < d; v++)
			for (int w = 0; w < d; w++)
				myout << m_gridDeformer.m_gradientMatrix[0](v + 1, w + 1)<<" ";
		myout.close();

		myout.open("Parameters.txt");
		myout<<m_collisionStiffness<<" collisionStiffness"<<std::endl;
		myout<<m_weightProportion<<" weightProportion"<<std::endl;
		myout << m_rangeMin << " rangeMin" << std::endl;
		myout << m_rangeMax << " rangeMax" << std::endl;
		myout << m_gridDeformer.m_uniformMu << " mu" << std::endl;
		myout.close();
	}
#endif
}

template<class T, int d>
void PDTetSolver<T, d>::updateCollisionConstraints()
{
	if (m_levelSet) {
		T threshold = (T)1e-5;
		for (int c = 0; c < m_gridDeformer.m_collisionConstraints.size(); c++) {
			auto &constraint = m_gridDeformer.m_collisionConstraints[c];
			VectorType pos = DiscretizationType::interpolateX(constraint.m_elementIndex, constraint.m_weights, m_gridDeformer.m_X);
			T phi = m_levelSet->Extended_Phi(pos);
			if (phi < -threshold) {
				// std::cout << "phi " << phi << std::endl;
				constraint.m_xT = pos - m_levelSet->Extended_Normal(pos)*phi;
				constraint.m_stiffness = m_collisionStiffness;
			}
			else {
				constraint.m_xT = pos;
				constraint.m_stiffness = 0;
			}
		}
	}

}

template<class T, int d>
void PDTetSolver<T, d>::updateCollisionSutures(const int length, const int* topI, const int* botI, const T* topW, const T* botW, const T* normal)
{
	T threshold = T(1e-6);
	for (auto& c : m_gridDeformer.m_collisionSutures)
		c.m_stiffness = 0;

	using Reshaped = const T (*)[d];
	const auto reshapedT = reinterpret_cast<Reshaped>(topW);
	const auto reshapedB = reinterpret_cast<Reshaped>(botW);
	const auto reshapedN = reinterpret_cast<Reshaped>(normal);


	// std::cout << "=====================================" << std::endl;
	for (int i = 0; i < length; i++) {
		auto& constraint = m_gridDeformer.m_collisionSutures[i];
		constraint.m_elementIndex1 = m_gridDeformer.m_elements[botI[i]];
		constraint.m_elementIndex2 = m_gridDeformer.m_elements[topI[i]];
		constraint.m_weights1[0] = T(1);
		constraint.m_weights2[0] = T(1);
		
		for (int v = 0; v < d; v++) {
			constraint.m_weights1[v + 1] = reshapedB[i][v];
			constraint.m_weights2[v + 1] = reshapedT[i][v];

			constraint.m_weights1[0] -= reshapedB[i][v];
			constraint.m_weights2[0] -= reshapedT[i][v];

			constraint.m_normal(v+1) = reshapedN[i][v];

			VectorType pos1 = DiscretizationType::interpolateX(constraint.m_elementIndex1, constraint.m_weights1, m_gridDeformer.m_X);
			VectorType pos2 = DiscretizationType::interpolateX(constraint.m_elementIndex2, constraint.m_weights2, m_gridDeformer.m_X);

			std::set<int> ind_set;
			bool dup = false;
			for (const auto j : constraint.m_elementIndex1)
				ind_set.insert(j);
			for (const auto j : constraint.m_elementIndex2)
				if (ind_set.find(j) != ind_set.end()) {
					dup = true;
					break;
				}


			if (VectorType::Dot_Product(pos2-pos1, constraint.m_normal.Normalized())<threshold /*&& VectorType::Dot_Product(pos2-pos1, constraint.m_normal.Normalized())>-.1*/ && !dup) {
				// std::cout << constraint.m_elementIndex1[0] << " " << constraint.m_elementIndex1[1] << " " << constraint.m_elementIndex1[2] << " " << constraint.m_elementIndex1[3] << "           " <<
				//	constraint.m_elementIndex2[0] << " " << constraint.m_elementIndex2[1] << " " << constraint.m_elementIndex2[2] << " " << constraint.m_elementIndex2[3] << std::endl;
				constraint.m_stiffness = m_selfCollisionStiffness;
			}
		}

	}
	// std::cout << "=====================================" << std::endl;


}



template<class T, int d>
inline void PDTetSolver<T, d>::initializeDeformer(const int(*elements)[d + 1], const T(*x)[d], const size_t nEls, const size_t nNodes)
{
	using namespace PhysBAM;
	m_gridDeformer.deallocateAuxiliaryStructures();

	m_gridDeformer.m_constraints.clear();
	m_gridDeformer.m_collisionConstraints.clear();
	m_gridDeformer.m_sutures.clear();
	m_gridDeformer.m_collisionSutures.clear();


	m_gridDeformer.m_X.resize(nNodes);
	for (int i = 0; i < nNodes; i++)
		for (int v = 0; v < d; v++)
			m_gridDeformer.m_X[i](v + 1) = x[i][v];
	
	m_gridDeformer.m_elements.resize(nEls);
	m_gridDeformer.m_muLow.resize(nEls, m_weightProportion * m_weightProportion * m_uniformMu);
	m_gridDeformer.m_muHigh.resize(nEls, m_uniformMu);
	m_gridDeformer.m_rangeMin.resize(nEls, m_rangeMin);
	m_gridDeformer.m_rangeMax.resize(nEls, m_rangeMax);

	for (int i = 0; i < nEls; i++)
		for (int v = 0; v < d + 1; v++)
			m_gridDeformer.m_elements[i][v] = elements[i][v];

	// assuming all nodes in x are active
	m_gridDeformer.m_nodeType.resize(nNodes);
	for (int i = 0; i < nNodes; i++)
		m_gridDeformer.m_nodeType[i] = NodeType::Active;

	m_gridDeformer.initializeDeformer();
	m_gridDeformer.initializeUndeformedState();
	 m_gridDeformer.initializeAuxiliaryStructures();

	std::cout << "Number of Elements = " << m_gridDeformer.m_elements.size() << std::endl;
	std::cout << "Number of Particles = " << m_gridDeformer.m_X.size() << std::endl;
}

template<class T, int d>
void PDTetSolver<T, d>::initializeDeformer(const int(*elements)[d + 1], const size_t nEls, const T gridSize)
{
	static_assert(d == 3, "this operation only support 3 dimension");
	using namespace PhysBAM;
	// m_gridDeformer.deallocate();
	m_gridDeformer.m_constraints.clear();
	m_gridDeformer.m_collisionConstraints.clear();
	m_gridDeformer.m_sutures.clear();
	m_gridDeformer.m_fakeSutures.clear();
	m_gridDeformer.m_collisionSutures.clear();


	int nNodes = 0;
	m_gridDeformer.m_elements.resize(nEls);
	m_gridDeformer.m_muLow.resize(nEls, m_weightProportion * m_weightProportion * m_uniformMu);
	m_gridDeformer.m_muHigh.resize(nEls, m_uniformMu);
	m_gridDeformer.m_rangeMin.resize(nEls, m_rangeMin);
	m_gridDeformer.m_rangeMax.resize(nEls, m_rangeMax);
	for (int i = 0; i < nEls; i++)
		for (int v = 0; v < d + 1; v++) {
			m_gridDeformer.m_elements[i][v] = elements[i][v];
			if (elements[i][v] > nNodes)
				nNodes = elements[i][v];
		}
	nNodes++;
	m_gridDeformer.m_X.resize(nNodes);
	m_gridDeformer.m_nodeType.resize(nNodes);
#if 0
	for (int i = 0; i < nEls; i++)
		for (int v = 0; v < d + 1; v++)
			m_gridDeformer.m_nodeType[elements[i][v]] = NodeType::Active;
#else
	for (int i = 0; i < nNodes; i++)
			m_gridDeformer.m_nodeType[i] = NodeType::Active;
#endif

#if 0
	// make sure all nodes in x are active
	for (int i = 0; i < nNodes; i++)
		if (m_gridDeformer.m_nodeType[i] != NodeType::Active || m_gridDeformer.m_nodeType[i] != NodeType::Collision)
			throw std::logic_error("node not active");
		// assert(m_gridDeformer.m_nodeType[i] == NodeType::Active);
#endif


	m_gridDeformer.initializeDeformer(); 
	{
		m_gridDeformer.m_gradientMatrix.resize(nEls);
		m_gridDeformer.m_elementRestVolume.resize(nEls);
		for (int i = 0; i < nEls; i++) {
			for (int v = 0; v < d; v++)
				for (int w = 0; w < d; w++)
					m_gridDeformer.m_gradientMatrix[i](v + 1, w + 1) = ::bccDmInv[(size_t)w * d + v] / T(gridSize);
			m_gridDeformer.m_elementRestVolume[i] = T(::vol) * gridSize * gridSize * gridSize;
		}
	}

}

template<class T, int d>
void PDTetSolver<T, d>::initializeDeformer_multires(const int(*elements)[d + 1], const uint8_t *tetSizeMultipliers, const size_t nEls, const T gridSize)  // COURT - new version for multires tets
{
	static_assert(d == 3, "this operation only support 3 dimension");
	using namespace PhysBAM;
	// m_gridDeformer.deallocate();
	m_gridDeformer.m_constraints.clear();
	m_gridDeformer.m_collisionConstraints.clear();
	m_gridDeformer.m_sutures.clear();
	m_gridDeformer.m_fakeSutures.clear();
	m_gridDeformer.m_collisionSutures.clear();


	int nNodes = 0;
	m_gridDeformer.m_elements.resize(nEls);
	m_gridDeformer.m_muLow.resize(nEls, m_weightProportion * m_weightProportion * m_uniformMu);
	m_gridDeformer.m_muHigh.resize(nEls, m_uniformMu);
	m_gridDeformer.m_rangeMin.resize(nEls, m_rangeMin);
	m_gridDeformer.m_rangeMax.resize(nEls, m_rangeMax);
	for (int i = 0; i < nEls; i++)
		for (int v = 0; v < d + 1; v++) {
			m_gridDeformer.m_elements[i][v] = elements[i][v];
			if (elements[i][v] > nNodes)
				nNodes = elements[i][v];
		}
	nNodes++;
	m_gridDeformer.m_X.resize(nNodes);
	m_gridDeformer.m_nodeType.resize(nNodes);
#if 0
	for (int i = 0; i < nEls; i++)
		for (int v = 0; v < d + 1; v++)
			m_gridDeformer.m_nodeType[elements[i][v]] = NodeType::Active;
#else
	for (int i = 0; i < nNodes; i++)
		m_gridDeformer.m_nodeType[i] = NodeType::Active;
#endif

#if 0
	// make sure all nodes in x are active
	for (int i = 0; i < nNodes; i++)
		if (m_gridDeformer.m_nodeType[i] != NodeType::Active || m_gridDeformer.m_nodeType[i] != NodeType::Collision)
			throw std::logic_error("node not active");
	// assert(m_gridDeformer.m_nodeType[i] == NodeType::Active);
#endif


	m_gridDeformer.initializeDeformer();
	{
		m_gridDeformer.m_gradientMatrix.resize(nEls);
		m_gridDeformer.m_elementRestVolume.resize(nEls);
		for (int i = 0; i < nEls; i++) {
			T sizeMult(tetSizeMultipliers[i]);
			for (int v = 0; v < d; v++)
				for (int w = 0; w < d; w++)
					m_gridDeformer.m_gradientMatrix[i](v + 1, w + 1) = ::bccDmInv[(size_t)w * d + v] / (T(gridSize) * sizeMult);
			m_gridDeformer.m_elementRestVolume[i] = T(::vol) * gridSize * gridSize * gridSize;  // *sizeMult* sizeMult* sizeMult;  // COURT same bug as Bouaziz?
		}
	}
}


template<class T, int d>
int PDTetSolver<T, d>::addInterNodeConstraint(const int microNode, int nMacros, const int *macroNodes, const T *macroWeights, const T stiffness) // COURT added 
{
	typename DeformerType::InternodeConstraint inC{};
	inC.m_microNodeNumber = microNode;
	inC.m_xT = { 0.0, 0.0, 0.0 };  // COURT unnecessary?
	inC.m_stiffness = stiffness;
	inC.m_macroNodes.assign(macroNodes, macroNodes + nMacros);
	inC.m_macroWeights.assign(macroWeights, macroWeights + nMacros);
	m_gridDeformer.m_InternodeConstraints.push_back(inC);
	return (int)m_gridDeformer.m_InternodeConstraints.size() - 1;
}


template<class T, int d>
int PDTetSolver<T, d>::addConstraint(const long tet, const T(&barycentricWeight)[d], const T(&hookPosition)[d], const T stiffness, T limit)
{
	return addConstraint(reinterpret_cast<const int(&)[d+1]>(m_gridDeformer.m_elements[tet]), barycentricWeight, hookPosition, stiffness, limit);
}

template<class T, int d>
int PDTetSolver<T, d>::addConstraint(const int(&index)[d + 1], const T(&barycentricWeight)[d], const T(&hookPosition)[d], const T stiffness, T limit)
{
	typename DeformerType::Constraint constraint{};
	constraint.m_weights[0] = T(1);
	for (int v = 0; v < d + 1; v++)
		constraint.m_elementIndex[v] = index[v];

	for (int v = 0; v < d; v++) {
		constraint.m_weights[0] -= barycentricWeight[v];
		constraint.m_weights[v + 1] = barycentricWeight[v];
	}
	constraint.m_stiffness = stiffness;
	constraint.m_stressLimit = limit;
	for (int v = 0; v < d; v++)
		constraint.m_xT(v + 1) = hookPosition[v];

	m_gridDeformer.m_constraints.push_back(constraint);

	if (index[0] > m_gridDeformer.m_X.size() || index[0] < 0)
		throw std::logic_error("index out of range");

	return (int)m_gridDeformer.m_constraints.size() - 1;
}

// template instantiation
template
class PDTetSolver<float, 3>;
