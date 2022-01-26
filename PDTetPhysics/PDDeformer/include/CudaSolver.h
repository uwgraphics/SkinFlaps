//#####################################################################
// Copyright (c) 2019, Eftychios Sifakis, Yutian Tao, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#pragma once

#include <mkl.h>
#include <array>

#include "MKLWrapper.h"
#include "PardisoWrapper.h"
#include "SimulationFlags.h"
#include "PDConstraints.h"
#include "CudaWrapper.h"

#include "Discretization.h"

// for debug
#include "dumper.h"



namespace PhysBAM {

    template <class Discretization, class IntType>
        struct CudaSolver {
            using DiscretizationType = Discretization;
            using StateVariableType = typename DiscretizationType::StateVariableType;
            using IteratorType = Iterator<StateVariableType>;

            using IndexType = typename IteratorType::IndexType;
            using VectorType = typename IteratorType::DataType;
            using T = typename VectorType::ELEMENT;
            static constexpr int d = VectorType::dimension;

            using GradientMatrixType = typename DiscretizationType::GradientMatrixType;
			static constexpr int elementNodes = Discretization::elementNodes;
			using ElementType = typename DiscretizationType::ElementType;
			using ElementIndexType = typename DiscretizationType::ElementIndexType;
#if 0
            static constexpr int elementNodes = d + 1;
            using ElementType = std::array<IndexType, elementNodes>;
            using ElementIndexType = std::array<IndexType, elementNodes>;
#endif

            using NodeArrayType = typename IteratorType::template ContainerType<typename PDSimulation::NodeType>;
            using NumberingArrayType = typename IteratorType::template ContainerType<IntType>;

            using Constraint = SoftConstraint<VectorType, elementNodes, IndexType>;
            using Suture = SutureConstraint<VectorType, elementNodes, IndexType>;
            using CollisionSuture = SlidingConstraint <VectorType, elementNodes, IndexType>;

            IntType schurSize = 0;
            IntType matrixSize = 0;
            NumberingArrayType m_numbering; // only number the active nodes, collisionNodes at the bottom
            std::vector<std::map<int, T>> m_tensor;
            T *m_Sigma1 = nullptr;
            T* m_A22 = nullptr;
            T* m_f2 = nullptr;
            T* m_schur = nullptr;

            T *rhs = nullptr;
            T *x = nullptr;

            mutable PardisoWrapper<T, IntType> m_pardiso;
            mutable CudaWrapper<T> m_cuda;
            int* m_rowPtrW_h = nullptr;
            int* m_colIdxPtrW_h = nullptr;
            T* m_valPtrW_h = nullptr;
            T* m_valPtrS_h = nullptr;
            T* m_valPtrD_h = nullptr;


            void initialize(const NodeArrayType &nodeType);

            template <int elementNodesN>
            void accumToTensor(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix,
                               const std::array<IndexType, elementNodesN> &elementIndex);

            template <int elementNodesN>
            void updateTensor(T* const result,
                              const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix,
                              const std::array<IndexType, elementNodesN> &elementIndex);

            // void updatePardiso(const std::vector<Constraint> &collisionConstraints);

            void updateCuda(const std::vector<Constraint> &collisionConstraints,
                            const std::vector<CollisionSuture> &collisionSutures);

            void computeE2Tensor (const std::vector<ElementType> &elements,
                                  const std::vector<PDSimulation::ElementFlag> &flags,
                                  const std::vector<GradientMatrixType> &gradients,
                                  const std::vector<T> &restVol, const T mu);

            inline void setTemp(int v) {
                for (int i=0; i<schurSize; i++) {
                    m_f2[i+v*schurSize] = x[i+matrixSize-schurSize];
                }
            }

            inline void updateTemp(int v) {
                CBLASPolicy<T>::mutiplyAdd(m_f2+schurSize*v, schurSize, T(-1), m_Sigma1, x+matrixSize-schurSize, T(1));
            }

            inline void updateForce(int v) {
                for (int i=0; i<matrixSize; i++)
                    x[i] = rhs[i];
                for (int i=0; i<schurSize; i++)
                    x[i+matrixSize-schurSize] += m_f2[i+v*schurSize];
            }

            void computeTensor(const std::vector<ElementType> &elements,
                               const std::vector<GradientMatrixType> &gradients,
                               const std::vector<T> &restVol, const T muLow, const T muHigh,
                               //const std::vector<Constraint> &constraints,
                               const std::vector<Suture> &sutures
			);

            void computeTensor(const std::vector<ElementType>& elements,
                const std::vector<GradientMatrixType>& gradients,
                const std::vector<T>& restVol, const std::vector<T>& muLow, const std::vector<T>& muHigh,
                //const std::vector<Constraint> &constraints,
                const std::vector<Suture>& sutures
            );



			inline void reInitializePardiso(const std::vector<Constraint> &constraints, const std::vector<Suture> &sutures, const std::vector<Constraint> &fakeSutures) { factPardiso(constraints, sutures, fakeSutures); }

            void reInitializeCuda(const std::vector<Constraint> &collisionConstraints,
                                  const std::vector<CollisionSuture> &collisionSutures);

            template <int elementNodesN>
            void accumToPardiso(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix,
                                const std::array<IndexType, elementNodesN> &elementIndex);


            void initializePardiso(const std::vector<Constraint> &constraints, const std::vector<Suture> &sutures, const std::vector<Constraint> &fakeSutures);

            void initializeCuda(const std::vector<Constraint> &collisionConstraints,
                                const std::vector<CollisionSuture> &collisionSutures);

            void copyIn(const StateVariableType &f, const int v) const {
                // copy in x
                for (Iterator<StateVariableType> iterator(f); !iterator.isEnd(); iterator.next()) {
                    const int number = iterator.value(m_numbering);
                    if (number >= 0)
                        rhs[number] = iterator.value(f)(v + 1);
                }
            }

            void copyOut(StateVariableType &f, const int v) const {
                // copy out x
                for (Iterator<StateVariableType> iterator(f); !iterator.isEnd(); iterator.next()) {
                    const int number = iterator.value(m_numbering);
                    if (number >= 0)
                        iterator.value(f)(v + 1) = x[number];
                }
            }

            inline void forwardSubstitution() const {m_pardiso.forwardSubstitution(rhs, x);}

            inline void diagSolve() const {
                //LOG::SCOPE scope("CudaSolver::diagSolve()");
                m_cuda.solve(&x[m_pardiso.n - schurSize], &rhs[m_pardiso.n - schurSize]);
                for (int i = 0; i < m_pardiso.n - schurSize; i++)
                    x[i] = rhs[i];
            }

            inline void backwardSubstitution() const { m_pardiso.backwardSubstitution(rhs, x);}

            void deallocate();

            void inline releasePardiso() {
                m_pardiso.releasePardisoInternal();
                m_pardiso.deallocate();
            }

            void inline releaseCuda() {
				m_cuda.release();
				if (m_rowPtrW_h) { delete[] m_rowPtrW_h; m_rowPtrW_h = nullptr; }
				if (m_colIdxPtrW_h) { delete[] m_colIdxPtrW_h; m_colIdxPtrW_h = nullptr; }
				if (m_valPtrW_h) { delete[] m_valPtrW_h; m_valPtrW_h = nullptr; }
				if (m_valPtrS_h) { delete[] m_valPtrS_h; m_valPtrS_h = nullptr; }
				if (m_valPtrD_h) { delete[] m_valPtrD_h; m_valPtrD_h = nullptr; }
			}

			void factPardiso(
				const std::vector<Constraint> &constraints,
				const std::vector<Suture> &sutures,
				const std::vector<Constraint> &fakeSutures
				);
        };

} // namespace PhysBAM
