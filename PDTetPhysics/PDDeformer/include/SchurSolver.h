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
#include <chrono>

namespace PhysBAM {

template <class Discretization, class IntType> struct SchurSolver {

    using DiscretizationType = Discretization;
    using StateVariableType = typename DiscretizationType::StateVariableType;
    using IteratorType = Iterator<StateVariableType>;

    using IndexType = typename IteratorType::IndexType;
    using VectorType = typename IteratorType::DataType;
    using T = typename VectorType::ELEMENT;
    static constexpr int d = VectorType::dimension;

    using GradientMatrixType = MATRIX<T, d>;

    static constexpr int elementNodes = d + 1;
    using ElementType = std::array<IndexType, elementNodes>;
    using ElementIndexType = std::array<IndexType, elementNodes>;
    using NodeArrayType = typename IteratorType::template ContainerType<typename PDSimulation::NodeType>;
    using NumberingArrayType = typename IteratorType::template ContainerType<IntType>;

    using Constraint = SoftConstraint<VectorType, elementNodes, IndexType>;
    using Suture = SutureConstraint<VectorType, elementNodes, IndexType>;

    IntType schurSize;
    NumberingArrayType m_numbering; // only number the active nodes, collisionNodes at the bottom
    std::vector<std::map<int, T>> m_tensor;
    T *m_originalValue = nullptr;
    T *m_schur = nullptr;
    T *m_x = nullptr;
    T *m_rhs = nullptr;
    mutable PardisoWrapper<T, IntType> m_pardiso;

    void initialize(const NodeArrayType& nodeType);

    template <int elementNodesN>
    void accumToTensor(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix,
        const std::array<IndexType, elementNodesN>& elementIndex);

    template <int elementNodesN>
    void updateTensor(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix,
        const std::array<IndexType, elementNodesN>& elementIndex);

    void updatePardiso(
        const std::vector<Constraint> &collisionConstraints
    );

    void computeTensor(const std::vector<ElementType> &elements,
                       const std::vector<GradientMatrixType> &gradients,
                       const std::vector<T> &restVol, const T mu,
                      // const std::vector<Constraint> &constraints,
        const std::vector<Suture>& sutures);

    inline void reInitializePardiso(const std::vector<Constraint>& constraints, const std::vector<Suture>& sutures, const std::vector<Constraint>& fakeSutures) { factPardiso(constraints, sutures, fakeSutures); }

    template <int elementNodesN>
    void accumToPardiso(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix,
        const std::array<IndexType, elementNodesN>& elementIndex);


    void initializePardiso(const std::vector<Constraint>& constraints, const std::vector<Suture>& sutures, const std::vector<Constraint>& fakeSutures);


    void copyIn(const StateVariableType &f, const int v) const {
        // copy in x
        for (Iterator<StateVariableType> iterator(f); !iterator.isEnd(); iterator.next()) {
            const int number = iterator.value(m_numbering);
            if (number >= 0)
                m_rhs[number] = iterator.value(f)(v + 1);
        }
    }

    void copyOut(StateVariableType &f, const int v) const {
        // copy out x
        for (Iterator<StateVariableType> iterator(f); !iterator.isEnd(); iterator.next()) {
            const int number = iterator.value(m_numbering);
            if (number >= 0)
                iterator.value(f)(v + 1) = m_x[number];
        }
    }

    void solve() const {
#if TIMING
        auto start1 = std::chrono::steady_clock::now();
#endif
        m_pardiso.forwardSubstitution(m_rhs, m_x);
#if TIMING
         auto end1 = std::chrono::steady_clock::now();
         std::chrono::duration<double> elapsed_seconds1 = end1 - start1;
         std::cout<<"Forward Substitution Time: "<<elapsed_seconds1.count()<<" s"<<std::endl;
#endif

        m_pardiso.diagSolve(m_x, m_rhs);

#if TIMING
        auto start2 = std::chrono::steady_clock::now();
#endif
        m_pardiso.backwardSubstitution(m_rhs, m_x);

#if TIMING
         auto end2 = std::chrono::steady_clock::now();
         std::chrono::duration<double> elapsed_seconds2 = end2 - start2;
         std::cout<<"Backward Substitution Time: "<<elapsed_seconds2.count()<<" s"<<std::endl;
#endif
    }

    void inline releasePardiso() {
        m_pardiso.releasePardisoInternal();
        m_pardiso.deallocate();
    }


    void inline deallocate() {
        if (m_originalValue) {
            delete[] m_originalValue;
            m_originalValue = NULL;
        }
        if (m_schur) {
            delete m_schur;
            m_schur = NULL;
            m_pardiso.schur = NULL;
        }
        if (m_x) {
            delete[] m_x;
            m_x = NULL;
        }
        if (m_rhs) {
            delete[] m_rhs;
            m_rhs = NULL;
        }
    }

    void factPardiso(
        const std::vector<Constraint>& constraints,
        const std::vector<Suture>& sutures,
        const std::vector<Constraint>& fakeSutures
    );

};


} // namespace PhysBAM
