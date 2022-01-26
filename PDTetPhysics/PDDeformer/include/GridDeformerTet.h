//#####################################################################
// Copyright (c) 2019, Eftychios Sifakis, Yutian Tao, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#pragma once

#include <functional>
//#include <map>
#include "Iterator.h"

#include "Discretization.h"
// #include "Algebra.h"
#include "SimulationFlags.h"
#include "PDConstraints.h"


#include <Common/KernelCommon.h>
#include "ReshapeDataStructure.h"

#include <PhysBAM_Tools/Math_Tools/FACTORIAL.h>



namespace PhysBAM {

    template <class StateVariable> struct GridDeformerTet;

    template <class dataType, int dim> struct GridDeformerTet<std::vector<VECTOR<dataType,dim>>>
    {

        using StateVariableType = std::vector<VECTOR<dataType, dim>>;

        using IteratorType = Iterator<StateVariableType>;
        using IndexType = typename IteratorType::IndexType;
        using VectorType = typename IteratorType::DataType;

        using T = dataType;
        static constexpr int d = VectorType::dimension;
        static constexpr int elementNodes = d + 1;

        using ElementType = std::array<IndexType, elementNodes>;
        using ElementIndexType = std::array<IndexType, elementNodes>;
        using ShapeMatrixType = MATRIX<T, d>;
        using GradientMatrixType = MATRIX<T, d>;
        using DiscretizationType = TetrahedralDiscretization<StateVariableType>;
        // using GeometryType = Geometry<d>;

        using WeightType = std::array<T, elementNodes>;

        //using MapType = Map<StateVariableType>;
        using Constraint = SoftConstraint<VectorType, elementNodes, IndexType>;
        using Suture = SutureConstraint<VectorType, elementNodes, IndexType>;
        using CollisionSuture = SlidingConstraint<VectorType, elementNodes, IndexType>;

        static_assert(d == 2 || d == 3, "only 2D/3D discretizations supported");

        using MatrixType = MATRIX<T, d>;
        using DiagonalMatrixType = DIAGONAL_MATRIX<T, d>;
        using NodeArrayType = typename IteratorType::template ContainerType<PDSimulation::NodeType>;

        using Tarch = typename SIMD_Numeric_Kernel::template SIMDArchitectureAVX2<T>;
        static constexpr int BlockWidth = 16;
        static constexpr int Alignment = 64;
        using BlockedShapeMatrixType = T (*) [d+1][d][BlockWidth];
        using BlockedMatrixType = T (*) [d*d][BlockWidth];
        using BlockedScalarType = T (*) [BlockWidth];
        using BlockedElementType = int (*) [d+1][BlockWidth];

        // T m_uniformMu;

        StateVariableType m_X;
        NodeArrayType m_nodeType;

        std::vector<ElementType> m_elements;
        std::vector<T> m_muLow;
        std::vector<T> m_muHigh;
        std::vector<T> m_rangeMin;
        std::vector<T> m_rangeMax;

        std::vector<PDSimulation::ElementFlag> m_elementFlags;
        std::vector<T> m_elementRestVolume;

        std::vector<Constraint> m_constraints;
		std::vector<Constraint> m_fakeSutures; // normal constraints that come in pairs so each 2 actually form a fake suture, x_T is never used for this struct
        std::vector<Constraint> m_collisionConstraints;
        std::vector<Suture> m_sutures;
        std::vector<CollisionSuture> m_collisionSutures;

        std::vector<GradientMatrixType> m_gradientMatrix;

        // reshaped data
        int m_nUncollisionBlocks;
        int m_nCollisionBlocks;
        BlockedShapeMatrixType m_reshapeUncollisionX;
        BlockedShapeMatrixType m_reshapeCollisionX;
        BlockedMatrixType m_reshapeUncollisionGradientMatrix;
        BlockedMatrixType m_reshapeCollisionGradientMatrix;
        BlockedScalarType m_reshapeUncollisionElementRestVolume;
        BlockedScalarType m_reshapeCollisionElementRestVolume;

        BlockedScalarType m_reshapeCollisionMuLow;
        BlockedScalarType m_reshapeCollisionMuHigh;
        BlockedScalarType m_reshapeCollisionRangeMin;
        BlockedScalarType m_reshapeCollisionRangeMax;

        BlockedScalarType m_reshapeUncollisionMuLow;
        BlockedScalarType m_reshapeUncollisionMuHigh;
        BlockedScalarType m_reshapeUncollisionRangeMin;
        BlockedScalarType m_reshapeUncollisionRangeMax;

        // auxilary structure
        std::vector<int> m_reshapeUncollisionIndicesOffsets;
        std::vector<int> m_reshapeCollisionIndicesOffsets;
        std::vector<int> m_reshapeUncollisionIndicesValues;
        std::vector<int> m_reshapeCollisionIndicesValues;
        BlockedElementType m_reshapeUncollisionElement;
        BlockedElementType m_reshapeCollisionElement;

        // std::function<void(const GeometryType &, const NodeArrayType &, StateVariableType &)> m_clearDirichlet;

		GridDeformerTet(T muIn = 100) /*:m_uniformMu(muIn)*/ {
			m_nCollisionBlocks = 0;
			m_nCollisionBlocks = 0;

			m_reshapeUncollisionX = nullptr;
			m_reshapeCollisionX = nullptr;
			m_reshapeUncollisionGradientMatrix = nullptr;
			m_reshapeCollisionGradientMatrix = nullptr;
			m_reshapeUncollisionElementRestVolume = nullptr;
			m_reshapeCollisionElementRestVolume = nullptr;

			m_reshapeUncollisionElement = nullptr;
			m_reshapeCollisionElement = nullptr;

			static_assert((BlockWidth * sizeof(T)) % Alignment == 0, "Blocks don't have requisite alignment");
		}

        void initializeDeformer();
        void initializeUndeformedState();
		void initializeAuxiliaryStructures();
        void updatePositionBasedState(const PDSimulation::ElementFlag flag, const T rangeMin = 1, const T rangeMax = 1);
        void addCollisionForce(StateVariableType &f) const;
        void addConstraintForce(StateVariableType &f) const;
        void addElasticForce(StateVariableType &SIMDf, const PDSimulation::ElementFlag flag, const T rangeMin, const T rangeMax, const T weightProportion) const;
        void clearDirichlet(StateVariableType &var) {
            /*
            for (RANGE_ITERATOR<d> nodeIterator(geometry.m_nodeRange); nodeIterator.Valid(); nodeIterator.Next()) {
                const auto &nodeIndex = nodeIterator.Index();
                if (IteratorType::at(nodeType, MapType::coordinateToIndex(nodeIndex, geometry.m_unpaddedCellRange)) == PDSimulation::DirichletNode) {
                    IteratorType::at(var, MapType::coordinateToIndex(nodeIndex, geometry.m_unpaddedCellRange)) =
                        VectorType();
                }
            }
            */
            for (IteratorType iterator(var); !iterator.isEnd(); iterator.next())
                if (iterator.value(m_nodeType) == PDSimulation::DirichletNode)
                    iterator.value(var) = VectorType();
        }

        void deallocateAuxiliaryStructures();
        void initializeCollisionElements();
    };

} // namespace PhysBAM
