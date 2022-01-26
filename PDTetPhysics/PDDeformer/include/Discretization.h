#pragma once

#include <map>

#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/DIAGONAL_MATRIX_3X3.h>

#include <PhysBAM_Tools/Matrices/MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>

#include <PhysBAM_Tools/Math_Tools/FACTORIAL.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <Common_Tools/Vectors/VectorSTLInterface.h>

#include "PDConstraints.h"
#include "Iterator.h"

namespace PhysBAM {

template <class StateVariable> struct TetrahedralDiscretization {
    using StateVariableType = StateVariable;
    using IteratorType = Iterator<StateVariableType>;
    using VectorType = typename IteratorType::DataType;
    using T = typename VectorType::ELEMENT;
    static constexpr int d = VectorType::dimension;
    using IndexType = typename IteratorType::IndexType;

    static constexpr int elementNodes = d + 1;
    using ElementType = std::array<IndexType, elementNodes>;
    using ElementIndexType = std::array<IndexType, elementNodes>;
    using WeightType = std::array<T, elementNodes>;
    using ConstraintType = SoftConstraint<VectorType, elementNodes, IndexType>;
    using SutureType = SutureConstraint<VectorType, elementNodes, IndexType>;


    using ShapeMatrixType = MATRIX<T, d>;
    using GradientMatrixType = MATRIX<T, d>;
    using MatrixType = MATRIX<T, d>;

    static inline ElementIndexType &getElementIndex(ElementType &element) { return element; }

    static inline const ElementIndexType &getElementIndex(const ElementType &element) { return element; }

    static inline VectorType interpolateX(const ElementIndexType &elementIndex, const WeightType &weights,
                                          const StateVariableType &x) {
        VectorType result;
        for (int v = 0; v < elementNodes; v++) {
            result += weights[v] * IteratorType::at(x, elementIndex[v]);
        }
        return result;
    }

    static inline void distributeForces(const VectorType force, const ElementIndexType &elementIndex, const WeightType weights,
                                        StateVariableType &f) {
        for (int v = 0; v < elementNodes; v++) {
            IteratorType::at(f, elementIndex[v]) += weights[v] * force;
        }
    }

    static bool isInsideElement(const ElementIndexType &elementIndex, const StateVariableType &x,
                                const VectorType &vector) {
        std::array<T, elementNodes> volumes;
        MatrixType M;
        for (int i = 0; i < elementNodes; i++) {
            int col = 1;
            for (int v = 0; v < elementNodes; v++) {
                if (v == i)
                    continue;
                M.Column(col) = IteratorType::at(x, elementIndex[v]) - vector;
                col++;
            }
            volumes[i] = std::abs(M.Determinant());
        }

        T sum = 0;
        for (int i = 0; i < elementNodes; i++)
            sum += volumes[i];

        for (int v = 0; v < elementNodes - 1; v++)
            M.Column(v + 1) = IteratorType::at(x, elementIndex[v + 1]) - IteratorType::at(x, elementIndex[0]);

        return std::abs(sum - std::abs(M.Determinant()))/sum < 1e-5;
    }

    static WeightType computeWeight(const ElementIndexType &elementIndex, const StateVariableType &x,
                                    const VectorType &vector) {
        std::array<T, elementNodes> weights;
        MatrixType M;
        for (int i = 0; i < elementNodes; i++) {
            int col = 1;
            for (int v = 0; v < elementNodes; v++) {
                if (v == i)
                    continue;
                M.Column(col) = IteratorType::at(x, elementIndex[v]) - vector;
                col++;
            }
            weights[i] = std::abs(M.Determinant());
        }

        T sum = 0;
        for (int i = 0; i < elementNodes; i++)
            sum += weights[i];
        for (int i = 0; i < elementNodes; i++)
            weights[i] /= sum;

        return weights;
    }

    #if 0
    static void initializeElements(std::vector<ElementType> &elements, const GeometryType &geometry) {
        LOG::SCOPE scope("TetrahedralDiscretization::initializeElements()");

        for (RANGE_ITERATOR<d> cellIterator(geometry.m_unpaddedCellRange); cellIterator.Valid(); cellIterator.Next()) {
            const auto &cellIndex = cellIterator.Index();
            if (geometry.m_cellType(cellIndex) == GeometryType::CellType::ActiveCell) {
                for (const auto &subelement : TesselationType::SubelementIndices) {
                    ElementType element;
                    for (int v = 0; v < d + 1; v++)
                        element[v] =
                            MapType::coordinateToIndex(cellIndex + STLInterface<VECTOR<int, d>>::Wrap(subelement[v]),
                                                       geometry.m_unpaddedCellRange);
                    elements.push_back(element);
                }
            }
        }

        LOG::cout << "size: " << elements.size() << std::endl;
    }

    static void addElements(std::vector<ElementType> &elements, const std::array<int, (1 << d)>& cube) {
        for (const auto &subelement : HypercubeIndices<d>::SubelementIndices) {
            ElementType element;
            for (int v = 0; v < d + 1; v++)
                element[v] = cube[subelement[v]];
            elements.push_back(element);
        }
    }
#endif

    static inline void computeShapeMatrix(ShapeMatrixType &shapeMatrix, const ElementType &element,
                                          const StateVariableType &x) {
        for (int v = 0; v < d; v++)
            shapeMatrix.Column(v + 1) = IteratorType::at(x, element[v + 1]) - IteratorType::at(x, element[0]);
    }

    static inline void computeGradientMatrixAndRestVolume(const ElementType &element, const StateVariableType &restX,
                                                          GradientMatrixType &gradientMatrix, T &elementRestVolume) {
        ShapeMatrixType restShapeMatrix;
        computeShapeMatrix(restShapeMatrix, element, restX);

        constexpr T scale = (T)1. / FACTORIAL<d>::value;
        gradientMatrix = restShapeMatrix.Inverse();
        elementRestVolume = scale * std::abs(restShapeMatrix.Determinant());
    }

    static inline MatrixType computeGradient(const ShapeMatrixType &shapeMatrix, const GradientMatrixType &gradientMatrix) {
        return shapeMatrix * gradientMatrix;
    }

    static inline ShapeMatrixType computeDivergence(const MatrixType &P, const GradientMatrixType &gradientMatrix,
                                                    const T scale) {
        return P.Times_Transpose(gradientMatrix) * scale;
    }

    static inline void distributeForces(const ShapeMatrixType &forceMatrix, const ElementType &element, StateVariableType &f) {
        for (int v = 0; v < d; v++) {
            IteratorType::at(f, element[v + 1]) += forceMatrix.Column(v + 1);
            IteratorType::at(f, element[0]) -= forceMatrix.Column(v + 1);
        }
    }
    static void computeElementTensor(MATRIX_MXN<T>& stiffnessMatrix,
                              const GradientMatrixType &gradientMatrix,
                              const T constant) {
        MATRIX_MXN<T> DmInverse(gradientMatrix);
        MATRIX_MXN<T> S(elementNodes, d);
        for (int i = 0; i < d; i++) {
            S(1, i + 1) = -1;
            S(i + 2, i + 1) = 1;
        }
        stiffnessMatrix = S * DmInverse * DmInverse.Transposed() * S.Transposed() * constant;
    }

    static void computeConstraintTensor(MATRIX_MXN<T>& stiffnessMatrix, const ConstraintType &constraint) {
        MATRIX_MXN<T> weightMatrix(1, elementNodes);
        for (int i = 0; i < elementNodes; i++) {
            weightMatrix(1, i + 1) = constraint.m_weights[i];
        }
        stiffnessMatrix = weightMatrix.Transpose_Times(weightMatrix) * -constraint.m_stiffness;
    }

    static void computeSutureTensor(MATRIX_MXN<T>& stiffnessMatrix, std::array<IndexType, elementNodes*2>& elementIndex, const SutureType &suture) {
        for  (int i=0; i<elementNodes; i++) {
            elementIndex[i] = suture.m_elementIndex1[i];
            elementIndex[i+elementNodes] = suture.m_elementIndex2[i];
        }
        MATRIX_MXN<T> weightMatrix(1, elementNodes*2);
        for (int i = 0; i < elementNodes; i++) {
            weightMatrix(1, i + 1) = suture.m_weights1[i];
            weightMatrix(1, elementNodes + i + 1) = -suture.m_weights2[i];
        }
        stiffnessMatrix = weightMatrix.Transpose_Times(weightMatrix) * -suture.m_stiffness;
    }
};

} // namespace PhysBAM
