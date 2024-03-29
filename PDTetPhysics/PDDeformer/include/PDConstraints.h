#pragma once

#include <array>
#include <vector>

namespace PhysBAM {
    template <class VectorType, int elementNodeNum, class IndexType> struct SoftConstraint {
        std::array<IndexType, elementNodeNum> m_elementIndex;
        VectorType m_xT;
        typename VectorType::ELEMENT m_stiffness;
        std::array<typename VectorType::ELEMENT, elementNodeNum> m_weights;
		typename VectorType::ELEMENT m_stressLimit;
    };

    template <class VectorType, int elementNodeNum, class IndexType> struct SutureConstraint {
        std::array<IndexType, elementNodeNum> m_elementIndex1, m_elementIndex2;
        typename VectorType::ELEMENT m_stiffness;
        std::array<typename VectorType::ELEMENT, elementNodeNum> m_weights1, m_weights2;
        typename VectorType::ELEMENT m_restLength;
        int m_elementNumber;
    };

    template <class VectorType, int elementNodeNum, class IndexType> struct SlidingConstraint {
        std::array<IndexType, elementNodeNum> m_elementIndex1{}, m_elementIndex2{};
        typename VectorType::ELEMENT m_stiffness;
        std::array<typename VectorType::ELEMENT, elementNodeNum> m_weights1, m_weights2;
        VectorType m_normal;
        int m_elementNumber;
    };

    template <class VectorType, class IndexType> struct NodeToNodesConstraint {  // COURT added
        std::array<IndexType,VectorType::dimension > m_macroNodes{}; // All the T junctions nodes are somehow imbedded on a face (3D) or edge (2D)
        typename VectorType::ELEMENT m_stiffness;
        std::array<typename VectorType::ELEMENT, VectorType::dimension> m_macroWeights;
        int m_microNodeNumber;
    };
}
