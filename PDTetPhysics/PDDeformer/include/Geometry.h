#pragma once

#include <vector>

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>

#include "SimulationFlags.h"
#include "Iterator.h"
#include "Map.h"

namespace PhysBAM {

template <int d> struct Geometry {
    using IndexType = VECTOR<int, d>;
    using RangeType = RANGE<IndexType>;
    using NodeType = PDSimulation::NodeType;

    IndexType m_gridSize;          // Domain size in cell counts
    RangeType m_paddedCellRange;   // One layer added, exclusively for marking cells as Dirichlet
    RangeType m_unpaddedCellRange; // The interior cell layer, which is allowed to include interior cells
    RangeType m_nodeRange;         // Node range spanning the unpadded cells

    enum CellType { InactiveCell = 0, ActiveCell = 1, DirichletCell = 2 };

    ARRAY<CellType, IndexType> m_cellType;

    Geometry(const IndexType gridSize) : m_gridSize(gridSize) { initializeArrays(); }

    void initializeArrays() {
        m_paddedCellRange = RangeType(-IndexType::All_Ones_Vector(), m_gridSize);
        m_unpaddedCellRange = RangeType(IndexType(), m_gridSize - 1);
        m_nodeRange = RangeType(IndexType(), m_gridSize);

        m_cellType.Resize(m_paddedCellRange);
    }

    template <class NodeArray> void initializeNodeType(NodeArray &nodeType) {
        using IteratorType = Iterator<NodeArray>;
        using MapType = Map<NodeArray>;

        MapType::resizeGeometry(nodeType, m_unpaddedCellRange);

        for (IteratorType iterator(nodeType); !iterator.isEnd(); iterator.next())
            iterator.value(nodeType) = PDSimulation::InactiveNode;

        for (RANGE_ITERATOR<d> cellIterator(m_unpaddedCellRange); cellIterator.Valid(); cellIterator.Next()) {
            const auto &cellIndex = cellIterator.Index();
            if (m_cellType(cellIndex) == ActiveCell || m_cellType(cellIndex) == DirichletCell)
                for (RANGE_ITERATOR<d> nodeIterator(RangeType(cellIndex, cellIndex + 1)); nodeIterator.Valid();
                     nodeIterator.Next())
                    IteratorType::at(nodeType, MapType::coordinateToIndex(nodeIterator.Index(), m_unpaddedCellRange)) =
                        PDSimulation::ActiveNode;
        }

        for (RANGE_ITERATOR<d> cellIterator(m_paddedCellRange); cellIterator.Valid(); cellIterator.Next()) {
            const auto &cellIndex = cellIterator.Index();
            if (m_cellType(cellIndex) == DirichletCell)
                for (RANGE_ITERATOR<d> nodeIterator(RangeType(cellIndex, cellIndex + 1)); nodeIterator.Valid();
                     nodeIterator.Next()) {
                    const auto &nodeIndex = nodeIterator.Index();
                    if (m_nodeRange.Lazy_Inside(nodeIndex) &&
                        IteratorType::at(nodeType, MapType::coordinateToIndex(nodeIndex, m_unpaddedCellRange)) ==
                        PDSimulation::ActiveNode)
                        IteratorType::at(nodeType, MapType::coordinateToIndex(nodeIndex, m_unpaddedCellRange)) =
                            static_cast<NodeType>(m_cellType(cellIndex));
                }
        }
    }
};
} // namespace PhysBAM
