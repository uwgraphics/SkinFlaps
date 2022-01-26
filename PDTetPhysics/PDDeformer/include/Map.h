#pragma once

#include <vector>

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>

namespace PhysBAM {

template <class StateVariableType> struct Map;

template <class DataType, int d> struct Map<ARRAY<DataType, VECTOR<int, d>>> {
    using IndexType = VECTOR<int, d>;

    static IndexType coordinateToIndex(VECTOR<int, d> coordinate, RANGE<VECTOR<int, d>> unpaddedCellRange) {
        return coordinate;
    }

    template <class ContainedData>
    static void resizeGeometry(ARRAY<ContainedData, VECTOR<int, d>> &v, const RANGE<VECTOR<int, d>> unpaddedCellRange) {
        v.Resize(RANGE<VECTOR<int, d>>(unpaddedCellRange.min_corner, unpaddedCellRange.max_corner + 1));
    }
};

template <class DataType> struct Map<std::vector<DataType>> {

    using IndexType = int;

    template <int d>
    static IndexType coordinateToIndex(VECTOR<int, d> coordinate, RANGE<VECTOR<int, d>> unpaddedCellRange) {
        VECTOR<int, d> gridSize = unpaddedCellRange.max_corner - unpaddedCellRange.min_corner + 1;
        VECTOR<int, d> node = coordinate - unpaddedCellRange.min_corner;

        int index = 0, base = 1;
        for (int i = d; i >= 1; i--) {
            index += base * node(i);
            base *= (gridSize(i) + 1);
        }

        return index;
    }

    template <class ContainedData, int d>
    static void resizeGeometry(std::vector<ContainedData> &v, const RANGE<VECTOR<int, d>> unpaddedCellRange) {
        v.clear();
        VECTOR<int, d> nodeSize = unpaddedCellRange.max_corner - unpaddedCellRange.min_corner + 2;
        v.resize(nodeSize.Product());
    }
};
} // namespace PhysBAM
