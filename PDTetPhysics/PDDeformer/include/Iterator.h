#pragma once

#include <vector>

#include <Common_Tools/Math_Tools/RANGE_ITERATOR.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>

namespace PhysBAM {

template <class StateVariableType> struct Iterator;

template <class T, int d> struct Iterator<ARRAY<T, VECTOR<int, d>>> {
    using DataType = T;
    using IndexType = VECTOR<int, d>;
    template <class ContainedData> using ContainerType = ARRAY<ContainedData, VECTOR<int, d>>;
    using ArrayType = ContainerType<DataType>;

    const RANGE<VECTOR<int, d>> m_range;
    VECTOR<int, d> index;

    Iterator(const ArrayType &v) : m_range(v.Domain_Indices()) { index = m_range.min_corner; }

    void begin() { index = m_range.min_corner; }

    bool isEnd() { return !m_range.Lazy_Inside(index); }

    void next() {
        for (int i = d; i >= 1; i--) {
            if (index(i) < m_range.max_corner(i) || i == 1) {
                index(i)++;
                break;
            } else
                index(i) = m_range.min_corner(i);
        }
    }

    template <class ContainedData> void resize(ContainerType<ContainedData> &v) const {
        v.Resize(m_range);
        v.Fill(ContainedData());
    }

    template <class ContainedData> ContainedData &value(ContainerType<ContainedData> &v) const { return v(index); }

    template <class ContainedData> const ContainedData &value(const ContainerType<ContainedData> &v) const {
        return v(index);
    }

    template <class ContainedData>
    static ContainedData &at(ContainerType<ContainedData> &v, VECTOR<int, d> nodeNumber) {
        return v(nodeNumber);
    }

    template <class ContainedData>
    static const ContainedData &at(const ContainerType<ContainedData> &v, VECTOR<int, d> nodeNumber) {
        return v(nodeNumber);
    }
};

template <class T> struct Iterator<std::vector<T>> {
    using DataType = T;
    using IndexType = int;
    template <class ContainedData> using ContainerType = std::vector<ContainedData>;
    using ArrayType = ContainerType<DataType>;

    const int m_range;
    IndexType index;

    Iterator(const ArrayType &v) : m_range((int) v.size()) { index = 0; }

    void begin() { index = 0; }

    bool isEnd() { return index < 0 || index >= m_range; }

    void next() { index++; }

    template <class ContainedData> void resize(ContainerType<ContainedData> &v) const {
        v.clear();
        v.resize(m_range);
    }

    template <class ContainedData> ContainedData &value(ContainerType<ContainedData> &v) const { return v[index]; }

    template <class ContainedData> const ContainedData &value(const ContainerType<ContainedData> &v) const {
        return v[index];
    }

    template <class ContainedData> static ContainedData &at(ContainerType<ContainedData> &v, IndexType nodeNumber) {
        return v[nodeNumber];
    }

    template <class ContainedData>
    static const ContainedData &at(const ContainerType<ContainedData> &v, IndexType nodeNumber) {
        return v[nodeNumber];
    }
};

} // namespace PhysBAM
