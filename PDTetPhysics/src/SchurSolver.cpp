#include "SchurSolver.h"

namespace PhysBAM {
    template<class Discretization, class IntType>
    inline void SchurSolver<Discretization, IntType>::initialize(const NodeArrayType& nodeType) {
        using IteratorType = Iterator<NodeArrayType>;
#ifndef _WIN32
        LOG::SCOPE scope("SchurSolver::initialize()");
#endif
        IteratorType iterator(nodeType);
        iterator.resize(m_numbering);

        int numOfActiveNodes = 0;
        schurSize = 0;
        for (iterator.begin(); !iterator.isEnd(); iterator.next())
            if (iterator.value(nodeType) == PDSimulation::ActiveNode) {
                numOfActiveNodes++;
            }
            else if (iterator.value(nodeType) == PDSimulation::CollisionNode) {
                numOfActiveNodes++;
                schurSize++;
            }
        LOG::cout << "    schursize   = " << schurSize << std::endl;
        LOG::cout << "    matrixsize  = " << numOfActiveNodes << std::endl;
        int activeIdx = 0;
        int collisionIdx = 0;
        for (iterator.begin(); !iterator.isEnd(); iterator.next())
            if (iterator.value(nodeType) == PDSimulation::ActiveNode)
                iterator.value(m_numbering) = activeIdx++;
            else if (iterator.value(nodeType) == PDSimulation::CollisionNode)
                iterator.value(m_numbering) = numOfActiveNodes - schurSize + collisionIdx++;
            else
                iterator.value(m_numbering) = -1;

        m_tensor.resize(numOfActiveNodes);
    }

    template<class Discretization, class IntType>
    inline void SchurSolver<Discretization, IntType>::updatePardiso(const std::vector<Constraint>& collisionConstraints) {
        //LOG::SCOPE scope("SchurSolver::updatePardiso");
        const IntType& n = m_pardiso.n;
        const IntType& nnz = m_pardiso.rowIndex[n];
        if (schurSize)
            for (int i = 0; i < schurSize * schurSize; i++)
                m_pardiso.schur[i] = m_originalValue[i];
        else
            for (int i = 0; i < nnz; i++)
                m_pardiso.value[i] = m_originalValue[i];

        for (int c = 0; c < collisionConstraints.size(); c++) {
            auto& constraint = collisionConstraints[c];
            if (constraint.m_stiffness != 0) {
                MATRIX_MXN<T> stiffnessMatrix;
                DiscretizationType::computeConstraintTensor(stiffnessMatrix, constraint);
                updateTensor<elementNodes>(stiffnessMatrix, constraint.m_elementIndex);
            }
        }
        if (schurSize) {
            m_pardiso.factSchur();
        }
        else {
            m_pardiso.factorize();
        }
    }

    template<class Discretization, class IntType>
    inline void SchurSolver<Discretization, IntType>::computeTensor(const std::vector<ElementType>& elements, const std::vector<GradientMatrixType>& gradients, const std::vector<T>& restVol, const T mu, const std::vector<Constraint>& constraints, const std::vector<Suture>& sutures) {

        for (int e = 0; e < elements.size(); e++) {
            MATRIX_MXN<T> stiffnessMatrix;
            DiscretizationType::computeElementTensor(stiffnessMatrix, gradients[e], -2 * mu * restVol[e]);
            accumToTensor<elementNodes>(stiffnessMatrix,
                DiscretizationType::getElementIndex(elements[e]));
        }

        for (int c = 0; c < constraints.size(); c++) {
            MATRIX_MXN<T> stiffnessMatrix;
            DiscretizationType::computeConstraintTensor(stiffnessMatrix, constraints[c]);

            accumToTensor<elementNodes>(stiffnessMatrix,
                constraints[c].m_elementIndex);
        }

        for (int c = 0; c < sutures.size(); c++) {
            MATRIX_MXN<T> stiffnessMatrix;
            std::array<IndexType, elementNodes * 2> elementIndex;
            DiscretizationType::computeSutureTensor(stiffnessMatrix, elementIndex, sutures[c]);
            accumToTensor<elementNodes * 2>(stiffnessMatrix,
                elementIndex);
        }
    }

    template<class Discretization, class IntType>
    inline void SchurSolver<Discretization, IntType>::initializePardiso() {

        IntType nnz = 0;
        for (int i = 0; i < m_tensor.size(); i++)
            nnz += (IntType)m_tensor[i].size();

        LOG::cout << "nnz = " << nnz << std::endl;
        m_pardiso.initialize((IntType)m_tensor.size(), nnz, schurSize);

        m_pardiso.rowIndex[0] = 0;
        for (int i = 0; i < m_pardiso.n; i++)
            m_pardiso.rowIndex[i + 1] = m_pardiso.rowIndex[i] + (IntType)m_tensor[i].size();

        const IntType& n = m_pardiso.n;

        size_t idx = 0;
        for (const auto& r : m_tensor)
            for (const auto& e : r) {
                m_pardiso.column[idx] = e.first;
                m_pardiso.value[idx] = -e.second;
                idx++;
            }
        m_rhs = new T[n];
        m_x = new T[n]();

        if (!schurSize) {
            m_originalValue = new T[nnz];
            for (IntType i = 0; i<nnz; i++)
                m_originalValue[i] = m_pardiso.value[i];
        }
        else {
            m_schur = new T[schurSize * schurSize];
            m_pardiso.schur = m_schur;
        }

        //dumper::writeCSRbyte(m_pardiso.n, m_pardiso.rowIndex, m_pardiso.column, m_pardiso.value);
        m_pardiso.factorize();

        if (schurSize) {
            m_originalValue = new T[schurSize * schurSize];
            for (IntType i = 0; i<schurSize * schurSize; i++)
                m_originalValue[i] = m_pardiso.schur[i];
            m_pardiso.factSchur();
        }

        m_tensor.resize(0);
    }

    template<class Discretization, class IntType>
    inline void SchurSolver<Discretization, IntType>::factPardiso(const std::vector<Constraint>& constraints, const std::vector<Suture>& sutures, const std::vector<Constraint>& fakeSutures) {
        size_t idx = 0;
        for (const auto& r : m_tensor)
            for (const auto& e : r) {
                m_pardiso.value[idx] = -e.second;
                idx++;
            }
        // dumper::writeCSRbyte(m_pardiso.n, m_pardiso.rowIndex, m_pardiso.column, m_pardiso.value, m_pardiso.n, "orig_i.txt", "orig_a.txt");

        for (int c = 0; c < constraints.size(); c++)
            if (constraints[c].m_stiffness != 0) {
                MATRIX_MXN<T> stiffnessMatrix;
                DiscretizationType::computeConstraintTensor(stiffnessMatrix, constraints[c]);
                accumToPardiso<elementNodes>(stiffnessMatrix,
                    constraints[c].m_elementIndex);
            }

        for (int c = 0; c < fakeSutures.size(); c++)
            if (fakeSutures[c].m_stiffness != 0) {
                MATRIX_MXN<T> stiffnessMatrix;
                DiscretizationType::computeConstraintTensor(stiffnessMatrix, fakeSutures[c]);
                accumToPardiso<elementNodes>(stiffnessMatrix,
                    fakeSutures[c].m_elementIndex);
            }

        for (int c = 0; c < sutures.size(); c++)
            if (sutures[c].m_stiffness != 0) {
                MATRIX_MXN<T> stiffnessMatrix;
                std::array<IndexType, elementNodes * 2> elementIndex;
                DiscretizationType::computeSutureTensor(stiffnessMatrix, elementIndex, sutures[c]);
                accumToPardiso<elementNodes * 2>(stiffnessMatrix,
                    elementIndex);
            }

        // dumper::writeCSRbyte(m_pardiso.n, m_pardiso.rowIndex, m_pardiso.column, m_pardiso.value, m_pardiso.n, "new_i.txt", "new_a.txt");

        m_pardiso.numericFact();

        if (schurSize) {
            for (IntType i = 0; i < schurSize * schurSize; i++)
                m_Sigma1[i] = m_pardiso.schur[i] - m_A22[i];
        }
    }

    template<class Discretization, class IntType>
    template<int elementNodesN>
    void SchurSolver<Discretization, IntType>::
        accumToTensor(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix, 
                const std::array<IndexType, elementNodesN>& elementIndex) {
        using IteratorType = Iterator<NodeArrayType>;
        for (int i = 0; i < elementNodes; i++) {
            int row = IteratorType::at(m_numbering, elementIndex[i]);
            if (row >= 0) {
                for (int j = 0; j < elementNodesN; j++) {
                    int col = IteratorType::at(m_numbering, elementIndex[j]);
                    if (col >= row) {
                        if (m_tensor[row].find(col) == m_tensor[row].end())
                            m_tensor[row].insert(std::pair<int, T>(col, stiffnessMatrix(i + 1, j + 1)));
                        else
                            m_tensor[row][col] += stiffnessMatrix(i + 1, j + 1);
                    }
                }
            }
        }
    }

    template<class Discretization, class IntType>
    template<int elementNodesN>
    inline void SchurSolver<Discretization, IntType>::
        updateTensor(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix, const std::array<IndexType, elementNodesN>& elementIndex) {
        // only cared about upper triangular in row major
        //LOG::SCOPE scope("SchurSolver::updateTensor");
        // TODO: think about sort the Indecies by numbering first; wether this will help with
        // branch prediction
        IntType& n = m_pardiso.n;
        using IteratorType = Iterator<NodeArrayType>;
        if (schurSize) {
            for (int i = 0; i < elementNodesN; i++) {
                int row = IteratorType::at(m_numbering, elementIndex[i]);
                if (row < 0)
                    continue;

                for (int j = 0; j < elementNodesN; j++) {
                    int col = IteratorType::at(m_numbering, elementIndex[j]);
                    if (col < row)
                        continue;

                    assert(row >= n - schurSize);
                    m_pardiso.schur[(row - n + schurSize) * schurSize + col - n + schurSize] -= stiffnessMatrix(i + 1, j + 1);
                }
            }
        }
        else {
            for (int i = 0; i < elementNodesN; i++) {
                int row = IteratorType::at(m_numbering, elementIndex[i]);
                if (row < 0)
                    continue;

                for (int j = 0; j < elementNodesN; j++) {
                    int col = IteratorType::at(m_numbering, elementIndex[j]);
                    if (col < row)
                        continue;

                    int index = -1;
                    for (int k = m_pardiso.rowIndex[row]; k < m_pardiso.rowIndex[row + 1]; k++) {
                        if (m_pardiso.column[k] == col) {
                            index = k;
                            break;
                        }
                    }

                    if (index == -1) {
                        std::cerr << "cannot find the entry" << std::endl;
                        exit(1);
                    }
                    else
                        // stiffnessMatrix is negative definite
                        m_pardiso.value[index] -= stiffnessMatrix(i + 1, j + 1);
                }
            }
        }
    }
}

namespace PhysBAM {
    template struct SchurSolver<TetrahedralDiscretization<std::vector<VECTOR<float, 3>>>, int>;
}