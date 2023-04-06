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
            if (iterator.value(nodeType) == NodeType::Active) {
                numOfActiveNodes++;
            }
            else if (iterator.value(nodeType) == NodeType::Collision) {
                numOfActiveNodes++;
                schurSize++;
               // LOG::cout << "Must not have collision nodes" << std::endl;
                //exit(1);
            }
        LOG::cout << "    schursize   = " << schurSize << std::endl;
        LOG::cout << "    matrixsize  = " << numOfActiveNodes << std::endl;
        int activeIdx = 0;
        int collisionIdx = 0;
        for (iterator.begin(); !iterator.isEnd(); iterator.next())
            if (iterator.value(nodeType) == NodeType::Active)
                iterator.value(m_numbering) = activeIdx++;
            else if (iterator.value(nodeType) == NodeType::Collision)
                iterator.value(m_numbering) = numOfActiveNodes - schurSize + collisionIdx++;
            else
                iterator.value(m_numbering) = -1;
        m_tensor.clear();
        m_tensor.resize(numOfActiveNodes);

        if (!schurSize) {
#if 0
            m_originalValue = new T[nnz];
            for (IntType i = 0; i < nnz; i++)
                m_originalValue[i] = m_pardiso.value[i];
#endif
        }
        else {
            m_originalValue = new T[schurSize * schurSize];
            m_schur = new T[schurSize * schurSize];
            m_pardiso.schur = m_schur;
        }

        m_rhs = new T[numOfActiveNodes];
        m_x = new T[numOfActiveNodes]();
    }

    template<class Discretization, class IntType>
    template<int elementNodesN>
    inline void SchurSolver<Discretization, IntType>::
        accumToTensor(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix, const std::array<IndexType, elementNodesN>& elementIndex) {
        using IteratorType = Iterator<NodeArrayType>;
        for (int i = 0; i < elementNodesN; i++) {
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


    template<class Discretization, class IntType>  // COURT added this debug subroutine to isolate Release only crash of physics thread
    template<int elementNodesN>
    inline void SchurSolver<Discretization, IntType>::
        accumToTensor_debug(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix, const std::array<IndexType, elementNodesN>& elementIndex) {
        using IteratorType = Iterator<NodeArrayType>;
        for (int i = 0; i < elementNodesN; i++) {
            int row = IteratorType::at(m_numbering, elementIndex[i]);
            if (row >= 0) {
                for (int j = 0; j < elementNodesN; j++) {
                    int col = IteratorType::at(m_numbering, elementIndex[j]);
                    if (col >= row) {
                        if (m_tensor[row].find(col) == m_tensor[row].end())
//                            ;
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
        updateTensor(
           const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix, 
          const std::array<IndexType, elementNodesN>& elementIndex) {
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

    template<class Discretization, class IntType>
    template<int elementNodesN>
    void SchurSolver<Discretization, IntType>::accumToPardiso(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix, const std::array<IndexType, elementNodesN>& elementIndex)
    {
        // need to check if repeated entry matters or not
        using IteratorType = Iterator<NodeArrayType>;
        for (int i = 0; i < elementNodesN; i++) {
            int row = IteratorType::at(m_numbering, elementIndex[i]);
            if (row >= 0) {
                //IntType idx = m_pardiso.rowIndex[row];
                for (int j = 0; j < elementNodesN; j++) {
                    int col = IteratorType::at(m_numbering, elementIndex[j]);
                    if (col >= row) {
                        bool found = false; // for debugging purposes
                        // maybe binary search here?
                        for (IntType jj = m_pardiso.rowIndex[row]; jj < m_pardiso.rowIndex[row + 1]; jj++)
                            if (col == m_pardiso.column[jj]) {
                                found = true;
                                m_pardiso.value[jj] -= stiffnessMatrix(i + 1, j + 1);
                            }
                        if (!found) {
                            std::cout << stiffnessMatrix << std::endl;
                            std::cout << "index = [ ";
                            for (int ii = 0; ii < elementNodesN; ii++)
                                std::cout << elementIndex[ii] << " ";
                            std::cout << "]" << std::endl;
                            std::cout << "row = " << row << " col = " << col << std::endl;
                            throw std::logic_error("entry (" + std::to_string(row) + " , " + std::to_string(col) + ") not found");
                        }
                        // if (m_tensor[row].find(col) == m_tensor[row].end())
                        //     m_tensor[row].insert(std::pair<int, T>(col, stiffnessMatrix(i + 1, j + 1)));
                        // else
                        //     m_tensor[row][col] += stiffnessMatrix(i + 1, j + 1);
                    }
                }
            }
        }

    }

#if 1
    template<class Discretization, class IntType>
    inline void SchurSolver<Discretization, IntType>::
        updatePardiso(const std::vector<Constraint>& collisionConstraints,
             const std::vector<CollisionSuture>& collisionSutures) {
#ifndef _WIN32
        LOG::SCOPE scope("SchurSolver::updatePardiso");
#endif
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

        for (int c = 0; c < collisionSutures.size(); c++) 
        if (collisionSutures[c].m_stiffness) {
            MATRIX_MXN<T> stiffnessMatrix;
            std::array<IndexType, elementNodes * 2> elementIndex;
            //CollisionSuture tmp = collisionSutures[c];
            //tmp.m_stiffness = 0;
            DiscretizationType::computeCollisionSutureTensor(stiffnessMatrix, elementIndex, collisionSutures[c]);
            updateTensor<elementNodes * 2>(stiffnessMatrix,
                elementIndex);
        }

        if (schurSize) {
            m_pardiso.factSchur();
        }
        else {
            m_pardiso.factorize();                      
        }
    }
#endif

    template<class Discretization, class IntType>
    inline void SchurSolver<Discretization, IntType>::
        computeTensor(
            const std::vector<ElementType>& elements,
            const std::vector<GradientMatrixType>& gradients,
            const std::vector<T>& restVol, const T mu,
            const std::vector<Suture>& sutures,
            const std::vector<InternodeConstraint>& microNodes
        ) {
#if TIMING
        auto startStamp = std::chrono::steady_clock::now();
        auto endStamp = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_second;
        double computeElTensorTime = 0;
        double accumElTensorTime = 0;
#endif

        for (int e = 0; e < elements.size(); e++) {
            MATRIX_MXN<T> stiffnessMatrix;
#if TIMING
            startStamp = std::chrono::steady_clock::now();
#endif
            DiscretizationType::computeElementTensor(stiffnessMatrix, gradients[e], -2 * mu * restVol[e]);
#if TIMING
            endStamp = std::chrono::steady_clock::now();
            elapsed_second = endStamp - startStamp;
            computeElTensorTime += elapsed_second.count();

            startStamp = std::chrono::steady_clock::now();
#endif
            accumToTensor<elementNodes>(stiffnessMatrix, DiscretizationType::getElementIndex(elements[e]));
#if TIMING
            endStamp = std::chrono::steady_clock::now();
            elapsed_second = endStamp - startStamp;
            accumElTensorTime += elapsed_second.count();
#endif
        }
#if TIMING
        LOG::cout << "        computeElTensor     Time : " << computeElTensorTime << std::endl;
        LOG::cout << "        accumElTensor       Time : " << accumElTensorTime << std::endl;
#endif
#if 0

        for (int c = 0; c < constraints.size(); c++) {
            MATRIX_MXN<T> stiffnessMatrix;
            DiscretizationType::computeConstraintTensor(stiffnessMatrix, constraints[c]);

            accumToTensor<elementNodes>(stiffnessMatrix,
                constraints[c].m_elementIndex);
        }
#endif

        for (int c = 0; c < sutures.size(); c++) {
            MATRIX_MXN<T> stiffnessMatrix;
            std::array<IndexType, elementNodes * 2> elementIndex;
            Suture tmp = sutures[c];
            tmp.m_stiffness = 0;
            DiscretizationType::computeSutureTensor(stiffnessMatrix, elementIndex, tmp);
            accumToTensor<elementNodes * 2>(stiffnessMatrix,
                elementIndex);
        }

        std::cout << "About to add microNode constraints in computeTensor()\n";  // COURT debug

        for (const auto& c : microNodes) {
            MATRIX_MXN<T> stiffnessMatrix;
            std::array<IndexType, d + 1> elementIndex;
            auto tmp = c;
            tmp.m_stiffness = 0;
            DiscretizationType::computeMicroNodeTensor(stiffnessMatrix, elementIndex, tmp);
            accumToTensor_debug<d + 1>(stiffnessMatrix,
                elementIndex);  // COURT created this routine to isolate Release only exception in accumToTensor()
        }

        std::cout << "Accumulated microNode constraints in computeTensor()\n";  // COURT debug

    }

    template<class Discretization, class IntType>
    void SchurSolver<Discretization, IntType>::computeTensor(const std::vector<ElementType>& elements, 
        const std::vector<GradientMatrixType>& gradients, 
        const std::vector<T>& restVol, 
        const std::vector<T>& muLow, 
        const std::vector<T>& muHigh, 
        const std::vector<Suture>& sutures,
        const std::vector<InternodeConstraint>& microNodes)
    {
        // only include things that will change the sparsity of stiffness matrix

#if TIMING
        auto startStamp = std::chrono::steady_clock::now();
        auto endStamp = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_second;
        double computeElTensorTime = 0;
        double accumElTensorTime = 0;
#endif

        for (int e = 0; e < elements.size(); e++) {
            MATRIX_MXN<T> stiffnessMatrix;

#if TIMING
            startStamp = std::chrono::steady_clock::now();
#endif
            DiscretizationType::computeElementTensor(stiffnessMatrix, gradients[e], -2 * (muLow[e] + muHigh[e]) * restVol[e]); // computeElementTensor
#if TIMING
            endStamp = std::chrono::steady_clock::now();
            elapsed_second = endStamp - startStamp;
            computeElTensorTime += elapsed_second.count();

            startStamp = std::chrono::steady_clock::now();
#endif
            accumToTensor<elementNodes>(stiffnessMatrix, DiscretizationType::getElementIndex(elements[e])); // accumToTensor

#if TIMING
            endStamp = std::chrono::steady_clock::now();
            elapsed_second = endStamp - startStamp;
            accumElTensorTime += elapsed_second.count();
#endif
        }

#if TIMING
        LOG::cout << "        computeElTensor     Time : " << computeElTensorTime << std::endl;
        LOG::cout << "        accumElTensor       Time : " << accumElTensorTime << std::endl;
#endif

        // It seems that I don't actually need to compute the matrix here, just need the sparsity
        for (int c = 0; c < sutures.size(); c++) {
            MATRIX_MXN<T> stiffnessMatrix;
            std::array<IndexType, elementNodes * 2> elementIndex;
            Suture tmp = sutures[c];
            tmp.m_stiffness = 0;
            DiscretizationType::computeSutureTensor(stiffnessMatrix, elementIndex, tmp);
            accumToTensor<elementNodes * 2>(stiffnessMatrix,
                elementIndex);
        }

        for (const auto& c:microNodes) {
            MATRIX_MXN<T> stiffnessMatrix;
            std::array<IndexType, d+1> elementIndex;
            auto tmp = c;
            tmp.m_stiffness = 0;
            DiscretizationType::computeMicroNodeTensor(stiffnessMatrix, elementIndex, tmp);
            accumToTensor_debug<d+1>(stiffnessMatrix,
                elementIndex);
        }

    }


    template<class Discretization, class IntType>
    inline void SchurSolver<Discretization, IntType>::initializePardiso(
        const std::vector<Constraint>& constraints,
        const std::vector<Suture>& sutures,
        const std::vector<Constraint>& fakeSutures, 
        const std::vector<InternodeConstraint>& microNodes
    ) {

        IntType nnz = 0;
        for (int i = 0; i < m_tensor.size(); i++)
            nnz += (IntType)m_tensor[i].size();

        LOG::cout << "nnz = " << nnz << std::endl;
        m_pardiso.initialize((IntType)m_tensor.size(), nnz, schurSize);

        m_pardiso.rowIndex[0] = 0;
        for (int i = 0; i < m_pardiso.n; i++)
            m_pardiso.rowIndex[i + 1] = m_pardiso.rowIndex[i] + (IntType)m_tensor[i].size();

        // const IntType& n = m_pardiso.n;

        size_t idx = 0;
        for (const auto& r : m_tensor)
            for (const auto& e : r) {
                m_pardiso.column[idx] = e.first;
                //m_pardiso.value[idx] = -e.second;
                idx++;
            }

        if (schurSize)
            m_pardiso.schur = m_schur;
        m_pardiso.symbolicFact();

        factPardiso(constraints, sutures, fakeSutures, microNodes);

        if (schurSize) {
            for (IntType i = 0; i < schurSize * schurSize; i++)
                m_originalValue[i] = m_pardiso.schur[i];
            m_pardiso.factSchur();
        }

        // m_tensor.resize(0);
    }

    template<class Discretization, class IntType>
    void SchurSolver<Discretization, IntType>::factPardiso(const std::vector<Constraint>& constraints, const std::vector<Suture>& sutures, const std::vector<Constraint>& fakeSutures, const std::vector<InternodeConstraint>& microNodes)
    {
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
#if 1
        for (const auto& c:microNodes)
            if (c.m_stiffness != 0) {
                MATRIX_MXN<T> stiffnessMatrix;
                std::array<IndexType, d+1> elementIndex;
                DiscretizationType::computeMicroNodeTensor(stiffnessMatrix, elementIndex, c);
                accumToPardiso<d+1>(stiffnessMatrix,
                    elementIndex);
            }
#endif
        // dumper::writeCSRbyte(m_pardiso.n, m_pardiso.rowIndex, m_pardiso.column, m_pardiso.value, m_pardiso.n, "new_i.txt", "new_a.txt");

        m_pardiso.numericFact();
        /*
        if (schurSize) {
            for (IntType i = 0; i < schurSize * schurSize; i++)
                m_Sigma1[i] = m_pardiso.schur[i] - m_A22[i];
        }
        */
    }

}

namespace PhysBAM {
    template struct SchurSolver<TetrahedralDiscretization<std::vector<VECTOR<float, 3>>>, int>;
}
