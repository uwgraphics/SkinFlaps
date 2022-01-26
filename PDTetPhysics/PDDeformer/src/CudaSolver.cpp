#include "CudaSolver.h"

#if TIMING
#include <chrono>
#endif

namespace PhysBAM {
    template <class Discretization, class IntType>
    void CudaSolver<Discretization, IntType>::initialize(const NodeArrayType& nodeType) {
        using IteratorType = Iterator<NodeArrayType>;
#ifndef _WIN32
        LOG::SCOPE scope("CudaSolver::initialize()");
#endif
        IteratorType iterator(nodeType);
        iterator.resize(m_numbering);

        matrixSize = 0;
        schurSize = 0;
        for (iterator.begin(); !iterator.isEnd(); iterator.next())
            if (iterator.value(nodeType) == PDSimulation::ActiveNode) {
                matrixSize++;
            }
            else if (iterator.value(nodeType) == PDSimulation::CollisionNode) {
                matrixSize++;
                schurSize++;
            }
        LOG::cout << "schursize  = " << schurSize << std::endl;
        LOG::cout << "matrixsize  = " << matrixSize << std::endl;
        if (schurSize == 0) {
            LOG::cout << "Must have collision nodes" << std::endl;
            exit(1);
        }

        int activeIdx = 0;
        int collisionIdx = 0;
        for (iterator.begin(); !iterator.isEnd(); iterator.next())
            if (iterator.value(nodeType) == PDSimulation::ActiveNode)
                iterator.value(m_numbering) = activeIdx++;
            else if (iterator.value(nodeType) == PDSimulation::CollisionNode)
                iterator.value(m_numbering) = matrixSize - schurSize + collisionIdx++;
            else
                iterator.value(m_numbering) = -1;

        m_tensor.clear();
        m_tensor.resize(matrixSize);
        if (schurSize) {
            m_A22 = new T[schurSize * schurSize]();
            m_f2 = new T[schurSize * d]();
            m_schur = new T[schurSize * schurSize];
            m_Sigma1 = new T[schurSize * schurSize];
        }

        rhs = new T[matrixSize];
        x = new T[matrixSize]();
    }

    template <class Discretization, class IntType>
    template<int elementNodesN>
    void CudaSolver<Discretization, IntType>::
        accumToTensor(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix,
            const std::array<IndexType, elementNodesN>& elementIndex) {
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


    template <class Discretization, class IntType>
    template <int elementNodesN>
    void CudaSolver<Discretization, IntType>::
        accumToPardiso(const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix,
            const std::array<IndexType, elementNodesN>& elementIndex) {
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

    template <class Discretization, class IntType>
    template <int elementNodesN>
    void CudaSolver<Discretization, IntType>::
        updateTensor(T* const result,
            const PhysBAM::MATRIX_MXN<T>& stiffnessMatrix,
            const std::array<IndexType, elementNodesN>& elementIndex) {
        // only cared about upper triangular in row major
        //LOG::SCOPE scope("CudaSolver::updateTensor");
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

                    assert(row >= matrixSize - schurSize);
                    result[(row - matrixSize + schurSize) * schurSize + col - matrixSize + schurSize] -= stiffnessMatrix(i + 1, j + 1);
                }
            }
        }
    }

    /*
    template <class StateVariable, class IntType>
    void CudaSolver<StateVariable, IntType>::updatePardiso(const std::vector<Constraint> &collisionConstraints)
    {
        // add A22 part to Sigma1 part
        if (schurSize)
            for (int i = 0; i < schurSize*schurSize; i++)
                m_pardiso.schur[i] = m_Sigma1[i]+m_A22[i];
        for (int c = 0; c < collisionConstraints.size(); c++) {
            auto &constraint = collisionConstraints[c];
            if (constraint.m_stiffness != 0) {
                MATRIX_MXN<T> stiffnessMatrix;
                computeConstraintTensor(stiffnessMatrix, constraint);
                updateTensor<elementNodes>(m_schur, stiffnessMatrix, constraint.m_elementIndex);
            }
        }

        if (schurSize) {
            m_pardiso.factSchur();
        }
    }
    */

    template <class Discretization, class IntType>
    void CudaSolver<Discretization, IntType>::
        updateCuda(const std::vector<Constraint>& collisionConstraints,
            const std::vector<CollisionSuture>& collisionSutures
        ) {

        const int pc = (int)collisionConstraints.size();
        const int ps = (int)collisionSutures.size();

        bool needUpdate = false;
        for (int i = 0; i < collisionConstraints.size(); i++) {
            if (m_valPtrD_h[i] != collisionConstraints[i].m_stiffness)
                needUpdate = true;
            m_valPtrD_h[i] = collisionConstraints[i].m_stiffness;
        }

        for (int i = 0; i < collisionSutures.size(); i++) {
            if (m_valPtrD_h[i + pc] != collisionSutures[i].m_stiffness) {
                needUpdate = true;
                m_valPtrD_h[i + pc] = collisionSutures[i].m_stiffness;
            }
        }
        {
#if TIMING
            auto stamp1 = std::chrono::steady_clock::now();
#endif

#pragma omp parallel for
            for (int i = 0; i < ps; i++) {
                const auto& suture = collisionSutures[i];
                std::pair<int, T> indexWeightPair[(d + 1) * 2];
                if (suture.m_stiffness) {
                    for (int j = 0; j < d + 1; j++) {
                        indexWeightPair[j].first = m_numbering[suture.m_elementIndex1[j]] - (m_pardiso.n - schurSize);
                        indexWeightPair[j].second = suture.m_weights1[j];
                    }

                    for (int j = 0; j < d + 1; j++) {
                        indexWeightPair[j + d + 1].first = m_numbering[suture.m_elementIndex2[j]] - (m_pardiso.n - schurSize);
                        indexWeightPair[j + d + 1].second = -suture.m_weights2[j];
                    }
                    std::sort(indexWeightPair, indexWeightPair + (d + 1) * 2);
                }
                else {
                    for (int j = 0; j < 2 * (d + 1); j++) {
                        indexWeightPair[j].first = j;
                        indexWeightPair[j].second = 0;
                    }
                }

                for (int j = 0; j < (d + 1) * 2; j++)
                    if (m_colIdxPtrW_h[pc * (d + 1) + i * 2 * (d + 1) + j] != indexWeightPair[j].first ||
                        m_valPtrW_h[pc * (d + 1) + i * 2 * (d + 1) + j] != indexWeightPair[j].second)
                    {
                        needUpdate = true;
                        m_colIdxPtrW_h[pc * (d + 1) + i * 2 * (d + 1) + j] = indexWeightPair[j].first;
                        m_valPtrW_h[pc * (d + 1) + i * 2 * (d + 1) + j] = indexWeightPair[j].second;
                    }

            }

#if TIMING
            auto stamp2 = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_second = stamp2 - stamp1;
            LOG::cout << "        Time for update W_h non const part " << elapsed_second.count() << std::endl;
#endif

        }


        if (needUpdate)
        {
            m_cuda.rankKUpdate(m_valPtrD_h, m_colIdxPtrW_h, m_valPtrW_h);
            m_cuda.factorize();
        }
    }

    template <class Discretization, class IntType>
    void CudaSolver<Discretization, IntType>::
        computeE2Tensor(const std::vector<ElementType>& elements,
            const std::vector<PDSimulation::ElementFlag>& flags,
            const std::vector<GradientMatrixType>& gradients,
            const std::vector<T>& restVol, const T mu
        )
    {
#ifndef _WIN32
        LOG::SCOPE scope("CudaSolver::computeE2Tensor()");
#endif
        for (int e = 0; e < elements.size(); e++)
            if (flags[e] == PDSimulation::CollisionEl)
            {
                MATRIX_MXN<T> stiffnessMatrix;
                DiscretizationType::computeElementTensor(stiffnessMatrix, gradients[e], -2 * mu * restVol[e]);
                updateTensor<elementNodes>(m_A22, stiffnessMatrix, DiscretizationType::getElementIndex(elements[e]));
            }
    }

    template <class Discretization, class IntType>
    void CudaSolver<Discretization, IntType>::
        computeTensor(
            const std::vector<ElementType>& elements,
            const std::vector<GradientMatrixType>& gradients,
            const std::vector<T>& restVol, const T muLow, const T muHigh,
            //const std::vector<Constraint>& constraints,
            const std::vector<Suture>& sutures
        ) {
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
            DiscretizationType::computeElementTensor(stiffnessMatrix, gradients[e], -2 * (muLow+muHigh) * restVol[e]); // computeElementTensor
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
        /*
         for (int c = 0; c < constraints.size(); c++) {
             MATRIX_MXN<T> stiffnessMatrix;
             computeConstraintTensor(stiffnessMatrix, constraints[c]);

            accumToTensor<elementNodes>(stiffnessMatrix,
              constraints[c].m_elementIndex);
         }
         */
        for (int c = 0; c < sutures.size(); c++) {
            MATRIX_MXN<T> stiffnessMatrix;
            std::array<IndexType, elementNodes * 2> elementIndex;
            Suture tmp = sutures[c];
            tmp.m_stiffness = 0;
            DiscretizationType::computeSutureTensor(stiffnessMatrix, elementIndex, tmp);
            accumToTensor<elementNodes * 2>(stiffnessMatrix,
                elementIndex);
        }

    }

    template <class Discretization, class IntType>
    void CudaSolver<Discretization, IntType>::
        computeTensor(
            const std::vector<ElementType>& elements,
            const std::vector<GradientMatrixType>& gradients,
            const std::vector<T>& restVol, const std::vector<T>& muLow, const std::vector<T>& muHigh,
            //const std::vector<Constraint>& constraints,
            const std::vector<Suture>& sutures
        ) {
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
        /*
         for (int c = 0; c < constraints.size(); c++) {
             MATRIX_MXN<T> stiffnessMatrix;
             computeConstraintTensor(stiffnessMatrix, constraints[c]);

            accumToTensor<elementNodes>(stiffnessMatrix,
              constraints[c].m_elementIndex);
         }
         */
        for (int c = 0; c < sutures.size(); c++) {
            MATRIX_MXN<T> stiffnessMatrix;
            std::array<IndexType, elementNodes * 2> elementIndex;
            Suture tmp = sutures[c];
            tmp.m_stiffness = 0;
            DiscretizationType::computeSutureTensor(stiffnessMatrix, elementIndex, tmp);
            accumToTensor<elementNodes * 2>(stiffnessMatrix,
                elementIndex);
        }

    }

    template <class Discretization, class IntType>
    void CudaSolver<Discretization, IntType>::
        initializePardiso(
            const std::vector<Constraint>& constraints,
            const std::vector<Suture>& sutures,
            const std::vector<Constraint>& fakeSutures
        ) {
#ifndef _WIN32
        LOG::SCOPE scope("CudaSolver::initializePardiso()");
#endif
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
                // m_pardiso.value[idx] = -e.second;
                idx++;
            }
        m_pardiso.schur = m_schur;
        m_pardiso.symbolicFact();

        factPardiso(constraints, sutures, fakeSutures);
        /*
        for (int c = 0; c < constraints.size(); c++)
            if (constraints[c].m_stiffness != 0) {
                MATRIX_MXN<T> stiffnessMatrix;
                DiscretizationType::computeConstraintTensor(stiffnessMatrix, constraints[c]);
                accumToPardiso<elementNodes>(stiffnessMatrix,
                    constraints[c].m_elementIndex);
            }

        for (int c = 0; c < sutures.size(); c++)
            if (sutures[c].m_stiffness != 0) {
                MATRIX_MXN<T> stiffnessMatrix;
                std::array<IndexType, elementNodes * 2> elementIndex;
                DiscretizationType::computeSutureTensor(stiffnessMatrix, elementIndex, sutures[c]);
                accumToPardiso<elementNodes * 2>(stiffnessMatrix,
                    elementIndex);
            }

        m_pardiso.factorize();
        if (schurSize) {
            for (IntType i=0; i<schurSize*schurSize; i++)
                m_Sigma1[i] = m_pardiso.schur[i] - m_A22[i];
        }
        */
        //m_tensor.resize(0);

    }

    template <class Discretization, class IntType>
    void CudaSolver<Discretization, IntType>::
        initializeCuda(
            const std::vector<Constraint>& collisionConstraints,
            const std::vector<CollisionSuture>& collisionSutures
        ) {
        if (schurSize == 0) throw std::logic_error("SchurSize can't be 0");

        int m = schurSize;
        int pc = (int)collisionConstraints.size();
        int ps = (int)collisionSutures.size();
        LOG::cout << "pc = " << pc << std::endl;
        LOG::cout << "ps = " << ps << std::endl;
        m_rowPtrW_h = new int[pc + ps + 1];
        m_colIdxPtrW_h = new int[pc * (d + 1) + ps * 2 * (d + 1)];
        m_valPtrW_h = new T[pc * (d + 1) + ps * 2 * (d + 1)];
        m_valPtrS_h = new T[m * m];
        m_valPtrD_h = new T[pc + ps];

        m_rowPtrW_h[0] = 0;
        for (int i = 0; i < pc; i++)
            m_rowPtrW_h[i + 1] = (i + 1) * (d + 1);
        for (int i = 0; i < ps; i++)
            m_rowPtrW_h[pc + i + 1] = pc * (d + 1) + (i + 1) * (d + 1) * 2;

        {
#ifndef _WIN32
            LOG::SCOPE scope("CudaSolver::initializeCuda - prepare W_h const part");
#endif

#if TIMING
            auto stamp1 = std::chrono::steady_clock::now();
#endif
            for (int i = 0; i < pc; i++) {
                const auto& constraint = collisionConstraints[i];
                std::pair<int, T> indexWeightPair[d + 1];
                for (int j = 0; j < d + 1; j++) {
                    indexWeightPair[j].first = m_numbering[constraint.m_elementIndex[j]] - (m_pardiso.n - m);
                    indexWeightPair[j].second = constraint.m_weights[j];
                }

                std::sort(indexWeightPair, indexWeightPair + d + 1);

                for (int j = 0; j < d + 1; j++) {
                    m_colIdxPtrW_h[i * (d + 1) + j] = indexWeightPair[j].first;
                    m_valPtrW_h[i * (d + 1) + j] = indexWeightPair[j].second;
                }
            }

            for (int i = 0; i < ps; i++) {
                const auto& suture = collisionSutures[i];
                std::pair<int, T> indexWeightPair[(d + 1) * 2];
                if (suture.m_stiffness) {
                    for (int j = 0; j < d + 1; j++) {
                        indexWeightPair[j].first = m_numbering[suture.m_elementIndex1[j]] - (m_pardiso.n - m);
                        indexWeightPair[j].second = suture.m_weights1[j];
                    }

                    for (int j = 0; j < d + 1; j++) {
                        indexWeightPair[j + d + 1].first = m_numbering[suture.m_elementIndex2[j]] - (m_pardiso.n - m);
                        indexWeightPair[j + d + 1].second = -suture.m_weights2[j];
                    }

                    std::sort(indexWeightPair, indexWeightPair + (d + 1) * 2);
                }
                else {
                    for (int j = 0; j < 2 * (d + 1); j++) {
                        indexWeightPair[j].first = j;
                        indexWeightPair[j].second = 0;
                    }
                }

                for (int j = 0; j < (d + 1) * 2; j++) {
                    m_colIdxPtrW_h[pc * (d + 1) + i * 2 * (d + 1) + j] = indexWeightPair[j].first;
                    m_valPtrW_h[pc * (d + 1) + i * 2 * (d + 1) + j] = indexWeightPair[j].second;
                }
            }

#if TIMING
            auto stamp2 = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed_second = stamp2 - stamp1;
            LOG::cout << "Time for prepare W const part " << elapsed_second.count() << std::endl;
#endif

        }

        for (int i = 0; i < m * m; i++)
            m_valPtrS_h[i] = m_pardiso.schur[i];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                if (i > j)
                    m_valPtrS_h[i * m + j] = m_valPtrS_h[j * m + i];

        for (int i = 0; i < pc; i++)
            m_valPtrD_h[i] = collisionConstraints[i].m_stiffness;
        for (int i = 0; i < ps; i++)
            m_valPtrD_h[i + pc] = collisionSutures[i].m_stiffness;

        m_cuda.initialize(m, pc, ps, m_rowPtrW_h, m_colIdxPtrW_h, m_valPtrW_h, m_valPtrS_h);
        {
            // LOG::SCOPE scope("RankKUpdate");
            m_cuda.rankKUpdate(m_valPtrD_h, m_colIdxPtrW_h, m_valPtrW_h);
        }
        m_cuda.factorize();
    }

    template <class Discretization, class IntType>
    void CudaSolver<Discretization, IntType>::
        reInitializeCuda(const std::vector<Constraint>& collisionConstraints,
            const std::vector<CollisionSuture>& collisionSutures) {

        for (int i = 0; i < schurSize * schurSize; i++)
            m_valPtrS_h[i] = m_pardiso.schur[i];
        for (int i = 0; i < schurSize; i++)
            for (int j = 0; j < schurSize; j++)
                if (i > j)
                    m_valPtrS_h[i * schurSize + j] = m_valPtrS_h[j * schurSize + i];
        for (int i = 0; i < collisionConstraints.size(); i++)
            m_valPtrD_h[i] = collisionConstraints[i].m_stiffness;
        for (int i = 0; i < collisionSutures.size(); i++)
            m_valPtrD_h[i + collisionConstraints.size()] = collisionSutures[i].m_stiffness;

        m_cuda.reInitialize(schurSize, m_valPtrS_h);
        {
            // LOG::SCOPE scope("RankKUpdate");
            m_cuda.rankKUpdate(m_valPtrD_h, m_colIdxPtrW_h, m_valPtrW_h);
        }
        m_cuda.factorize();
    }

    template <class Discretization, class IntType>
    void CudaSolver<Discretization, IntType>::deallocate() {
        if (m_A22) {
            delete[] m_A22;
            m_A22 = nullptr;
        }
        if (m_f2) {
            delete[] m_f2;
            m_f2 = nullptr;
        }
        if (m_schur) {
            delete[] m_schur;
            m_schur = nullptr;
            m_pardiso.schur = nullptr;
        }
        if (m_Sigma1) {
            delete[] m_Sigma1;
            m_Sigma1 = nullptr;
        }
        if (rhs) {
            delete[] rhs;
            rhs = nullptr;
        }
        if (x) {
            delete[] x;
            x = nullptr;
        }
    }

    template<class Discretization, class IntType>
    void CudaSolver<Discretization, IntType>::factPardiso(
        const std::vector<Constraint>& constraints,
        const std::vector<Suture>& sutures,
        const std::vector<Constraint>& fakeSutures
    )
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

        // dumper::writeCSRbyte(m_pardiso.n, m_pardiso.rowIndex, m_pardiso.column, m_pardiso.value, m_pardiso.n, "new_i.txt", "new_a.txt");

        m_pardiso.numericFact();

        if (schurSize) {
            for (IntType i = 0; i < schurSize * schurSize; i++)
                m_Sigma1[i] = m_pardiso.schur[i] - m_A22[i];
        }
    }

}

namespace PhysBAM {
    template struct CudaSolver<TetrahedralDiscretization<std::vector<VECTOR<float, 3>>>, int>;
}
