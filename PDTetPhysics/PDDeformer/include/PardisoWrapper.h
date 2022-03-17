//#####################################################################
// Copyright (c) 2019, Eftychios Sifakis, Yutian Tao, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################
#pragma once
#include <iostream>


template <class T, class IntType_> struct PardisoWrapper {
    using IntType = IntType_;

    IntType  n = 0; // dimension of the matrix
    int      m = 0;// number of nodes in the schur complement part
    IntType *schurNodes = nullptr;

    IntType *rowIndex = nullptr;
    IntType *column = nullptr;
    T       *value = nullptr;
    T       *schur = nullptr;

    IntType mtype = 0;
    int     nrhs = 0;
    void   *pt[64] = {
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr,
        nullptr, nullptr, nullptr, nullptr
    };      // Internal solver memory pointer pt
    IntType iparm[64]{}; // Pardiso control parameters.
    IntType maxfct=0, mnum=0, msglvl=0;

    void initialize(const IntType _n, const IntType _nnz, const IntType _m = 0);

    void  factSchur();

    void factorize() {
        symbolicFact(); // symFact
        numericFact(); // numFact
    }

    void symbolicFact();
    void numericFact();

    void releasePardisoInternal();
    void deallocate();

    void forwardSubstitution(T* const _rhs, T* const _x);
    void diagSolve(T* const _rhs, T* const _x);
    void backwardSubstitution(T* const _rhs, T* const _x);

    void printSchur () {
        std::cout<<std::endl;
        for (int i=0; i<m; i++) {
            for (int j=0; j<m; j++)
                std::cout<<schur[i*m+j]<<" ";
            std::cout<<std::endl;
        }
    }

};
