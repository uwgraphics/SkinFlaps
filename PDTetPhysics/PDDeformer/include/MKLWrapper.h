#pragma once

#include "mkl_types.h"

template<class T, class IntType> struct PardisoPolicy;
    template<class T> struct PardisoPolicy<T, int> {
        using IntType = int;
        static inline IntType exec(void** pt, const IntType maxfct, const IntType mnum, const IntType mtype, const IntType phase, const IntType n, T* a, IntType* ia, IntType* ja, IntType* perm, const IntType nrhs, IntType* iparm, const IntType msglvl, T* b, T* x) {
            IntType error;
            pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);
            return error;
        }
    };

    template<class T> struct PardisoPolicy<T, long long int> {
        using IntType = long long int;
        static inline IntType exec(void** pt, const IntType maxfct, const IntType mnum, const IntType mtype, const IntType phase, const IntType n, T* a, IntType* ia, IntType* ja, IntType* perm, const IntType nrhs, IntType* iparm, const IntType msglvl, T* b, T* x) {
            IntType error;
            pardiso_64(pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, perm, &nrhs, iparm, &msglvl, b, x, &error);
            return error;
        }
    };



template<class T> struct LAPACKPolicy;
template<> struct LAPACKPolicy<double> {
        static constexpr int matrix_order=LAPACK_ROW_MAJOR;
        static constexpr char uplo = 'U';

        using T = double;


        static inline int fact(const int m, T* a) {
            return (int)LAPACKE_dpotrf(matrix_order,uplo,m,a,m);
        }

        static inline int solve(const int m, const int nrhs, const T* a, T* b) {
            return (int)LAPACKE_dpotrs(matrix_order,uplo,m,nrhs,a,m,b,nrhs);
        }
    };

    template<> struct LAPACKPolicy<float> {
        static constexpr int matrix_order=LAPACK_ROW_MAJOR;
        static constexpr char uplo = 'U';

        using T = float;

        static inline int fact(const int m, T* a) {
            return (int)LAPACKE_spotrf(matrix_order,uplo,m,a,m);
        }

        static inline int solve(const int m, const int nrhs, const T* a, T* b) {
            return (int)LAPACKE_spotrs(matrix_order,uplo,m,nrhs,a,m,b,nrhs);
        }
    };

template<class T> struct CBLASPolicy;
template<> struct CBLASPolicy<double> {
    static constexpr CBLAS_LAYOUT matrix_order = CblasRowMajor;
    static constexpr CBLAS_UPLO uplo = CblasUpper;
    using T = double;
    static inline void mutiplyAdd(T* result, const int n, const T alpha, const T* a, const T* x, const T beta) {
        cblas_dsymv (matrix_order, uplo, n, alpha, a, n, x, 1, beta, result, 1);
    }
};

template<> struct CBLASPolicy<float> {
    static constexpr CBLAS_LAYOUT matrix_order = CblasRowMajor;
    static constexpr CBLAS_UPLO uplo = CblasUpper;
    using T = float;
    static inline void mutiplyAdd(T* result, const int n, const T alpha, const T* a, const T* x, const T beta) {
        cblas_ssymv (matrix_order, uplo, n, alpha, a, n, x, 1, beta, result, 1);
    }
};
