//#####################################################################
// Copyright (c) 2019, Eftychios Sifakis, Yutian Tao, Qisi Wang
// Distributed under the FreeBSD license (see license.txt)
//#####################################################################

#include <mkl.h>
#include "PardisoWrapper.h"
#include "MKLWrapper.h"

#if TIMING
#include "chrono" // for timing only
#endif


template<class T, class IntType>
void PardisoWrapper<T, IntType>::initialize(const IntType _n, const IntType _nnz, const IntType _m) {
    //PhysBAM::LOG::SCOPE scope("PardisoWrapper::initialize()");

        n = _n;
        m = (int)_m;

        // allocate spaces
        rowIndex = new IntType[n+1];
        column = new IntType[_nnz];
        value = new T[_nnz];
        // initialize schur and schurNodes
        if (m) {
            schurNodes = new IntType[n];
        }

        mtype = 2; /* Real symmetric positive definite matrix */
                   // IntType mtype = -2;       /* Real symmetric (maybe indefinite) matrix */
        nrhs = 1; /* Number of right hand sides. */

        //      Auxiliary variables.

        for (int i = 0; i < 64; i++) {
            iparm[i] = 0;
        }

        iparm[0] = 1; /* No solver default */
        iparm[1] = 2; /* Use omp */
        //iparm[3] = 0; /* No iterative-direct algorithm */
        //iparm[4] = 0; /* No user fill-in reducing permutation, return permutaion in perm (could change this, but will ignore iparm[2]) */
          //iparm[5] = 0;  /* write result to x*/
//iparm[6] = 0;
          //      iparm[7] = 0;  /* Max numbers of iterative refinement steps */
//iparm[8] = 0;  /* Not in use */
        iparm[9] = 8;  /* not needed for mtype = 2 */
        iparm[10] = 0; /* not needed for mtype = 2 */
//iparm[11] = 0; /* Not transpose or conjugate A */
iparm[12] = 0;  /* Maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm[12] = 1
                           in case of inappropriate accuracy */
iparm[13] = 0; /* Output: Number of perturbed pivots (not important for mtype = 2) */
//iparm[14] = 0; /* Output: Peak memory on symbolic fact */
//      iparm[15] = 0; /* Output: Permanent memory on symbolic fact and solve */
//iparm[16] = 0; /* Output: Total mem consumption */
        // total peak mem = max(iparm[14],iparm[15])+iparm[16]
        iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
        iparm[18] = -1; /* Output: Mflops for LU factorization */
//iparm[19] = 0;  /* Output: Numbers of CG Iterations */
//      iparm[23] = 1;  /* Two-level factorization*/
//      iparm[24] = 2;  /*parallel solve*/
//      iparm[26] = 1;  /* Check matrix for errors */
        iparm[27] = std::is_same<T, float>::value;
        /* float or double precision */
        iparm[34] = 1; /*1 for 0-based index*/
        iparm[35] = m?1:0;        /* Use Schur complement */ // using 2 somehow stuck on the numerical fact phase for some cases
        // iparm[36]=1; /* Use Blocked CSR format for input, need to have iparm[12] = 0 */
        maxfct = 1; /* Maximum number of numerical factorizations. */
        mnum = 1;   /* Number of matrix */
        msglvl = 0; /* 0:No statistical information print*/

        /* -------------------------------------------------------------------- */
        /* .. Initialize the internal solver memory pointer. This is only */
        /* necessary for the FIRST call of the PARDISO solver. */
        /* -------------------------------------------------------------------- */
        for (int i = 0; i < 64; i++) {
            pt[i] = nullptr;
        }
    }

template<class T, class IntType>
void PardisoWrapper<T, IntType>::factSchur() {
        //PhysBAM::LOG::SCOPE scope("PardisoWrapper::factSchur()");
        if (m) {
            IntType info = LAPACKPolicy<T>::fact(m, schur);
            if(info != 0) {
                std::cerr<<"info after LAPACKE_dspotrf = "<<info<<std::endl;
                exit(22);
            }
        }
    }

template<class T, class IntType>
void PardisoWrapper<T, IntType>::symbolicFact() {
#if TIMING
    auto startStamp = std::chrono::steady_clock::now();
    auto endStamp = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_second;
    startStamp = std::chrono::steady_clock::now();
#endif

    IntType error;
    T ddum;       /* Scalar dummy */
    IntType idum; /* Integer dummy. */
    IntType phase = 11;

    if (m) {
        for (IntType i=0; i<n-m; i++)
            schurNodes[i] = 0;
        for (IntType i=n-m; i<n; i++)
            schurNodes[i] = 1;

        error = PardisoPolicy<T, IntType>::exec(pt, maxfct, mnum, mtype, phase,
                                                n, value, rowIndex, column, schurNodes, nrhs,
                                                iparm, msglvl, &ddum, &ddum);
    } else {
        error = PardisoPolicy<T, IntType>::exec(pt, maxfct, mnum, mtype, phase, n, value, rowIndex, column, &idum, nrhs, iparm, msglvl, &ddum, &ddum);
    }

    if ( error != 0 ) {
        std::cerr<<"ERROR during symbolic factorization: "<<error<<std::endl;
        exit ((int)phase);
    }
#if TIMING
    endStamp = std::chrono::steady_clock::now();
    elapsed_second = endStamp - startStamp;
    std::cout<<std::endl<<"            Symbolic Fact           Time : "<<elapsed_second.count()<<std::endl;
#endif
}

template<class T, class IntType>
void PardisoWrapper<T, IntType>::numericFact() {
#if TIMING
    auto startStamp = std::chrono::steady_clock::now();
    auto endStamp = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_second;
    startStamp = std::chrono::steady_clock::now();
#endif

    IntType error;
    T ddum;       /* Scalar dummy */
    IntType idum; /* Integer dummy. */
    IntType phase = 22;

    if (m) {
        error = PardisoPolicy<T, IntType>::exec(pt, maxfct, mnum, mtype, phase,
                                                n, value, rowIndex, column, schurNodes, nrhs,
                                                iparm, msglvl, &ddum, schur);
    } else {
#if TIMING
        auto start = std::chrono::steady_clock::now();
#endif
        error = PardisoPolicy<T, IntType>::exec(pt, maxfct, mnum, mtype, phase, n, value, rowIndex, column, &idum, nrhs, iparm, msglvl, &ddum, &ddum);
#if TIMING
        auto end = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout<<"Numeric Fact Time:  "<<elapsed_seconds.count()<<" s"<<std::endl;
#endif
    }

    if ( error != 0 ) {
        std::cerr<<"ERROR during numerical factorization: "<<error<<std::endl;
        exit ((int)phase);
    }
#if TIMING
    endStamp = std::chrono::steady_clock::now();
    elapsed_second = endStamp - startStamp;

    std::cout<<std::endl<<"            Numeric Fact            Time : "<<elapsed_second.count()<<std::endl;
#endif
}

template<class T, class IntType>
void PardisoWrapper<T, IntType>::releasePardisoInternal() {
        // Termination and release of memory.
        IntType phase = -1; /* Release internal memory. */
        IntType error;
        IntType idum;
        T ddum;
        error = PardisoPolicy<T, IntType>::exec(pt, maxfct, mnum, mtype, phase, n, &ddum, rowIndex, column, &idum, nrhs, iparm, msglvl, &ddum, &ddum);
        if ( error != 0 ) {
            std::cerr<<"ERROR during release phase "<<phase<<": "<<error<<std::endl;
            exit ((int)phase);
        }
    }

template<class T, class IntType>
void PardisoWrapper<T, IntType>::deallocate() {
        if (value) {
            delete[] value;
            value = nullptr;
        }
        if (column) {
            delete[] column;
            column = nullptr;
        }
        if (rowIndex) {
            delete[] rowIndex;
            rowIndex = nullptr;
        }

		if (m && schurNodes) {
			delete[] schurNodes;
			schurNodes = nullptr;
		}
    }

template<class T, class IntType>
void PardisoWrapper<T, IntType>::forwardSubstitution(T* const _rhs, T* const _x) {
    IntType error;
    const IntType phase = 331;
    iparm[7] = 0; /* Max numbers of iterative refinement steps. */

    IntType idum;

    if (m) {
        error = PardisoPolicy<T, IntType>::exec(pt, maxfct, mnum, mtype, phase, n, value, rowIndex, column, schurNodes, nrhs, iparm, msglvl, _rhs, _x);

        if ( error != 0 ) {
            std::cerr<<"ERROR during solution phase "<<phase<<": "<<error<<std::endl;
            exit (phase);
        }

    } else {

        error = PardisoPolicy<T, IntType>::exec(pt, maxfct, mnum, mtype, phase, n, value, rowIndex, column, &idum, nrhs, iparm, msglvl, _rhs, _x);
        if (error != 0) {
            std::cerr<<"ERROR during solution: "<<error<<std::endl;
            exit(331);
        }
    }
}

template<class T, class IntType>
void PardisoWrapper<T, IntType>::diagSolve(T* const _rhs, T* const _x) {
        // this may not support 64 bit ints
        const IntType phase = 332;
        if (m) {

            IntType info = LAPACKPolicy<T>::solve(m,nrhs,schur,&_rhs[n-m]);
            if (info != 0)
            {
                std::cerr<<"info after LAPACKE_dspotrs = "<<info<<std::endl;
                exit(phase);
            }
            for (IntType i = 0; i < n; i++ ) {
                _x[i] = _rhs[i];
            }

        } else {
            IntType error;
            iparm[7] = 0; /* Max numbers of iterative refinement steps. */
            IntType idum;
            // T ddum;

            error = PardisoPolicy<T, IntType>::exec(pt, maxfct, mnum, mtype, phase, n, value, rowIndex, column, &idum, nrhs, iparm, msglvl, _rhs,
                    _x);
            if (error != 0) {
                std::cerr<<"ERROR during solution: "<<error<<std::endl;
                exit(33);
            }

        }
    }


template<class T, class IntType>
void PardisoWrapper<T, IntType>::backwardSubstitution(T* const _rhs, T* const _x) {
        const IntType phase = 333;
        IntType error;
        IntType idum;
        iparm[7] = 0; /* Max numbers of iterative refinement steps. */

        if (m) {
            error = PardisoPolicy<T, IntType>::exec(pt, maxfct, mnum, mtype, phase,
                    n, value, rowIndex, column, schurNodes, nrhs,
                    iparm, msglvl, _rhs, _x);
            if ( error != 0 )
            {
                std::cerr<<"ERROR during solution phase "<<phase<<": "<<error<<std::endl;
                exit (phase);
            }
        } else {
            error = PardisoPolicy<T, IntType>::exec(pt, maxfct, mnum, mtype, phase, n, value, rowIndex, column, &idum, nrhs, iparm, msglvl, _rhs,
                    _x);

            if (error != 0) {
                std::cerr<<"ERROR during solution: "<<error<<std::endl;
                exit(33);
            }
        }
    }

 template struct PardisoWrapper<double, int>;
template struct PardisoWrapper<float, int>;

 template struct PardisoWrapper<double, long long int>;
 template struct PardisoWrapper<float, long long int>;