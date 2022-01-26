//#####################################################################
// Copyright 2005-2013, Avi Robinson-Mosher, Craig Schroeder, Andrew Selle, Rahul Sheth, Tamar Shinar, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class SYMMETRIC_MATRIX_NXN
//#####################################################################
#ifndef __SYMMETRIC_MATRIX_NXN__
#define __SYMMETRIC_MATRIX_NXN__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <cassert>
#include <iostream>
namespace PhysBAM{
template<class T,class GENERATOR> class RANDOM_NUMBERS;
template<class T> class VECTOR_ND;
template<class T> class MATRIX_MXN;
template<class T,int d> class VECTOR;

template<class T>
class SYMMETRIC_MATRIX_NXN
{
public:
    int n; // size of the n by n matrix
    T *x; // pointer to the one dimensional data
    int size; // number of elements in the matrix: (n*n+n)/2

    SYMMETRIC_MATRIX_NXN()
        :n(0),x(0),size(0)
    {}

    SYMMETRIC_MATRIX_NXN(const int n_input);
    SYMMETRIC_MATRIX_NXN(const SYMMETRIC_MATRIX_NXN<T>& matrix_input);

    ~SYMMETRIC_MATRIX_NXN();

    T& operator()(int i,int j)
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    const T& operator()(int i,int j) const
    {return i<j?Element_Upper(i,j):Element_Lower(i,j);}

    T& Element_Upper(int i,int j)
    {return Element_Lower(j,i);}

    const T& Element_Upper(int i,int j) const
    {return Element_Lower(j,i);}

    T& Element_Lower(int i,int j)
    {assert(i<=n && j>=1 && j<=i);return x[((2*n-j)*(j-1)>>1)+i-1];}

    const T& Element_Lower(int i,int j) const
    {assert(i<=n && j>=1 && j<=i);return x[((2*n-j)*(j-1)>>1)+i-1];}

    inline const int& Rows() const{return n;}
    inline const int& Columns() const{return n;}

//#####################################################################
    static SYMMETRIC_MATRIX_NXN<T> Outer_Product(const VECTOR_ND<T>& u);
    SYMMETRIC_MATRIX_NXN<T> Sqr() const;
    void Givens_Conjugate(const int i,const int j,const T c,const T s);
    void Jacobi_Solve_Eigenproblem(ARRAY<VECTOR<int,2> >& givens_pairs,ARRAY<VECTOR<T,2> >& givens_coefficients,const T tolerance=(T)1e-5,
        const int max_iterations=1000000);
    template<class GENERATOR>
    void Maximum_Eigenvalue_Eigenvector_Pair(T& max_eigenvalue,VECTOR_ND<T>& max_eigenvector,RANDOM_NUMBERS<T,GENERATOR>* random_numbers=0,const T tolerance=(T)1e-5,
        const T randomization_decay_factor=(T)0.9,const int max_iterations=1000000);
    void In_Place_Cholesky_Factorization(MATRIX_MXN<T>& L);
    SYMMETRIC_MATRIX_NXN<T>& operator=(const SYMMETRIC_MATRIX_NXN<T>& A);
    SYMMETRIC_MATRIX_NXN<T>& operator+=(const SYMMETRIC_MATRIX_NXN<T>& A);
    SYMMETRIC_MATRIX_NXN<T>& operator-=(const SYMMETRIC_MATRIX_NXN<T>& A);
    SYMMETRIC_MATRIX_NXN<T>& operator*=(const T a);
    SYMMETRIC_MATRIX_NXN<T> operator+(const SYMMETRIC_MATRIX_NXN<T>& A) const;
    SYMMETRIC_MATRIX_NXN<T> operator-(const SYMMETRIC_MATRIX_NXN<T>& A) const;
    SYMMETRIC_MATRIX_NXN<T> operator*(const T a) const;
    VECTOR_ND<T> operator*(const VECTOR_ND<T>& y) const;
    void Set_Identity_Matrix();
    void Set_Zero_Matrix();
    static SYMMETRIC_MATRIX_NXN<T> Identity_Matrix(const int n);
    T Trace() const;
    void Resize(const int n_new);
    template<class T_MATRIX> void Similarity_Transform(SYMMETRIC_MATRIX_NXN<T>& destination,const MATRIX_BASE<T,T_MATRIX>& basis);
//#####################################################################
};
template<class T>
inline SYMMETRIC_MATRIX_NXN<T> operator*(const T a,const SYMMETRIC_MATRIX_NXN<T>& A)
{return A*a;}
template<class T> std::ostream& operator<<(std::ostream& output_stream,const PhysBAM::SYMMETRIC_MATRIX_NXN<T>& A) {
    output_stream<<"[";for(int i=1;i<=A.Rows();i++){for(int j=1;j<=A.Columns();j++){output_stream<<A(i,j);if(j<A.Columns()) output_stream<<" ";}if(i<A.Rows()) output_stream<<"; ";}output_stream<<"]";return output_stream;}
}
#endif
