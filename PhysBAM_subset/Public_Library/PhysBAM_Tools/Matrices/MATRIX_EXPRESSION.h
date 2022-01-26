//#####################################################################
// Copyright 2013, Rahul Sheth.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_EXPRESSION
//#####################################################################
#ifndef __MATRIX_EXPRESSION__
#define __MATRIX_EXPRESSION__

#include <PhysBAM_Tools/Vectors/VECTOR_EXPRESSION.h>
#include <PhysBAM_Tools/Matrices/MATRIX_FORWARD.h>
#include <cassert>
namespace PhysBAM{

template<class T_MATRIX,class T_VECTOR>
class LAZY_MATRIX_VECTOR_PRODUCT:public VECTOR_EXPRESSION<typename PRODUCT<typename T_MATRIX::SCALAR,typename T_VECTOR::SCALAR>::TYPE,LAZY_MATRIX_VECTOR_PRODUCT<T_MATRIX,T_VECTOR> >
{
  public:
    typedef typename T_MATRIX::SCALAR T1; typedef typename T_VECTOR::SCALAR T2;
    typedef typename PRODUCT<T1,T2>::TYPE SCALAR;
    typedef VECTOR_EXPRESSION<SCALAR,LAZY_MATRIX_VECTOR_PRODUCT<T_MATRIX,T_VECTOR> > BASE;

    typedef SCALAR ELEMENT;

    const T_MATRIX& matrix;
    const T_VECTOR& vector;

    LAZY_MATRIX_VECTOR_PRODUCT(const T_MATRIX& matrix_in,const T_VECTOR& vector_in)
        :matrix(matrix_in),vector(vector_in)
    {}

    LAZY_MATRIX_VECTOR_PRODUCT(const LAZY_MATRIX_VECTOR_PRODUCT& v_in)
        :BASE(),matrix(v_in.matrix),vector(v_in.vector)
    {}

    int Size() const
    {assert(matrix.Columns()==vector.Size());return matrix.Rows();}

    const SCALAR operator()(const int i) const
    {assert(i<=matrix.Rows() && i>=1);SCALAR sum=(SCALAR)0;        
        for(int j=1;j<=matrix.Columns();j++){sum+=matrix(i,j)*vector(j);}
        return sum;}

//#####################################################################
};

}
#endif
