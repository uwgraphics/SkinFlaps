//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Matrix_Times_Transpose
#include <assert.h>
#include "KernelCommon.h"
#else
namespace {
#endif
    using namespace SIMD_Numeric_Kernel;

template<class Tw,class T_DATA>
#ifdef SUBROUTINE_Matrix_Times_Transpose
__forceinline
#endif
void Matrix_Times_Transpose(const T_DATA (&A)[9], const T_DATA (&B)[9], T_DATA (&C)[9])
{
    typedef enum { x11=0,x21,x31,x12,x22,x32,x13,x23,x33} Matrix_Entry;

    typedef SIMD_Numeric_Kernel::Number<Tw> Tn; // icc says that there is ambiguity here, wondering why this is the case
    typedef typename std::remove_extent<T_DATA>::type T_Base;
    //typedef typename std::remove_extent<I_DATA>::type I_Base;
    //static_assert( std::rank<T_DATA>::value == std::rank<I_DATA>::value, "Error: Base datatypes mismatch in width." );
    static_assert( std::rank<T_DATA>::value == 1 ||
                   std::rank<T_DATA>::value == 0 , "Error: Base datatypes must be a rank 1 array or a scalar." );

    const int T_SIZE = std::rank<T_DATA>::value == 0 ? 1 : std::extent<T_DATA>::value;
    typedef const T_Base (&VEC3_TYPE)[3][T_SIZE];

    Tn rB;

    Vector3<Tn> vA;
    Vector3<Tn> vC1;
    Vector3<Tn> vC2;
    Vector3<Tn> vC3;


    // Use 1st column of A & B

    vA.Load_Aligned(A[x11],A[x21],A[x31]);

    rB.Load_Aligned(B[x11]);
    vC1=vA*rB;

    rB.Load_Aligned(B[x21]);
    vC2=vA*rB;

    rB.Load_Aligned(B[x31]);
    vC3=vA*rB;

    // Use 2nd column of A & B

    vA.Load_Aligned(A[x12],A[x22],A[x32]);

    rB.Load_Aligned(B[x12]);
    vC1=vC1+vA*rB;

    rB.Load_Aligned(B[x22]);
    vC2=vC2+vA*rB;

    rB.Load_Aligned(B[x32]);
    vC3=vC3+vA*rB;

    // Use 3rd column of A & B

    vA.Load_Aligned(A[x13],A[x23],A[x33]);

    rB.Load_Aligned(B[x13]);
    vC1=vC1+vA*rB;

    rB.Load_Aligned(B[x23]);
    vC2=vC2+vA*rB;

    rB.Load_Aligned(B[x33]);
    vC3=vC3+vA*rB;

    // Write result

    vC1.Store( C[x11], C[x21], C[x31] );
    vC2.Store( C[x12], C[x22], C[x32] );
    vC3.Store( C[x13], C[x23], C[x33] );
}

#ifndef SUBROUTINE_Matrix_Times_Transpose
#define INSTANCE_KERNEL_Matrix_Times_Transpose(WIDTH,TYPE) const WIDETYPE(TYPE,WIDTH) (&A)[9], const WIDETYPE(TYPE,WIDTH) (&B)[9], WIDETYPE(TYPE,WIDTH) (&C)[9]
INSTANCE_KERNEL(Matrix_Times_Transpose);
#undef INSTANCE_KERNEL_Matrix_Times_Transpose
#else
}
#endif
