//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################


#ifndef SUBROUTINE_Singular_Value_Decomposition
#include "KernelCommon.h"
#else
namespace {
#endif

#include "Singular_Value_Decomposition.h"
    using namespace SIMD_Numeric_Kernel;

#define COMPUTE_V_AS_MATRIX
#define COMPUTE_U_AS_MATRIX

namespace{
    //float rsqrt(const float f)
    //{__m128 val=_mm_load_ss(&f);val=_mm_rsqrt_ss(val);float result;_mm_store_ss(&result,val);return result;}

typedef enum { x11=0, x21, x31, x12, x22, x32, x13, x23, x33 } Entry;
}

  template<class Tarch, class T_DATA>
#ifdef SUBROUTINE_Singular_Value_Decomposition
__forceinline
#endif
  void Singular_Value_Decomposition(const T_DATA (&A)[9], T_DATA (&U)[9], T_DATA (&S)[3], T_DATA (&V)[9])
{

#include "Singular_Value_Decomposition_Kernel_Declarations.hpp"

    Va11.Load_Aligned(A[x11]);
    Va21.Load_Aligned(A[x21]);
    Va31.Load_Aligned(A[x31]);
    Va12.Load_Aligned(A[x12]);
    Va22.Load_Aligned(A[x22]);
    Va32.Load_Aligned(A[x32]);
    Va13.Load_Aligned(A[x13]);
    Va23.Load_Aligned(A[x23]);
    Va33.Load_Aligned(A[x33]);


#include "Singular_Value_Decomposition_Main_Kernel_Body.hpp"


    Store(U[x11],Vu11);
    Store(U[x21],Vu21);
    Store(U[x31],Vu31);
    Store(U[x12],Vu12);
    Store(U[x22],Vu22);
    Store(U[x32],Vu32);
    Store(U[x13],Vu13);
    Store(U[x23],Vu23);
    Store(U[x33],Vu33);

    Store(S[0],Va11);
    Store(S[1],Va22);
    Store(S[2],Va33);

    Store(V[x11],Vv11);
    Store(V[x21],Vv21);
    Store(V[x31],Vv31);
    Store(V[x12],Vv12);
    Store(V[x22],Vv22);
    Store(V[x32],Vv32);
    Store(V[x13],Vv13);
    Store(V[x23],Vv23);
    Store(V[x33],Vv33);

}

#ifndef SUBROUTINE_Singular_Value_Decomposition
#define INSTANCE_KERNEL_Singular_Value_Decomposition(WIDTH,TYPE) const WIDETYPE(TYPE,WIDTH) (&A)[9], WIDETYPE(TYPE,WIDTH) (&U)[9], WIDETYPE(TYPE,WIDTH) (&S)[3], WIDETYPE(TYPE,WIDTH) (&V)[9]
INSTANCE_KERNEL(Singular_Value_Decomposition);
#undef INSTANCE_KERNEL_Singular_Value_Decomposition
#else
}
#endif
