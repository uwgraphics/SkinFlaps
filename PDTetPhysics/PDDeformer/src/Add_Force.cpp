// #pragma once
#include <Common/KernelCommon.h>

#ifdef FORCE_INLINE
#include <Kernels/Matrix_Times_Matrix/Matrix_Times_Matrix.h>
#include <Kernels/Matrix_Times_Transpose/Matrix_Times_Transpose.h>
#include <Kernels/Singular_Value_Decomposition/Singular_Value_Decomposition.h>
#else

#define SUBROUTINE_Matrix_Times_Transpose
#include <Kernels/Matrix_Times_Transpose/Matrix_Times_Transpose.cpp>
#undef SUBROUTINE_Matrix_Times_Transpose

#define SUBROUTINE_Singular_Value_Decomposition
#include <Kernels/Singular_Value_Decomposition/Singular_Value_Decomposition.cpp>
#undef SUBROUTINE_Singular_Value_Decomposition

#define SUBROUTINE_Matrix_Times_Matrix
#include <Kernels/Matrix_Times_Matrix/Matrix_Times_Matrix.cpp>
#undef SUBROUTINE_Matrix_Times_Matrix
#endif

template<class Tarch,class T_DATA>
void Add_Force(const T_DATA (&x_Blocked)[4][3],
               const T_DATA (&DmInverse_Blocked)[9],
               const T_DATA &restVolume,
               const T_DATA &muLow,
               const T_DATA &muHigh,
               const T_DATA &strainMin,
               const T_DATA &strainMax,
               T_DATA (&f_Blocked)[4][3])
{
    using namespace SIMD_Numeric_Kernel;
    constexpr int d = 3;

    using WideNumberType = Number<Tarch>;
    using WideVectorType = Vector3<WideNumberType>;

    using T = typename Tarch::Scalar;
    T TWO[Tarch::Width]{};
    for (int i=0; i<Tarch::Width; i++) TWO[i] = 2;

    alignas(sizeof(T_DATA)) T_DATA  F_Blocked[d * d]{};
    alignas(sizeof(T_DATA)) T_DATA  R_Blocked[d * d]{};
   // alignas(sizeof(T_DATA)) T_DATA  P_Blocked[d*d];
    alignas(sizeof(T_DATA)) T_DATA  U_Blocked[d * d]{};
    alignas(sizeof(T_DATA)) T_DATA  V_Blocked[d * d]{};
    alignas(sizeof(T_DATA)) T_DATA  S_Blocked[d]{};

    WideVectorType v0, v1, v2, v3, v4, v5, v6;
    WideNumberType s0, s1;

    v0.Load_Aligned(x_Blocked[0]);
    v1.Load_Aligned(x_Blocked[1]);
    v2.Load_Aligned(x_Blocked[2]);
    v3.Load_Aligned(x_Blocked[3]);
    v1 = v1-v0;
    v2 = v2-v0;
    v3 = v3-v0;

    v1.Store(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[0][0]));
    v2.Store(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[3][0]));
    v3.Store(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[6][0]));

    Matrix_Times_Matrix<Tarch, T_DATA>(F_Blocked,
                                       DmInverse_Blocked,
                                       F_Blocked);

    Singular_Value_Decomposition<Tarch, T_DATA>(F_Blocked, U_Blocked, S_Blocked, V_Blocked);

    //Sigma[v] = std::min( std::max( Sigma[v], strainMin[eee] ), strainMax[eee] );
    v0.Load_Aligned(S_Blocked);

    s0.Load_Aligned(strainMin);
    s1.Load_Aligned(strainMax);

    v0.x = min(max(v0.x, s0), s1);
    v0.y = min(max(v0.y, s0), s1);
    v0.z = min(max(v0.z, s0), s1);

    //Sigma[v] = muLow[eee] + muHigh[eee] * Sigma[v];
    s0.Load_Aligned(muLow);
    s1.Load_Aligned(muHigh);

    v0 *= s1;
    v0 += s0;

    v1.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(V_Blocked[0][0]));
    v2.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(V_Blocked[3][0]));
    v3.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(V_Blocked[6][0]));

    v1 *= v0.x;
    v2 *= v0.y;
    v3 *= v0.z;

    // R = U * Sigma.asDiagonal() * V.transpose();
    v1.Store(reinterpret_cast<T_DATA(&)[3]>(V_Blocked[0][0]));
    v2.Store(reinterpret_cast<T_DATA(&)[3]>(V_Blocked[3][0]));
    v3.Store(reinterpret_cast<T_DATA(&)[3]>(V_Blocked[6][0]));

    Matrix_Times_Transpose<Tarch, T_DATA>(U_Blocked,
                                          V_Blocked,
                                          R_Blocked);

// MatrixType P =  -2. * ((muHigh[eee] + muLow[eee]) * F - R);
    v0.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[0][0]));
    v1.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[3][0]));
    v2.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[6][0]));

    v3.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(R_Blocked[0][0]));
    v4.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(R_Blocked[3][0]));
    v5.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(R_Blocked[6][0]));

    s0 = s0 + s1;

    v0 *= s0;
    v1 *= s0;
    v2 *= s0;

    // v0 = v0 - v3;
    // v1 = v1 - v4;
    // v2 = v2 - v5;

    v0 = v3 - v0;
    v1 = v4 - v1;
    v2 = v5 - v2;

    s0.Load_Aligned(restVolume);
    s1.Load_Aligned(TWO);

    v0 *= s0;
    v1 *= s0;
    v2 *= s0;

    v0 *= s1;
    v1 *= s1;
    v2 *= s1;

    //v0.Store(S_Blocked);
    v0.Store(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[0][0]));
    v1.Store(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[3][0]));
    v2.Store(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[6][0]));

    //MatrixType H = -restVolume[eee] * P * DmInverse.transpose();
    Matrix_Times_Transpose<Tarch, T_DATA>(F_Blocked,
                                          DmInverse_Blocked,
                                          F_Blocked);

    v0.Load_Aligned(f_Blocked[0]);
    v1.Load_Aligned(f_Blocked[1]);
    v2.Load_Aligned(f_Blocked[2]);
    v3.Load_Aligned(f_Blocked[3]);

    v4.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[0][0]));
    v5.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[3][0]));
    v6.Load_Aligned(reinterpret_cast<T_DATA(&)[3]>(F_Blocked[6][0]));

    v0 = v0 - v4;
    v0 = v0 - v5;
    v0 = v0 - v6;

    v1 = v1 + v4;
    v2 = v2 + v5;
    v3 = v3 + v6;

    v0.Store(f_Blocked[0]);
    v1.Store(f_Blocked[1]);
    v2.Store(f_Blocked[2]);
    v3.Store(f_Blocked[3]);
}

#define INSTANCE_KERNEL_Add_Force(WIDTH,TYPE)               \
    const WIDETYPE(TYPE,WIDTH) (&x_Blocked)[4][3],          \
        const WIDETYPE(TYPE,WIDTH) (&DmInverse_Blocked)[9], \
        const WIDETYPE(TYPE,WIDTH) &restVolume,             \
        const WIDETYPE(TYPE,WIDTH) &muLow,                  \
        const WIDETYPE(TYPE,WIDTH) &muHigh,                 \
        const WIDETYPE(TYPE,WIDTH) &strainMin,              \
        const WIDETYPE(TYPE,WIDTH) &strainMax,              \
        WIDETYPE(TYPE,WIDTH) (&f_Blocked)[4][3]

INSTANCE_KERNEL_SIMD_AVX_FLOAT( Add_Force, 16)
INSTANCE_KERNEL_SIMD_MIC_FLOAT( Add_Force, 16)
#undef INSTANCE_KERNEL_Add_Force
