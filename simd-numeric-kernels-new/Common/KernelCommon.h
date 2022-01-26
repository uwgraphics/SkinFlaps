//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################
#pragma once

#include <type_traits>

//
//   Don't define here. Use FIB=1 on the makefile command line.
//
//#define FORCE_IDENTICAL_BEHAVIOR

//
//   Enable use of c++ std io to print Number's
//
#define ENABLE_IO_SUPPORT

//
//   Enable INTEL VECTOR ARCHETECTURES
//
//#define ENABLE_SSE_INSTRUCTION_SET
//#define ENABLE_AVX_INSTRUCTION_SET
//#define ENABLE_MIC_INSTRUCTION_SET

//
//   Enable ARM VECTOR ARCHETECTURES
//
//#define ENABLE_NEON_INSTRUCTION_SET




#if defined(ENABLE_AVX_INSTRUCTION_SET)
#include <immintrin.h>
#ifndef __INTEL_COMPILER
#include <pmmintrin.h>
#endif
#endif

#if defined(ENABLE_MIC_INSTRUCTION_SET)
#include <immintrin.h>
#include <zmmintrin.h>
#endif

//#include "NumberPolicy.h"
#include "Mask.h"
#include "Mask.Scalar.h"
#include "Number.h"
#include "Number.Scalar.h"
//#include "Discrete.h"
#include "Vector3.h"
#include "Constants.h"


// Load Architecture Specific Versions of Number

template <typename Vtype>
struct VTYPE_POLICY{
    const static int V_WIDTH=0;
};

template <>
struct VTYPE_POLICY<float>{
    const static int V_WIDTH=1;
};

#if defined(ENABLE_AVX_INSTRUCTION_SET)
template <>
struct VTYPE_POLICY<__m256>{
    const static int V_WIDTH=8;
};
//#include "arch/x86_64/Mask.AVX2.h"
#include "arch/x86_64/Number.AVX2.h"
//#include "arch/x86_64/Discrete.AVX2.h"
#endif


#if defined(ENABLE_MIC_INSTRUCTION_SET)
template <>
struct VTYPE_POLICY<__m512>{
    const static int V_WIDTH=16;
};
//#include "arch/x86_64/Mask.AVX512.h"
#include "arch/x86_64/Number.AVX512.h"
//#include "arch/x86_64/Discrete.AVX512.h"
#endif


// Define Architecture Specific Helper Macros and Miscellanious

#if defined(ENABLE_AVX_INSTRUCTION_SET)
#define KERNEL_MEM_ALIGN __declspec(align(32))
#else
#define KERNEL_MEM_ALIGN
#endif


// Define this to switch to non-augmented material formulas
#define USE_NONMIXED_FORMULAS

struct COROTATED_TAG;
struct NEOHOOKEAN_TAG;
struct BIPHASIC_TAG;

template <class T, int width>
struct widetype{
   typedef T type[width];
};

template <class T>
struct widetype<T,1>{
   typedef T type;
};

#define WIDETYPE(TYPE,WIDTH) widetype<TYPE, WIDTH>::type
#define VWIDTH(TYPE) VTYPE_POLICY<TYPE>::V_WIDTH

#define INSTANCE_KERNEL_SCALAR_FLOAT(KERNEL,DATA_WIDTH) template void KERNEL<SIMDArchitectureScalar<float>,WIDETYPE(float,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH, float) );
#define INSTANCE_KERNEL_SIMD_FLOAT(KERNEL,DATA_WIDTH) template void KERNEL<SIMDArchitectureScalar<float>,WIDETYPE(float,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH, float) );

#define INSTANCE_KERNEL_SCALAR_DOUBLE(KERNEL,DATA_WIDTH) template void KERNEL<SIMDArchitectureScalar<double>,WIDETYPE(double,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH, double) );
#define INSTANCE_KERNEL_SIMD_DOUBLE(KERNEL,DATA_WIDTH) template void KERNEL<SIMDArchitectureScalar<double>,WIDETYPE(double,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH, double) );


#ifdef ENABLE_AVX_INSTRUCTION_SET
#define INSTANCE_KERNEL_SIMD_AVX_FLOAT(KERNEL,DATA_WIDTH) template void KERNEL<SIMDArchitectureAVX2<float>,WIDETYPE(float,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH, float) );
#define INSTANCE_KERNEL_SIMD_AVX_DOUBLE(KERNEL,DATA_WIDTH) template void KERNEL<SIMDArchitectureAVX2<double>,WIDETYPE(double,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH, double) );
#else
#define INSTANCE_KERNEL_SIMD_AVX_FLOAT(KERNEL,DATA_WIDTH)
#define INSTANCE_KERNEL_SIMD_AVX_DOUBLE(KERNEL,DATA_WIDTH)
#endif
#ifdef ENABLE_MIC_INSTRUCTION_SET
#define INSTANCE_KERNEL_SIMD_MIC_FLOAT(KERNEL,DATA_WIDTH) template void KERNEL<SIMDArchitectureAVX512<float>,WIDETYPE(float,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH, float) );
#define INSTANCE_KERNEL_SIMD_MIC_DOUBLE(KERNEL,DATA_WIDTH) template void KERNEL<SIMDArchitectureAVX512<double>,WIDETYPE(double,DATA_WIDTH)>( INSTANCE_KERNEL_ ## KERNEL(DATA_WIDTH, double) );
#else
#define INSTANCE_KERNEL_SIMD_MIC_FLOAT(KERNEL,DATA_WIDTH)
#define INSTANCE_KERNEL_SIMD_MIC_DOUBLE(KERNEL,DATA_WIDTH)
#endif
#define INSTANCE_KERNEL( KERNEL )                \
    INSTANCE_KERNEL_SCALAR_FLOAT( KERNEL, 1 )    \
    INSTANCE_KERNEL_SIMD_FLOAT( KERNEL, 4 )      \
    INSTANCE_KERNEL_SIMD_FLOAT( KERNEL, 8 )      \
    INSTANCE_KERNEL_SIMD_FLOAT( KERNEL, 16 )     \
    INSTANCE_KERNEL_SIMD_AVX_FLOAT( KERNEL, 8 )  \
    INSTANCE_KERNEL_SIMD_AVX_FLOAT( KERNEL, 16 ) \
    INSTANCE_KERNEL_SIMD_MIC_FLOAT( KERNEL, 16 ) \
    INSTANCE_KERNEL_SCALAR_DOUBLE( KERNEL, 1 )   \
    INSTANCE_KERNEL_SIMD_DOUBLE( KERNEL, 2 )     \
    INSTANCE_KERNEL_SIMD_DOUBLE( KERNEL, 4 )     \
    INSTANCE_KERNEL_SIMD_DOUBLE( KERNEL, 8 )     \
    INSTANCE_KERNEL_SIMD_DOUBLE( KERNEL, 16 )    \
    INSTANCE_KERNEL_SIMD_AVX_DOUBLE( KERNEL, 4 ) \
    INSTANCE_KERNEL_SIMD_AVX_DOUBLE( KERNEL, 8 ) \
    INSTANCE_KERNEL_SIMD_AVX_DOUBLE( KERNEL, 16 )\
    INSTANCE_KERNEL_SIMD_MIC_DOUBLE( KERNEL, 8 ) \
    INSTANCE_KERNEL_SIMD_MIC_DOUBLE( KERNEL, 16 ) \

// end of file
