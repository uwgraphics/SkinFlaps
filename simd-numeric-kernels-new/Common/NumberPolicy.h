//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#pragma once

template<class Tw> class Number;
template<class Tw> class Discrete;
template<class Twi> class Mask;


//======================================================
//
//      NUMBER POLICY: Links Number to a Mask
//
//======================================================


struct ERROR_NO_MASK_TYPE;
struct ERROR_NO_DISCRETE_TYPE;
struct ERROR_NO_NUMBER_TYPE;

template<class NumberType> struct NumberPolicy{
    typedef Mask<ERROR_NO_MASK_TYPE> MASK_TYPE;
    typedef Discrete<ERROR_NO_DISCRETE_TYPE> DISCRETE_TYPE;
    typedef Number<ERROR_NO_NUMBER_TYPE> NUMBER_TYPE;
    static const int width=0;
};

template<> struct NumberPolicy<Number<float> >{
    typedef Mask<bool> MASK_TYPE;
    typedef Discrete<int> DISCRETE_TYPE;
    typedef Number<float> NUMBER_TYPE;
    static const int width=1;
};

#if defined(ENABLE_AVX_INSTRUCTION_SET)
template<> struct NumberPolicy<Number<__m256> >{
    typedef Mask<__m256> MASK_TYPE;
    typedef Discrete<__m256i> DISCRETE_TYPE;
    typedef Number<__m256> NUMBER_TYPE;
    static const int width=8;
};
#endif

#if defined(ENABLE_MIC_INSTRUCTION_SET)
template<> struct NumberPolicy<Number<__m512> >{
    typedef Mask<__mmask16> MASK_TYPE;
    typedef Discrete<__m512i> DISCRETE_TYPE;
    typedef Number<__m512> NUMBER_TYPE;
    static const int width=16;
};
#endif

//======================================================
//
//      DISCRETE POLICY: Links Discrete to a Mask
//
//======================================================

template<> struct NumberPolicy<Discrete<int> >{
    typedef Mask<bool> MASK_TYPE;
    typedef Discrete<int> DISCRETE_TYPE;
    typedef Number<float> NUMBER_TYPE;
    static const int width=1;
};

#if defined(ENABLE_AVX_INSTRUCTION_SET)
template<> struct NumberPolicy<Discrete<__m256i> >{
    typedef Mask<__m256> MASK_TYPE;
    typedef Discrete<__m256i> DISCRETE_TYPE;
    typedef Number<__m256> NUMBER_TYPE;
    static const int width=8;
};
#endif

#if defined(ENABLE_MIC_INSTRUCTION_SET)
template<> struct NumberPolicy<Discrete<__m512i> >{
    typedef Mask<__mmask16> MASK_TYPE;
    typedef Discrete<__m512i> DISCRETE_TYPE;
    typedef Number<__m512> NUMBER_TYPE;
    static const int width=16;
};
#endif

//======================================================
//
//      MASK POLICY: Links Mask to a common type
//
//======================================================

struct ERROR_NO_COMMON_TYPE;

template<class MaskType> struct MaskPolicy{
    typedef ERROR_NO_COMMON_TYPE MASK_EXTERNAL_TYPE;
    static const int width=0;
};

template<> struct MaskPolicy<Mask<bool> >{
    typedef float MASK_EXTERNAL_TYPE;
    static const int width=1;
};

#if defined(ENABLE_AVX_INSTRUCTION_SET)
template<> struct MaskPolicy<Mask<__m256> >{
    typedef float MASK_EXTERNAL_TYPE;
    static const int width=8;
};
#endif

#if defined(ENABLE_MIC_INSTRUCTION_SET)
template<> struct MaskPolicy<Mask<__mmask16> >{
    typedef int MASK_EXTERNAL_TYPE;
    static const int width=16;
};
#endif
