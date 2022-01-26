//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#pragma once

#include "arch/SIMDArchitectureScalar.h"

namespace {
    typedef union {
        int i;
        float f;
    } floatConverter;

    typedef union {
        int i[2];
        double d;
    } doubleConverter;
}

namespace SIMD_Numeric_Kernel {

    template<class T>
        class Number<SIMDArchitectureScalar<T>>
    {
        using Tarch = SIMDArchitectureScalar<T>;
        using Tw = typename Tarch::ScalarRegister;
        Tw value;
    public:
        Number();
        Number operator+(const Number& other) const;
        Number operator*(const Number& other) const;
        Number operator-(const Number& other) const;
        Mask<Tarch> operator<(const Number& other) const;
        Mask<Tarch> operator<=(const Number& other) const;
        Mask<Tarch> operator>=(const Number& other) const;
        Number operator^(const Number& other) const;
        Number sqrt() const;
        Number rsqrt() const;
        friend Number max<>(const Number& A, const Number& B);
        friend Number min<>(const Number& A, const Number& B);
        friend Number blend<>(const Mask<Tarch>& mask, const Number& A, const Number& B);
        Number mask(const Mask<Tarch>& mask) const;
        void Load_Aligned(const T* data);
        void Load_Aligned(const T& data);
        friend void Store<>(T* data,const Number& number);
        friend void Store<>(T& data,const Number& number);
    };

// float implementation
//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

    template<class T> inline
        Number<SIMDArchitectureScalar<T>>::Number()
    {value=0.;}

#if 0
    template<> inline
        Number<float>::Number(const Mask& mask)
    {
        if(mask.value){
            floatConverter ALL_ONES;
            ALL_ONES.i = (0xFFFFFFFF);
            value = ALL_ONES.f;
        }
        else
            value = 0;
    }

//==============================================================//
//==============================================================//

//==============================================================//
//                                                              //
//                      BASIC OPERATIONS                        //
//                                                              //
//==============================================================//

#endif
//------------------------------------//
//             ADDITION               //
//------------------------------------//

    template<> inline
        Number<SIMDArchitectureScalar<float>> Number<SIMDArchitectureScalar<float>>::operator+(const Number& other) const
    {

        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value=value+other.value;
#else
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&other.value);
        __m128 val3=_mm_add_ps(val1,val2);
        _mm_store_ss(&result.value,val3);
#endif
        return result;
    }

    template<> inline
        Number<SIMDArchitectureScalar<double>> Number<SIMDArchitectureScalar<double>>::operator+(const Number& other) const
    {

        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value=value+other.value;
#else
        // TODO: figure this out
#endif
        return result;
    }

//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//


    template<class T> inline
        Number<SIMDArchitectureScalar<T>> Number<SIMDArchitectureScalar<T>>::operator*(const Number& other) const
    {Number result;result.value=value*other.value;return result;}


//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//


    template<class T> inline
        Number<SIMDArchitectureScalar<T>> Number<SIMDArchitectureScalar<T>>::operator-(const Number& other) const
    {Number result;result.value=value-other.value;return result;}


#if 0
//------------------------------------//
//            DIVISION                //
//------------------------------------//

    template<> inline
        Number<float> Number<float>::operator/(const Number& other) const
    {Number result;result.value=value/other.value;return result;}


//==============================================================//
//                                                              //
//                      COMPARISON OPERATIONS                   //
//                                                              //
//==============================================================//
#endif
//------------------------------------//
//             LESS THAN              //
//------------------------------------//


    template<> inline
        Mask<SIMDArchitectureScalar<float>> Number<SIMDArchitectureScalar<float>>::operator<(const Number& other) const
    {
        Mask<Tarch> result;
        floatConverter ALL_ONES;
        ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value = value < other.value ? true : false;
#else
        float f_result;
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&other.value);
        __m128 val3=_mm_cmp_ss(val1,val2,_CMP_LT_OS);
        _mm_store_ss(&f_result,val3);
        ALL_ONES.f=f_result;
        result.value = ALL_ONES.i==0 ? false : true;
#endif
        return result;
    }


    template<> inline
        Mask<SIMDArchitectureScalar<double>> Number<SIMDArchitectureScalar<double>>::operator<(const Number& other) const
    {
        Mask<Tarch> result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value = value < other.value ? true : false;
#else
        // TODO: figure this out
#endif
        return result;
    }


#if 0
//------------------------------------//
//             GREATER THAN              //
//------------------------------------//

    template<> inline
        Number<float>::Mask Number<float>::operator>(const Number& other) const
    {
        Mask result;
        floatConverter ALL_ONES;
        ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value = value > other.value ? true : false;
#else
        float f_result;
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&other.value);
        __m128 val3=_mm_cmp_ss(val1,val2,_CMP_GT_OS);
        _mm_store_ss(&f_result,val3);
        ALL_ONES.f=f_result;
        result.value = ALL_ONES.i==0 ? false : true;
#endif
        return result;
    }
#endif

//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//


    template<> inline
        Mask<SIMDArchitectureScalar<float>> Number<SIMDArchitectureScalar<float>>::operator<=(const Number& other) const
    {
        Mask<Tarch> result;
        floatConverter ALL_ONES;
        ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value = value <= other.value ? true : false;
#else
        float f_result;
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&other.value);
        __m128 val3=_mm_cmp_ss(val1,val2,_CMP_LE_OS);
        _mm_store_ss(&f_result,val3);
        ALL_ONES.f=f_result;
        result.value = ALL_ONES.i==0 ? false : true;
#endif
        return result;
    }


    template<> inline
        Mask<SIMDArchitectureScalar<double>> Number<SIMDArchitectureScalar<double>>::operator<=(const Number& other) const
    {
        Mask<Tarch> result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value = value <= other.value ? true : false;
#else
        // TODO: figure this out
#endif
        return result;
    }


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//


    template<> inline
        Mask<SIMDArchitectureScalar<float>> Number<SIMDArchitectureScalar<float>>::operator>=(const Number& other) const
    {
        Mask<Tarch> result;
        floatConverter ALL_ONES;
        ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value = value >= other.value ? true : false;
#else
        float f_result;
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&other.value);
        __m128 val3=_mm_cmp_ss(val1,val2,_CMP_GE_OS);
        _mm_store_ss(&f_result,val3);
        ALL_ONES.f=f_result;
        result.value = ALL_ONES.i==0 ? false : true;
#endif
        return result;
    }


    template<> inline
        Mask<SIMDArchitectureScalar<double>> Number<SIMDArchitectureScalar<double>>::operator>=(const Number& other) const
    {
        Mask<Tarch> result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value = value >= other.value ? true : false;
#else
        // TODO: figure this out
#endif
        return result;
    }


#if 0
//------------------------------------//
//               EQUALS               //
//------------------------------------//

    template<> inline
        Number<float>::Mask Number<float>::operator==(const Number& other) const
    {
        Mask result;
        floatConverter ALL_ONES;
        ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value = value == other.value ? true : false;
#else
        float f_result;
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&other.value);
        __m128 val3=_mm_cmp_ss(val1,val2,_CMP_EQ_OS);
        _mm_store_ss(&f_result,val3);
        ALL_ONES.f=f_result;
        result.value = ALL_ONES.i==0 ? false : true;
#endif
        return result;
    }

//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//

    template<> inline
        Number<float> Number<float>::operator&(const Number& other) const
    {
        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        floatConverter VAL1;
        floatConverter VAL2;
        floatConverter RESULT;
        VAL1.f = value;
        VAL2.f = other.value;
        RESULT.i = VAL1.i & VAL2.i;
        result.value = RESULT.f;
#else
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&other.value);
        __m128 val3=_mm_and_ps(val1,val2);
        _mm_store_ss(&result.value,val3);
#endif
        return result;
    }


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

    template<> inline
        Number<float> Number<float>::operator|(const Number& other) const
    {
        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        floatConverter VAL1;
        floatConverter VAL2;
        floatConverter RESULT;
        VAL1.f = value;
        VAL2.f = other.value;
        RESULT.i = VAL1.i | VAL2.i;
        result.value = RESULT.f;
#else
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&other.value);
        __m128 val3=_mm_or_ps(val1,val2);
        _mm_store_ss(&result.value,val3);
#endif
        return result;
    }

#endif
//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//


    template<> inline
        Number<SIMDArchitectureScalar<float>> Number<SIMDArchitectureScalar<float>>::operator^(const Number& other) const
    {
        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        floatConverter VAL1;
        floatConverter VAL2;
        floatConverter RESULT;
        VAL1.f = value;
        VAL2.f = other.value;
        RESULT.i = VAL1.i ^ VAL2.i;
        result.value = RESULT.f;
#else
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&other.value);
        __m128 val3=_mm_xor_ps(val1,val2);
        _mm_store_ss(&result.value,val3);
#endif
        return result;
    }


    template<> inline
        Number<SIMDArchitectureScalar<double>> Number<SIMDArchitectureScalar<double>>::operator^(const Number& other) const
    {
        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        doubleConverter VAL1;
        doubleConverter VAL2;
        doubleConverter RESULT;
        VAL1.d = value;
        VAL2.d = other.value;
        RESULT.i[0] = VAL1.i[0] ^ VAL2.i[0];
        RESULT.i[1] = VAL1.i[1] ^ VAL2.i[1];
        result.value = RESULT.d;
#else
        // TODO: figure this out
#endif
        return result;
    }


#if 0
//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

    template<> inline
        Number<float> Number<float>::operator~() const
    {
        Number result;
        floatConverter ALL_ONES;
        ALL_ONES.i = (0xFFFFFFFF);
#ifndef FORCE_IDENTICAL_BEHAVIOR
        floatConverter VAL1;
        floatConverter RESULT;
        VAL1.f = value;
        RESULT.i = ~VAL1.i & ALL_ONES.i;
        result.value = RESULT.f;
#else
        const float __one = ALL_ONES.f;
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&__one);
        __m128 val3=_mm_andnot_ps(val1,val2);
        _mm_store_ss(&result.value,val3);
#endif
        return result;
    }

//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

    template<> inline
        Number<float> Number<float>::andnot(const Number& other) const
    {
        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        floatConverter VAL1;
        floatConverter VAL2;
        floatConverter RESULT;
        VAL1.f = value;
        VAL2.f = other.value;
        RESULT.i = ~VAL1.i & VAL2.i;
        result.value = RESULT.f;
#else
        const float __one = 1.0f;
        __m128 val1=_mm_load_ss(&value);
        __m128 val2=_mm_load_ss(&other.value);
        __m128 val3=_mm_andnot_ps(val1,val2);
        _mm_store_ss(&result.value,val3);
#endif
        return result;
    }


//==============================================================//
//==============================================================//

//==============================================================//
//                                                              //
//                  ADVANCED OPERATIONS                         //
//                                                              //
//==============================================================//
#endif
//------------------------------------//
//              SQUARE ROOT           //
//------------------------------------//


    template<class T> inline
        Number<SIMDArchitectureScalar<T>> Number<SIMDArchitectureScalar<T>>::sqrt() const
    {Number result;result.value=::sqrt(value);return result;}


//------------------------------------//
//       RECIPROCAL SQUARE ROOT       //
//------------------------------------//


    template<> inline
        Number<SIMDArchitectureScalar<float>> Number<SIMDArchitectureScalar<float>>::rsqrt() const
    {
        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value=1.0f/ (float)::sqrt(value);
#else
        __m128 val1=_mm_broadcast_ss(&value);
        val1=_mm_rsqrt_ss(val1);
        _mm_store_ss(&result.value,val1);
#endif
        return result;
    }


    template<> inline
        Number<SIMDArchitectureScalar<double>> Number<SIMDArchitectureScalar<double>>::rsqrt() const
    {
        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value=1.0/::sqrt(value);
#else
        // TODO: figure this out
#endif
        return result;
    }


#if 0
//------------------------------------//
//                 LOG                //
//------------------------------------//

    template<> inline
        Number<float> Number<float>::log() const
    {
        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value=::log(value);
#else
        __m128 val1=_mm_broadcast_ss(&value);
        val1=_mm_log_ps(val1);
        _mm_store_ss(&result.value,val1);
#endif
        return result;
    }



//------------------------------------//
//                 EXP                //
//------------------------------------//

    template<> inline
        Number<float> Number<float>::exp() const
    {
        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value=::exp(value);
#else
        __m128 val1=_mm_broadcast_ss(&value);
        val1=_mm_exp_ps(val1);
        _mm_store_ss(&result.value,val1);
#endif
        return result;
    }


//------------------------------------//
//        INVERSE (RECIPROCAL)        //
//------------------------------------//

    template<> inline
        Number<float> Number<float>::inverse() const
    {
        Number result;
#ifndef FORCE_IDENTICAL_BEHAVIOR
        result.value=1.f/value;
#else
        // 22-bit precision -- See Intel IA32 Architecture Optimization Manual
        __m128 val1=_mm_broadcast_ss(&value);
        __m128 val3=_mm_rcp_ps(val1);
        __m128 val2=_mm_add_ss(val3,val3);
        val3=_mm_mul_ps(val3,val3);
        val3=_mm_mul_ps(val3,val1);
        val2=_mm_sub_ps(val2,val3);
        _mm_store_ss(&result.value,val2);
#endif
        return result;
    }



//------------------------------------//
//        Absolulte Value             //
//------------------------------------//

    template<> inline
        Number<float> Number<float>::abs() const
    {
        Number result;
        result.value=(value>=0.f)?value:-value;
        return result;
    }

//------------------------------------//
//               SIGN                 //
//------------------------------------//

    template<> inline
        Number<float> Number<float>::sign() const
    {
        Number result;
        result.value=(value<0.f)?-1.0f:1.0f;
        return result;
    }
#endif

//------------------------------------//
//               MINIMUM              //
//------------------------------------//


    template<> inline
        Number<SIMDArchitectureScalar<float>> min(const Number<SIMDArchitectureScalar<float>>& A, const Number<SIMDArchitectureScalar<float>>& B)
    {Number<SIMDArchitectureScalar<float>> result;result.value=A.value < B.value ? A.value : B.value;return result;}


    template<> inline
        Number<SIMDArchitectureScalar<double>> min(const Number<SIMDArchitectureScalar<double>>& A, const Number<SIMDArchitectureScalar<double>>& B)
    {Number<SIMDArchitectureScalar<double>> result;result.value=A.value < B.value ? A.value : B.value;return result;}


//------------------------------------//
//               MAXIMUM              //
//------------------------------------//


    template<> inline
        Number<SIMDArchitectureScalar<float>> max(const Number<SIMDArchitectureScalar<float>>& A, const Number<SIMDArchitectureScalar<float>>& B)
    {Number<SIMDArchitectureScalar<float>> result;result.value=A.value > B.value ? A.value : B.value;return result;}


    template<> inline
        Number<SIMDArchitectureScalar<double>> max(const Number<SIMDArchitectureScalar<double>>& A, const Number<SIMDArchitectureScalar<double>>& B)
    {Number<SIMDArchitectureScalar<double>> result;result.value=A.value > B.value ? A.value : B.value;return result;}


//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//


    template<> inline
        Number<SIMDArchitectureScalar<float>> blend(const Mask<SIMDArchitectureScalar<float>>& mask, const Number<SIMDArchitectureScalar<float>>& A, const Number<SIMDArchitectureScalar<float>>& B)
    {
#ifndef FORCE_IDENTICAL_BEHAVIOR
        Number<SIMDArchitectureScalar<float>> result;
        result.value = mask.value ? B.value : A.value;
        return result;
#else
        Number<float> result;
        floatConverter ALL_ONES,ALL_ZEROS;
        ALL_ONES.i = (0xFFFFFFFF);
        ALL_ZEROS.i = (0x0);
        __m128 val1=_mm_load_ss(&A.value);
        __m128 val2=_mm_load_ss(&B.value);
        __m128 val3;
        if(mask.value)
            val3 =_mm_load_ss(&ALL_ONES.f);
        else
            val3=_mm_load_ss(&ALL_ZEROS.f);
        val1=_mm_blendv_ps(val1,val2,val3);
        _mm_store_ss(&result.value,val1);
        return result;
#endif
    }


    template<> inline
        Number<SIMDArchitectureScalar<double>> blend(const Mask<SIMDArchitectureScalar<double>>& mask, const Number<SIMDArchitectureScalar<double>>& A, const Number<SIMDArchitectureScalar<double>>& B)
    {
#ifndef FORCE_IDENTICAL_BEHAVIOR
        Number<SIMDArchitectureScalar<double>> result;
        result.value = mask.value ? B.value : A.value;
        return result;
#else
        // TODO: figure this out
#endif
    }


//------------------------------------//
//      Masked assignment (mask )     //
//------------------------------------//
 
 
    template<> inline
        Number<SIMDArchitectureScalar<float>> Number<SIMDArchitectureScalar<float>>::mask(const Mask<SIMDArchitectureScalar<float>>& mask) const
    {
#ifndef FORCE_IDENTICAL_BEHAVIOR
        Number<Tarch> result;
        result.value = mask.value ? value : (float)(0.0);
        return result;
#else
        Number<Tarch> result;
        floatConverter ALL_ONES,ALL_ZEROS;
        ALL_ONES.i = (0xFFFFFFFF);
        ALL_ZEROS.i = (0x0);
        __m128 val1=_mm_load_ss(&value);
        __m128 val2;
        if(mask.value)
            val2 =_mm_load_ss(&ALL_ONES.f);
        else
            val2=_mm_load_ss(&ALL_ZEROS.f);
        val1=_mm_and_ps(val1,val2);
        _mm_store_ss(&result.value,val1);
        return result;
#endif
    }


    template<> inline
        Number<SIMDArchitectureScalar<double>> Number<SIMDArchitectureScalar<double>>::mask(const Mask<SIMDArchitectureScalar<double>>& mask) const
    {
#ifndef FORCE_IDENTICAL_BEHAVIOR
        Number<Tarch> result;
        result.value = mask.value ? value : (double)(0.0);
        return result;
#else
        // TODO: figure this out
#endif
    }


#if 0
//==============================================================//
//==============================================================//

#ifdef ENABLE_IO_SUPPORT
//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//

    template<> inline
        std:: ostream & operator<<( std:: ostream & os, const Number<float>& number)
    {
        os << "[ " << number.value << " ]";
        return os;
    }

//==============================================================//
//==============================================================//
#endif


//==============================================================//
//                                                              //
//                  LOADS AND STORES                            //
//                                                              //
//==============================================================//

//------------------------------------//
//              LOADS                 //
//------------------------------------//

    template<> inline
        void Number<float>::Load(const float* data)
    {value=*data;}


    template<> inline
        void Number<float>::Load(const float& data)
    {value=data;}


//------------------------------------//
//              LOADS                 //
//------------------------------------//

    template<> inline
        void Number<float>::Gather(const float* data,const int* offsets)
    {value=data[*offsets];}

    template<> inline
        void Number<float>::Gather(const float* data,const int& offsets)
    {value=data[offsets];}

#endif
//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

    template<class T> inline
        void Number<SIMDArchitectureScalar<T>>::Load_Aligned(const T* data)
    {value=*data;}


    template<class T> inline
        void Number<SIMDArchitectureScalar<T>>::Load_Aligned(const T& data)
    {value=data;}

//------------------------------------//
//             STORES                 //
//------------------------------------//


    template<> inline
        void Store(float* data,const Number<SIMDArchitectureScalar<float>>& number)
    {*data=number.value;}


    template<> inline
        void Store(double* data,const Number<SIMDArchitectureScalar<double>>& number)
    {*data=number.value;}


    template<> inline
        void Store(float& data,const Number<SIMDArchitectureScalar<float>>& number)
    {data=number.value;}


    template<> inline
        void Store(double& data,const Number<SIMDArchitectureScalar<double>>& number)
    {data=number.value;}


#if 0
//==============================================================//
//                                                              //
//                  HORIZONTAL OPERATIONS                       //
//                                                              //
//==============================================================//

    template<> inline
        Number<float> Number<float>::Spread(int i){
        return *this;
    };


    template<> inline
        Number<float> Number<float>::Distribute(int i){
        return *this;
    };


    template<> inline
        Number<float> Number<float>::SwizzleAdd(int i, const Number& other){
        return *this;
    };

    template<> inline
        Number<float> Number<float>::Horizontal_Quad_Add(){
        return *this;
    };

    template<> inline
        Number<float> Number<float>::Quad_Mask(int i){
        return *this;
    };

    template<> inline
        void StoreQuadIn16(float& data, const Number<float>& number, int quad){
        data=number.value;
    }

    template<> inline
        void StoreQuadIn16(float* data, const Number<float>& number, int quad){
        *data=number.value;
    }

//==============================================================//
//==============================================================//
#endif
}
