//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#pragma once

#include "arch/x86_64/SIMDArchitectureAVX512.h"
namespace SIMD_Numeric_Kernel {
//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//


template<> inline
Number<SIMDArchitectureAVX512<float>>::Number()
{value=_mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(value),_mm512_castps_si512(value)));}


template<> inline
Number<SIMDArchitectureAVX512<double>>::Number()
{value=_mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(value),_mm512_castpd_si512(value)));}


//==============================================================//
//==============================================================//

//==============================================================//
//                                                              //
//                      BASIC OPERATIONS                        //
//                                                              //
//==============================================================//


//------------------------------------//
//             ADDITION               //
//------------------------------------//


template<> inline
Number<SIMDArchitectureAVX512<float>> Number<SIMDArchitectureAVX512<float>>::operator+(const Number& other) const
{Number<SIMDArchitectureAVX512<float>> result;result.value=_mm512_add_ps(value,other.value);return result;}


template<> inline
Number<SIMDArchitectureAVX512<double>> Number<SIMDArchitectureAVX512<double>>::operator+(const Number& other) const
{Number<SIMDArchitectureAVX512<double>> result;result.value=_mm512_add_pd(value,other.value);return result;}


//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//


template<> inline
Number<SIMDArchitectureAVX512<float>> Number<SIMDArchitectureAVX512<float>>::operator*(const Number& other) const
{Number<SIMDArchitectureAVX512<float>> result;result.value=_mm512_mul_ps(value,other.value);return result;}


template<> inline
Number<SIMDArchitectureAVX512<double>> Number<SIMDArchitectureAVX512<double>>::operator*(const Number& other) const
{Number<SIMDArchitectureAVX512<double>> result;result.value=_mm512_mul_pd(value,other.value);return result;}


//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//


template<> inline
Number<SIMDArchitectureAVX512<float>> Number<SIMDArchitectureAVX512<float>>::operator-(const Number& other) const
{Number<SIMDArchitectureAVX512<float>> result;result.value=_mm512_sub_ps(value,other.value);return result;}


template<> inline
Number<SIMDArchitectureAVX512<double>> Number<SIMDArchitectureAVX512<double>>::operator-(const Number& other) const
{Number<SIMDArchitectureAVX512<double>> result;result.value=_mm512_sub_pd(value,other.value);return result;}


#if 0
//------------------------------------//
//            DIVISION                //
//------------------------------------//


template<> inline
Number<__m512> Number<__m512>::operator/(const Number& other) const
{Number<__m512> result;result.value=_mm512_div_ps(value,other.value);return result;}


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
Mask<SIMDArchitectureAVX512<float>> Number<SIMDArchitectureAVX512<float>>::operator<(const Number& other) const
{
    Mask<Tarch> result;
    result.value=_mm512_cmp_ps_mask(value,other.value,_MM_CMPINT_LT);
    return result;
}


template<> inline
Mask<SIMDArchitectureAVX512<double>> Number<SIMDArchitectureAVX512<double>>::operator<(const Number& other) const
{
    Mask<Tarch> result;
    result.value=_mm512_cmp_pd_mask(value,other.value,_MM_CMPINT_LT);
    return result;
}


#if 0
//------------------------------------//
//             GREATER THAN              //
//------------------------------------//


template<> inline
typename Number<__m512>::Mask Number<__m512>::operator>(const Number& other) const
{
    Mask result;
    result.value=_mm512_cmp_ps_mask(value,other.value,_MM_CMPINT_GT);
    return result;
}

#endif
//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//


template<> inline
Mask<SIMDArchitectureAVX512<float>> Number<SIMDArchitectureAVX512<float>>::operator<=(const Number& other) const
{
    Mask<Tarch> result;
    result.value=_mm512_cmp_ps_mask(value,other.value,_MM_CMPINT_LE);
    return result;
}


template<> inline
Mask<SIMDArchitectureAVX512<double>> Number<SIMDArchitectureAVX512<double>>::operator<=(const Number& other) const
{
    Mask<Tarch> result;
    result.value=_mm512_cmp_pd_mask(value,other.value,_MM_CMPINT_LE);
    return result;
}


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//



template<> inline
Mask<SIMDArchitectureAVX512<float>> Number<SIMDArchitectureAVX512<float>>::operator>=(const Number& other) const
{
    Mask<Tarch> result;
    result.value=_mm512_cmp_ps_mask(value,other.value,_MM_CMPINT_GE);
    return result;
}


template<> inline
Mask<SIMDArchitectureAVX512<double>> Number<SIMDArchitectureAVX512<double>>::operator>=(const Number& other) const
{
    Mask<Tarch> result;
    result.value=_mm512_cmp_pd_mask(value,other.value,_MM_CMPINT_GE);
    return result;
}


#if 0
//------------------------------------//
//               EQUALS               //
//------------------------------------//


template<> inline
typename Number<__m512>::Mask Number<__m512>::operator==(const Number& other) const
{
    Mask result;
    result.value=_mm512_cmp_ps_mask(value,other.value,_MM_CMPINT_EQ);
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
Number<__m512> Number<__m512>::operator&(const Number& other) const
{Number<__m512> result;result.value=_mm512_castsi512_ps(_mm512_and_epi32(_mm512_castps_si512(value),_mm512_castps_si512(other.value)));return result;}


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//


template<> inline
Number<__m512> Number<__m512>::operator|(const Number& other) const
{Number<__m512> result;result.value=_mm512_castsi512_ps(_mm512_or_epi32(_mm512_castps_si512(value),_mm512_castps_si512(other.value)));return result;}

#endif
//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//


template<> inline
Number<SIMDArchitectureAVX512<float>> Number<SIMDArchitectureAVX512<float>>::operator^(const Number& other) const
{Number<Tarch> result;result.value=_mm512_castsi512_ps(_mm512_xor_epi32(_mm512_castps_si512(value),_mm512_castps_si512(other.value)));return result;}


template<> inline
Number<SIMDArchitectureAVX512<double>> Number<SIMDArchitectureAVX512<double>>::operator^(const Number& other) const
{Number<Tarch> result;result.value=_mm512_castsi512_pd(_mm512_xor_epi32(_mm512_castpd_si512(value),_mm512_castpd_si512(other.value)));return result;}


#if 0
//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//


template<> inline
Number<__m512> Number<__m512>::operator~() const
{
    Number<__m512> result;
    floatConverter X;
    X.i = 0xFFFFFFFF;
    BUILD_CONSTANT_16(__one, X.f);
    __m512 val2=_mm512_load_ps(__one);
    result.value=_mm512_castsi512_ps(_mm512_andnot_epi32(_mm512_castps_si512(value),_mm512_castps_si512(val2)));
    return result;
}


//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//



template<> inline
Number<__m512> Number<__m512>::andnot(const Number& other) const
{
    Number<__m512> result;
    __m512i A, B;
    A = _mm512_castps_si512(value);
    B = _mm512_castps_si512(other.value);
    result.value=_mm512_castsi512_ps(_mm512_and_epi32( B, _mm512_xor_epi32( B, A )));
    return result;
}



//==============================================================//
//==============================================================//

//==============================================================//
//                                                              //
//                  ADVANCED OPERATIONS                         //
//                                                              //
//==============================================================//

//------------------------------------//
//              SQUARE ROOT           //
//------------------------------------//


template<> inline
Number<__m512> Number<__m512>::sqrt() const
{Number<__m512> result;result.value=_mm512_sqrt_ps(value);return result;}

#endif
//------------------------------------//
//       RECIPROCAL SQUARE ROOT       //
//------------------------------------//


template<> inline
Number<SIMDArchitectureAVX512<float>> Number<SIMDArchitectureAVX512<float>>::rsqrt() const
{
    Number<Tarch> result;
    result.value=_mm512_rsqrt14_ps(value); // original 23, illegal 28
    return result;
}


template<> inline
Number<SIMDArchitectureAVX512<double>> Number<SIMDArchitectureAVX512<double>>::rsqrt() const
{
    Number<Tarch> result;
    result.value=_mm512_rsqrt14_pd(value); // original 23, illegal 28
    return result;
}


#if 0
//------------------------------------//
//                 LOG                //
//------------------------------------//


template<> inline
Number<__m512> Number<__m512>::log() const
{
    Number<__m512> result;
    result.value = _mm512_log_ps(value);
return result;
}



//------------------------------------//
//                 EXP                //
//------------------------------------//


template<> inline
Number<__m512> Number<__m512>::exp() const
{
    Number<__m512> result;
    result.value = _mm512_exp_ps(value);
    return result;
}



//------------------------------------//
//        INVERSE (RECIPROCAL)        //
//------------------------------------//


template<> inline
Number<__m512> Number<__m512>::inverse() const
{
   Number<__m512> result;
   result.value = _mm512_rcp23_ps( value );
   return result;
}



//------------------------------------//
//        Absolulte Value             //
//------------------------------------//


template<> inline
Number<__m512> Number<__m512>::abs() const
{
    Number<__m512> neg,result;
    neg.value=_mm512_sub_ps(neg.value,value);
    result.value=_mm512_max_ps(value,neg.value);
    return result;
}



//------------------------------------//
//               SIGN                 //
//------------------------------------//


template<> inline
Number<__m512> Number<__m512>::sign() const
{
    Number<__m512> zero,poss,result;
    Mask negs;
    BUILD_CONSTANT_16(__one,1.0f);
    __m512 none=_mm512_load_ps(__one);
    none = _mm512_sub_ps(zero.value,none);
    __m512 one=_mm512_load_ps(__one);
    negs.value = _mm512_cmp_ps_mask(value,zero.value,_MM_CMPINT_GE);
    result.value = _mm512_mask_blend_ps(negs.value,none,one);

    return result;
}


#endif

//------------------------------------//
//               MINIMUM              //
//------------------------------------//


template<> inline
    Number<SIMDArchitectureAVX512<float>> min(const Number<SIMDArchitectureAVX512<float>>& A, const Number<SIMDArchitectureAVX512<float>>& B)
{Number<SIMDArchitectureAVX512<float>> result;result.value=_mm512_min_ps(A.value, B.value);return result;}


template<> inline
    Number<SIMDArchitectureAVX512<double>> min(const Number<SIMDArchitectureAVX512<double>>& A, const Number<SIMDArchitectureAVX512<double>>& B)
{Number<SIMDArchitectureAVX512<double>> result;result.value=_mm512_min_pd(A.value, B.value);return result;}


//------------------------------------//
//               MAXIMUM              //
//------------------------------------//


template<> inline
Number<SIMDArchitectureAVX512<float>> max(const Number<SIMDArchitectureAVX512<float>>& A, const Number<SIMDArchitectureAVX512<float>>& B)
{Number<SIMDArchitectureAVX512<float>> result;result.value=_mm512_max_ps(A.value, B.value);return result;}


template<> inline
Number<SIMDArchitectureAVX512<double>> max(const Number<SIMDArchitectureAVX512<double>>& A, const Number<SIMDArchitectureAVX512<double>>& B)
{Number<SIMDArchitectureAVX512<double>> result;result.value=_mm512_max_pd(A.value, B.value);return result;}


//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//


template<> inline
Number<SIMDArchitectureAVX512<float>> blend(const Mask<SIMDArchitectureAVX512<float>>& mask, const Number<SIMDArchitectureAVX512<float>>& A, const Number<SIMDArchitectureAVX512<float>>& B)
{Number<SIMDArchitectureAVX512<float>> result;result.value=_mm512_mask_blend_ps(mask.value,A.value,B.value);return result;}


template<> inline
Number<SIMDArchitectureAVX512<double>> blend(const Mask<SIMDArchitectureAVX512<double>>& mask, const Number<SIMDArchitectureAVX512<double>>& A, const Number<SIMDArchitectureAVX512<double>>& B)
{Number<SIMDArchitectureAVX512<double>> result;result.value=_mm512_mask_blend_pd(mask.value,A.value,B.value);return result;}


//------------------------------------//
//      Masked assignment (mask)     //
//------------------------------------//


template<> inline
Number<SIMDArchitectureAVX512<float>> Number<SIMDArchitectureAVX512<float>>::mask(const Mask<Tarch>& mask) const
{
    Number<Tarch> result;
/*
    result.value=_mm512_castsi512_ps(_mm512_mask_xor_epi32(_mm512_castps_si512(result.value),
                                                           mask.value,
                                                           _mm512_castps_si512(value),
                                                           _mm512_castps_si512(value)
                                                           )
                                                           );
*/
    result.value = _mm512_mask_mov_ps( result.value,  mask.value,  value );

    return result;
}


template<> inline
Number<SIMDArchitectureAVX512<double>> Number<SIMDArchitectureAVX512<double>>::mask(const Mask<Tarch>& mask) const
{
    Number<Tarch> result;
/*
    result.value=_mm512_castsi512_ps(_mm512_mask_xor_epi32(_mm512_castps_si512(result.value),
                                                           mask.value,
                                                           _mm512_castps_si512(value),
                                                           _mm512_castps_si512(value)
                                                           )
                                                           );
*/
    result.value = _mm512_mask_mov_pd( result.value,  mask.value,  value );

    return result;
}


//==============================================================//
//==============================================================//

#if 0
#ifdef ENABLE_IO_SUPPORT
//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//


template<> inline
std:: ostream & operator<<( std:: ostream & os, const Number<__m512>& number)
{
    float *fp = (float*)&number.value;
    os << "[ "
       << *(fp+15)
       << " " << *(fp+14)
       << " " << *(fp+13)
       << " " << *(fp+12)
       << " " << *(fp+11)
       << " " << *(fp+10)
       << " " << *(fp+9)
       << " " << *(fp+8)
       << " " << *(fp+7)
       << " " << *(fp+6)
       << " " << *(fp+5)
       << " " << *(fp+4)
       << " " << *(fp+3)
       << " " << *(fp+2)
       << " " << *(fp+1)
       << " " << *(fp) << " ]";
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
void Number<__m512>::Load(const float* data)
{value=_mm512_load_ps((void*)data);}




template<> inline
void Number<__m512>::Load(const float& data)
{value=_mm512_load_ps((void*)&data);}

//------------------------------------//
//              LOADS                 //
//------------------------------------//



template<> inline
void Number<__m512>::Gather(const float* data,const int* offsets)
{
    __m512i index = _mm512_load_epi32((void*)offsets);
    value = _mm512_i32gather_ps(index, (void*)data, 4 );
}




template<> inline
void Number<__m512>::Gather(const float* data,const int& offsets)
{
    __m512i index = _mm512_load_epi32((void*)&offsets);
    value = _mm512_i32gather_ps(index, (void*)data, 4 );
}


#endif
//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//


template<> inline
void Number<SIMDArchitectureAVX512<float>>::Load_Aligned(const float* data)
{value=_mm512_load_ps((void*)data);}


template<> inline
void Number<SIMDArchitectureAVX512<double>>::Load_Aligned(const double* data)
{value=_mm512_load_pd((void*)data);}


template<> inline
void Number<SIMDArchitectureAVX512<float>>::Load_Aligned(const float& data)
{value=_mm512_load_ps((void*)&data);}


template<> inline
void Number<SIMDArchitectureAVX512<double>>::Load_Aligned(const double& data)
{value=_mm512_load_pd((void*)&data);}


//------------------------------------//
//             STORES                 //
//------------------------------------//


template<> inline
void Store(float* data,const Number<SIMDArchitectureAVX512<float>>& number)
{_mm512_store_ps(data,number.value);}


template<> inline
void Store(double* data,const Number<SIMDArchitectureAVX512<double>>& number)
{_mm512_store_pd(data,number.value);}


template<> inline
void Store(float& data,const Number<SIMDArchitectureAVX512<float>>& number)
{_mm512_store_ps(&data,number.value);}


template<> inline
void Store(double& data,const Number<SIMDArchitectureAVX512<double>>& number)
{_mm512_store_pd(&data,number.value);}


#if 0
//==============================================================//
//                                                              //
//                  HORIZONTAL OPERATIONS                       //
//                                                              //
//==============================================================//

template<> inline
Number<__m512> Number<__m512>::Spread(int i){
    Number<__m512> result;
    __mmask16 mask;
    mask = _mm512_int2mask( 0xFFFF );
    switch(i){
    case 0:
        result.value = _mm512_mask_permute4f128_ps( value, mask, value, (_MM_PERM_ENUM)0x0 );
        break;
    case 1:
        result.value = _mm512_mask_permute4f128_ps( value, mask, value, (_MM_PERM_ENUM)0x55 );
        break;
    case 2:
        result.value = _mm512_mask_permute4f128_ps( value, mask, value, (_MM_PERM_ENUM)0xAA );
        break;
    case 3:
        result.value = _mm512_mask_permute4f128_ps( value, mask, value, (_MM_PERM_ENUM)0xFF );
        break;
    }
    return result;
};


template<> inline
Number<__m512> Number<__m512>::Distribute(int i){
    Number<__m512> result, tmp;
    tmp.value = _mm512_swizzle_ps(value, _MM_SWIZ_REG_BADC );
    __mmask16 mask1, mask2;
    mask1 = _mm512_int2mask( 0x33CC );
    result.value = _mm512_mask_blend_ps( mask1, value, tmp.value );
    tmp.value = _mm512_swizzle_ps(result.value, _MM_SWIZ_REG_CDAB );
    mask2 = _mm512_int2mask( 0x5A5A );
    result.value = _mm512_mask_blend_ps( mask2, result.value, tmp.value );
    return result;
};


template<> inline
Number<__m512> Number<__m512>::SwizzleAdd(int i, const Number& other){
    Number<__m512> result;
    __mmask16 mask;
    mask = _mm512_int2mask( 0xFFFF );
    switch ( i ){
    case 1:
        result.value = _mm512_mask_permute4f128_ps( value, mask, value, (_MM_PERM_ENUM)0xB1 ); // ???
        result.value = _mm512_add_ps( result.value, value );
        break;
    case 2:
        result.value = _mm512_mask_permute4f128_ps( value, mask, value, (_MM_PERM_ENUM)0x4E ); // ???
        result.value = _mm512_add_ps( result.value, value );
        break;
    };
    return result;
};

template<> inline
Number<__m512> Number<__m512>::Horizontal_Quad_Add(){
    Number<__m512> result;
    Number<__m512> temp;
    temp.value = _mm512_swizzle_ps(value, _MM_SWIZ_REG_CDAB);
    result.value = _mm512_add_ps( value, temp.value);
    temp.value = _mm512_swizzle_ps(result.value, _MM_SWIZ_REG_BADC);
    result.value = _mm512_add_ps( temp.value, result.value);
    return result;
};

template<> inline
Number<__m512> Number<__m512>::Quad_Mask(int i){
    Number<__m512> result;
    Number<__m512> zero;
    __mmask16 mask;
    mask = _mm512_int2mask( 0x8421 );
    result.value = _mm512_mask_blend_ps( mask, zero.value, value );
    return result;
};

template<> inline
void StoreQuadIn16(float& data, const Number<__m512>& number, int quad){
    __mmask16 mask;
    switch(quad){
    case 0:
        mask = _mm512_int2mask( 0x000F );
        break;
    case 1:
        mask = _mm512_int2mask( 0x00F0 );
        break;
    case 2:
        mask = _mm512_int2mask( 0x0F00 );
        break;
    case 3:
        mask = _mm512_int2mask( 0xF000 );
        break;
    }
    _mm512_mask_store_ps(&data,mask,number.value);
}

template<> inline
void StoreQuadIn16(float* data, const Number<__m512>& number, int quad){
    __mmask16 mask;
    switch(quad){
    case 0:
        mask = _mm512_int2mask( 0x000F );
        break;
    case 1:
        mask = _mm512_int2mask( 0x00F0 );
        break;
    case 2:
        mask = _mm512_int2mask( 0x0F00 );
        break;
    case 3:
        mask = _mm512_int2mask( 0xF000 );
        break;
    }
    _mm512_mask_store_ps(data,mask,number.value);
}

//==============================================================//
//==============================================================//



#endif
}
