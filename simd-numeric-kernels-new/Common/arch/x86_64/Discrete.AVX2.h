//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_DISCRETE_AVX_H__
#define __KERNEL_DISCRETE_AVX_H__


//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Discrete<__m256i>::Discrete()
{value=_mm256_xor_si256(value,value);}

template<> inline
Discrete<__m256i>::Discrete(const Mask& mask)
{
    Discrete A; A.value = _mm256_set1_epi32( 0 );
    Discrete B; B.value = _mm256_set1_epi32( 0xFFFFFFFF );
    value = _mm256_castps_si256( _mm256_blendv_ps( _mm256_castsi256_ps(A.value), _mm256_castsi256_ps(B.value), mask.value ));
}

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
Discrete<__m256i> Discrete<__m256i>::operator+(const Discrete& other) const
{Discrete<__m256i> result;result.value=_mm256_add_epi32(value,other.value);return result;}

//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//

template<> inline
Discrete<__m256i> Discrete<__m256i>::operator*(const Discrete& other) const
{Discrete<__m256i> result;result.value=_mm256_mul_epi32(value,other.value);return result;}



//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//

template<> inline
Discrete<__m256i> Discrete<__m256i>::operator-(const Discrete& other) const
{Discrete<__m256i> result;result.value=_mm256_sub_epi32(value,other.value);return result;}


//------------------------------------//
//            DIVISION                //
//------------------------------------//

template<> inline
Discrete<__m256i> Discrete<__m256i>::operator/(const Discrete& other) const
{
    Discrete<__m256i> result;
#ifdef __INTEL_COMPILER
    result.value=_mm256_div_epi32(value,other.value);
#else
    __m256 this_ps = _mm256_cvtepi32_ps(value);
    __m256 other_ps = _mm256_cvtepi32_ps(other.value);
    __m256 result_ps = _mm256_div_ps(this_ps,other_ps);
    result.value = _mm256_cvttps_epi32(result_ps);    
#endif
    return result;
}



//==============================================================//
//                                                              //
//                      COMPARISON OPERATIONS                   //
//                                                              //
//==============================================================//

//------------------------------------//
//               EQUALS               //
//------------------------------------//

template<> inline
Discrete<__m256i>::Mask Discrete<__m256i>::operator==(const Discrete& other) const
{
    Mask result;
    result.value=_mm256_castsi256_ps(_mm256_cmpeq_epi32(value,other.value));
    return result;
}

//------------------------------------//
//             GREATER THAN              //
//------------------------------------//

template<> inline
Discrete<__m256i>::Mask Discrete<__m256i>::operator>(const Discrete& other) const
{
    Mask result;
    result.value=_mm256_castsi256_ps(_mm256_cmpgt_epi32(value,other.value));
    return result;
}

//------------------------------------//
//             LESS THAN              //
//------------------------------------//

template<> inline
Discrete<__m256i>::Mask Discrete<__m256i>::operator<(const Discrete& other) const
{
    Mask result;
    result = ~(((*this) > other) | ((*this) == other ));
    return result;
}

//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//

template<> inline
Discrete<__m256i>::Mask Discrete<__m256i>::operator<=(const Discrete& other) const
{
    Mask result;
    result = ~(((*this) > other));    
    return result;
}


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//

template<> inline
Discrete<__m256i>::Mask Discrete<__m256i>::operator>=(const Discrete& other) const
{
    Mask result;
    result = (((*this) > other) | ((*this) == other ));
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
Discrete<__m256i> Discrete<__m256i>::operator&(const Discrete& other) const
{
    Discrete result;
    result.value = _mm256_and_si256( value, other.value );
    return result;
}


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Discrete<__m256i> Discrete<__m256i>::operator|(const Discrete& other) const
{
    Discrete result;
    result.value = _mm256_or_si256( value, other.value );
    return result;
}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Discrete<__m256i> Discrete<__m256i>::operator^(const Discrete& other) const
{
    Discrete result;
    result.value = _mm256_xor_si256( value, other.value );
    return result;
}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Discrete<__m256i> Discrete<__m256i>::operator~() const
{
    Discrete result;
    BUILD_ICONSTANT_8(__one, 0xFFFFFFFF);
    __m256i val2=_mm256_load_si256((const __m256i*)(__one));
    result.value=_mm256_xor_si256(value,val2);
    return result;
}

//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Discrete<__m256i> Discrete<__m256i>::andnot(const Discrete& other) const
{
    Discrete result;
    result.value = _mm256_andnot_si256( value, other.value );
    return result;
}


//==============================================================//
//==============================================================//

//------------------------------------//
//               MINIMUM              //
//------------------------------------//

template<> inline
Discrete<__m256i> min(const Discrete<__m256i>& A, const Discrete<__m256i>& B) 
{Discrete<__m256i> result;result.value=_mm256_min_epi32(A.value, B.value);return result;}

//------------------------------------//
//               MAXIMUM              //
//------------------------------------//

template<> inline
Discrete<__m256i> max(const Discrete<__m256i>& A, const Discrete<__m256i>& B)
{Discrete<__m256i> result;result.value=_mm256_min_epi32(A.value, B.value);return result;}


//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//

template<> inline
Discrete<__m256i> blend(const NumberPolicy<Discrete<__m256i> >::MASK_TYPE& mask, const Discrete<__m256i>& A, const Discrete<__m256i>& B)
{
    Discrete<__m256i> result;
    result.value = _mm256_castps_si256( _mm256_blendv_ps( _mm256_castsi256_ps(A.value), _mm256_castsi256_ps(B.value), mask.value ));
    return result;
}

//------------------------------------//
//      Masked assignment (mask )     //
//------------------------------------//
template<> inline
Discrete<__m256i> Discrete<__m256i>::mask(const Mask& mask) const
{
    Discrete<__m256i> result;
    Discrete<__m256i> temp;
    result.value = _mm256_castps_si256( _mm256_blendv_ps( _mm256_castsi256_ps(temp.value), _mm256_castsi256_ps(value), mask.value ));
    return result;
}



//==============================================================//
//==============================================================//

#ifdef ENABLE_IO_SUPPORT
//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//

template<> inline
std:: ostream & operator<<( std:: ostream & os, const Discrete<__m256i>& discrete)
{
    int *fp = (int*)&discrete.value;
    os << "[ " << *(fp+7)
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
void Discrete<__m256i>::Load(const int* data)
{value=_mm256_loadu_si256((const __m256i*)(data));}


template<> inline
void Discrete<__m256i>::Load(const int& data)
{value=_mm256_loadu_si256((const __m256i*)(&data));}


//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Discrete<__m256i>::Gather(const int* data,const int* offsets)
{
    // Cheap hacks here: treat int* as float* for loading and bit movement
    __m128 val1=_mm_load_ss(((float*)data)+offsets[0]);
    __m128 val2=_mm_load_ss(((float*)data)+offsets[1]);
    __m128 val3=_mm_load_ss(((float*)data)+offsets[2]);
    __m128 val4=_mm_load_ss(((float*)data)+offsets[3]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    val1=_mm_movelh_ps(val1,val3);
    value=_mm256_castps_si256(_mm256_insertf128_ps(_mm256_castsi256_ps(value),val1,0));
    val1=_mm_load_ss(((float*)data)+offsets[4]);
    val2=_mm_load_ss(((float*)data)+offsets[5]);
    val3=_mm_load_ss(((float*)data)+offsets[6]);
    val4=_mm_load_ss(((float*)data)+offsets[7]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    val1=_mm_movelh_ps(val1,val3);
    value=_mm256_castps_si256(_mm256_insertf128_ps(_mm256_castsi256_ps(value),val1,1));
}

template<> inline
void Discrete<__m256i>::Gather(const int* data,const int& offsets)
{
    // Cheap hacks here: treat int* as float* for loading and bit movement
    __m128 val1=_mm_load_ss(((float*)data)+(&offsets)[0]);
    __m128 val2=_mm_load_ss(((float*)data)+(&offsets)[1]);
    __m128 val3=_mm_load_ss(((float*)data)+(&offsets)[2]);
    __m128 val4=_mm_load_ss(((float*)data)+(&offsets)[3]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    val1=_mm_movelh_ps(val1,val3);
    value=_mm256_castps_si256(_mm256_insertf128_ps(_mm256_castsi256_ps(value),val1,0));
    val1=_mm_load_ss(((float*)data)+(&offsets)[4]);
    val2=_mm_load_ss(((float*)data)+(&offsets)[5]);
    val3=_mm_load_ss(((float*)data)+(&offsets)[6]);
    val4=_mm_load_ss(((float*)data)+(&offsets)[7]);
    val1=_mm_insert_ps(val1,val2,0x10);
    val3=_mm_insert_ps(val3,val4,0x10);
    val1=_mm_movelh_ps(val1,val3);
    value=_mm256_castps_si256(_mm256_insertf128_ps(_mm256_castsi256_ps(value),val1,1));
}


//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Discrete<__m256i>::Load_Aligned(const int* data)
{value=_mm256_load_si256((const __m256i*)(data));}

template<> inline
void Discrete<__m256i>::Load_Aligned(const int& data)
{value=_mm256_load_si256((const __m256i*)(&data));}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(int* data,const Discrete<__m256i>& discrete)
{_mm256_store_si256((__m256i*)(data),discrete.value);}

template<> inline
void Store(int& data,const Discrete<__m256i>& discrete)
{_mm256_store_si256((__m256i*)(&data),discrete.value);}


//==============================================================//
//==============================================================//
#if 0

#endif

#endif
