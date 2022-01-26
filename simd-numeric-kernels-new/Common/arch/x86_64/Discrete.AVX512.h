//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_DISCRETE_MIC_H__
#define __KERNEL_DISCRETE_MIC_H__


//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Discrete<__m512i>::Discrete()
{value=_mm512_xor_epi32(value,value);}

template<> inline
Discrete<__m512i>::Discrete(const Mask& mask)
{
    Discrete A;
    BUILD_ICONSTANT_8(__one, 0xFFFFFFFF);
    Discrete B; B.value = _mm512_load_epi32(__one);
    value = _mm512_mask_blend_epi32( mask.value, A.value, B.value);
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
Discrete<__m512i> Discrete<__m512i>::operator+(const Discrete& other) const
{Discrete<__m512i> result;result.value=_mm512_add_epi32(value,other.value);return result;}

//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//

template<> inline
Discrete<__m512i> Discrete<__m512i>::operator*(const Discrete& other) const
{Discrete<__m512i> result;result.value=_mm512_mul_epi32(value,other.value);return result;}



//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//

template<> inline
Discrete<__m512i> Discrete<__m512i>::operator-(const Discrete& other) const
{Discrete<__m512i> result;result.value=_mm512_sub_epi32(value,other.value);return result;}


//------------------------------------//
//            DIVISION                //
//------------------------------------//

template<> inline
Discrete<__m512i> Discrete<__m512i>::operator/(const Discrete& other) const
{Discrete<__m512i> result;result.value=_mm512_div_epi32(value,other.value);return result;}



//==============================================================//
//                                                              //
//                      COMPARISON OPERATIONS                   //
//                                                              //
//==============================================================//

//------------------------------------//
//               EQUALS               //
//------------------------------------//

template<> inline
Discrete<__m512i>::Mask Discrete<__m512i>::operator==(const Discrete& other) const
{
    Mask result;
    result.value=_mm512_cmpeq_epi32_mask(value,other.value);
    return result;
}

//------------------------------------//
//             GREATER THAN              //
//------------------------------------//

template<> inline
Discrete<__m512i>::Mask Discrete<__m512i>::operator>(const Discrete& other) const
{
    Mask result;
    result.value=_mm512_cmpgt_epi32_mask(value,other.value);
    return result;
}

//------------------------------------//
//             LESS THAN              //
//------------------------------------//

template<> inline
Discrete<__m512i>::Mask Discrete<__m512i>::operator<(const Discrete& other) const
{
    Mask result;
    result = ~(((*this) > other) | ((*this) == other ));
    return result;
}

//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//

template<> inline
Discrete<__m512i>::Mask Discrete<__m512i>::operator<=(const Discrete& other) const
{
    Mask result;
    result = ~(((*this) > other));    
    return result;
}


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//

template<> inline
Discrete<__m512i>::Mask Discrete<__m512i>::operator>=(const Discrete& other) const
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
Discrete<__m512i> Discrete<__m512i>::operator&(const Discrete& other) const
{
    Discrete result;
    result.value = _mm512_and_epi32( value, other.value );
    return result;
}


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Discrete<__m512i> Discrete<__m512i>::operator|(const Discrete& other) const
{
    Discrete result;
    result.value = _mm512_or_epi32( value, other.value );
    return result;
}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Discrete<__m512i> Discrete<__m512i>::operator^(const Discrete& other) const
{
    Discrete result;
    result.value = _mm512_xor_epi32( value, other.value );
    return result;
}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Discrete<__m512i> Discrete<__m512i>::operator~() const
{
    Discrete result;
    BUILD_ICONSTANT_8(__one, 0xFFFFFFFF);
    __m512i val2=_mm512_load_epi32(__one);
    result.value=_mm512_xor_epi32(value,val2);
    return result;
}

//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Discrete<__m512i> Discrete<__m512i>::andnot(const Discrete& other) const
{
    Discrete result;
    result.value = _mm512_andnot_epi32( value, other.value );
    return result;
}


//==============================================================//
//==============================================================//

//------------------------------------//
//               MINIMUM              //
//------------------------------------//

template<> inline
Discrete<__m512i> min(const Discrete<__m512i>& A, const Discrete<__m512i>& B) 
{Discrete<__m512i> result;result.value=_mm512_min_epi32(A.value, B.value);return result;}

//------------------------------------//
//               MAXIMUM              //
//------------------------------------//

template<> inline
Discrete<__m512i> max(const Discrete<__m512i>& A, const Discrete<__m512i>& B)
{Discrete<__m512i> result;result.value=_mm512_min_epi32(A.value, B.value);return result;}


//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//

template<> inline
Discrete<__m512i> blend(const NumberPolicy<Discrete<__m512i> >::MASK_TYPE& mask, const Discrete<__m512i>& A, const Discrete<__m512i>& B)
{
    Discrete<__m512i> result;
    result.value = _mm512_mask_blend_epi32( mask.value, A.value, B.value);
    return result;
}

//------------------------------------//
//      Masked assignment (mask )     //
//------------------------------------//
template<> inline
Discrete<__m512i> Discrete<__m512i>::mask(const Mask& mask) const
{
    Discrete<__m512i> result;
    Discrete<__m512i> temp;
    result.value = _mm512_mask_blend_epi32( mask.value, temp.value, value);
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
std:: ostream & operator<<( std:: ostream & os, const Discrete<__m512i>& discrete)
{
    int *fp = (int*)&discrete.value;
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
void Discrete<__m512i>::Load(const int* data)
{value=_mm512_load_si512(data);}


template<> inline
void Discrete<__m512i>::Load(const int& data)
{value=_mm512_load_si512(&data);}


//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Discrete<__m512i>::Gather(const int* data,const int* offsets)
{
    __m512i index = _mm512_load_epi32((void*)offsets);
    value = _mm512_i32gather_epi32(index, (void*)data, 4 );
}

template<> inline
void Discrete<__m512i>::Gather(const int* data,const int& offsets)
{
    __m512i index = _mm512_load_epi32((void*)&offsets);
    value = _mm512_i32gather_epi32(index, (void*)data, 4 );
}


//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Discrete<__m512i>::Load_Aligned(const int* data)
{value=_mm512_load_si512(data);}

template<> inline
void Discrete<__m512i>::Load_Aligned(const int& data)
{value=_mm512_load_si512(&data);}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(int* data,const Discrete<__m512i>& discrete)
{_mm512_storeu_si512(data,discrete.value);}

template<> inline
void Store(int& data,const Discrete<__m512i>& discrete)
{_mm512_storeu_si512(&data,discrete.value);}


//==============================================================//
//==============================================================//
#if 0

#endif

#endif
