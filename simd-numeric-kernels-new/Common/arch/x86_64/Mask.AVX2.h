//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_MASK_AVX_H__
#define __KERNEL_MASK_AVX_H__

//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Mask<__m256>::Mask()
{value=_mm256_xor_ps(value,value);}

//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//

template<> inline
Mask<__m256> Mask<__m256>::operator&(const Mask& other) const
{Mask<__m256> result;result.value=_mm256_and_ps(value,other.value);return result;}




//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Mask<__m256> Mask<__m256>::operator|(const Mask& other) const
{Mask<__m256> result;result.value=_mm256_or_ps(value,other.value);return result;}




//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//


template<> inline
Mask<__m256> Mask<__m256>::operator^(const Mask& other) const
{Mask<__m256> result;result.value=_mm256_xor_ps(value,other.value);return result;}



//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//


template<> inline
Mask<__m256> Mask<__m256>::operator~() const
{
    Mask<__m256> result;
    const float __one[8] = {1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f,1.0f};
    __m256 val2=_mm256_load_ps(__one);
    result.value=_mm256_andnot_ps(value,val2);
     return result;
}


//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//


template<> inline
Mask<__m256> Mask<__m256>::andnot(const Mask& other) const
{
    Mask<__m256> result;
    result.value=_mm256_andnot_ps(value,other.value);
    return result;
}



//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//

template<> inline
std:: ostream & operator<<( std:: ostream & os, const Mask<__m256>& number)
{
    float *fp = (float*)&number.value;
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
//                                                              //
//                  LOADS AND STORES                            //
//                                                              //
//==============================================================//

//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Mask<__m256>::Load(const float* data)
{value=_mm256_loadu_ps(data);}

template<> inline
void Mask<__m256>::Load(const float& data)
{value=_mm256_loadu_ps(&data);}

//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Mask<__m256>::Load_Aligned(const float* data)
{value=_mm256_load_ps(data);}

template<> inline
void Mask<__m256>::Load_Aligned(const float& data)
{value=_mm256_load_ps(&data);}


//------------------------------------//
//             STORES                 //
//------------------------------------//
template<> inline
void Store(float* data,const Mask<__m256>& number)
{_mm256_store_ps(data,number.value);}

template<> inline
void Store(float& data,const Mask<__m256>& number)
{_mm256_store_ps(&data,number.value);}




#endif
