//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_MASK_MIC_H__
#define __KERNEL_MASK_MIC_H__

#include <cmath>
#include <iostream>

//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Mask<__mmask16>::Mask()
{value=_mm512_kxor(value, value);}

template<> inline
Mask<__mmask16> Mask<__mmask16>::True()
{Mask result; result.value = _mm512_knot(result.value); return result;}

template<> inline
Mask<__mmask16> Mask<__mmask16>::False()
{Mask result; return result;}

//==============================================================//
//                                                              //
//                      BITWISE  OPERATIONS                     //
//                                                              //
//==============================================================//

//------------------------------------//
//             BITWISE AND            //
//------------------------------------//

template<> inline
Mask<__mmask16> Mask<__mmask16>::operator&(const Mask& other) const
{
    Mask<__mmask16> result;
    result.value = _mm512_kand(value,other.value);
    return result;
}

//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Mask<__mmask16> Mask<__mmask16>::operator|(const Mask& other) const
{
    Mask<__mmask16> result;
    result.value = _mm512_kor(value, other.value);
    return result;
}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Mask<__mmask16> Mask<__mmask16>::operator^(const Mask& other) const
{
    Mask<__mmask16> result;
    result.value  = _mm512_kxor(value, other.value);
    return result;
}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Mask<__mmask16> Mask<__mmask16>::operator~() const
{
    Mask<__mmask16> result;
    result.value = _mm512_knot(value);
    return result;
}


//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Mask<__mmask16> Mask<__mmask16>::andnot(const Mask& other) const
{
    Mask<__mmask16> result;
    result.value = _mm512_kandn(value, other.value);
    return result;
}

//==============================================================//
//                                                              //
//                  DEBUG PRINT SUPPORT                         //
//                                                              //
//==============================================================//

template<> inline
std:: ostream & operator<<( std:: ostream & os, const Mask<__mmask16>& number)
{
    int mask = _mm512_mask2int( number.value );
    os << "[ ";
    os << (int)(mask & 0x1) << " ";
    os << (int)(mask & 0x2) << " ";
    os << (int)(mask & 0x4) << " ";
    os << (int)(mask & 0x8) << " ";
    os << (int)(mask & 0x10) << " ";
    os << (int)(mask & 0x20) << " ";
    os << (int)(mask & 0x40) << " ";
    os << (int)(mask & 0x80) << " ";
    os << (int)(mask & 0x100) << " ";
    os << (int)(mask & 0x200) << " ";
    os << (int)(mask & 0x400) << " ";
    os << (int)(mask & 0x800) << " ";
    os << (int)(mask & 0x1000) << " ";
    os << (int)(mask & 0x2000) << " ";
    os << (int)(mask & 0x4000) << " ";
    os << (int)(mask & 0x8000) << " ";
    os << " ]";
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
void Mask<__mmask16>::Load(const Ext_Type* data)
{value=_mm512_int2mask(*data);}

template<> inline
void Mask<__mmask16>::Load(const Ext_Type& data)
{value=_mm512_int2mask(data);}

//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Mask<__mmask16>::Load_Aligned(const Ext_Type* data)
{value=_mm512_int2mask(*data);}

template<> inline
void Mask<__mmask16>::Load_Aligned(const Ext_Type& data)
{value=_mm512_int2mask(data);}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(typename Mask<__mmask16>::Ext_Type* data,const Mask<__mmask16>& number)
{*data=_mm512_mask2int(number.value);}

template<> inline
void Store(typename Mask<__mmask16>::Ext_Type& data,const Mask<__mmask16>& number)
{data=_mm512_mask2int(number.value);}

#endif
