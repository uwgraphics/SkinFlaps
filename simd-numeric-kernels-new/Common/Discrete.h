//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifndef __KERNEL_DISCRETE_H__
#define __KERNEL_DISCRETE_H__

#include <cmath>

#ifdef ENABLE_IO_SUPPORT
#include <iostream>
#endif


template<class Tw> class Discrete;
template<class Tw> Discrete<Tw> min(const Discrete<Tw>& A, const Discrete<Tw>& B);
template<class Tw> Discrete<Tw> max(const Discrete<Tw>& A, const Discrete<Tw>& B);
template<class Tw> Discrete<Tw> blend(const typename NumberPolicy<Discrete<Tw> >::MASK_TYPE& mask, const Discrete<Tw>& A, const Discrete<Tw>& B);
template<class Tw> void Store(int* data,const Discrete<Tw>& discrete);
template<class Tw> void Store(int& data,const Discrete<Tw>& discrete);

#ifdef ENABLE_IO_SUPPORT
template<class Tw> std::ostream& operator<<( std::ostream& os, const Discrete<Tw>& discrete);
#endif

template<class Tw>
class Discrete
{
    Tw value;
public:
    typedef typename NumberPolicy<Discrete<Tw> >::MASK_TYPE Mask;
    Discrete();
    explicit Discrete(const Mask& mask);

    Discrete operator+(const Discrete& other) const;
    Discrete operator*(const Discrete& other) const;
    Discrete operator-(const Discrete& other) const;
    Discrete operator/(const Discrete& other) const;

    Mask operator<(const Discrete& other) const;
    Mask operator>(const Discrete& other) const;
    Mask operator<=(const Discrete& other) const;
    Mask operator>=(const Discrete& other) const;
    Mask operator==(const Discrete& other) const;

    Discrete operator&(const Discrete& other) const;
    Discrete operator|(const Discrete& other) const;
    Discrete operator^(const Discrete& other) const;
    Discrete andnot(const Discrete& other) const;
    Discrete operator~() const;
    
    friend Discrete min<>(const Discrete& A, const Discrete& B);
    friend Discrete max<>(const Discrete& A, const Discrete& B);
    friend Discrete blend<>(const Mask& mask, const Discrete& A, const Discrete& B);
    Discrete mask(const Mask& mask) const;

    void Load_Aligned(const int* data);
    void Load_Aligned(const int& data);

    void Load(const int* data);
    void Load(const int& data);

    void Gather(const int* data,const int* offsets);
    void Gather(const int* data,const int& offsets);

    friend void Store<>(int* data,const Discrete& discrete);
    friend void Store<>(int& data,const Discrete& discrete);

#ifdef ENABLE_IO_SUPPORT
    friend std::ostream& operator<< <>( std:: ostream & os, const Discrete& discrete);
#endif

};


//==============================================================//
//                                                              //
//                      CONSTRUCTORS                            //
//                                                              //
//==============================================================//

template<> inline
Discrete<int>::Discrete()
{value=0.;}

template<> inline
Discrete<int>::Discrete(const Mask& mask)
{
    if(mask.value)
        value = 1;   
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


//------------------------------------//
//             ADDITION               //
//------------------------------------//

template<> inline
Discrete<int> Discrete<int>::operator+(const Discrete& other) const
{Discrete result; result.value=value+other.value; return result;}

//------------------------------------//
//           MULTIPLICATION           //
//------------------------------------//

template<> inline
Discrete<int> Discrete<int>::operator*(const Discrete& other) const
{Discrete result;result.value=value*other.value;return result;}


//------------------------------------//
//           SUBTRACTION              //
//------------------------------------//

template<> inline
Discrete<int> Discrete<int>::operator-(const Discrete& other) const
{Discrete result;result.value=value-other.value;return result;}

//------------------------------------//
//            DIVISION                //
//------------------------------------//

template<> inline
Discrete<int> Discrete<int>::operator/(const Discrete& other) const
{Discrete result;result.value=value/other.value;return result;}


//==============================================================//
//                                                              //
//                      COMPARISON OPERATIONS                   //
//                                                              //
//==============================================================//

//------------------------------------//
//             LESS THAN              //
//------------------------------------//

template<> inline
Discrete<int>::Mask Discrete<int>::operator<(const Discrete& other) const
{
    Mask result;
    result.value = value < other.value ? true : false;
    return result;
}


//------------------------------------//
//             GREATER THAN              //
//------------------------------------//

template<> inline
Discrete<int>::Mask Discrete<int>::operator>(const Discrete& other) const
{
    Mask result;
    result.value = value > other.value ? true : false;
    return result;
}

//------------------------------------//
//             LESS EQUALS              //
//------------------------------------//

template<> inline
Discrete<int>::Mask Discrete<int>::operator<=(const Discrete& other) const
{
    Mask result;
    result.value = value <= other.value ? true : false;
    return result;
}


//------------------------------------//
//             GREATER EQUALS              //
//------------------------------------//

template<> inline
Discrete<int>::Mask Discrete<int>::operator>=(const Discrete& other) const
{
    Mask result;
    result.value = value >= other.value ? true : false;
    return result;
}


//------------------------------------//
//               EQUALS               //
//------------------------------------//

template<> inline
Discrete<int>::Mask Discrete<int>::operator==(const Discrete& other) const
{
    Mask result;
    result.value = value == other.value ? true : false;
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
Discrete<int> Discrete<int>::operator&(const Discrete& other) const
{
    Discrete result;
    result.value = value & other.value;
    return result;
}


//------------------------------------//
//             BITWISE OR             //
//------------------------------------//

template<> inline
Discrete<int> Discrete<int>::operator|(const Discrete& other) const
{
    Discrete result;
    result.value = value | other.value;
    return result;
}


//------------------------------------//
//             BITWISE XOR            //
//------------------------------------//

template<> inline
Discrete<int> Discrete<int>::operator^(const Discrete& other) const
{
    Discrete result;
    result.value = value ^ other.value;
    return result;
}


//------------------------------------//
//             BITWISE NOT            //
//------------------------------------//

template<> inline
Discrete<int> Discrete<int>::operator~() const
{
    Discrete result;
    result.value = ~value;
    return result;
}

//------------------------------------//
//             BITWISE ANDNOT         //
//------------------------------------//

template<> inline
Discrete<int> Discrete<int>::andnot(const Discrete& other) const
{
    Discrete result;
    result.value = ~value & other.value;
    return result;
}


//==============================================================//
//==============================================================//

//------------------------------------//
//               MINIMUM              //
//------------------------------------//

template<> inline
Discrete<int> min(const Discrete<int>& A, const Discrete<int>& B) 
{Discrete<int> result;result.value=A.value < B.value ? A.value : B.value;return result;}

//------------------------------------//
//               MAXIMUM              //
//------------------------------------//

template<> inline
Discrete<int> max(const Discrete<int>& A, const Discrete<int>& B)
{Discrete<int> result;result.value=A.value > B.value ? A.value : B.value;return result;}


//------------------------------------//
//      Masked assignment (blend)     //
//------------------------------------//

template<> inline
Discrete<int> blend(const NumberPolicy<Discrete<int> >::MASK_TYPE& mask, const Discrete<int>& A, const Discrete<int>& B)
{
    Discrete<int> result;
    result.value = mask.value ? B.value : A.value;
    return result;
}

//------------------------------------//
//      Masked assignment (mask )     //
//------------------------------------//
template<> inline
Discrete<int> Discrete<int>::mask(const Mask& mask) const
{
    Discrete<int> result;
    result.value = mask.value ? value : (int)(0);
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
std:: ostream & operator<<( std:: ostream & os, const Discrete<int>& discrete)
{
    os << "[ " << discrete.value << " ]";
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
void Discrete<int>::Load(const int* data)
{value=*data;}


template<> inline
void Discrete<int>::Load(const int& data)
{value=data;}


//------------------------------------//
//              LOADS                 //
//------------------------------------//

template<> inline
void Discrete<int>::Gather(const int* data,const int* offsets)
{value=data[*offsets];}

template<> inline
void Discrete<int>::Gather(const int* data,const int& offsets)
{value=data[offsets];}


//------------------------------------//
//           ALIGNED LOADS            //
//------------------------------------//

template<> inline
void Discrete<int>::Load_Aligned(const int* data)
{value=*data;}

template<> inline
void Discrete<int>::Load_Aligned(const int& data)
{value=data;}

//------------------------------------//
//             STORES                 //
//------------------------------------//

template<> inline
void Store(int* data,const Discrete<int>& discrete)
{*data=discrete.value;}

template<> inline
void Store(int& data,const Discrete<int>& discrete)
{data=discrete.value;}


//==============================================================//
//==============================================================//
#if 0

#endif

#endif
