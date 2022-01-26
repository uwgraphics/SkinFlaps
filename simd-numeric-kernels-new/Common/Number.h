//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#pragma once


#include <cmath>

#ifdef ENABLE_IO_SUPPORT
#include <iostream>
#endif

namespace SIMD_Numeric_Kernel {
template<class Tw> class Number;

template<class Tarch> Number<Tarch> min(const Number<Tarch>& A, const Number<Tarch>& B);
template<class Tarch> Number<Tarch> max(const Number<Tarch>& A, const Number<Tarch>& B);
template<class Tarch> Number<Tarch> blend(const Mask<Tarch>& mask, const Number<Tarch>& A, const Number<Tarch>& B);

template<class Tarch> void Store(typename Tarch::Scalar* data,const Number<Tarch>& number);
template<class Tarch> void Store(typename Tarch::Scalar& data,const Number<Tarch>& number);
#if 0
template<class Tw> void StoreQuadIn16(float* data,const Number<Tw>& number,int quad);
template<class Tw> void StoreQuadIn16(float& data,const Number<Tw>& number,int quad);

#ifdef ENABLE_IO_SUPPORT
template<class Tw> std::ostream& operator<<( std::ostream& os, const Number<Tw>& number);
#endif

#endif
template<class ArchType>
class Number
{
    using Tarch = ArchType;
    using Tw = typename Tarch::ScalarRegister;
    Tw value;
public:
#if 0
    typedef typename NumberPolicy<Number<Tw> >::MASK_TYPE Mask;
#endif
    Number();
#if 0
    explicit Number(const Mask& mask);

    #endif
    Number operator+(const Number& other) const;
    Number operator*(const Number& other) const;
    Number operator-(const Number& other) const;
#if 0
    Number operator/(const Number& other) const;
#endif
    Mask<Tarch> operator<(const Number& other) const;
#if 0
    Mask operator>(const Number& other) const;
#endif
    Mask<Tarch> operator<=(const Number& other) const;
    Mask<Tarch> operator>=(const Number& other) const;
#if 0
    Mask operator==(const Number& other) const;

    Number operator&(const Number& other) const;
    Number operator|(const Number& other) const;
#endif
    Number operator^(const Number& other) const;
#if 0
    Number andnot(const Number& other) const;
    Number operator~() const;
#endif
    Number sqrt() const;
    Number rsqrt() const;
#if 0
    Number log() const;
    Number exp() const;
    Number inverse() const;
    Number abs() const;
    Number sign() const;
#endif
    friend Number min<>(const Number& A, const Number& B);
    friend Number max<>(const Number& A, const Number& B);
    friend Number blend<>(const Mask<Tarch>& mask, const Number& A, const Number& B);

    Number mask(const Mask<Tarch>& mask) const;
    void Load_Aligned(const typename Tarch::Scalar* data);
    void Load_Aligned(const typename Tarch::Scalar& data);
#if 0
    void Load(const float* data);
    void Load(const float& data);

    void Gather(const float* data,const int* offsets);
    void Gather(const float* data,const int& offsets);
#endif
    friend void Store<>(typename Tarch::Scalar* data, const Number& number);
    friend void Store<>(typename Tarch::Scalar& data, const Number& number);
#if 0
    friend void StoreQuadIn16<>(float* data,const Number& number,int quad);
    friend void StoreQuadIn16<>(float& data,const Number& number,int quad);


    Number Spread(int i);
    Number Distribute(int i );
    Number SwizzleAdd(int i, const Number& other);
    Number Horizontal_Quad_Add();
    Number Quad_Mask(int i);

#ifdef ENABLE_IO_SUPPORT
    friend std::ostream& operator<< <>( std:: ostream & os, const Number& number);
#endif

    enum NUMBER_CONSTANTS  {NEG_TEN=0,NEG_NINE,NEG_EIGHT,NEG_SEVEN,NEG_SIX,NEG_FIVE,NEG_FOUR,NEG_THREE,
                            NEG_TWO,NEG_ONE,NEG_ONE_OVER_FOUR,NEG_ONE_HALF,ZERO,ONE_HALF,ONE_OVER_FOUR,
                            ONE,TWO,THREE,FOUR,FIVE,SIX,SEVEN,EIGHT,NINE,TEN};


#endif
};
}
