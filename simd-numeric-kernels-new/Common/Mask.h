//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################


#pragma once

#include <cmath>
#include <iostream>
namespace SIMD_Numeric_Kernel {
    template<class Tw> class Mask;
#if 0
    template<class Tw> void Store(typename MaskPolicy<Mask<Tw> >::MASK_EXTERNAL_TYPE* data,const Mask<Tw>& number);
    template<class Tw> void Store(typename MaskPolicy<Mask<Tw> >::MASK_EXTERNAL_TYPE& data,const Mask<Tw>& number);
    template<class Tw> std::ostream& operator<<( std::ostream& os, const Mask<Tw>& number);

    namespace MASK_HELPERS{
        typedef union {
            int i;
            float f;
        } floatConverter;
    }
#endif

    template<class Tarch>
        class Mask
    {
        using Tw = typename Tarch::MaskRegister;
    public:
        Tw value;
#if 0
        typedef typename MaskPolicy<Mask<Tw> >::MASK_EXTERNAL_TYPE Ext_Type;

        Mask();

        static Mask True();
        static Mask False();

        Mask operator&(const Mask& other) const;
        Mask operator|(const Mask& other) const;
        Mask operator^(const Mask& other) const;
        Mask andnot(const Mask& other) const;
        Mask operator~() const;

        void Load_Aligned(const Ext_Type* data);
        void Load_Aligned(const Ext_Type& data);

        void Load(const Ext_Type* data);
        void Load(const Ext_Type& data);

        friend void Store<>(Ext_Type* data,const Mask& number);
        friend void Store<>(Ext_Type& data,const Mask& number);

        friend std::ostream& operator<< <>( std:: ostream & os, const Mask& number);
#endif
    };

}
