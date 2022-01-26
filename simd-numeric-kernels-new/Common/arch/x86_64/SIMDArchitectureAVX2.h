//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#pragma once

namespace SIMD_Numeric_Kernel {
    template<class T>
        struct SIMDArchitectureAVX2;

    template<>
        struct SIMDArchitectureAVX2<float>
    {
        static constexpr int Width = 8;
        using Scalar = float;
        using ScalarRegister = __m256;
        using MaskRegister = __m256;
        static_assert(sizeof(ScalarRegister)/sizeof(Scalar) == Width, "size not matching");
        static_assert(sizeof(ScalarRegister)%sizeof(Scalar) == 0, "size needs to be multiples");
    };
    
    template<>
        struct SIMDArchitectureAVX2<double>
    {
        static constexpr int Width = 4;
        using Scalar = double;
        using ScalarRegister = __m256d;
        using MaskRegister = __m256d;
        static_assert(sizeof(ScalarRegister)/sizeof(Scalar) == Width, "size not matching");
        static_assert(sizeof(ScalarRegister)%sizeof(Scalar) == 0, "size needs to be multiples");
    };
}
