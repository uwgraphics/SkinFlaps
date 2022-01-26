//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#pragma once

namespace SIMD_Numeric_Kernel {
    template<class T>
    struct SIMDArchitectureAVX512;

    template<>
    struct SIMDArchitectureAVX512<float>
    {
        static constexpr int Width = 16;
        using Scalar = float;
        using ScalarRegister = __m512;
        using MaskRegister = __mmask16;
        static_assert(sizeof(ScalarRegister)/sizeof(Scalar) == Width, "size not matching");
        static_assert(sizeof(ScalarRegister)%sizeof(Scalar) == 0, "size needs to be multiples");
    };
    
    template<>
    struct SIMDArchitectureAVX512<double>
    {
        static constexpr int Width = 8;
        using Scalar = double;
        using ScalarRegister = __m512d;
        using MaskRegister = __mmask8;
        static_assert(sizeof(ScalarRegister)/sizeof(Scalar) == Width, "size not matching");
        static_assert(sizeof(ScalarRegister)%sizeof(Scalar) == 0, "size needs to be multiples");
    };
}
