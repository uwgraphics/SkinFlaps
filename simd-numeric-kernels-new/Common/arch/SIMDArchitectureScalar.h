//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#pragma once

namespace SIMD_Numeric_Kernel {
template<class T>
struct SIMDArchitectureScalar
{
    static constexpr int Width = 1;
    using Scalar = T;
    using ScalarRegister = T;
    using MaskRegister = bool;
};
}
