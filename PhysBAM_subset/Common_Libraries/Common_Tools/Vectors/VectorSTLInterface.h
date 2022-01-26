//#####################################################################
// Copyright 2019, Eftychios Sifakis.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file
// PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Casting helpers between PhysBAM VECTOR classes and STL arrays
//#####################################################################
#pragma once

#include <array>
#include <PhysBAM_Tools/Vectors/VECTOR.h>

namespace PhysBAM {

template <class PhysBAMType> struct STLInterface;

template <class T, int d> struct STLInterface<VECTOR<T, d>> {
    using PhysBAMType = VECTOR<T, d>;
    using STLType = std::array<T, d>;

    static PhysBAMType &Wrap(STLType &STLv) { return reinterpret_cast<PhysBAMType &>(STLv); }
    static const PhysBAMType &Wrap(const STLType &STLv) { return reinterpret_cast<const PhysBAMType &>(STLv); }
    static const PhysBAMType &WrapConst(const STLType &STLv) { return reinterpret_cast<const PhysBAMType &>(STLv); }
};

} // namespace PhysBAM
