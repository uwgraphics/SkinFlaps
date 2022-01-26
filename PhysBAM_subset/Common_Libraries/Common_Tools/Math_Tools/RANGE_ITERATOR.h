//#####################################################################
// Copyright 2009, Eftychios Sifakis, Yongning Zhu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class RANGE_ITERATOR
//#####################################################################
#ifndef __RANGE_ITERATOR__
#define __RANGE_ITERATOR__

#include <array>
#include <type_traits>
#include <PhysBAM_Tools/Math_Tools/RANGE.h>

namespace PhysBAM{

template<int d,class TV_INT=VECTOR<int,d> > 
class RANGE_ITERATOR
{
    typedef std::array<typename TV_INT::value_type,d> T_INDEX;

    union {const T_INDEX Rmin;const TV_INT min_corner;};
    union {const T_INDEX Rmax;const TV_INT max_corner;};
    union {T_INDEX Idx;TV_INT index;};

public:
    template<class T_RANGE>
    RANGE_ITERATOR(const T_RANGE& range, typename std::enable_if<std::is_same<typename T_RANGE::VECTOR_T,TV_INT>::value>::type* = nullptr)
        :min_corner(range.min_corner),max_corner(range.max_corner)
    {
        Reset();
    }

    RANGE_ITERATOR(const TV_INT& min_corner_input, const TV_INT& max_corner_input)
        :min_corner(min_corner_input),max_corner(max_corner_input)
    {
        Reset();
    }

    ~RANGE_ITERATOR() {}

    void Reset()
    {Idx=Rmin;}

    bool Valid() const
    {return Idx[0]<=Rmax[0];}

    void Next()
    {
        for(int i=d-1;i>=0;i--) if(Idx[i]<Rmax[i] || i==0){Idx[i]++;return;} else Idx[i]=Rmin[i];
    }

    const TV_INT& Index()
    {return index;}
};

//#####################################################################
}
#endif
