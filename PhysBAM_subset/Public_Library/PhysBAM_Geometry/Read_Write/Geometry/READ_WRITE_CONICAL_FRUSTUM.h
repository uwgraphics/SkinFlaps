//#####################################################################
// Copyright 2013, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_CONICAL_FRUSTUM
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_CONICAL_FRUSTUM__
#define __READ_WRITE_CONICAL_FRUSTUM__

#include <PhysBAM_Tools/Read_Write/Utilities/READ_WRITE.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Geometry/Basic_Geometry/CONICAL_FRUSTUM.h>
namespace PhysBAM{

template<class RW,class T,int d>
class Read_Write<CONICAL_FRUSTUM<VECTOR<T,d> >,RW>
{
public:
    static void Read(std::istream& input,CONICAL_FRUSTUM<VECTOR<T,d> >& object)
    {Read_Binary<RW>(input,object.plane1.x1,object.plane2.x1,object.radius1,object.radius2);
    object.Set_Endpoints(object.plane1.x1,object.plane2.x1);}

    static void Write(std::ostream& output,const CONICAL_FRUSTUM<VECTOR<T,d> >& object)
    {Write_Binary<RW>(output,object.plane1.x1,object.plane2.x1,object.radius1,object.radius2);}
};
}

#endif
#endif
