//#####################################################################
// Copyright 2009, Avi Robinson-Mosher.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class READ_WRITE_IMPLICIT_OBJECT_COMBINED
//#####################################################################
#ifndef COMPILE_WITHOUT_READ_WRITE_SUPPORT
#ifndef __READ_WRITE_IMPLICIT_OBJECT_COMBINED__
#define __READ_WRITE_IMPLICIT_OBJECT_COMBINED__

#include <PhysBAM_Geometry/Implicit_Objects_Uniform/IMPLICIT_OBJECT_COMBINED.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
namespace PhysBAM{

template<class RW,class TV>
class Read_Write<IMPLICIT_OBJECT_COMBINED<TV>,RW>:public Read_Write<IMPLICIT_OBJECT<TV>,RW>
{
public:
    static void Read_Helper(std::istream& input,STRUCTURE<TV>& structure_object)
    {IMPLICIT_OBJECT_COMBINED<TV>& object=dynamic_cast<IMPLICIT_OBJECT_COMBINED<TV>&>(structure_object);
    Read_Binary<RW>(input,object.implicit_object2);object.implicit_object1=object.implicit_object2;object.Update_Box();object.Update_Minimum_Cell_Size();}

    static void Write_Helper(std::ostream& output,const STRUCTURE<TV>& structure_object)
    {typedef typename TV::SCALAR T;typedef VECTOR<int,TV::dimension> TV_INT;
    const IMPLICIT_OBJECT_COMBINED<TV>& object=dynamic_cast<const IMPLICIT_OBJECT_COMBINED<TV>&>(structure_object);
    GRID<TV> grid(object.Minimum_Cell_Size(),object.box);
    ARRAY<T,TV_INT> phi=(grid.Domain_Indices());
    LEVELSET_IMPLICIT_OBJECT<TV> levelset(grid,phi);
    for(typename GRID<TV>::CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){phi(iterator.Cell_Index())=object.Extended_Phi(grid.X(iterator.Cell_Index()));}
    Write_Binary<RW>(output,levelset);}
};
}
#endif
#endif
