//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __EXTRAPOLATION_PDE_UNIFORM__
#define __EXTRAPOLATION_PDE_UNIFORM__

#include <PhysBAM_Tools/Arrays/ARRAYS_FORWARD.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_1D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_2D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_3D.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/LEVELSET_POLICY_UNIFORM.h>
namespace PhysBAM{

template<class T_GRID,class T2>
class EXTRAPOLATION_PDE_UNIFORM:public NONCOPYABLE
{
    typedef typename T_GRID::VECTOR_T TV;typedef typename TV::SCALAR T;
    typedef typename T_GRID::VECTOR_INT TV_INT;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_BASE T_ARRAYS_BASE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR T_ARRAYS_SCALAR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<TV>::TYPE T_ARRAYS_VECTOR;
    typedef typename T_ARRAYS_SCALAR::template REBIND<T2>::TYPE T_ARRAYS_T2;
    typedef typename T_ARRAYS_BASE::template REBIND<T2>::TYPE T_ARRAYS_T2_BASE;
    typedef typename GRID_ARRAYS_POLICY<T_GRID>::ARRAYS_SCALAR::template REBIND<bool>::TYPE T_ARRAYS_BOOL;
    typedef typename T_ARRAYS_BASE::template REBIND<bool>::TYPE T_ARRAYS_BOOL_BASE;
    typedef typename T_ARRAYS_BASE::template REBIND<VECTOR<bool,T_GRID::dimension> >::TYPE T_ARRAYS_BOOL_DIMENSION_BASE;typedef typename T_GRID::NODE_ITERATOR NODE_ITERATOR;
    typedef typename BOUNDARY_POLICY<T_GRID>::BOUNDARY_SCALAR T_BOUNDARY_SCALAR;
    typedef typename REBIND<T_BOUNDARY_SCALAR,T2>::TYPE T_BOUNDARY_T2;
    typedef typename LEVELSET_POLICY<GRID<TV> >::LEVELSET T_LEVELSET;
public:
    const T_LEVELSET& levelset_ghost;
    T_ARRAYS_T2_BASE& u_ghost;
    T band_width;
    int degree;
    int number_of_ghost_cells;

    T_GRID node_grid;
    ARRAY<T_ARRAYS_BOOL> initial_fixed_indices;
    ARRAY<T_ARRAYS_T2> initial_u_values;

    EXTRAPOLATION_PDE_UNIFORM(const T_GRID& grid,const T_LEVELSET& levelset_ghost_input,const int degree_input,T_ARRAYS_T2_BASE& u_ghost_input,const T_ARRAYS_BOOL_BASE& fixed_indices_ghost_input);
    ~EXTRAPOLATION_PDE_UNIFORM();

    void Set_Band_Width(const T number_of_cells=(T)5)
    {band_width=number_of_cells*node_grid.dX.Max();}

//#####################################################################
    void Extrapolate();
private:
    void Calculate_Initial_Conditions();
    void Extrapolate_Forward_Euler(const T distance,const T dt);
    void Extrapolate_Second_Order_Runge_Kutta(const T distance,const T dt);
//#####################################################################
}; 
}  
#endif
