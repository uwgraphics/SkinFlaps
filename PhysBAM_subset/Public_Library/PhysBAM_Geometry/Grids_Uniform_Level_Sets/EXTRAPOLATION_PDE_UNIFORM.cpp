//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_PDE_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> EXTRAPOLATION_PDE_UNIFORM<T_GRID,T2>::
EXTRAPOLATION_PDE_UNIFORM(const T_GRID& grid,const T_LEVELSET& levelset_ghost_input,const int degree_input,T_ARRAYS_T2_BASE& u_ghost_input,const T_ARRAYS_BOOL_BASE& fixed_indices_ghost_input):
    levelset_ghost(levelset_ghost_input),u_ghost(u_ghost_input),degree(degree_input),number_of_ghost_cells(2*degree+2)
{
    node_grid=grid.Is_MAC_Grid()?grid.Get_Regular_Grid_At_MAC_Positions():grid;
    Set_Band_Width();
    initial_u_values.Resize(degree+1,false);
    initial_fixed_indices.Resize(degree+1,false);
    for(int i=1;i<=initial_u_values.m;i++){
        initial_u_values(i).Resize(node_grid.Domain_Indices(number_of_ghost_cells+1-i),false);
        initial_fixed_indices(i).Resize(node_grid.Domain_Indices(number_of_ghost_cells+1-i),false);}

    T_ARRAYS_BOOL::Get(initial_fixed_indices(1),fixed_indices_ghost_input);
    T_ARRAYS_T2::Get(initial_u_values(1),u_ghost_input);
    for(int i=2;i<=initial_u_values.m;i++){
        initial_u_values(i).Fill((T)0.);
        initial_fixed_indices(i).Fill(true);}
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> EXTRAPOLATION_PDE_UNIFORM<T_GRID,T2>::
~EXTRAPOLATION_PDE_UNIFORM()
{
}
//#####################################################################
// Extrapolate
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_PDE_UNIFORM<T_GRID,T2>::
Extrapolate()
{
    //TODO: Improve normal computations
    Calculate_Initial_Conditions();
    Extrapolate_Forward_Euler(band_width,node_grid.dX.Max()/(TV::dimension+1));
}
//#####################################################################
// Calculate_Initial_Conditions
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_PDE_UNIFORM<T_GRID,T2>::
Calculate_Initial_Conditions()
{
    for(typename T_GRID::NODE_ITERATOR iterator(node_grid,number_of_ghost_cells);iterator.Valid();iterator.Next()){TV_INT node_index=iterator.Node_Index();
        if(levelset_ghost.Extended_Phi(iterator.Location())>band_width){
            initial_fixed_indices(1)(node_index)=true;
            initial_u_values(1)(node_index)=(T)0.;}}

    for(int i=2;i<=initial_fixed_indices.m;i++){
        for(typename T_GRID::NODE_ITERATOR iterator(node_grid,number_of_ghost_cells+1-i);iterator.Valid();iterator.Next()){TV_INT node_index=iterator.Node_Index();
            for(int axis=1;axis<=TV::dimension;axis++){
                TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                initial_fixed_indices(i)(node_index)&=initial_fixed_indices(i-1)(node_index+axis_vector)&&initial_fixed_indices(i-1)(node_index-axis_vector);}}}

    for(int i=2;i<=initial_u_values.m;i++){
        for(typename T_GRID::NODE_ITERATOR iterator(node_grid,number_of_ghost_cells+1-i);iterator.Valid();iterator.Next()){
            TV_INT node_index=iterator.Node_Index();
            TV normal=levelset_ghost.Extended_Normal(iterator.Location());
            if(initial_fixed_indices(i)(node_index)){
                for(int axis=1;axis<=TV::dimension;axis++){
                    TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                    initial_u_values(i)(node_index)+=normal(axis)*(initial_u_values(i-1)(node_index+axis_vector)-initial_u_values(i-1)(node_index-axis_vector))*((T)0.5)*node_grid.one_over_dX(axis);}}
            else initial_u_values(i)(node_index)=(T)0.0;}}
}
//#####################################################################
// Extrapolate_Forward_Euler
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_PDE_UNIFORM<T_GRID,T2>::
Extrapolate_Forward_Euler(const T distance,const T dt)
{
    for(int i=initial_u_values.m;i>=1;i--){
        for(int iteration=1;iteration<=2*ceil(distance/dt);iteration++){
            T_ARRAYS_T2 advection_term;advection_term.Resize(node_grid.Domain_Indices(number_of_ghost_cells),false);advection_term.Fill((T)0.);
            for(typename GRID<TV>::NODE_ITERATOR iterator(node_grid,i-1);iterator.Valid();iterator.Next()){TV_INT node_index=iterator.Node_Index();
                if(!initial_fixed_indices(i)(node_index)){
                    TV node_normal=levelset_ghost.Extended_Normal(iterator.Location());
                    for(int axis=1;axis<=TV::dimension;axis++){
                        T du_plus_three_halves=(initial_u_values(i)(node_index+TV_INT::Axis_Vector(axis)*2)-initial_u_values(i)(node_index+TV_INT::Axis_Vector(axis)))*node_grid.one_over_dX(axis);
                        T du_plus_one_half=(initial_u_values(i)(node_index+TV_INT::Axis_Vector(axis))-initial_u_values(i)(node_index))*node_grid.one_over_dX(axis);
                        T du_minus_one_half=(initial_u_values(i)(node_index)-initial_u_values(i)(node_index-TV_INT::Axis_Vector(axis)))*node_grid.one_over_dX(axis);
                        T du_minus_three_halves=(initial_u_values(i)(node_index-TV_INT::Axis_Vector(axis))-initial_u_values(i)(node_index-TV_INT::Axis_Vector(axis)*2))*node_grid.one_over_dX(axis);
                        if(node_normal(axis)<0){
                            advection_term(node_index)+=node_normal(axis)*du_plus_one_half;
                            T d2u_dx=du_plus_one_half-du_minus_one_half;
                            T d2u_plus_one_dx=du_plus_three_halves-du_plus_one_half;
                            if(abs(d2u_dx)>abs(d2u_plus_one_dx)) advection_term(node_index)-=node_normal(axis)*d2u_plus_one_dx/2;
                            else advection_term(node_index)-=node_normal(axis)*d2u_dx/2;
                        }
                        else{
                            advection_term(node_index)+=node_normal(axis)*du_minus_one_half;
                            T d2u_dx=du_plus_one_half-du_minus_one_half;
                            T d2u_minus_one_dx=du_minus_one_half-du_minus_three_halves;
                            if(abs(d2u_dx)>abs(d2u_minus_one_dx)) advection_term(node_index)+=node_normal(axis)*d2u_minus_one_dx/2;
                            else advection_term(node_index)+=node_normal(axis)*d2u_dx/2;
                        }
                    }
                }
            }
            for(typename GRID<TV>::NODE_ITERATOR iterator(node_grid,i-1);iterator.Valid();iterator.Next()){TV_INT node_index=iterator.Node_Index();
                if(!initial_fixed_indices(i)(node_index)){
                    initial_u_values(i)(node_index)+=dt*(-advection_term(node_index));
                    if(i<initial_u_values.m) initial_u_values(i)(node_index)+=dt*initial_u_values(i+1)(node_index);
                }
            }
        }
    }
    T_ARRAYS_T2::Get(u_ghost,initial_u_values(1));
}
//#####################################################################
// Extrapolate_Second_Order_Runge_Kutta
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_PDE_UNIFORM<T_GRID,T2>::
Extrapolate_Second_Order_Runge_Kutta(const T distance,const T dt)
{
}
namespace PhysBAM{
template class EXTRAPOLATION_PDE_UNIFORM<GRID<VECTOR<float,1> >,float>;
template class EXTRAPOLATION_PDE_UNIFORM<GRID<VECTOR<float,2> >,float>;
template class EXTRAPOLATION_PDE_UNIFORM<GRID<VECTOR<float,3> >,float>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EXTRAPOLATION_PDE_UNIFORM<GRID<VECTOR<double,1> >,double>;
template class EXTRAPOLATION_PDE_UNIFORM<GRID<VECTOR<double,2> >,double>;
template class EXTRAPOLATION_PDE_UNIFORM<GRID<VECTOR<double,3> >,double>;
#endif
}
