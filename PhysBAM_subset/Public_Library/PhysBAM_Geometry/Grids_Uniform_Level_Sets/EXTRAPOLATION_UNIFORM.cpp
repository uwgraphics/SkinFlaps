//####################################################################
// Copyright 2002-2005, Ronald Fedkiw, Eran Guendelman, Geoffrey Irving, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EXTRAPOLATION_UNIFORM  
//##################################################################### 
#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Grids_Uniform_Boundaries/BOUNDARY_UNIFORM.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_2X2.h>
#include <PhysBAM_Tools/Matrices/SYMMETRIC_MATRIX_3X3.h>
#include <PhysBAM_Geometry/Grids_Uniform_Level_Sets/EXTRAPOLATION_UNIFORM.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_GRID,class T2> EXTRAPOLATION_UNIFORM<T_GRID,T2>::
EXTRAPOLATION_UNIFORM(const T_GRID& grid,const T_ARRAYS_BASE& phi_input,T_ARRAYS_T2_BASE& u_input,const int ghost_cells_input)
    :u(u_input),phi(phi_input),seed_indices(0),seed_done(0),collision_aware_extrapolation(false),neighbors_visible(0),weights(0),ghost_cells(ghost_cells_input)
{
    Set_Band_Width();
    Set_Isobaric_Fix_Width();
    Set_Small_Number();
    node_grid=grid.Is_MAC_Grid()?grid.Get_Regular_Grid_At_MAC_Positions():grid;
    TV DX2=node_grid.dX*node_grid.dX;
    if(T_GRID::dimension>=2){optimization_scale[2]=DX2[2]/DX2[1]; // dy^2/dx^2
        if(T_GRID::dimension==3){optimization_scale[0]=DX2[3]/DX2[2];optimization_scale[1]=DX2[3]/DX2[1];}} // dz^2/dy^2 and dz^2/dx^2 
    RANGE<TV_INT> domain_indices=node_grid.Domain_Indices().Thickened(ghost_cells);dimension_start=domain_indices.Minimum_Corner();dimension_end=domain_indices.Maximum_Corner();
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_GRID,class T2> EXTRAPOLATION_UNIFORM<T_GRID,T2>::
~EXTRAPOLATION_UNIFORM()
{} 
//#####################################################################
// Function Extrapolate
//#####################################################################
template<class T_GRID,class T2> void EXTRAPOLATION_UNIFORM<T_GRID,T2>::
Extrapolate(const T time,const bool fill_ghost_cells,const bool interpolate_high_low)
{
    T_ARRAYS_SCALAR phi_ghost(node_grid.Domain_Indices(ghost_cells),false);T_ARRAYS_T2 u_ghost(node_grid.Domain_Indices(ghost_cells),false); 
    if(fill_ghost_cells){
        T_BOUNDARY_SCALAR phi_boundary;phi_boundary.Fill_Ghost_Cells(node_grid,phi,phi_ghost,0,time,ghost_cells);
        boundary->Fill_Ghost_Cells(node_grid,u,u_ghost,0,time,ghost_cells);}
    else{
        phi_ghost=phi;
        u_ghost=u;}
    if(isobaric_fix_width) phi_ghost+=isobaric_fix_width;   

    int heap_length=0;
    T_ARRAYS_BOOL done(node_grid.Domain_Indices(ghost_cells+1)); // extra cell for Update_Close_Point
    T_ARRAYS_BOOL close(node_grid.Domain_Indices(ghost_cells));
    ARRAY<TV_INT> heap(close.array.Size());
    Initialize(phi_ghost,done,close,heap,heap_length);
    if(weights){weights->Clean_Memory();weights->Resize(node_grid.Domain_Indices(ghost_cells),true,false);}

    while(heap_length && phi_ghost(heap(1)) <= band_width+isobaric_fix_width){
        TV_INT index=Remove_Root_From_Heap(phi_ghost,heap,heap_length,close);
        done(index)=true;
        Update_Close_Point(u_ghost,phi_ghost,done,index,interpolate_high_low);
    
        if(collision_aware_extrapolation){
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && !close(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector))
                    Add_To_Heap(phi_ghost,heap,heap_length,close,index-axis_vector);
                if(index[axis] != dimension_end[axis] && !done(index+axis_vector) && !close(index+axis_vector) && Neighbor_Visible(axis,index))
                    Add_To_Heap(phi_ghost,heap,heap_length,close,index+axis_vector);}}
        else{
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && !close(index-axis_vector))
                    Add_To_Heap(phi_ghost,heap,heap_length,close,index-axis_vector);
                if(index[axis] != dimension_end[axis] && !done(index+axis_vector) && !close(index+axis_vector))
                    Add_To_Heap(phi_ghost,heap,heap_length,close,index+axis_vector);}}}

    T_ARRAYS_T2_BASE::Get(u,u_ghost);
    boundary->Apply_Boundary_Condition(node_grid,u,time);
}
//#####################################################################
// Function Initialize
//#####################################################################
// pass heap_length by reference
template<class T_GRID,class T2> void EXTRAPOLATION_UNIFORM<T_GRID,T2>::
Initialize(const T_ARRAYS_BASE& phi,T_ARRAYS_BOOL_BASE& done,T_ARRAYS_BOOL_BASE& close,ARRAY<TV_INT>& heap,int& heap_length)
{  
    assert(!seed_indices||!seed_done);
    if(seed_indices) for(int i=1;i<=seed_indices->m;i++) done((*seed_indices)(i))=true;
    else if(seed_done) T_ARRAYS_BOOL_BASE::Put(*seed_done,done);
    else for(NODE_ITERATOR iterator(node_grid,ghost_cells);iterator.Valid();iterator.Next()) if(phi(iterator.Node_Index()) <= 0) done(iterator.Node_Index())=true;

    // find neighbors of done nodes which have positive phi
    if(collision_aware_extrapolation){
        for(NODE_ITERATOR iterator(node_grid,ghost_cells);iterator.Valid();iterator.Next()) if(done(iterator.Node_Index())){TV_INT index=iterator.Node_Index();
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && !close(index-axis_vector) && Neighbor_Visible(axis,index-axis_vector))
                    Add_To_Heap(phi,heap,heap_length,close,index-axis_vector);
                if(index[axis] != dimension_end[axis] && !done(index+axis_vector) && !close(index+axis_vector) && Neighbor_Visible(axis,index))
                    Add_To_Heap(phi,heap,heap_length,close,index+axis_vector);}}}
    else{
        for(NODE_ITERATOR iterator(node_grid,ghost_cells);iterator.Valid();iterator.Next()) if(done(iterator.Node_Index())){TV_INT index=iterator.Node_Index();
            for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis);
                if(index[axis] != dimension_start[axis] && !done(index-axis_vector) && !close(index-axis_vector))
                    Add_To_Heap(phi,heap,heap_length,close,index-axis_vector);
                if(index[axis] != dimension_end[axis] && !done(index+axis_vector) && !close(index+axis_vector))
                    Add_To_Heap(phi,heap,heap_length,close,index+axis_vector);}}}
}
template<class T,int d>
void Add_Weight(ARRAY<ARRAY<PAIR<VECTOR<int,d>,T> >,VECTOR<int,d> >* weights,const VECTOR<int,d>& index,const VECTOR<int,d>& closest_node_index,T weight=(T)1)
{
    if(weights)weights->operator()(index).Append(Tuple(closest_node_index,weight));
}
template<class T,int d>
void Scale_Weight(ARRAY<ARRAY<PAIR<VECTOR<int,d>,T> >,VECTOR<int,d> >* weights,const VECTOR<int,d>& index,const T& fraction)
{
    typedef VECTOR<int,d> TV_INT;
    if(weights){ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& w=*weights;
        int node_index=1;bool node_half_done=false;
        for(int k=1;k<=w(index).m;++k){
            bool full_node=(w(index)(k).y>(T).75); //==1
            assert(node_index<=2);assert(!(full_node&&node_half_done));
            T scaling(node_index==1?fraction:1-fraction);
            w(index)(k).y*=scaling;
            ++node_index;
            if(!full_node){
                node_half_done=!node_half_done;
                if(node_half_done)--node_index;}}
    }
}
template<class T,int d>
void Scale_Weight(ARRAY<ARRAY<PAIR<VECTOR<int,d>,T> >,VECTOR<int,d> >* weights,const VECTOR<int,d>& index,const T& fraction_1,const T& fraction_2)
{
    typedef VECTOR<int,d> TV_INT;
    if(weights){ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& w=*weights;
        int node_index=1;bool node_half_done=false;
        for(int k=1;k<=w(index).m;++k){
            bool full_node=(w(index)(k).y>(T).75); //==1
            assert(node_index<=3);assert(!(full_node&&node_half_done));
            T scaling(node_index==3?((T)1-fraction_1-fraction_2):(node_index==1?fraction_1:fraction_2));
            w(index)(k).y*=scaling;
            ++node_index;
            if(!full_node){
                node_half_done=!node_half_done;
                if(node_half_done)--node_index;}}
    }
}
//#####################################################################
// Function Update_Close_Point
//##################################################################### 
// needs done=0 around the outside of the domain
// note that sqrt(phix^2+phiy^2+phiz^2)=1 if it's a distance function
template<class T_GRID,class T2> void EXTRAPOLATION_UNIFORM<T_GRID,T2>::
Update_Close_Point(T_ARRAYS_T2_BASE& u,const T_ARRAYS_BASE& phi,const T_ARRAYS_BOOL_BASE& done,const TV_INT& index,const bool interpolate_high_low)
{
    //T2 value[T_GRID::dimension]={T2()}; // the value to use in the given direction
    T2 value[T_GRID::dimension]; // copy initialization not allowed for explicit VECTOR() constructor
    for(int i=0;i<T_GRID::dimension;i++) value[i]=T2();
    T phix_dx[T_GRID::dimension]={0}; // the difference in phi value for the direction
    int number_of_axis=0; // the number of axis that we want to use later
    int missing_axis=3; // used in number_of_axis==2 case only, so it gives you which axis is missing (==3 for 2d)

    // check each principal axis
    for(int axis=1;axis<=T_GRID::dimension;axis++){TV_INT axis_vector=TV_INT::Axis_Vector(axis),low=index-axis_vector,high=index+axis_vector;
        bool check_low=done(low),check_high=done(high);
        if(collision_aware_extrapolation){
            if(check_low && !Neighbor_Visible(axis,low)) check_low=false;
            if(check_high && !Neighbor_Visible(axis,index)) check_high=false;}
        if(!check_low){
            if(check_high){value[number_of_axis]=u(high);phix_dx[number_of_axis]=phi(index)-phi(high);Add_Weight(weights,index,high);}
            else{missing_axis=axis;number_of_axis--;}}
        else if(!check_high){value[number_of_axis]=u(low);phix_dx[number_of_axis]=phi(index)-phi(low);Add_Weight(weights,index,low);}
        else{
            if(phi(low)<=phi(high)){if(interpolate_high_low){value[number_of_axis]=((T).5)*(u(low)+u(high));Add_Weight(weights,index,low,(T).5);Add_Weight(weights,index,high,(T).5);}
                else{value[number_of_axis]=u(low);Add_Weight(weights,index,low);}
                phix_dx[number_of_axis]=phi(index)-phi(low);}
            else{if(interpolate_high_low){value[number_of_axis]=((T).5)*(u(high)+u(low));Add_Weight(weights,index,high,(T).5);Add_Weight(weights,index,low,(T).5);}
                else{value[number_of_axis]=u(high);Add_Weight(weights,index,high);}
                phix_dx[number_of_axis]=phi(index)-phi(high);}}
        number_of_axis++;}
       
    assert(number_of_axis);
    if(number_of_axis==1) u(index)=value[0];
    else if(T_GRID::dimension==2 || number_of_axis==2){
        T a=phix_dx[0]*optimization_scale[missing_axis-1],b=phix_dx[1],denominator=a+b,fraction=(T).5;
        if(denominator > small_number) fraction=clamp(a/denominator,(T)0,(T)1);
        Scale_Weight(weights,index,fraction);
        u(index)=fraction*value[0]+(1-fraction)*value[1];}
    else{PHYSBAM_ASSERT(T_GRID::dimension==3); // should only get here in 3D
        T a=phix_dx[0]*optimization_scale[1],b=phix_dx[1]*optimization_scale[0],c=phix_dx[2],denominator=a+b+c;
        T fraction_1=(T)one_third,fraction_2=(T)one_third;
        if(denominator > small_number){fraction_1=clamp(a/denominator,(T)0,(T)1);fraction_2=clamp(b/denominator,(T)0,1-fraction_1);}
        Scale_Weight(weights,index,fraction_1,fraction_2);
        u(index)=fraction_1*value[0]+fraction_2*value[1]+((T)1-fraction_1-fraction_2)*value[2];}
    
    if(weights){ARRAY<ARRAY<PAIR<TV_INT,T> >,TV_INT>& w=*weights;
        const int w_index_size_pre_insert=w(index).m;
        ARRAY<bool> node_to_delete(w_index_size_pre_insert);node_to_delete.Fill(false);
        for(int j=1;j<=w_index_size_pre_insert;++j){const TV_INT jndex=w(index)(j).x;
            assert(jndex!=index);
            if(w(index)(j).y<1e-8) node_to_delete(j)=true;
            else if(w(jndex).m){
                node_to_delete(j)=true;
                for(int k=1;k<=w(jndex).m;++k){const TV_INT kndex=w(jndex)(k).x;
                    const T weight_product=w(index)(j).y*w(jndex)(k).y;
                    if(weight_product>=1e-8)w(index).Append(Tuple(kndex,weight_product));}}}
        // delete the marked compound jndex from w(index) and merge duplicated kndex's in w(index)
        for(int j=w_index_size_pre_insert;j>=1;--j)if(node_to_delete(j))w(index).Remove_Index_Lazy(j);
        for(int k=1;k<=w(index).m;++k){const TV_INT kndex=w(index)(k).x;
            for(int h=w(index).m;h>k;--h){const TV_INT hndex=w(index)(h).x;
                if(hndex==kndex){
                    w(index)(k).y+=w(index)(h).y;
                    w(index).Remove_Index_Lazy(h);}}}
        T weight_sum(0);for(int k=1;k<=w(index).m;k++)weight_sum+=w(index)(k).y;
        for(int k=1;k<=w(index).m;++k)w(index)(k).y/=weight_sum;}
}
//#####################################################################
namespace PhysBAM {
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,1> >,float>;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,2> >,float>;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,3> >,float>;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,1> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,1> >,VECTOR<float,3> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,2> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,2> >,VECTOR<float,4> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,3> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,3> >,VECTOR<float,5> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,2> >,SYMMETRIC_MATRIX<float,2> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<float,3> >,SYMMETRIC_MATRIX<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,1> >,double>;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,2> >,double>;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,3> >,double>;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,1> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,1> >,VECTOR<double,3> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,2> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,2> >,VECTOR<double,4> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,3> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,3> >,VECTOR<double,5> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,2> >,SYMMETRIC_MATRIX<double,2> >;
template class EXTRAPOLATION_UNIFORM<GRID<VECTOR<double,3> >,SYMMETRIC_MATRIX<double,3> >;
#endif
}
