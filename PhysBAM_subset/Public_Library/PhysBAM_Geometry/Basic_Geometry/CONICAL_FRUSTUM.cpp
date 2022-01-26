//#####################################################################
// Copyright 2002-2013, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Jonathan Su, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONICAL_FRUSTUM
//##################################################################### 
#include <PhysBAM_Geometry/Basic_Geometry/CONICAL_FRUSTUM.h>
#include <PhysBAM_Geometry/Basic_Geometry/ORIENTED_BOX.h>
using namespace PhysBAM;
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV CONICAL_FRUSTUM<TV>::
Normal(const TV& location,const int aggregate) const
{
    assert(aggregate >= 1 && aggregate <= 3);
    if(aggregate == 1){
        TV radial_direction=(location-plane1.x1).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Normalized();
        return (plane1.normal*(radius2-radius1)+radial_direction*height).Normalized();}
    else if(aggregate == 2) return plane1.Normal();
    else return plane2.Normal();

}
//#####################################################################
// Function Inside
//#####################################################################
template<class TV> bool CONICAL_FRUSTUM<TV>::
Inside(const TV& location,const T thickness_over_two) const //TODO the tolerence breaks down when the frustum is very steep (almost planar). Should use real cylinder distance instead to fix that.
{
    if(!plane1.Inside(location,thickness_over_two) || !plane2.Inside(location,thickness_over_two)) return false;
    T alpha=TV::Dot_Product(location-plane1.x1,-plane1.normal)/height;
    T radius=(1-alpha)*radius1+alpha*radius2;
    return (location-plane1.x1).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() <= sqr(radius-thickness_over_two);
}
//#####################################################################
// Function Lazy_Inside
//#####################################################################
template<class TV> bool CONICAL_FRUSTUM<TV>::
Lazy_Inside(const TV& location) const 
{
    if(!plane1.Lazy_Inside(location) || !plane2.Lazy_Inside(location)) return false;
    T alpha=TV::Dot_Product(location-plane1.x1,-plane1.normal)/height;
    T radius=(1-alpha)*radius1+alpha*radius2;
    return (location-plane1.x1).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() <= sqr(radius);
}
//#####################################################################
// Function Outside
//#####################################################################
template<class TV> bool CONICAL_FRUSTUM<TV>::
Outside(const TV& location,const T thickness_over_two) const  
{
    if(plane1.Outside(location,thickness_over_two) || plane2.Outside(location,thickness_over_two)) return true;
    T alpha=TV::Dot_Product(location-plane1.x1,-plane1.normal)/height;
    T radius=(1-alpha)*radius1+alpha*radius2;
    return (location-plane1.x1).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() >= sqr(radius+thickness_over_two);
}
//#####################################################################
// Function Lazy_Outside
//#####################################################################
template<class TV> bool CONICAL_FRUSTUM<TV>::
Lazy_Outside(const TV& location) const  
{
    if(plane1.Lazy_Outside(location) || plane2.Lazy_Outside(location)) return true;
    T alpha=TV::Dot_Product(location-plane1.x1,-plane1.normal)/height;
    T radius=(1-alpha)*radius1+alpha*radius2;
    return (location-plane1.x1).Projected_Orthogonal_To_Unit_Direction(plane1.normal).Magnitude_Squared() >= sqr(radius);
}
//#####################################################################
// Function Boundary
//#####################################################################
template<class TV> bool CONICAL_FRUSTUM<TV>::
Boundary(const TV& location,const T thickness_over_two) const
{
    return !Inside(location,thickness_over_two) && !Outside(location,thickness_over_two);
}
//#####################################################################
// Function Surface
//#####################################################################
template<class TV> TV CONICAL_FRUSTUM<TV>::
Surface(const TV& location) const
{
    TV v=location-plane1.x1;T axial_distance=-TV::Dot_Product(v,plane1.normal);
    TV radial_direction=v.Projected_Orthogonal_To_Unit_Direction(plane1.normal).Normalized();
    TV normal_direction=(plane1.normal*(radius2-radius1)+radial_direction*height).Normalized();
    T radial_distance_minus_radius=TV::Dot_Product(v-radial_direction*radius1,normal_direction);
    if(radial_distance_minus_radius>0 || (radial_distance_minus_radius>-axial_distance && radial_distance_minus_radius>axial_distance-height)) // closest point is on infinite cylinder
        return plane1.x1+(radial_direction*radius1).Projected_On_Unit_Direction(normal_direction);
    if(axial_distance < height-axial_distance) return location+axial_distance*plane1.normal; // closest point is on plane1
    else return location-(height-axial_distance)*plane1.normal; // closest point is on plane2
}
//#####################################################################
// Function Signed_Distance
//#####################################################################
template<class TV> typename TV::SCALAR CONICAL_FRUSTUM<TV>::
Signed_Distance(const TV& location) const
{  
    TV v=location-plane1.x1;
    T plane1_distance=TV::Dot_Product(v,plane1.normal),plane_distance=max(plane1_distance,-height-plane1_distance);
    TV radial_direction=v.Projected_Orthogonal_To_Unit_Direction(plane1.normal).Normalized();
    TV normal_direction=(plane1.normal*(radius2-radius1)+radial_direction*height).Normalized();
    T cylinder_distance=TV::Dot_Product(v-radial_direction*radius1,normal_direction);
    return cylinder_distance>0 && plane_distance>0?sqrt(sqr(cylinder_distance)+sqr(plane_distance)):max(cylinder_distance,plane_distance);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV CONICAL_FRUSTUM<TV>::
Normal(const TV& location) const
{
    TV v=location-plane1.x1;
    T plane1_distance=TV::Dot_Product(v,plane1.normal),plane2_distance=-height-plane1_distance;
    TV radial_direction=v.Projected_Orthogonal_To_Unit_Direction(plane1.normal).Normalized();
    TV normal_direction=(plane1.normal*(radius2-radius1)+radial_direction*height).Normalized();
    TV infinite_cylinder_normal=v.Projected_On_Unit_Direction(normal_direction);
    T cylinder_distance=TV::Dot_Product(v-radial_direction*radius1,normal_direction);
    if(plane1_distance>=plane2_distance){
        if(cylinder_distance>0 && plane1_distance>0){
            T magnitude=sqrt(sqr(cylinder_distance)+sqr(plane1_distance));
            return cylinder_distance/magnitude*infinite_cylinder_normal+plane1_distance/magnitude*plane1.normal;}
        else if(cylinder_distance>plane1_distance) return infinite_cylinder_normal;
        return plane1.normal;}
    if(cylinder_distance>0 && plane2_distance>0){
        T magnitude=sqrt(sqr(cylinder_distance)+sqr(plane2_distance));
        return cylinder_distance/magnitude*infinite_cylinder_normal+plane2_distance/magnitude*plane2.normal;}
    else if(cylinder_distance>plane2_distance) return infinite_cylinder_normal;
    return plane2.normal;
}
//#####################################################################
// Function Bounding_Box
//#####################################################################
template<class TV> RANGE<TV> CONICAL_FRUSTUM<TV>::
Bounding_Box() const
{
    TV edge1=(T)2*max(radius1,radius2)*plane1.normal.Unit_Orthogonal_Vector();
    typename TV::SPIN edge2_spin=TV::Cross_Product(plane1.normal,edge1);
    TV edge2;if(TV::m==3)for(int d=1;d<=3;++d)edge2(d)=edge2_spin(d);
    TV corner=plane1.x1-(T).5*(edge1+edge2);
    MATRIX<T,TV::m> edges;edges.Column(1)=-height*plane1.normal;edges.Column(2)=edge1;if(TV::m==3)edges.Column(3)=edge2;
    return ORIENTED_BOX<TV>(corner,edges).Axis_Aligned_Bounding_Box();
}
//#####################################################################
namespace PhysBAM {
template class CONICAL_FRUSTUM<VECTOR<float,2> >;
template class CONICAL_FRUSTUM<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class CONICAL_FRUSTUM<VECTOR<double,2> >;
template class CONICAL_FRUSTUM<VECTOR<double,3> >;
#endif
}
