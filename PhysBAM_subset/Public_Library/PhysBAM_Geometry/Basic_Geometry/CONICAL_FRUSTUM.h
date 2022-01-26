//#####################################################################
// Copyright 2002-2013, Ronald Fedkiw, Geoffrey Irving, Neil Molino, Jonathan Su, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CONICAL_FRUSTUM
//##################################################################### 
#ifndef __CONICAL_FRUSTUM__
#define __CONICAL_FRUSTUM__

#include <PhysBAM_Geometry/Basic_Geometry/PLANE.h>
namespace PhysBAM{


template<class T>
class PLANE_2D
{
    typedef VECTOR<T,2> TV;
public:
    typedef TV VECTOR_T;
    TV normal;
    TV x1;
    PLANE_2D() :normal(0,1),x1(0,0) {}
    PLANE_2D(const TV& normal_input,const TV& x1_input) :normal(normal_input),x1(x1_input) {}
    TV Normal() const {return normal;}
    T Signed_Distance(const TV& location) const {return TV::Dot_Product(normal,location-x1);}
    // inside is the half space behind the normal
    bool Inside(const TV& location,const T thickness_over_two) const {return Signed_Distance(location)<=-thickness_over_two;}
    bool Lazy_Inside(const TV& location) const {return Signed_Distance(location)<=0;}
    bool Outside(const TV& location,const T thickness_over_two) const {return !Inside(location,-thickness_over_two);}
    bool Lazy_Outside(const TV& location) const {return !Lazy_Inside(location);}
};

template<class TV>
class CONICAL_FRUSTUM
{
    typedef typename TV::SCALAR T;
public:
    typedef TV VECTOR_T;

    T height,radius1,radius2;
    typename IF<TV::m==3,PLANE<T>,PLANE_2D<T> >::TYPE plane1,plane2; // plane2 is height units behind circle

    CONICAL_FRUSTUM()
        :height(1),radius1(1),radius2(.5),plane2(plane1.x1-height*plane1.normal,-plane1.normal)
    {}

    CONICAL_FRUSTUM(const TV& point1,const TV& point2,const T radius1,const T radius2)
        :radius1(radius1),radius2(radius2)
    {
        Set_Endpoints(point1,point2);
    }

    void Set_Height(const T height_input)
    {height=height_input;plane2.x1=plane1.x1-height*plane1.normal;}

    void Set_Endpoints(const TV& point1,const TV& point2)
    {TV height_vector=point1-point2;height=height_vector.Normalize();
    plane1.x1=point1;plane1.normal=height_vector;
    plane2.x1=point2;plane2.normal=-plane1.normal;}

    VECTOR<T,TV::m-1> Principal_Curvatures(const TV& X) const
    {VECTOR<T,TV::m-1> result;
    if(TV::m<3) return result;
    TV v=X-plane1.x1;
    T plane1_distance=TV::Dot_Product(v,plane1.normal),plane_distance=max(plane1_distance,-height-plane1_distance); //TODO fix this
    TV radial_direction=v.Projected_Orthogonal_To_Unit_Direction(plane1.normal).Normalized();
    TV normal_direction=(plane1.normal*(radius2-radius1)+radial_direction*height).Normalized();
    T cylinder_distance=TV::Dot_Product(v-radial_direction*radius1,normal_direction);
    if(abs(plane_distance)<abs(cylinder_distance)) return result;
    T a=abs((radius1-radius2)/height);
    T u=TV::Dot_Product(radial_direction,(radial_direction*radius1).Projected_On_Unit_Direction(normal_direction))/a;
    result(1)=-1/(a*sqrt(1+a*a)*u);
    return result;}

    static std::string Name() 
    {return "CONICAL_FRUSTUM<T>";}

//#####################################################################
    TV Normal(const TV& location) const;
    TV Normal(const TV& location,const int aggregate) const;
    bool Inside(const TV& location,const T thickness_over_two) const;
    bool Lazy_Inside(const TV& location) const;
    bool Outside(const TV& location,const T thickness_over_two) const;
    bool Lazy_Outside(const TV& location) const;
    bool Boundary(const TV& location,const T thickness_over_two) const;
    TV Surface(const TV& location) const;
    T Signed_Distance(const TV& location) const;
    RANGE<TV> Bounding_Box() const;
//#####################################################################
};   
}
#endif
