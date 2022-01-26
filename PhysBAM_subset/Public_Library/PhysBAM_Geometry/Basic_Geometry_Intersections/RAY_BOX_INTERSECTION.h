//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace INTERSECTION
//##################################################################### 
#ifndef __RAY_BOX_INTERSECTION__
#define __RAY_BOX_INTERSECTION__
#include <PhysBAM_Geometry/Basic_Geometry/BASIC_GEOMETRY_FORWARD.h>
namespace PhysBAM{
namespace INTERSECTION{
//#####################################################################
template<class T> bool Intersects(RAY<VECTOR<T,1> >& ray,const RANGE<VECTOR<T,1> >& box, const T thickness_over_two=(T)0);
template<class T> bool Intersects(RAY<VECTOR<T,2> >& ray,const RANGE<VECTOR<T,2> >& box, const T thickness_over_two=(T)0,const T segment_intersect_epsilon=(T)0);
template<class T> bool Intersects(RAY<VECTOR<T,3> >& ray,const RANGE<VECTOR<T,3> >& box, const T thickness_over_two=(T)0);
template<class T> bool Lazy_Outside(RAY<VECTOR<T,3> >& ray,const RANGE<VECTOR<T,3> >& box);
template<class T> bool Get_Intersection_Range(const RAY<VECTOR<T,3> >& ray,const BOX<VECTOR<T,3> >& box,T& start_t,T& end_t);
template<class TV> bool Get_Intersection_Range(const RAY<TV>& ray,const RANGE<TV>& box,typename TV::SCALAR& start_t,typename TV::SCALAR& end_t);
template<class T> bool Lazy_Intersects(RAY<VECTOR<T,3> >& ray,const RANGE<VECTOR<T,3> >& box, const T thickness_over_two=(T)0);
//#####################################################################
};
};
#endif
