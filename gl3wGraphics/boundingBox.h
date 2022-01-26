//#####################################################################
// Author: Court Cutting
// Date: 12/14/07
// Purpose: 3D bounding box utility class
//#####################################################################

#ifndef __boundingBox__
#define __boundingBox__

#include <math.h>
#include <limits>
#include <algorithm>
#include <float.h>

template<class T> class boundingBox
{
public:
	union {
		struct
		{
			T xmin, xmax, ymin, ymax, zmin, zmax;
		};
		T val[6];
	};

    boundingBox()
        :xmin(1),xmax(-1),ymin(1),ymax(-1),zmin(1),zmax(-1)
    {}

    boundingBox(const T xmin_input,const T xmax_input,const T ymin_input,const T ymax_input,const T zmin_input,const T zmax_input)
        :xmin(xmin_input),xmax(xmax_input),ymin(ymin_input),ymax(ymax_input),zmin(zmin_input),zmax(zmax_input)
    {}

    boundingBox(const T (&minimum_corner)[3],const T (&maximum_corner)[3])
        :xmin(minimum_corner[0]),xmax(maximum_corner[0]),ymin(minimum_corner[1]),ymax(maximum_corner[1]),zmin(minimum_corner[2]),zmax(maximum_corner[2])
    {}

    template<class T2> explicit boundingBox(const boundingBox<T2>& box)
       :xmin((T)box.xmin),xmax((T)box.xmax),ymin((T)box.ymin),ymax((T)box.ymax),zmin((T)box.zmin),zmax((T)box.zmax)
    {}

    void Empty_Box()
	{xmin=FLT_MAX; xmax=-FLT_MAX; ymin=FLT_MAX; ymax=-FLT_MAX; zmin=FLT_MAX; zmax=-FLT_MAX;}

    bool IsEmpty() const
    {return xmax<xmin || ymax<ymin || zmax<zmin;}

    bool operator==(const boundingBox<T>& b) const
    {return xmin==b.xmin&&xmax==b.xmax&&ymin==b.ymin&&ymax==b.ymax&&zmin==b.zmin&&zmax==b.zmax;}

    bool operator!=(const boundingBox<T>& b) const
    {return !(*this==b);}

    boundingBox<T>& operator+=(const boundingBox<T>& r)
    {xmin+=r.xmin;xmax+=r.xmax;ymin+=r.ymin;ymax+=r.ymax;zmin+=r.zmin;zmax+=r.zmax;return *this;}

    boundingBox<T> operator+(const boundingBox<T>& r) const
    {return boundingBox<T>(xmin+r.xmin,xmax+r.xmax,ymin+r.ymin,ymax+r.ymax,zmin+r.zmin,zmax+r.zmax);}

    boundingBox<T> operator*(const T a) const
    {return a>=0?boundingBox<T>(xmin*a,xmax*a,ymin*a,ymax*a,zmin*a,zmax*a):boundingBox<T>(xmax*a,xmin*a,ymax*a,ymin*a,zmax*a,zmin*a);}

    boundingBox<T>& operator*=(const T a)
    {return *this=*this*a;}

    boundingBox<T> operator/(const T a) const
    {assert(a!=0);return *this*(1/a);}

    boundingBox<T>& operator/=(const T a)
    {return *this=*this/a;}

    void Edge_Lengths(T (&lengths)[3]) const
    {lengths[0]=xmax-xmin; lengths[1]=ymax-ymin; lengths[2]=zmax-zmin;}

    void Center(T (&center)[3]) const
    {center[0]=(T).5*(xmin+xmax); center[1]=(T).5*(ymin+ymax); center[2]=(T).5*(zmin+zmax);}
    
    void Minimum_Corner(T (&minimum_corner)[3]) const
    {minimum_corner[0]=xmin; minimum_corner[1]=ymin; minimum_corner[2]=zmin;}

    void Maximum_Corner(T (&maximum_corner)[3]) const
    {maximum_corner[0]=xmax; maximum_corner[1]=ymax; maximum_corner[2]=zmax;}

    T Volume() const
    {return (xmax-xmin)*(ymax-ymin)*(zmax-zmin);}

    void Enlarge_To_Include_Point(const T (&point)[3])
	{
#ifdef NOMINMAX
		xmin = std::min(xmin, point[0]);
		xmax = std::max(xmax, point[0]); ymin = std::min(ymin, point[1]); ymax = std::max(ymax, point[1]); zmin = std::min(zmin, point[2]); zmax = std::max(zmax, point[2]);
#else
		xmin = min(xmin, point[0]); xmax = max(xmax, point[0]); ymin = min(ymin, point[1]); ymax = max(ymax, point[1]); zmin = min(zmin, point[2]); zmax = max(zmax, point[2]);
#endif
	}

    void Enlarge_To_Include_Box(const boundingBox<T>& box)
    {
#ifdef NOMINMAX
		xmin = std::min(xmin, box.xmin); xmax = std::max(xmax, box.xmax); ymin = std::min(ymin, box.ymin); ymax = std::max(ymax, box.ymax); zmin = std::min(zmin, box.zmin); zmax = std::max(zmax, box.zmax);
#else
		xmin = min(xmin, box.xmin); xmax = max(xmax, box.xmax); ymin = min(ymin, box.ymin); ymax = max(ymax, box.ymax); zmin = min(zmin, box.zmin); zmax = max(zmax, box.zmax);
#endif
	}

    void Scale_About_Center(const T factor)
    {T x_center=(T).5*(xmin+xmax),y_center=(T).5*(ymin+ymax),z_center=(T).5*(zmin+zmax),x_length_over_two=factor*(T).5*(xmax-xmin),y_length_over_two=factor*(T).5*(ymax-ymin),
        z_length_over_two=factor*(T).5*(zmax-zmin);
    xmin=x_center-x_length_over_two;xmax=x_center+x_length_over_two;ymin=y_center-y_length_over_two;ymax=y_center+y_length_over_two;
    zmin=z_center-z_length_over_two;zmax=z_center+z_length_over_two;}

    void Scale_About_Center(const T x_factor,const T y_factor,const T z_factor)
    {T x_center=(T).5*(xmin+xmax),y_center=(T).5*(ymin+ymax),z_center=(T).5*(zmin+zmax),x_length_over_two=x_factor*(T).5*(xmax-xmin),y_length_over_two=y_factor*(T).5*(ymax-ymin),
        z_length_over_two=z_factor*(T).5*(zmax-zmin);
    xmin=x_center-x_length_over_two;xmax=x_center+x_length_over_two;ymin=y_center-y_length_over_two;ymax=y_center+y_length_over_two;
    zmin=z_center-z_length_over_two;zmax=z_center+z_length_over_two;}
	
    bool Inside(const T (&location)[3]) const
    {return location[0] >= xmin && location[0] <= xmax && location[1] >= ymin && location[1] <= ymax && location[2] >= zmin && location[2] <= zmax;}

    bool Outside(const T (&location)[3]) const
    {return location[0] < xmin || location[0] > xmax || location[1] < ymin || location[1] > ymax || location[2] < zmin || location[2] > zmax;}

    bool Intersection(const boundingBox<T>& box) const
    {return !(xmin > box.xmax || xmax < box.xmin || ymin > box.ymax || ymax < box.ymin || zmin > box.zmax || zmax < box.zmin);}

};
#endif
