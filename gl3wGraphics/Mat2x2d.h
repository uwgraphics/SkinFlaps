#ifndef __MAT2X2D__
#define __MAT2X2D__

#include <assert.h>
#include "Vec2d.h"

class Mat2x2d
{
public:
    double x[4];

    Mat2x2d()
    {
        for(int i=0;i<4;i++) x[i]=0;
    }

    Mat2x2d(const double x11,const double x21,const double x12,const double x22)
    {	x[0] = x11; x[1] = x21; x[2] = x12; x[3] = x22;    }

    Mat2x2d(const Vec2d& column1,const Vec2d& column2)
    {
		x[0]=column1.xy[0];x[1]=column1.xy[1];x[2]=column2.xy[0];x[3]=column2.xy[1];
    }

    void Initialize_With_Column_Vectors(const Vec2d& column1,const Vec2d& column2)
	{
		x[0] = column1.xy[0]; x[1] = column1.xy[1]; x[2] = column2.xy[0]; x[3] = column2.xy[1];
	}

    double& operator()(const int i,const int j)
    {assert(i > -1 && i < 2);assert(j > -1 && j < 2);return x[i+(j<<1)];}

    const double& operator()(const int i,const int j) const
	{
		assert(i > -1 && i < 2); assert(j > -1 && j < 2); return x[i + (j << 1)];
	}

    Mat2x2d operator-() const
    {return Mat2x2d(-x[0],-x[1],-x[2],-x[3]);}
    
    Mat2x2d& operator+=(const Mat2x2d& A)
    {for(int i=0;i<4;i++) x[i]+=A.x[i];return *this;}
        
    Mat2x2d& operator+=(const double& a)
    {x[0]+=a;x[3]+=a;return *this;}
        
    Mat2x2d& operator-=(const Mat2x2d& A)
    {for(int i=0;i<4;i++) x[i]-=A.x[i];return *this;}
    
    Mat2x2d& operator-=(const double& a)
    {x[0]-=a;x[3]-=a;return *this;}
        
    Mat2x2d& operator*=(const Mat2x2d& A)
    {return *this=*this*A;}
    
    Mat2x2d& operator*=(const double a)
    {for(int i=0;i<4;i++) x[i]*=a;return *this;}
    
    Mat2x2d& operator/=(const double a)
    {assert(a!=0);double s=1/a;for(int i=0;i<4;i++) x[i]*=s;return *this;}
    
    Mat2x2d operator+(const Mat2x2d& A) const
    {return Mat2x2d(x[0]+A.x[0],x[1]+A.x[1],x[2]+A.x[2],x[3]+A.x[3]);}

    Mat2x2d operator+(const double a) const
    {return Mat2x2d(x[0]+a,x[1],x[2],x[3]+a);}

    Mat2x2d operator-(const Mat2x2d& A) const
    {return Mat2x2d(x[0]-A.x[0],x[1]-A.x[1],x[2]-A.x[2],x[3]-A.x[3]);}

    Mat2x2d operator-(const double a) const
    {return Mat2x2d(x[0]-a,x[1],x[2],x[3]-a);}

	Mat2x2d operator*(const Mat2x2d& A) const
	{
		return Mat2x2d(x[0] * A.x[0] + x[2] * A.x[1], x[1] * A.x[0] + x[3] * A.x[1],
			x[0] * A.x[2] + x[2] * A.x[3], x[1] * A.x[2] + x[3] * A.x[3]);
	}
      
    Mat2x2d operator*(const double a) const
    {return Mat2x2d(a*x[0],a*x[1],a*x[2],a*x[3]);}
    
    Mat2x2d operator/(const double a) const
    {assert(a!=0);double s=1/a;return Mat2x2d(s*x[0],s*x[1],s*x[2],s*x[3]);}
    
    Vec2d operator*(const Vec2d& v) const
    {return Vec2d(x[0]*v.xy[0]+x[2]*v.xy[1],x[1]*v.xy[0]+x[3]*v.xy[1]);}

    double Determinant() const
    {return x[0]*x[3]-x[1]*x[2];}
    
    void Invert()
    {*this=Inverse();}

    Mat2x2d Inverse() const
	{
	double determinant = x[0] * x[3] - x[1] * x[2];
	assert(determinant != 0); double s = 1 / determinant;
    return Mat2x2d(x[3],-x[1],-x[2],x[0])*s;}

    Mat2x2d Inverse_Transpose() const
    {return Inverse().Transposed();}

    bool Solve_Linear_System(const Vec2d b, Vec2d &result) const
    {
	double determinant = x[0] * x[3] - x[1] * x[2];
	if (fabs(determinant)<1e-8f) 	return false;
	double s = 1 / determinant;
	result.xy[0] = (b.xy[0] * x[3] - b.xy[1] * x[2])*s;
	result.xy[1] = (x[0] * b.xy[1] - x[1] * b.xy[0])*s;
	return true; }

	Vec2d Robust_Solve_Linear_System(const Vec2d b, const double tolerance = 1e-8) const
	{
		double determinant = x[0] * x[3] - x[1] * x[2];
		if (fabs(determinant) < tolerance) determinant = determinant >= 0 ? tolerance : -tolerance;
		double s = 1 / determinant;
		return Vec2d((b.xy[0] * x[3] - b.xy[1] * x[2])*s, (x[0] * b.xy[1] - x[1] * b.xy[0])*s);
	}
    
    void Transpose()
    {double ex=x[1];x[1]=x[2];x[2]=ex;}

    Mat2x2d Transposed() const
    {return Mat2x2d(x[0],x[2],x[1],x[3]);}
        
    double Trace() const
    {return x[0]+x[3];}
    
    static Mat2x2d Transpose(const Mat2x2d& A)
    {return Mat2x2d(A.x[0],A.x[2],A.x[1],A.x[3]);}

    static Mat2x2d Identity_Matrix()
    {return Mat2x2d(1,0,0,1);}

    static Mat2x2d Diagonal_Matrix(const double d)
    {return Mat2x2d(d,0,0,d);}
    
    static Mat2x2d Rotation_Matrix(const double radians)
    {double c=cos(radians),s=sin(radians);return Mat2x2d(c,s,-s,c);}

    static Mat2x2d Outer_Product(const Vec2d& u,const Vec2d& v)
    {return Mat2x2d(u.xy[0]*v.xy[0],u.xy[1]*v.xy[0],u.xy[0]*v.xy[1],u.xy[1]*v.xy[1]);}

    static double Inner_Product(const Mat2x2d& A,const Mat2x2d& B)
    {return A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2]+A.x[3]*B.x[3];}

//#####################################################################
};

// global functions 
inline Mat2x2d operator+(const double a,const Mat2x2d& A)
{return A+a;}

inline Mat2x2d operator*(const double a,const Mat2x2d& A)
{return A*a;}

inline Vec2d operator*(const Vec2d& v,const Mat2x2d& A)
{return Vec2d(v.xy[0]*A.x[0]+v.xy[1]*A.x[1],v.xy[0]*A.x[2]+v.xy[1]*A.x[3]);}

//#####################################################################

#endif  // __Mat2x2d__
