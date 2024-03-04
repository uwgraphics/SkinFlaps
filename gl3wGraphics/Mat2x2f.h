#ifndef __MAT2X2F__
#define __MAT2X2F__

#include <assert.h>
#include "Vec2f.h"

class Mat2x2f
{
public:
    float x[4];

    Mat2x2f()
    {
        for(int i=0;i<4;i++) x[i]=0;
    }

    Mat2x2f(const float x11,const float x21,const float x12,const float x22)
    {	x[0] = x11; x[1] = x21; x[2] = x12; x[3] = x22;    }

    Mat2x2f(const Vec2f& column1,const Vec2f& column2)
    {
		x[0]=column1.xy[0];x[1]=column1.xy[1];x[2]=column2.xy[0];x[3]=column2.xy[1];
    }

    void Initialize_With_Column_Vectors(const Vec2f& column1,const Vec2f& column2)
	{
		x[0] = column1.xy[0]; x[1] = column1.xy[1]; x[2] = column2.xy[0]; x[3] = column2.xy[1];
	}

    float& operator()(const int i,const int j)
    {assert(i > -1 && i < 2);assert(j > -1 && j < 2);return x[i+(j<<1)];}

    const float& operator()(const int i,const int j) const
	{
		assert(i > -1 && i < 2); assert(j > -1 && j < 2); return x[i + (j << 1)];
	}

    Mat2x2f operator-() const
    {return Mat2x2f(-x[0],-x[1],-x[2],-x[3]);}
    
    Mat2x2f& operator+=(const Mat2x2f& A)
    {for(int i=0;i<4;i++) x[i]+=A.x[i];return *this;}
        
    Mat2x2f& operator+=(const float& a)
    {x[0]+=a;x[3]+=a;return *this;}
        
    Mat2x2f& operator-=(const Mat2x2f& A)
    {for(int i=0;i<4;i++) x[i]-=A.x[i];return *this;}
    
    Mat2x2f& operator-=(const float& a)
    {x[0]-=a;x[3]-=a;return *this;}
        
    Mat2x2f& operator*=(const Mat2x2f& A)
    {return *this=*this*A;}
    
    Mat2x2f& operator*=(const float a)
    {for(int i=0;i<4;i++) x[i]*=a;return *this;}
    
    Mat2x2f& operator/=(const float a)
    {assert(a!=0);float s=1/a;for(int i=0;i<4;i++) x[i]*=s;return *this;}
    
    Mat2x2f operator+(const Mat2x2f& A) const
    {return Mat2x2f(x[0]+A.x[0],x[1]+A.x[1],x[2]+A.x[2],x[3]+A.x[3]);}

    Mat2x2f operator+(const float a) const
    {return Mat2x2f(x[0]+a,x[1],x[2],x[3]+a);}

    Mat2x2f operator-(const Mat2x2f& A) const
    {return Mat2x2f(x[0]-A.x[0],x[1]-A.x[1],x[2]-A.x[2],x[3]-A.x[3]);}

    Mat2x2f operator-(const float a) const
    {return Mat2x2f(x[0]-a,x[1],x[2],x[3]-a);}

	Mat2x2f operator*(const Mat2x2f& A) const
	{
		return Mat2x2f(x[0] * A.x[0] + x[2] * A.x[1], x[1] * A.x[0] + x[3] * A.x[1],
			x[0] * A.x[2] + x[2] * A.x[3], x[1] * A.x[2] + x[3] * A.x[3]);
	}
      
    Mat2x2f operator*(const float a) const
    {return Mat2x2f(a*x[0],a*x[1],a*x[2],a*x[3]);}
    
    Mat2x2f operator/(const float a) const
    {assert(a!=0);float s=1/a;return Mat2x2f(s*x[0],s*x[1],s*x[2],s*x[3]);}
    
    Vec2f operator*(const Vec2f& v) const
    {return Vec2f(x[0]*v.xy[0]+x[2]*v.xy[1],x[1]*v.xy[0]+x[3]*v.xy[1]);}

    float Determinant() const
    {return x[0]*x[3]-x[1]*x[2];}
    
    void Invert()
    {*this=Inverse();}

    Mat2x2f Inverse() const
	{
	float determinant = x[0] * x[3] - x[1] * x[2];
	assert(determinant != 0); float s = 1 / determinant;
    return Mat2x2f(x[3],-x[1],-x[2],x[0])*s;}

    Mat2x2f Inverse_Transpose() const
    {return Inverse().Transposed();}

    bool Solve_Linear_System(const Vec2f b, Vec2f &result) const
    {
	float determinant = x[0] * x[3] - x[1] * x[2];
	if (fabs(determinant)<1e-8f) 	return false;
	float s = 1 / determinant;
	result.xy[0] = (b.xy[0] * x[3] - b.xy[1] * x[2])*s;
	result.xy[1] = (x[0] * b.xy[1] - x[1] * b.xy[0])*s;
	return true; }

	Vec2f Robust_Solve_Linear_System(const Vec2f b, const float tolerance = 1e-8) const
	{
		float determinant = x[0] * x[3] - x[1] * x[2];
		if (fabs(determinant) < tolerance) determinant = determinant >= 0 ? tolerance : -tolerance;
		float s = 1 / determinant;
		return Vec2f((b.xy[0] * x[3] - b.xy[1] * x[2])*s, (x[0] * b.xy[1] - x[1] * b.xy[0])*s);
	}
    
    void Transpose()
    {float ex=x[1];x[1]=x[2];x[2]=ex;}

    Mat2x2f Transposed() const
    {return Mat2x2f(x[0],x[2],x[1],x[3]);}
        
    float Trace() const
    {return x[0]+x[3];}
    
    static Mat2x2f Transpose(const Mat2x2f& A)
    {return Mat2x2f(A.x[0],A.x[2],A.x[1],A.x[3]);}

    static Mat2x2f Identity_Matrix()
    {return Mat2x2f(1,0,0,1);}

    static Mat2x2f Diagonal_Matrix(const float d)
    {return Mat2x2f(d,0,0,d);}
    
    static Mat2x2f Rotation_Matrix(const float radians)
    {float c=cos(radians),s=sin(radians);return Mat2x2f(c,s,-s,c);}

    static Mat2x2f Outer_Product(const Vec2f& u,const Vec2f& v)
    {return Mat2x2f(u.xy[0]*v.xy[0],u.xy[1]*v.xy[0],u.xy[0]*v.xy[1],u.xy[1]*v.xy[1]);}

    static float Inner_Product(const Mat2x2f& A,const Mat2x2f& B)
    {return A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2]+A.x[3]*B.x[3];}

//#####################################################################
};

// global functions 
inline Mat2x2f operator+(const float a,const Mat2x2f& A)
{return A+a;}

inline Mat2x2f operator*(const float a,const Mat2x2f& A)
{return A*a;}

inline Vec2f operator*(const Vec2f& v,const Mat2x2f& A)
{return Vec2f(v.xy[0]*A.x[0]+v.xy[1]*A.x[1],v.xy[0]*A.x[2]+v.xy[1]*A.x[3]);}

//#####################################################################

#endif  // __Mat2x2f__
