#ifndef __MAT3X3D__
#define __MAT3X3D__

#include <assert.h>
#include "Vec3d.h"

class Mat3x3d
{
public:
    double x[9];

    Mat3x3d()
    {
        for(int i=0;i<9;i++) x[i]=0;
    }

    Mat3x3d(const double x11,const double x21,const double x31,const double x12,const double x22,const double x32,const double x13,const double x23,const double x33)
    {
        x[0]=x11;x[1]=x21;x[2]=x31;x[3]=x12;x[4]=x22;x[5]=x32;x[6]=x13;x[7]=x23;x[8]=x33;
    }

    Mat3x3d(const Vec3d& column1,const Vec3d& column2,const Vec3d& column3)
    {
		x[0]=column1.xyz[0];x[1]=column1.xyz[1];x[2]=column1.xyz[2];x[3]=column2.xyz[0];x[4]=column2.xyz[1];x[5]=column2.xyz[2];x[6]=column3.xyz[0];x[7]=column3.xyz[1];x[8]=column3.xyz[2];
    }

    void Initialize_With_Column_Vectors(const Vec3d& column1,const Vec3d& column2,const Vec3d& column3)
    {x[0]=column1.xyz[0];x[1]=column1.xyz[1];x[2]=column1.xyz[2];x[3]=column2.xyz[0];x[4]=column2.xyz[1];x[5]=column2.xyz[2];x[6]=column3.xyz[0];x[7]=column3.xyz[1];x[8]=column3.xyz[2];}

    double& operator()(const int i,const int j)
    {assert(i > -1 && i < 3);assert(j > -1 && j < 3);return x[i+(j<<1)+j];}

    const double& operator()(const int i,const int j) const
	{
		assert(i > -1 && i < 3); assert(j > -1 && j < 3); return x[i + (j << 1) + j];
	}

/*
    VECTOR_3D<T>& Column(const int j)
    {assert(1<=j && j<=3);return *(VECTOR_3D<T>*)(x+3*(j-1));}

    const VECTOR_3D<T>& Column(const int j) const
    {assert(1<=j && j<=3);return *(VECTOR_3D<T>*)(x+3*(j-1));} */

    Mat3x3d operator-() const
    {return Mat3x3d(-x[0],-x[1],-x[2],-x[3],-x[4],-x[5],-x[6],-x[7],-x[8]);}
    
    Mat3x3d& operator+=(const Mat3x3d& A)
    {for(int i=0;i<9;i++) x[i]+=A.x[i];return *this;}
        
    Mat3x3d& operator+=(const double& a)
    {x[0]+=a;x[4]+=a;x[8]+=a;return *this;}
        
    Mat3x3d& operator-=(const Mat3x3d& A)
    {for(int i=0;i<9;i++) x[i]-=A.x[i];return *this;}
    
    Mat3x3d& operator-=(const double& a)
    {x[0]-=a;x[4]-=a;x[8]-=a;return *this;}
        
    Mat3x3d& operator*=(const Mat3x3d& A)
    {return *this=*this*A;}
    
    Mat3x3d& operator*=(const double a)
    {for(int i=0;i<9;i++) x[i]*=a;return *this;}
    
    Mat3x3d& operator/=(const double a)
    {assert(a!=0);double s=1/a;for(int i=0;i<9;i++) x[i]*=s;return *this;}
    
    Mat3x3d operator+(const Mat3x3d& A) const
    {return Mat3x3d(x[0]+A.x[0],x[1]+A.x[1],x[2]+A.x[2],x[3]+A.x[3],x[4]+A.x[4],x[5]+A.x[5],x[6]+A.x[6],x[7]+A.x[7],x[8]+A.x[8]);}

    Mat3x3d operator+(const double a) const
    {return Mat3x3d(x[0]+a,x[1],x[2],x[3],x[4]+a,x[5],x[6],x[7],x[8]+a);}

    Mat3x3d operator-(const Mat3x3d& A) const
    {return Mat3x3d(x[0]-A.x[0],x[1]-A.x[1],x[2]-A.x[2],x[3]-A.x[3],x[4]-A.x[4],x[5]-A.x[5],x[6]-A.x[6],x[7]-A.x[7],x[8]-A.x[8]);}

    Mat3x3d operator-(const double a) const
    {return Mat3x3d(x[0]-a,x[1],x[2],x[3],x[4]-a,x[5],x[6],x[7],x[8]-a);}

    Mat3x3d operator*(const Mat3x3d& A) const // 27 mults, 18 adds
    {return Mat3x3d(x[0]*A.x[0]+x[3]*A.x[1]+x[6]*A.x[2],x[1]*A.x[0]+x[4]*A.x[1]+x[7]*A.x[2],x[2]*A.x[0]+x[5]*A.x[1]+x[8]*A.x[2],
        x[0]*A.x[3]+x[3]*A.x[4]+x[6]*A.x[5],x[1]*A.x[3]+x[4]*A.x[4]+x[7]*A.x[5],x[2]*A.x[3]+x[5]*A.x[4]+x[8]*A.x[5],
        x[0]*A.x[6]+x[3]*A.x[7]+x[6]*A.x[8],x[1]*A.x[6]+x[4]*A.x[7]+x[7]*A.x[8],x[2]*A.x[6]+x[5]*A.x[7]+x[8]*A.x[8]);}
      
    Mat3x3d operator*(const double a) const
    {return Mat3x3d(a*x[0],a*x[1],a*x[2],a*x[3],a*x[4],a*x[5],a*x[6],a*x[7],a*x[8]);}
    
    Mat3x3d operator/(const double a) const
    {assert(a!=0);double s=1/a;return Mat3x3d(s*x[0],s*x[1],s*x[2],s*x[3],s*x[4],s*x[5],s*x[6],s*x[7],s*x[8]);}
    
    Vec3d operator*(const Vec3d& v) const // 9 mults, 6 adds
    {return Vec3d(x[0]*v.xyz[0]+x[3]*v.xyz[1]+x[6]*v.xyz[2],x[1]*v.xyz[0]+x[4]*v.xyz[1]+x[7]*v.xyz[2],x[2]*v.xyz[0]+x[5]*v.xyz[1]+x[8]*v.xyz[2]);}

/*    VECTOR_2D<T> operator*(const VECTOR_2D<T>& v) const // assumes w=1 is the 3rd coordinate of v
    {T w=x[2]*v.x+x[5]*v.y+x[8];assert(w!=0);
    if(w==1) return VECTOR_2D<T>(x[0]*v.x+x[3]*v.y+x[6],x[1]*v.x+x[4]*v.y+x[7]);
    else{T s=1/w;return VECTOR_2D<T>(s*(x[0]*v.x+x[3]*v.y+x[6]),s*(x[1]*v.x+x[4]*v.y+x[7]));}} // rescale so w=1

    VECTOR_2D<T> Transform_2X2(const VECTOR_2D<T>& v) const // multiplies vector by upper 2x2 of matrix only
    {return VECTOR_2D<T>(x[0]*v.x+x[3]*v.y,x[1]*v.x+x[4]*v.y);} */

    double Determinant() const
    {return x[0]*(x[4]*x[8]-x[7]*x[5])+x[3]*(x[7]*x[2]-x[1]*x[8])+x[6]*(x[1]*x[5]-x[4]*x[2]);}
    
    void Invert()
    {*this=Inverse();}

    Mat3x3d Inverse() const
    {double cofactor11=x[4]*x[8]-x[7]*x[5],cofactor12=x[7]*x[2]-x[1]*x[8],cofactor13=x[1]*x[5]-x[4]*x[2];
    double determinant=x[0]*cofactor11+x[3]*cofactor12+x[6]*cofactor13;assert(determinant!=0);double s=1/determinant;
    return Mat3x3d(cofactor11,cofactor12,cofactor13,x[6]*x[5]-x[3]*x[8],x[0]*x[8]-x[6]*x[2],x[3]*x[2]-x[0]*x[5],x[3]*x[7]-x[6]*x[4],x[6]*x[1]-x[0]*x[7],x[0]*x[4]-x[3]*x[1])*s;}

    Mat3x3d Inverse_Transpose() const
    {return Inverse().Transposed();}

//    Mat3x3d Deviatoric() const
//    {double one_third_trace=(double)one_third*(x[0]+x[4]+x[8]);return Mat3x3d(x[0]-one_third_trace,x[1],x[2],x[3],x[4]-one_third_trace,x[5],x[6],x[7],x[8]-one_third_trace);}

//    Mat3x3d Dilational() const
//    {double one_third_trace=(double)one_third*(x[0]+x[4]+x[8]);return Mat3x3d(one_third_trace,0,0,0,one_third_trace,0,0,0,one_third_trace);}

    bool Solve_Linear_System(const Vec3d b, Vec3d &result) const // 33 mults, 17 adds, 1 div
    {double cofactor11=x[4]*x[8]-x[7]*x[5],cofactor12=x[7]*x[2]-x[1]*x[8],cofactor13=x[1]*x[5]-x[4]*x[2];
    double determinant=x[0]*cofactor11+x[3]*cofactor12+x[6]*cofactor13;
	if(fabs(determinant)<1e-16) 	return false;
    result = Mat3x3d(cofactor11,cofactor12,cofactor13,(x[6]*x[5]-x[3]*x[8]),(x[0]*x[8]-x[6]*x[2]),(x[3]*x[2]-x[0]*x[5]),(x[3]*x[7]-x[6]*x[4]),(x[6]*x[1]-x[0]*x[7]),(x[0]*x[4]-x[3]*x[1]))*b/determinant;
	return true; }

    Vec3d Robust_Solve_Linear_System(const Vec3d b,const double tolerance=1e-16) const // 33 mults, 17 adds, 1 div
    {double cofactor11=x[4]*x[8]-x[7]*x[5],cofactor12=x[7]*x[2]-x[1]*x[8],cofactor13=x[1]*x[5]-x[4]*x[2];
    double determinant=x[0]*cofactor11+x[3]*cofactor12+x[6]*cofactor13;if(fabs(determinant)<tolerance) determinant=determinant>=0?tolerance:-tolerance;
    return Mat3x3d(cofactor11,cofactor12,cofactor13,(x[6]*x[5]-x[3]*x[8]),(x[0]*x[8]-x[6]*x[2]),(x[3]*x[2]-x[0]*x[5]),(x[3]*x[7]-x[6]*x[4]),(x[6]*x[1]-x[0]*x[7]),(x[0]*x[4]-x[3]*x[1]))*b/determinant;}
    
    void Transpose()
    {double ex=x[1];x[1]=x[3];x[3]=ex; ex=x[2];x[2]=x[6];x[6]=ex; ex=x[5];x[5]=x[7];x[7]=ex;}

    Mat3x3d Transposed() const
    {return Mat3x3d(x[0],x[3],x[6],x[1],x[4],x[7],x[2],x[5],x[8]);}
        
    double Trace() const
    {return x[0]+x[4]+x[8];}
    
    Mat3x3d Q_From_Gram_Schmidt_QR_Factorization() const 
    {int k;Mat3x3d Q=*this;
    double one_over_r11=1/sqrt(Q.x[0]*Q.x[0]+Q.x[1]*Q.x[1]+Q.x[2]*Q.x[2]);for(k=0;k<=2;k++) Q.x[k]=one_over_r11*Q.x[k];
    double r12=Q.x[0]*Q.x[3]+Q.x[1]*Q.x[4]+Q.x[2]*Q.x[5];Q.x[3]-=r12*Q.x[0];Q.x[4]-=r12*Q.x[1];Q.x[5]-=r12*Q.x[2];
    double r13=Q.x[0]*Q.x[6]+Q.x[1]*Q.x[7]+Q.x[2]*Q.x[8];Q.x[6]-=r13*Q.x[0];Q.x[7]-=r13*Q.x[1];Q.x[8]-=r13*Q.x[2];
    double one_over_r22=1/sqrt(Q.x[3]*Q.x[3]+Q.x[4]*Q.x[4]+Q.x[5]*Q.x[5]);for(k=3;k<=5;k++) Q.x[k]=one_over_r22*Q.x[k];
    double r23=Q.x[3]*Q.x[6]+Q.x[4]*Q.x[7]+Q.x[5]*Q.x[8];Q.x[6]-=r23*Q.x[3];Q.x[7]-=r23*Q.x[4];Q.x[8]-=r23*Q.x[5];
    double one_over_r33=1/sqrt(Q.x[6]*Q.x[6]+Q.x[7]*Q.x[7]+Q.x[8]*Q.x[8]);for(k=6;k<=8;k++) Q.x[k]=one_over_r33*Q.x[k];
    return Q;}

/*    UPPER_TRIANGULAR_Mat3x3d R_From_Gram_Schmidt_QR_Factorization() const 
    {int k;Mat3x3d Q=*this;UPPER_TRIANGULAR_Mat3x3d R;
    R.x11=sqrt((sqr(Q.x[0])+sqr(Q.x[1])+sqr(Q.x[2])));double one_over_r11=1/R.x11;for(k=0;k<=2;k++) Q.x[k]=one_over_r11*Q.x[k];
    R.x12=Q.x[0]*Q.x[3]+Q.x[1]*Q.x[4]+Q.x[2]*Q.x[5];Q.x[3]-=R.x12*Q.x[0];Q.x[4]-=R.x12*Q.x[1];Q.x[5]-=R.x12*Q.x[2];
    R.x13=Q.x[0]*Q.x[6]+Q.x[1]*Q.x[7]+Q.x[2]*Q.x[8];Q.x[6]-=R.x13*Q.x[0];Q.x[7]-=R.x13*Q.x[1];Q.x[8]-=R.x13*Q.x[2];
    R.x22=sqrt((sqr(Q.x[3])+sqr(Q.x[4])+sqr(Q.x[5])));double one_over_r22=1/R.x22;for(k=3;k<=5;k++) Q.x[k]=one_over_r22*Q.x[k];
    R.x23=Q.x[3]*Q.x[6]+Q.x[4]*Q.x[7]+Q.x[5]*Q.x[8];Q.x[6]-=R.x23*Q.x[3];Q.x[7]-=R.x23*Q.x[4];Q.x[8]-=R.x23*Q.x[5];
    R.x33=sqrt((sqr(Q.x[6])+sqr(Q.x[7])+sqr(Q.x[8])));
    return R;} */
    
//    Mat3x3d Higham_Iterate(const double tolerance=1e-5,const int max_iterations=20,const bool exit_on_max_iterations=false) const
//    {Mat3x3d X=*this,Y;int iterations=0;
//    for(;;){
//        Y=(double).5*(X+X.Inverse_Transpose());if((X-Y).Max_Abs_Element()<tolerance) return Y;X=Y;
//        if(++iterations >= max_iterations){if(exit_on_max_iterations){std::cerr<<"Failed in Higham iteration!"<<std::endl;exit(1);}else return X;}}}
    
/*    Mat3x3d<T> Cofactor_Matrix() const
    {return Mat3x3d<T>(x[4]*x[8]-x[5]*x[7],-x[3]*x[8]+x[5]*x[6],x[3]*x[7]-x[4]*x[6],
                                          -x[1]*x[8]+x[2]*x[7],x[0]*x[8]-x[2]*x[6],-x[0]*x[7]+x[1]*x[6],
                                          x[1]*x[5]-x[2]*x[4],-x[0]*x[5]+x[2]*x[3],x[0]*x[4]-x[1]*x[3]);}

    SYMMETRIC_Mat3x3d<T> Outer_Product_Matrix() const
    {return SYMMETRIC_Mat3x3d<T>(x[0]*x[0]+x[3]*x[3]+x[6]*x[6],x[1]*x[0]+x[4]*x[3]+x[7]*x[6],x[2]*x[0]+x[5]*x[3]+x[8]*x[6],
                                    x[1]*x[1]+x[4]*x[4]+x[7]*x[7],x[2]*x[1]+x[5]*x[4]+x[8]*x[7],x[2]*x[2]+x[5]*x[5]+x[8]*x[8]);}

    SYMMETRIC_Mat3x3d<T> Normal_Equations_Matrix() const // 18 mults, 12 adds
    {return SYMMETRIC_Mat3x3d<T>(x[0]*x[0]+x[1]*x[1]+x[2]*x[2],x[3]*x[0]+x[4]*x[1]+x[5]*x[2],x[6]*x[0]+x[7]*x[1]+x[8]*x[2],
                                    x[3]*x[3]+x[4]*x[4]+x[5]*x[5],x[6]*x[3]+x[7]*x[4]+x[8]*x[5],x[6]*x[6]+x[7]*x[7]+x[8]*x[8]);}

    SYMMETRIC_Mat3x3d<T> Symmetric_Part() const
    {return SYMMETRIC_Mat3x3d<T>(x[0],(T).5*(x[1]+x[3]),(T).5*(x[2]+x[6]),x[4],(T).5*(x[5]+x[7]),x[8]);}
    
    SYMMETRIC_Mat3x3d<T> Twice_Symmetric_Part() const // 3 mults, 3 adds
    {return SYMMETRIC_Mat3x3d<T>(2*x[0],x[1]+x[3],x[2]+x[6],2*x[4],x[5]+x[7],2*x[8]);}

    DIAGONAL_Mat3x3d<T> Diagonal_Part() const
    {return DIAGONAL_Mat3x3d<T>(x[0],x[4],x[8]);}

    void Normalize_Columns()
    {T magnitude=sqrt(sqr(x[0])+sqr(x[1])+sqr(x[2]));assert(magnitude != 0);T s=1/magnitude;x[0]*=s;x[1]*=s;x[2]*=s;
    magnitude=sqrt(sqr(x[3])+sqr(x[4])+sqr(x[5]));assert(magnitude != 0);s=1/magnitude;x[3]*=s;x[4]*=s;x[5]*=s;
    magnitude=sqrt(sqr(x[6])+sqr(x[7])+sqr(x[8]));assert(magnitude != 0);s=1/magnitude;x[6]*=s;x[7]*=s;x[8]*=s;}
    
    VECTOR_3D<T> Largest_Normalized_Column() const
    {T scale1=sqr(x[0])+sqr(x[1])+sqr(x[2]),scale2=sqr(x[3])+sqr(x[4])+sqr(x[5]),scale3=sqr(x[6])+sqr(x[7])+sqr(x[8]);
    if(scale1>scale2){if(scale1>scale3)return VECTOR_3D<T>(x[0],x[1],x[2])/sqrt(scale1);}
    else if(scale2>scale3)return VECTOR_3D<T>(x[3],x[4],x[5])/sqrt(scale2);
    return VECTOR_3D<T>(x[6],x[7],x[8])/sqrt(scale3);} */

    static Mat3x3d Transpose(const Mat3x3d& A)
    {return Mat3x3d(A.x[0],A.x[3],A.x[6],A.x[1],A.x[4],A.x[7],A.x[2],A.x[5],A.x[8]);}

/*    static Mat3x3d<T> Translation_Matrix(const VECTOR_2D<T>& translation) // treating the 3x3 matrix as a homogeneous transformation on 2d vectors
    {return Mat3x3d<T>(1,0,0,0,1,0,translation.x,translation.y,1);} */

    static Mat3x3d Identity_Matrix()
    {return Mat3x3d(1,0,0,0,1,0,0,0,1);}

    static Mat3x3d Diagonal_Matrix(const double d)
    {return Mat3x3d(d,0,0,0,d,0,0,0,d);}
    
    static Mat3x3d Rotation_Matrix_X_Axis(const double radians)
    {double c=cos(radians),s=sin(radians);return Mat3x3d(1,0,0,0,c,s,0,-s,c);}

    static Mat3x3d Rotation_Matrix_Y_Axis(const double radians)
    {double c=cos(radians),s=sin(radians);return Mat3x3d(c,0,-s,0,1,0,s,0,c);}

    static Mat3x3d Rotation_Matrix_Z_Axis(const double radians)
    {double c=cos(radians),s=sin(radians);return Mat3x3d(c,s,0,-s,c,0,0,0,1);}

    static Mat3x3d Rotation_Matrix(const Vec3d& axis,const double radians)
    {double c=cos(radians),s=sin(radians);
	return Mat3x3d(axis.xyz[0]*axis.xyz[0]+(1-axis.xyz[0]*axis.xyz[0])*c,axis.xyz[0]*axis.xyz[1]*(1-c)+axis.xyz[2]*s,axis.xyz[0]*axis.xyz[2]*(1-c)-axis.xyz[1]*s,
                                         axis.xyz[0]*axis.xyz[1]*(1-c)-axis.xyz[2]*s,axis.xyz[1]*axis.xyz[1]+(1-axis.xyz[1]*axis.xyz[1])*c,axis.xyz[1]*axis.xyz[2]*(1-c)+axis.xyz[0]*s,
                                         axis.xyz[0]*axis.xyz[2]*(1-c)+axis.xyz[1]*s,axis.xyz[1]*axis.xyz[2]*(1-c)-axis.xyz[0]*s,axis.xyz[2]*axis.xyz[2]+(1-axis.xyz[2]*axis.xyz[2])*c);}

/*    static Mat3x3d<T> Rotation_Matrix(const VECTOR_3D<T>& rotation)
    {T angle=rotation.Magnitude();return angle?Rotation_Matrix(rotation/angle,angle):Identity_Matrix();}

    static Mat3x3d<T> Rotation_Matrix(const VECTOR_3D<T>& x_final,const VECTOR_3D<T>& y_final,const VECTOR_3D<T>& z_final)
    {return Mat3x3d<T>(x_final.x,y_final.x,z_final.x,x_final.y,y_final.y,z_final.y,x_final.z,y_final.z,z_final.z);}
    
    static Mat3x3d<T> Rotation_Matrix(const VECTOR_3D<T>& initial_vector,const VECTOR_3D<T>& final_vector)
    {VECTOR_3D<T> initial_unit=initial_vector/initial_vector.Magnitude(),final_unit=final_vector/final_vector.Magnitude();
    T cos_theta=VECTOR_3D<T>::Dot_Product(initial_unit,final_unit);
    if(cos_theta > 1-1e-14) return Identity_Matrix();
    if(cos_theta < -1+1e-14) return Mat3x3d<T>(-1,0,0,0,-1,0,0,0,-1); // note-this is actually a reflection
    VECTOR_3D<T> axis=VECTOR_3D<T>::Cross_Product(initial_unit,final_unit);axis.Normalize();
    return Rotation_Matrix(axis,acos(cos_theta));} */

    static Mat3x3d Outer_Product(const Vec3d& u,const Vec3d& v)
    {return Mat3x3d(u.xyz[0]*v.xyz[0],u.xyz[1]*v.xyz[0],u.xyz[2]*v.xyz[0],u.xyz[0]*v.xyz[1],u.xyz[1]*v.xyz[1],u.xyz[2]*v.xyz[1],u.xyz[0]*v.xyz[2],u.xyz[1]*v.xyz[2],u.xyz[2]*v.xyz[2]);}

    static double Inner_Product(const Mat3x3d& A,const Mat3x3d& B)
    {return A.x[0]*B.x[0]+A.x[1]*B.x[1]+A.x[2]*B.x[2]+A.x[3]*B.x[3]+A.x[4]*B.x[4]+A.x[5]*B.x[5]+A.x[6]*B.x[6]+A.x[7]*B.x[7]+A.x[8]*B.x[8];}

/*    static Mat3x3d<T> Cross_Product_Matrix(const VECTOR_3D<T>& v)
    {return Mat3x3d<T>(0,v.z,-v.y,-v.z,0,v.x,v.y,-v.x,0);}
    
    VECTOR_3D<T> Antisymmetric_Part_Cross_Product_Vector() const
    {return (T).5*VECTOR_3D<T>(x[5]-x[7],x[6]-x[2],x[1]-x[3]);}

    static SYMMETRIC_Mat3x3d<T> Right_Multiply_With_Symmetric_Result(const Mat3x3d<T>& A, const DIAGONAL_Mat3x3d<T>& B)
    {return SYMMETRIC_Mat3x3d<T>(B.x11*A.x[0],B.x11*A.x[1],B.x11*A.x[2],B.x22*A.x[4],B.x22*A.x[5],B.x33*A.x[8]);}

    T Max_Abs_Element() const
    {return max(fabs(x[0]),fabs(x[1]),fabs(x[2]),fabs(x[3]),fabs(x[4]),fabs(x[5]),fabs(x[6]),fabs(x[7]),fabs(x[8]));}
    
    T Frobenius_Norm() const
    {return sqrt(sqr(x[0])+sqr(x[1])+sqr(x[2])+sqr(x[3])+sqr(x[4])+sqr(x[5])+sqr(x[6])+sqr(x[7])+sqr(x[8]));}
    
    Mat3x3d<T> operator*(const DIAGONAL_Mat3x3d<T>& A) const // 9 mults
    {return Mat3x3d<T>(x[0]*A.x11,x[1]*A.x11,x[2]*A.x11,x[3]*A.x22,x[4]*A.x22,x[5]*A.x22,x[6]*A.x33,x[7]*A.x33,x[8]*A.x33);}
    
    Mat3x3d<T> operator*(const UPPER_TRIANGULAR_Mat3x3d<T>& A) const // 18 mults, 8 adds
    {return Mat3x3d<T>(x[0]*A.x11,x[1]*A.x11,x[2]*A.x11,x[0]*A.x12+x[3]*A.x22,x[1]*A.x12+x[4]*A.x22,x[2]*A.x12+x[5]*A.x22,
                          x[0]*A.x13+x[3]*A.x23+x[6]*A.x33,x[1]*A.x13+x[4]*A.x23+x[7]*A.x33,x[2]*A.x13+x[5]*A.x23+x[8]*A.x33);}
    
    Mat3x3d<T> operator*(const SYMMETRIC_Mat3x3d<T>& A) const // 27 mults, 18 adds
    {return Mat3x3d<T>(x[0]*A.x11+x[3]*A.x21+x[6]*A.x31,x[1]*A.x11+x[4]*A.x21+x[7]*A.x31,x[2]*A.x11+x[5]*A.x21+x[8]*A.x31,
                          x[0]*A.x21+x[3]*A.x22+x[6]*A.x32,x[1]*A.x21+x[4]*A.x22+x[7]*A.x32,x[2]*A.x21+x[5]*A.x22+x[8]*A.x32,
                          x[0]*A.x31+x[3]*A.x32+x[6]*A.x33,x[1]*A.x31+x[4]*A.x32+x[7]*A.x33,x[2]*A.x31+x[5]*A.x32+x[8]*A.x33);}
      
    Mat3x3d<T> Multiply_With_Transpose(const UPPER_TRIANGULAR_Mat3x3d<T>& A) const
    {return Mat3x3d<T>(x[0]*A.x11+x[3]*A.x12+x[6]*A.x13,x[1]*A.x11+x[4]*A.x12+x[7]*A.x13,x[2]*A.x11+x[5]*A.x12+x[8]*A.x13,
                          x[3]*A.x22+x[6]*A.x23,x[4]*A.x22+x[7]*A.x23,x[5]*A.x22+x[8]*A.x23,x[6]*A.x33,x[7]*A.x33,x[8]*A.x33);}
      
    Mat3x3d<T> Multiply_With_Transpose(const Mat3x3d<T>& A) const
    {return Mat3x3d<T>(x[0]*A.x[0]+x[3]*A.x[3]+x[6]*A.x[6],x[1]*A.x[0]+x[4]*A.x[3]+x[7]*A.x[6],x[2]*A.x[0]+x[5]*A.x[3]+x[8]*A.x[6],
                          x[0]*A.x[1]+x[3]*A.x[4]+x[6]*A.x[7],x[1]*A.x[1]+x[4]*A.x[4]+x[7]*A.x[7],x[2]*A.x[1]+x[5]*A.x[4]+x[8]*A.x[7],
                          x[0]*A.x[2]+x[3]*A.x[5]+x[6]*A.x[8],x[1]*A.x[2]+x[4]*A.x[5]+x[7]*A.x[8],x[2]*A.x[2]+x[5]*A.x[5]+x[8]*A.x[8]);}
      
    Mat3x3d<T> Transpose_Times(const Mat3x3d<T>& A) const
    {return Mat3x3d<T>(x[0]*A.x[0]+x[1]*A.x[1]+x[2]*A.x[2],x[3]*A.x[0]+x[4]*A.x[1]+x[5]*A.x[2],x[6]*A.x[0]+x[7]*A.x[1]+x[8]*A.x[2],
                          x[0]*A.x[3]+x[1]*A.x[4]+x[2]*A.x[5],x[3]*A.x[3]+x[4]*A.x[4]+x[5]*A.x[5],x[6]*A.x[3]+x[7]*A.x[4]+x[8]*A.x[5],
                          x[0]*A.x[6]+x[1]*A.x[7]+x[2]*A.x[8],x[3]*A.x[6]+x[4]*A.x[7]+x[5]*A.x[8],x[6]*A.x[6]+x[7]*A.x[7]+x[8]*A.x[8]);}

    VECTOR_3D<T> Transpose_Times(const VECTOR_3D<T>& v) const
    {return VECTOR_3D<T>(x[0]*v.x+x[1]*v.y+x[2]*v.z,x[3]*v.x+x[4]*v.y+x[5]*v.z,x[6]*v.x+x[7]*v.y+x[8]*v.z);}

    static Mat3x3d<T> Left_Procrustes_Rotation(const Mat3x3d<T>& A,const Mat3x3d<T>& B)
    {Mat3x3d<T> U,V;DIAGONAL_Mat3x3d<T> D;A.Multiply_With_Transpose(B).Fast_Singular_Value_Decomposition(U,D,V);return U.Multiply_With_Transpose(V);} */

//    void Fast_Singular_Value_Decomposition(Mat3x3d<T>& U,DIAGONAL_Mat3x3d<T>& singular_values,Mat3x3d<T>& V) const
//    {Mat3x3d<double> U_double,V_double;DIAGONAL_Mat3x3d<double> singular_values_double;
//    Fast_Singular_Value_Decomposition_Double(U_double,singular_values_double,V_double);U=U_double;singular_values=singular_values_double;V=V_double;}

//    void Fast_Singular_Value_Decomposition_Double(Mat3x3d<double>& U,DIAGONAL_Mat3x3d<double>& singular_values,Mat3x3d<double>& V,const double tolerance=1e-7) const;
//    T Tetrahedron_Minimum_Altitude() const;
//#####################################################################
};      
// global functions 
inline Mat3x3d operator+(const double a,const Mat3x3d& A)
{return A+a;}

inline Mat3x3d operator*(const double a,const Mat3x3d& A)
{return A*a;}

inline Vec3d operator*(const Vec3d& v,const Mat3x3d& A)
{return Vec3d(v.xyz[0]*A.x[0]+v.xyz[1]*A.x[1]+v.xyz[2]*A.x[2],v.xyz[0]*A.x[3]+v.xyz[1]*A.x[4]+v.xyz[2]*A.x[5],v.xyz[0]*A.x[6]+v.xyz[1]*A.x[7]+v.xyz[2]*A.x[8]);}

//#####################################################################

#endif  // __MAT3X3D__
