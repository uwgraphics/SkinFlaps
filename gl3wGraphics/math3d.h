// file: math3d.h
// Author: Court Cutting
// Purpose: Bunch of useful 3D math routines

#ifndef __MATH3D_H__
#define __MATH3D_H__

#include <math.h>

inline bool m3dSolveLinearSystem3D(const float *m, const float *b, float (&result)[3])
{	// Input 3x3 matrix m, and 3 element vector b. Outputs result if true. If returns false determinant==0
	float cofactor11=m[4]*m[8]-m[7]*m[5],cofactor12=m[7]*m[2]-m[1]*m[8],cofactor13=m[1]*m[5]-m[4]*m[2];
	float determinant=m[0]*cofactor11+m[3]*cofactor12+m[6]*cofactor13;
	if(fabs(determinant)<1e-15f)
		return false;
	float n[9]= {cofactor11,cofactor12,cofactor13,m[6]*m[5]-m[3]*m[8],m[0]*m[8]-m[6]*m[2],m[3]*m[2]-m[0]*m[5],m[3]*m[7]-m[6]*m[4],m[6]*m[1]-m[0]*m[7],m[0]*m[4]-m[3]*m[1]};
	result[0] = b[0]*n[0]+b[1]*n[1]+b[2]*n[2];
	result[1] = b[0]*n[3]+b[1]*n[4]+b[2]*n[5];
	result[2] = b[0]*n[6]+b[1]*n[7]+b[2]*n[8];
	result[0]/=determinant; result[1]/=determinant; result[2]/=determinant;
	return true;
}

inline bool m3dRayTriangleIntersection(const float *rayOrigin, const float *rayDirection, const float *t0, const float *t1, const float *t2, float &rayParameter, float (&intersect)[3])
{
	float b[3],m[9],r[3],U[3],V[3];
	for(int i=0; i<3; ++i) {
		b[i]=rayOrigin[i]-t0[i];
		m[(i<<1)+i] = -rayDirection[i];
		U[i] = t1[i]-t0[i];
		m[(i<<1)+i+1] = U[i];
		V[i] = t2[i]-t0[i];
		m[(i<<1)+i+2] = V[i];
	}
	if(!m3dSolveLinearSystem3D(m,b,r))
		return false;
	if(r[1]<-1e-4f || r[2]<-1e-4f || r[1]+r[2]>1.0001f) // allows for roundoff error
		return false;
	rayParameter = r[0];
	for(int i=0; i<3; ++i)
		intersect[i] = t0[i] + r[1]*U[i] + r[2]*V[i];
	return true;
}


#endif // __MATH3D_H__
