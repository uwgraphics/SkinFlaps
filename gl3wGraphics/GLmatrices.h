// File: GLmatrices.h
// Author: Court Cutting
// Date: February 3, 2012
// Purpose: Rudimentary matrix handling for openGL purposes.
//	Much of the code is modified from Richard Wright's excellent
//	OpenGL Superbible fifth edition who should get much of the credit for this.

#ifndef __GLMATRICES_H__
#define __GLMATRICES_H__

#include <math.h>
#include <GL/gl3w.h>

class GLmatrices
{
public:
	void getDragVector(float dScreenX, float dScreenY, float (&position)[3], float (&dragVector)[3]);
	void getViewVector(float (&view)[3]) {view[0]=_mFR[2]; view[1]=_mFR[6]; view[2]=_mFR[10];}		// unscaled
	const GLfloat* getProjectionMatrix() {return &_mProj[0];}
	const GLfloat* getFrameAndRotationMatrix() {return &_mFR[0];}
	void setFrameAndRotation(GLfloat *rotMatrix);
	void setView(float angleRadians, float screenAspect);
	void resetPerspective();	// must be called whenever scene changes
	float getSceneRadius() {return _radius;}
	void getSceneCenter(float (&center)[3]) {center[0]=_center[0]; center[1]=_center[1]; center[2]=_center[2];}
	void changeZoom(float dZoom);
	void setPan(float frX, float frY);
	void setSceneRadius(float radius) {_radius=radius;}
	void setSceneCenter(float (&center)[3]) {_center[0]=center[0]; _center[1]=center[1]; _center[2]=center[2];}
	void getCameraData(float &zCenter, float &sceneRadius, float &verticalAngle, float &screenAspect)
		{
		zCenter=_zCenter; sceneRadius=_radius; verticalAngle=_angleRadians; screenAspect=_screenAspect;
	}
	GLmatrices();
	~GLmatrices();

private:
	void MakePerspectiveMatrix(float fFov, float fAspect, float zMin, float zMax);
	GLfloat _mProj[16],_mFR[16],_rotNow[16];
	GLfloat _center[3];
	GLfloat _zCenter,_radius;
	float _zmin,_angleRadians,_screenAspect;
};

// public functions
	inline void loadIdentity4x4(GLfloat *mat)
	{
		mat[0]=1.0f; mat[1]=0.0f; mat[2]=0.0f; mat[3]=0.0f; 
		mat[4]=0.0f; mat[5]=1.0f; mat[6]=0.0f; mat[7]=0.0f; 
		mat[8]=0.0f; mat[9]=0.0f; mat[10]=1.0f; mat[11]=0.0f; 
		mat[12]=0.0f; mat[13]=0.0f; mat[14]=0.0f; mat[15]=1.0f; 
	}

	inline void translateMatrix4x4(GLfloat *m, float x, float y, float z) {m[12]+=x; m[13]+=y; m[14]+=z;} // assume no perspective. Last row 0,0,0,1

	inline void scaleMatrix4x4(GLfloat *m, float x, float y, float z) {m[0]*=x; m[4]*=x; m[8]*=x; m[1]*=y; m[5]*=y; m[9]*=y; m[2]*=z; m[6]*=z; m[10]*=z;} // assume no perspective. Last row 0,0,0,1

	inline void rotateMatrix4x4(GLfloat *m, char axis, float r) { // rotates matrix m r radians about axis. Assume no perspective. Last row 0,0,0,1
		int c1=0,c2=1;
		float s=(float)sin(r),c= (float)cos(r);
		if(axis=='x' || axis=='X') {c1=1; c2=2;}
		else if(axis=='y' || axis=='Y') {c1=0; c2=2; s=-s;}
		else ;
		for(int i=0; i<4; ++i) {
			float u=m[(i<<2)+c1]*c-m[(i<<2)+c2]*s,v=m[(i<<2)+c1]*s+m[(i<<2)+c2]*c;
			m[(i<<2)+c1]=u;	m[(i<<2)+c2]=v;
		}
	}

	inline void axisAngleRotateMatrix4x4(GLfloat *m, float (&axis)[3], float angle) { // rotates matrix m angle radians about axis. Assume no perspective. Last row 0,0,0,1
		float x=axis[0],y=axis[1],z=axis[2],w;
		w = x*x + y*y + z*z;
		if(w<1e-10f)	return;
		if(fabs(1.0f-w)>0.00001f) {
			w = (float)sqrt(w);
			x/=w; y/=w; z/=w;
		}
		w=angle*0.5f;
		float sA = (float)sin(w);
		x*=sA; y*=sA; z*=sA;
		w = (float)cos(w);	// have quat
		float x2=x*x,y2=y*y,z2=z*z,xy=x*y,xz=x*z,yz=y*z,wx=w*x,wy=w*y,wz=w*z;
		float r[][3]={1.0f-2.0f*(y2+z2), 2.0f*(xy-wz), 2.0f*(xz+wy),
					2.0f*(xy+wz), 1.0f-2.0f*(x2+z2), 2.0f*(yz-wx),
					2.0f*(xz-wy), 2.0f*(yz+wx), 1.0f-2.0f*(x2+y2)};
		for(int i=0; i<4; ++i)	{
			x=*(m+(i<<2)); y=*(m+(i<<2)+1); z=*(m+(i<<2)+2);
			*(m+(i<<2))=x*r[0][0]+y*r[0][1]+z*r[0][2];
			*(m+(i<<2)+1)=x*r[1][0]+y*r[1][1]+z*r[1][2];
			*(m+(i<<2)+2)=x*r[2][0]+y*r[2][1]+z*r[2][2];
		}
	}

	static float glMatricesDetIJ(const GLfloat *m, const int i, const int j)
	{ // 3x3 determinant
		int x, y, ii, jj;
		float ret, mat[3][3];
		x=0;
		for(ii=0; ii<4; ii++)
		{
			if (ii == i) continue;
			y = 0;
			for(jj=0; jj<4; jj++)
			{
				if (jj == j) continue;
				mat[x][y] = m[(ii*4)+jj];
				y++;
			}
			x++;
		}
		ret =  mat[0][0]*(mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2]);
		ret -= mat[0][1]*(mat[1][0]*mat[2][2]-mat[2][0]*mat[1][2]);
		ret += mat[0][2]*(mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1]);
		return ret;
	}

	inline void invertMatrix4x4(const GLfloat *m, GLfloat *mInverse)
	{
		int i, j;
		float det, detij;
		// calculate 4x4 determinant
		det = 0.0f;
		for (i=0; i<4; i++)
		{
			det += (i & 0x1) ? (-m[i] * glMatricesDetIJ(m, 0, i)) : (m[i] * glMatricesDetIJ(m, 0,i));
	}
		det = 1.0f / det;
		// calculate inverse
		for(i=0; i<4; i++)
		{
			for(j=0; j<4; j++)
			{
				detij = glMatricesDetIJ(m, j, i);
				mInverse[(i*4)+j] = ((i+j) & 0x1) ? (-detij * det) : (detij *det); 
			}
		}
	}

	inline void transformVector4(const float (&v)[4], const GLfloat *m, float (&vOut)[4])
	{
		vOut[0] = m[0] * v[0] + m[4] * v[1] + m[8] *  v[2] + m[12] * v[3]; 
		vOut[1] = m[1] * v[0] + m[5] * v[1] + m[9] *  v[2] + m[13] * v[3];	
		vOut[2] = m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14] * v[3];	
		vOut[3] = m[3] * v[0] + m[7] * v[1] + m[11] * v[2] + m[15] * v[3];
    }

	inline void transformVector3(const float (&v)[3], const GLfloat *m, float (&vOut)[3])
	{ // assumes no perspective in m and v[3]==1
		vOut[0] = m[0] * v[0] + m[4] * v[1] + m[8] *  v[2] + m[12]; 
		vOut[1] = m[1] * v[0] + m[5] * v[1] + m[9] *  v[2] + m[13];	
		vOut[2] = m[2] * v[0] + m[6] * v[1] + m[10] * v[2] + m[14];	
    }

#endif	// __GLMATRICES_H__
