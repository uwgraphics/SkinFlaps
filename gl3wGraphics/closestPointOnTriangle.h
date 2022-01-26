// does what it says

#ifndef __CLOSE_TRIANGLE_POINT__
#define __CLOSE_TRIANGLE_POINT__

#include "Vec3f.h"

class closestPointOnTriangle
{
private:
	Vec3f _U, _V, _vert0, _t0minusp;
	float s, t;

public:
	closestPointOnTriangle()
	{
	}

	void getClosePoint(const Vec3f &point, const Vec3f triangle[3])
	{
		_U -= triangle[1] - triangle[0];
		_V -= triangle[2] - triangle[0];
		_vert0 = triangle[0];
		_t0minusp = triangle[0] - point;
		float a = _U*_U, b = _U*_V, c = _V*_V, d = _U*_t0minusp, e = _V*_t0minusp;
		float det = a*c - b*b;
		s = b*e - c*d;
		t = b*d - a*e;
		if (s + t < det){
			if (s < 0.0f){
				if (t < 0.0f){
					if (d < 0.0f){
						s = -d / a;
						if (s < 0.0f) s = 0.0f;
						else if (s > 1.0f) s = 1.0f;
						else;
						t = 0.0f;
					}
					else{
						s = 0.0f;
						t = -e / c;
						if (t < 0.0f) t = 0.0f;
						else if (t > 1.0f) t = 1.0f;
						else;
					}
				}
				else{
					s = 0.0f;
					t = -e / c;
					if (t < 0.0f) t = 0.0f;
					else if (t > 1.0f) t = 1.0f;
					else;
				}
			}
			else if (t < 0.0f){
				s = -d / a;
				if (s < 0.0f) s = 0.0f;
				else if (s > 1.0f) s = 1.0f;
				else;
				t = 0.0f;
			}
			else{
				float invDet = 1.0f / det;
				s *= invDet;
				t *= invDet;
			}
		}
		else{
			if (s < 0.0f){
				float tmp0 = b + d, tmp1 = c + e;
				if (tmp1 > tmp0){
					float numer = tmp1 - tmp0, denom = a - 2 * b + c;
					s = numer / denom;
					if (s < 0.0f) s = 0.0f;
					else if (s > 1.0f) s = 1.0f;
					else;
					t = 1.0f - s;
				}
				else{
					t = -e / c;
					if (t < 0.0f) t = 0.0f;
					else if (t > 1.0f) t = 1.0f;
					else;
					s = 0.0f;
				}
			}
			else if (t < 0.0f){
				if (a + d > b + e){
					float numer = c + e - b - d, denom = a - 2 * b + c;
					s = numer / denom;
					if (s < 0.0f) s = 0.0f;
					else if (s > 1.0f) s = 1.0f;
					else;
					t = 1.0f - s;
				}
				else{
					s = -e / c;
					if (s < 0.0f) s = 0.0f;
					else if (s > 1.0f) s = 1.0f;
					else;
					t = 0.0f;
				}
			}
			else{
				float numer = c + e - b - d, denom = a - 2 * b + c;
				s = numer / denom;
				if (s < 0.0f) s = 0.0f;
				else if (s > 1.0f) s = 1.0f;
				else;
				t = 1.0f - s;
			}
		}
	}

	inline void returnTriangleParameters(float(&params)[2]){ params[0] = s; params[1] = t; }

	inline void returnCloseTrianglePoint(Vec3f &closePoint){ closePoint = _vert0 + _U*s + _V*t; }

	inline float returnSignedDistanceSquared(){  // sign reflects which side of triangle the imput point is on
		float len2 = (_U*s + _V*t + _t0minusp).length2();
		return (_t0minusp * (_U^_V) < 0.0f) ? -len2 : len2;  // _t0minusp direction reversed
	}
};

#endif  // __CLOSE_TRIANGLE_POINT__