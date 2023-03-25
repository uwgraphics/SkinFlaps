/* -*-c++-*- OpenSceneGraph - Copyright (C) 1998-2006 Robert Osfield
*
* This library is open source and may be redistributed and/or modified under
* the terms of the OpenSceneGraph Public License (OSGPL) version 0.0 or
* (at your option) any later version.  The full license is in LICENSE file
* included with this distribution, and on the openscenegraph.org website.
*
* This library is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* OpenSceneGraph Public License for more details.
*/

#ifndef __VEC2D__
#define __VEC2D__

#include <math.h>

/** General purpose double pair. Uses include representation of
* texture coordinates.
* No support yet added for double * Vec2d - is it necessary?
* Need to define a non-member non-friend operator* etc.
* BTW: Vec2d * double is okay
*/

class Vec2d
{
public:
	/** Type of Vec class.*/
	typedef double value_type;
	/** Number of vector components. */
	enum { num_components = 2 };
	/** Vec member varaible. */
	union {
		struct {
			value_type X, Y;
		};
		value_type xy[2];
	};


	Vec2d() { xy[0] = 0.0; xy[1] = 0.0; }
	Vec2d(value_type x, value_type y) { xy[0] = x; xy[1] = y; }

	inline bool operator == (const Vec2d& v) const { return xy[0] == v.xy[0] && xy[1] == v.xy[1]; }

	inline bool operator != (const Vec2d& v) const { return xy[0] != v.xy[0] || xy[1] != v.xy[1]; }

	inline bool operator <  (const Vec2d& v) const
	{
		if (xy[0]<v.xy[0]) return true;
		else if (xy[0]>v.xy[0]) return false;
		else return (xy[1]<v.xy[1]);
	}

	inline value_type * ptr() { return xy; }
	inline const value_type * ptr() const { return xy; }

	inline void set(value_type x, value_type y) { xy[0] = x; xy[1] = y; }

	inline void set(const float x, const float y) { xy[0] = (double)x; xy[1] = (double)y; }

	inline void set(const float tx[2]) { xy[0] = (double)tx[0]; xy[1] = (double)tx[1]; }

	inline value_type & operator [] (int i) { return xy[i]; }
	inline value_type operator [] (int i) const { return xy[i]; }

	inline value_type & x() { return xy[0]; }
	inline value_type & y() { return xy[1]; }

	inline value_type x() const { return xy[0]; }
	inline value_type y() const { return xy[1]; }

	/** Dot product. */
	inline value_type operator * (const Vec2d& rhs) const
	{
		return xy[0] * rhs.xy[0] + xy[1] * rhs.xy[1];
	}

	/** Multiply by scalar. */
	inline const Vec2d operator * (value_type rhs) const
	{
		return Vec2d(xy[0] * rhs, xy[1] * rhs);
	}

	/** Unary multiply by scalar. */
	inline Vec2d& operator *= (value_type rhs)
	{
		xy[0] *= rhs;
		xy[1] *= rhs;
		return *this;
	}

	/** Divide by scalar. */
	inline const Vec2d operator / (value_type rhs) const
	{
		return Vec2d(xy[0] / rhs, xy[1] / rhs);
	}

	/** Unary divide by scalar. */
	inline Vec2d& operator /= (value_type rhs)
	{
		xy[0] /= rhs;
		xy[1] /= rhs;
		return *this;
	}

	/** Binary vector add. */
	inline const Vec2d operator + (const Vec2d& rhs) const
	{
		return Vec2d(xy[0] + rhs.xy[0], xy[1] + rhs.xy[1]);
	}

	// Unary vector add. Slightly more efficient because no temporary intermediate object.
	inline Vec2d& operator += (const Vec2d& rhs)
	{
		xy[0] += rhs.xy[0];
		xy[1] += rhs.xy[1];
		return *this;
	}

	// Binary vector subtract.
	inline const Vec2d operator - (const Vec2d& rhs) const
	{
		return Vec2d(xy[0] - rhs.xy[0], xy[1] - rhs.xy[1]);
	}

	/** Unary vector subtract. */
	inline Vec2d& operator -= (const Vec2d& rhs)
	{
		xy[0] -= rhs.xy[0];
		xy[1] -= rhs.xy[1];
		return *this;
	}

	/** Negation operator. Returns the negative of the Vec2d. */
	inline const Vec2d operator - () const
	{
		return Vec2d(-xy[0], -xy[1]);
	}

	/** Length of the vector = sqrt( vec . vec ) */
	inline value_type length() const
	{
		return sqrt(xy[0] * xy[0] + xy[1] * xy[1]);
	}

	/** Length squared of the vector = vec . vec */
	inline value_type length2(void) const
	{
		return xy[0] * xy[0] + xy[1] * xy[1];
	}

	// Normalize the vector so that it has length unity.
	// Returns the previous length of the vector.
	inline value_type normalize()
	{
		value_type norm = Vec2d::length();
		if (norm>0.0)
		{
			value_type inv = 1.0f / norm;
			xy[0] *= inv;
			xy[1] *= inv;
		}
		return(norm);
	}

};    // end of class Vec2d

#endif

