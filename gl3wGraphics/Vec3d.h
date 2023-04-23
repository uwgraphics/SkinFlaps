/** This class is from Robert Osfield's OpenSceneGraph and is repeated here
  * for convenience. License for its use is the same as stated in the
  * GNU pulic license for the rest of that software suite. It can be found
  * at http://www.openscenegraph.org
  * General purpose float triple for use as vertices, vectors and normals.
  * Provides general math operations from addition through to cross products.
  * No support yet added for float * Vec3d - is it necessary?
  * Need to define a non-member non-friend operator*  etc.
  * Vec3d * float is okay
*/

#ifndef _Vec3d_
#define _Vec3d_

#include "Vec3f.h"
#include <math.h>

class Vec3d
{
    public:

        /** Type of Vec class.*/
        typedef double value_type;

        /** Number of vector components. */
        enum { num_components = 3 };

		/** Vec member variable. */
		union {
			struct {
				value_type X, Y, Z;
			};
			value_type xyz[3];
		};

        Vec3d() { X=0.0; Y=0.0; Z=0.0;}
        Vec3d(value_type x,value_type y,value_type z) { X=x; Y=y; Z=z; }
        Vec3d(value_type (&v)[3]) { X=v[0]; Y=v[1]; Z=v[2]; }
		Vec3d(const float(&v)[3]) { X = (double)v[0]; Y = (double)v[1]; Z = (double)v[2]; }
		Vec3d(const short(&v)[3]) { X = (double)v[0]; Y = (double)v[1]; Z = (double)v[2]; }
		Vec3d(const unsigned short(&v)[3]) { X = (double)v[0]; Y = (double)v[1]; Z = (double)v[2]; }
		Vec3d(const Vec3f &v) { X = (double)v.X; Y = (double)v.Y; Z = (double)v.Z; }

        inline bool operator == (const Vec3d& v) const { return X==v.X && Y==v.Y && Z==v.Z; }
        
        inline bool operator != (const Vec3d& v) const { return X!=v.X || Y!=v.Y || Z!=v.Z; }

        inline value_type* ptr() { return xyz; }
        inline const value_type* ptr() const { return xyz; }

        inline void set( value_type x, value_type y, value_type z)
        {
            X=x; Y=y; Z=z;
        }

        inline void set( const Vec3d& rhs)
        {
            X=rhs.X; Y=rhs.Y; Z=rhs.Z;
        }

		inline void set(const short(&rhs)[3])
		{
			X = rhs[0]; Y = rhs[1]; Z = rhs[2];
		}

		inline void set(const unsigned short(&rhs)[3])
		{
			X = rhs[0]; Y = rhs[1]; Z = rhs[2];
		}

		inline void set(const double(&rhs)[3])
        {
            X=rhs[0]; Y=rhs[1]; Z=rhs[2];
        }

		inline void set(const float(&rhs)[3])
		{
			X = (double)rhs[0]; Y = (double)rhs[1]; Z = (double)rhs[2];
		}

		inline void set(const Vec3f& rhs)
		{
			X = rhs.X; Y = rhs.Y; Z = rhs.Z;
		}

		inline value_type& operator [] (int i) { return xyz[i]; }
        inline value_type operator [] (int i) const { return xyz[i]; }

        inline value_type& x() { return X; }
        inline value_type& y() { return Y; }
        inline value_type& z() { return Z; }

        inline value_type x() const { return X; }
        inline value_type y() const { return Y; }
        inline value_type z() const { return Z; }

        inline bool valid() const { return !isNaN(); }
		inline bool isNaN() const { return std::isnan(X) || std::isnan(Y) || std::isnan(Z); }

        /** Dot product. */
        inline value_type operator * (const Vec3d& rhs) const
        {
            return X*rhs.X+Y*rhs.Y+Z*rhs.Z;
        }

		inline value_type operator * (const Vec3f& rhs) const
		{
			return X * rhs.X + Y * rhs.Y + Z * rhs.Z;
		}

		/** Cross product. */
        inline const Vec3d operator ^ (const Vec3d& rhs) const
        {
            return Vec3d(Y*rhs.Z-Z*rhs.Y,
                         Z*rhs.X-X*rhs.Z,
                         X*rhs.Y-Y*rhs.X);
        }

		/** Cross product. */
		inline const Vec3d operator ^ (const Vec3f& rhs) const
		{
			return Vec3d(Y * rhs.Z - Z * rhs.Y,
				Z * rhs.X - X * rhs.Z,
				X * rhs.Y - Y * rhs.X);
		}

		/** Multiply by scalar. */
        inline const Vec3d operator * (value_type rhs) const
        {
            return Vec3d(X*rhs, Y*rhs, Z*rhs);
        }

        /** Unary multiply by scalar. */
        inline Vec3d& operator *= (value_type rhs)
        {
            X*=rhs;
            Y*=rhs;
            Z*=rhs;
            return *this;
        }

        /** Divide by scalar. */
        inline const Vec3d operator / (value_type rhs) const
        {
            return Vec3d(X/rhs, Y/rhs, Z/rhs);
        }

        /** Unary divide by scalar. */
        inline Vec3d& operator /= (value_type rhs)
        {
            X/=rhs;
            Y/=rhs;
            Z/=rhs;
            return *this;
        }

        /** Binary vector add. */
        inline const Vec3d operator + (const Vec3d& rhs) const
        {
            return Vec3d(X+rhs.X, Y+rhs.Y, Z+rhs.Z);
        }

        /** Unary vector add. Slightly more efficient because no temporary
          * intermediate object.
        */
        inline Vec3d& operator += (const Vec3d& rhs)
        {
            X += rhs.X;
            Y += rhs.Y;
            Z += rhs.Z;
            return *this;
        }

        /** Binary vector subtract. */
        inline const Vec3d operator - (const Vec3d& rhs) const
        {
            return Vec3d(X-rhs.X, Y-rhs.Y, Z-rhs.Z);
        }

        /** Unary vector subtract. */
        inline Vec3d& operator -= (const Vec3d& rhs)
        {
            X -=rhs.X;
            Y -=rhs.Y;
            Z -= rhs.Z;
            return *this;
        }

        /** Negation operator. Returns the negative of the Vec3d. */
        inline const Vec3d operator - () const
        {
            return Vec3d (-X, -Y, -Z);
        }

        /** Length of the vector = sqrt( vec . vec ) */
        inline value_type length() const
        {
            return sqrt( X*X + Y*Y + Z*Z );
        }

        /** Length squared of the vector = vec . vec */
        inline value_type length2() const
        {
            return X * X + Y * Y + Z * Z;
        }

        /** Normalize the vector so that it has length unity.
          * Returns the previous length of the vector.
        */
        inline value_type normalize()
        {
            value_type norm = Vec3d::length();
            if (norm>0.0)
            {
                value_type inv = 1.0f/norm;
                X *= inv;
                Y *= inv;
                Z *= inv;
            }                
            return( norm );
        }

        inline void q_normalize() {  // does very fast, very approximate normalization (Quake algorithm).  Court added.  Better done in assembly code.
            float l2 = length2();
            union {
                float    f;
                uint32_t i;
            } conv;
            conv.f = l2;
            conv.i = 0x5f3759df - (conv.i >> 1);
            conv.f *= 1.5F - (l2 * 0.5F * conv.f * conv.f);
            X *= conv.f;
            Y *= conv.f;
            Z *= conv.f;
        }

		inline const Vec3d floor() const
		{
			return Vec3d(std::floor(X), std::floor(Y), std::floor(Z));
		}

		inline const Vec3d fraction() const
		{
			return Vec3d(*this - this->floor());
		}

};    // end of class Vec3d
#endif	// ifndef _Vec3d_
