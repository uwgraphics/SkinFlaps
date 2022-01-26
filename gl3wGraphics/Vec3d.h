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
			value_type _v[3];
		};

        Vec3d() { _v[0]=0.0; _v[1]=0.0; _v[2]=0.0;}
        Vec3d(value_type x,value_type y,value_type z) { _v[0]=x; _v[1]=y; _v[2]=z; }
        Vec3d(value_type (&v)[3]) { _v[0]=v[0]; _v[1]=v[1]; _v[2]=v[2]; }
		Vec3d(const float(&v)[3]) { _v[0] = (double)v[0]; _v[1] = (double)v[1]; _v[2] = (double)v[2]; }
		Vec3d(const short(&v)[3]) { _v[0] = (double)v[0]; _v[1] = (double)v[1]; _v[2] = (double)v[2]; }
		Vec3d(const unsigned short(&v)[3]) { _v[0] = (double)v[0]; _v[1] = (double)v[1]; _v[2] = (double)v[2]; }
		Vec3d(const Vec3f &v) { X = (double)v.X; Y = (double)v.Y; Z = (double)v.Z; }

        inline bool operator == (const Vec3d& v) const { return _v[0]==v._v[0] && _v[1]==v._v[1] && _v[2]==v._v[2]; }
        
        inline bool operator != (const Vec3d& v) const { return _v[0]!=v._v[0] || _v[1]!=v._v[1] || _v[2]!=v._v[2]; }

        inline value_type* ptr() { return _v; }
        inline const value_type* ptr() const { return _v; }

        inline void set( value_type x, value_type y, value_type z)
        {
            _v[0]=x; _v[1]=y; _v[2]=z;
        }

        inline void set( const Vec3d& rhs)
        {
            _v[0]=rhs._v[0]; _v[1]=rhs._v[1]; _v[2]=rhs._v[2];
        }

		inline void set(const short(&rhs)[3])
		{
			_v[0] = rhs[0]; _v[1] = rhs[1]; _v[2] = rhs[2];
		}

		inline void set(const unsigned short(&rhs)[3])
		{
			_v[0] = rhs[0]; _v[1] = rhs[1]; _v[2] = rhs[2];
		}

		inline void set(const double(&rhs)[3])
        {
            _v[0]=rhs[0]; _v[1]=rhs[1]; _v[2]=rhs[2];
        }

		inline void set(const float(&rhs)[3])
		{
			_v[0] = (double)rhs[0]; _v[1] = (double)rhs[1]; _v[2] = (double)rhs[2];
		}

		inline void set(const Vec3f& rhs)
		{
			_v[0] = rhs._v[0]; _v[1] = rhs._v[1]; _v[2] = rhs._v[2];
		}

		inline value_type& operator [] (int i) { return _v[i]; }
        inline value_type operator [] (int i) const { return _v[i]; }

        inline value_type& x() { return _v[0]; }
        inline value_type& y() { return _v[1]; }
        inline value_type& z() { return _v[2]; }

        inline value_type x() const { return _v[0]; }
        inline value_type y() const { return _v[1]; }
        inline value_type z() const { return _v[2]; }

        inline bool valid() const { return !isNaN(); }
		inline bool isNaN() const { return std::isnan(_v[0]) || std::isnan(_v[1]) || std::isnan(_v[2]); }

        /** Dot product. */
        inline value_type operator * (const Vec3d& rhs) const
        {
            return _v[0]*rhs._v[0]+_v[1]*rhs._v[1]+_v[2]*rhs._v[2];
        }

		inline value_type operator * (const Vec3f& rhs) const
		{
			return _v[0] * rhs._v[0] + _v[1] * rhs._v[1] + _v[2] * rhs._v[2];
		}

		/** Cross product. */
        inline const Vec3d operator ^ (const Vec3d& rhs) const
        {
            return Vec3d(_v[1]*rhs._v[2]-_v[2]*rhs._v[1],
                         _v[2]*rhs._v[0]-_v[0]*rhs._v[2] ,
                         _v[0]*rhs._v[1]-_v[1]*rhs._v[0]);
        }

		/** Cross product. */
		inline const Vec3d operator ^ (const Vec3f& rhs) const
		{
			return Vec3d(_v[1] * rhs._v[2] - _v[2] * rhs._v[1],
				_v[2] * rhs._v[0] - _v[0] * rhs._v[2],
				_v[0] * rhs._v[1] - _v[1] * rhs._v[0]);
		}

		/** Multiply by scalar. */
        inline const Vec3d operator * (value_type rhs) const
        {
            return Vec3d(_v[0]*rhs, _v[1]*rhs, _v[2]*rhs);
        }

        /** Unary multiply by scalar. */
        inline Vec3d& operator *= (value_type rhs)
        {
            _v[0]*=rhs;
            _v[1]*=rhs;
            _v[2]*=rhs;
            return *this;
        }

        /** Divide by scalar. */
        inline const Vec3d operator / (value_type rhs) const
        {
            return Vec3d(_v[0]/rhs, _v[1]/rhs, _v[2]/rhs);
        }

        /** Unary divide by scalar. */
        inline Vec3d& operator /= (value_type rhs)
        {
            _v[0]/=rhs;
            _v[1]/=rhs;
            _v[2]/=rhs;
            return *this;
        }

        /** Binary vector add. */
        inline const Vec3d operator + (const Vec3d& rhs) const
        {
            return Vec3d(_v[0]+rhs._v[0], _v[1]+rhs._v[1], _v[2]+rhs._v[2]);
        }

        /** Unary vector add. Slightly more efficient because no temporary
          * intermediate object.
        */
        inline Vec3d& operator += (const Vec3d& rhs)
        {
            _v[0] += rhs._v[0];
            _v[1] += rhs._v[1];
            _v[2] += rhs._v[2];
            return *this;
        }

        /** Binary vector subtract. */
        inline const Vec3d operator - (const Vec3d& rhs) const
        {
            return Vec3d(_v[0]-rhs._v[0], _v[1]-rhs._v[1], _v[2]-rhs._v[2]);
        }

        /** Unary vector subtract. */
        inline Vec3d& operator -= (const Vec3d& rhs)
        {
            _v[0]-=rhs._v[0];
            _v[1]-=rhs._v[1];
            _v[2]-=rhs._v[2];
            return *this;
        }

        /** Negation operator. Returns the negative of the Vec3d. */
        inline const Vec3d operator - () const
        {
            return Vec3d (-_v[0], -_v[1], -_v[2]);
        }

        /** Length of the vector = sqrt( vec . vec ) */
        inline value_type length() const
        {
            return sqrt( _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] );
        }

        /** Length squared of the vector = vec . vec */
        inline value_type length2() const
        {
            return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2];
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
                _v[0] *= inv;
                _v[1] *= inv;
                _v[2] *= inv;
            }                
            return( norm );
        }

		inline const Vec3d floor() const
		{
			return Vec3d(std::floor(_v[0]), std::floor(_v[1]), std::floor(_v[2]));
		}

		inline const Vec3d fraction() const
		{
			return Vec3d(*this - this->floor());
		}

};    // end of class Vec3d
#endif	// ifndef _Vec3d_
