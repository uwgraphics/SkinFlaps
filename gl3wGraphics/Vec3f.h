/** This class is from Robert Osfield's OpenSceneGraph and is repeated here
  * for convenience. License for its use is the same as stated in the
  * GNU pulic license for the rest of that software suite. It can be found
  * at http://www.openscenegraph.org
  * General purpose float triple for use as vertices, vectors and normals.
  * Provides general math operations from addition through to cross products.
  * No support yet added for float * Vec3f - is it necessary?
  * Need to define a non-member non-friend operator*  etc.
  * Vec3f * float is okay
*/

#ifndef _Vec3f_
#define _Vec3f_

#include <cmath>
#include <cstdint>

class Vec3f
{
    public:

        /** Type of Vec class.*/
        typedef float value_type;

        /** Number of vector components. */
        enum { num_components = 3 };
        
        /** Vec member variable. */
		union {
			struct {
				value_type X, Y, Z;
			};
			value_type xyz[3];
		};

        

        Vec3f() { xyz[0]=0.0f; xyz[1]=0.0f; xyz[2]=0.0f;}
        Vec3f(value_type x,value_type y,value_type z) { xyz[0]=x; xyz[1]=y; xyz[2]=z; }
        Vec3f(value_type (&v)[3]) { xyz[0]=v[0]; xyz[1]=v[1]; xyz[2]=v[2]; }
		Vec3f(const short(&v)[3]) { xyz[0] = (value_type)v[0]; xyz[1] = (value_type)v[1]; xyz[2] = (value_type)v[2]; }
		Vec3f(const unsigned short(&v)[3]) { xyz[0] = (value_type)v[0]; xyz[1] = (value_type)v[1]; xyz[2] = (value_type)v[2]; }
		Vec3f(int(&v)[3]) { xyz[0] = (value_type)v[0]; xyz[1] = (value_type)v[1]; xyz[2] = (value_type)v[2]; }
		Vec3f(const double(&v)[3]) { xyz[0] = (value_type)v[0]; xyz[1] = (value_type)v[1]; xyz[2] = (value_type)v[2]; }

        inline bool operator == (const Vec3f& v) const { return xyz[0]==v.xyz[0] && xyz[1]==v.xyz[1] && xyz[2]==v.xyz[2]; }
        
        inline bool operator != (const Vec3f& v) const { return xyz[0]!=v.xyz[0] || xyz[1]!=v.xyz[1] || xyz[2]!=v.xyz[2]; }

        inline bool operator < (const Vec3f& v) const {  // for std ordered associative containers
            if (xyz[0] < v.xyz[0]) return true; if (xyz[0] > v.xyz[0]) return false;
            if (xyz[1] < v.xyz[1]) return true; if (xyz[1] > v.xyz[1]) return false;
            return xyz[2] < v.xyz[2]; }

        inline value_type* ptr() { return xyz; }
        inline const value_type* ptr() const { return xyz; }

        inline void set( value_type x, value_type y, value_type z)
        {
            xyz[0]=x; xyz[1]=y; xyz[2]=z;
        }

        inline void set( const Vec3f& rhs)
        {
            xyz[0]=rhs.xyz[0]; xyz[1]=rhs.xyz[1]; xyz[2]=rhs.xyz[2];
        }

        inline void set( const float (&rhs)[3])
        {
            xyz[0]=rhs[0]; xyz[1]=rhs[1]; xyz[2]=rhs[2];
        }

		inline void set(const double(&rhs)[3])
		{
			xyz[0] = (float)rhs[0]; xyz[1] = (float)rhs[1]; xyz[2] = (float)rhs[2];
		}

        inline void set(const short(&rhs)[3])
        {
            xyz[0] = (float)rhs[0]; xyz[1] = (float)rhs[1]; xyz[2] = (float)rhs[2];
        }

        inline void set(const int (&rhs)[3])
        {
            xyz[0] = (float)rhs[0]; xyz[1] = (float)rhs[1]; xyz[2] = (float)rhs[2];
        }

        inline value_type& operator [] (int i) { return xyz[i]; }
        inline value_type operator [] (int i) const { return xyz[i]; }

        inline bool valid() const { return !isNaN(); }
		inline bool isNaN() const { return std::isnan(xyz[0]) || std::isnan(xyz[1]) || std::isnan(xyz[2]); }

        /** Dot product. */
        inline value_type operator * (const Vec3f& rhs) const
        {
            return xyz[0]*rhs.xyz[0]+xyz[1]*rhs.xyz[1]+xyz[2]*rhs.xyz[2];
        }

        /** Cross product. */
        inline const Vec3f operator ^ (const Vec3f& rhs) const
        {
            return Vec3f(xyz[1]*rhs.xyz[2]-xyz[2]*rhs.xyz[1],
                         xyz[2]*rhs.xyz[0]-xyz[0]*rhs.xyz[2] ,
                         xyz[0]*rhs.xyz[1]-xyz[1]*rhs.xyz[0]);
        }

        /** Multiply by scalar. */
        inline const Vec3f operator * (value_type rhs) const
        {
            return Vec3f(xyz[0]*rhs, xyz[1]*rhs, xyz[2]*rhs);
        }

        /** Unary multiply by scalar. */
        inline Vec3f& operator *= (value_type rhs)
        {
            xyz[0]*=rhs;
            xyz[1]*=rhs;
            xyz[2]*=rhs;
            return *this;
        }

        /** Divide by scalar. */
        inline const Vec3f operator / (value_type rhs) const
        {
            return Vec3f(xyz[0]/rhs, xyz[1]/rhs, xyz[2]/rhs);
        }

        /** Unary divide by scalar. */
        inline Vec3f& operator /= (value_type rhs)
        {
            xyz[0]/=rhs;
            xyz[1]/=rhs;
            xyz[2]/=rhs;
            return *this;
        }

        /** Binary vector add. */
        inline const Vec3f operator + (const Vec3f& rhs) const
        {
            return Vec3f(xyz[0]+rhs.xyz[0], xyz[1]+rhs.xyz[1], xyz[2]+rhs.xyz[2]);
        }

        /** Unary vector add. Slightly more efficient because no temporary
          * intermediate object.
        */
        inline Vec3f& operator += (const Vec3f& rhs)
        {
            xyz[0] += rhs.xyz[0];
            xyz[1] += rhs.xyz[1];
            xyz[2] += rhs.xyz[2];
            return *this;
        }

        /** Binary vector subtract. */
        inline const Vec3f operator - (const Vec3f& rhs) const
        {
            return Vec3f(xyz[0]-rhs.xyz[0], xyz[1]-rhs.xyz[1], xyz[2]-rhs.xyz[2]);
        }

        /** Unary vector subtract. */
        inline Vec3f& operator -= (const Vec3f& rhs)
        {
            xyz[0]-=rhs.xyz[0];
            xyz[1]-=rhs.xyz[1];
            xyz[2]-=rhs.xyz[2];
            return *this;
        }

        /** Negation operator. Returns the negative of the Vec3f. */
        inline const Vec3f operator - () const
        {
            return Vec3f (-xyz[0], -xyz[1], -xyz[2]);
        }

        /** Length of the vector = sqrt( vec . vec ) */
        inline value_type length() const
        {
            return sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
        }

        /** Length squared of the vector = vec . vec */
        inline value_type length2() const
        {
            return xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2];
        }

        /** Normalize the vector so that it has length unity.
          * Returns the previous length of the vector.
        */
        inline value_type normalize()
        {
            value_type norm = Vec3f::length();
            if (norm>0.0)
            {
                value_type inv = 1.0f/norm;
                xyz[0] *= inv;
                xyz[1] *= inv;
                xyz[2] *= inv;
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

		inline const Vec3f floor() const
		{
			return Vec3f(std::floor(xyz[0]), std::floor(xyz[1]), std::floor(xyz[2]));
		}

		inline const Vec3f fraction () const
		{
			return Vec3f(*this - this->floor());
		}

};    // end of class Vec3f
#endif	// ifndef _Vec3f_
