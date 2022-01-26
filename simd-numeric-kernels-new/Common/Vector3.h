//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################

#pragma once

#include "stdlib.h"
#include <type_traits>
namespace SIMD_Numeric_Kernel {

    template<class Tn> struct Vector3 {
        typedef Tn TnType;
        typedef typename std::decay<Tn>::type TnBase;
        Tn x, y, z;

    Vector3() :
        x(), y(), z()
            {}

    Vector3(const Vector3& v) :
        x(v.x), y(v.y), z(v.z)
            {}

        template<class VectorClass>
        explicit Vector3(VectorClass& v) :
        x(v.x), y(v.y), z(v.z)
            {}

        explicit Vector3(const Tn& x_in, const Tn& y_in, const Tn& z_in) :
        x(x_in), y(y_in), z(z_in)
            {}


        Vector3& operator=(const Vector3 &v){
            static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                           "Implicit assignment operator (Vector3) cannot be used for const-types.");
            x = v.x;
            y = v.y;
            z = v.z;
            return *this;
        }

            template<class U>
            Vector3& operator=(const Vector3<U> &v){
                static_assert( std::is_same<typename std::decay<Tn>::type, typename std::decay<U>::type >::value,
                               "Assignment (Vector3) must operate on same base type." );
                static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                               "Implicit assignment operator (Vector3) cannot be used for const-types." );
                x = v.x;
                y = v.y;
                z = v.z;
                return *this;
            }

                template<class T_DATA>
                void Load_Aligned(const T_DATA (&A)[3]){
                    static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                                   "Load Aligned (Vector3) cannot be used for const-types." );
                    x.Load_Aligned(A[0]);
                    y.Load_Aligned(A[1]);
                    z.Load_Aligned(A[2]);
                }

        template<class T_DATA>
        void Load_Aligned(const T_DATA (&A1), const T_DATA (&A2), const T_DATA (&A3)){
            static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                           "Load Aligned (Vector3) cannot be used for const-types." );
            x.Load_Aligned(A1);
            y.Load_Aligned(A2);
            z.Load_Aligned(A3);
        }

        template<class T_DATA>
        void Store(T_DATA (&A)[3]) const {
            SIMD_Numeric_Kernel::Store(A[0],x);
            SIMD_Numeric_Kernel::Store(A[1],y);
            SIMD_Numeric_Kernel::Store(A[2],z);
        }

        template<class T_DATA>
        void Store(T_DATA (&A1),T_DATA (&A2),T_DATA (&A3) ) const {
            SIMD_Numeric_Kernel::Store(A1,x);
            SIMD_Numeric_Kernel::Store(A2,y);
            SIMD_Numeric_Kernel::Store(A3,z);
        }

        Tn& operator()(int i){
            if( i<1 || i>3)
                exit(1);
            if(i==1)
                return x;
            if(i==2)
                return y;
            if(i==3)
                return z;
        }

        const Tn& operator()(int i) const {
            if( i<1 || i>3)
                exit(1);
            if(i==1)
                return x;
            if(i==2)
                return y;
            if(i==3)
                return z;
        }

        // Scalar Math

        Vector3& operator+=( const Tn& n) {
            static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                           "Scalar operator (Vector3) cannot be used for const-types." );
            x = x+n;
            y = y+n;
            z = z+n;
            return *this;
        }

        Vector3& operator-=( const Tn& n) {
            static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                           "Scalar operator (Vector3) cannot be used for const-types." );
            x = x-n;
            y = y-n;
            z = z-n;
            return *this;
        }

        Vector3& operator*=( const Tn& n) {
            static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                           "Scalar operator (Vector3) cannot be used for const-types." );
            x = x*n;
            y = y*n;
            z = z*n;
            return *this;
        }

        Vector3& operator/=( const Tn& n) {
            static_assert( !std::is_const<std::remove_reference<Tn> >::value,
                           "Scalar operator (Vector3) cannot be used for const-types." );
            x = x/n;
            y = y/n;
            z = z/n;
            return *this;
        }



        Vector3 operator+( const Tn& n) const {
            Vector3 c;
            c.x = x+n;
            c.y = y+n;
            c.z = z+n;
            return c;
        }

        Vector3 operator-( const Tn& n) const {
            Vector3 c;
            c.x = x-n;
            c.y = y-n;
            c.z = z-n;
            return c;
        }

        Vector3 operator*( const Tn& n) const {
            Vector3 c;
            c.x = x*n;
            c.y = y*n;
            c.z = z*n;
            return c;
        }

        Vector3 operator/( const Tn& n) const {
            Vector3 c;
            c.x = x/n;
            c.y = y/n;
            c.z = z/n;
            return c;
        }

        template<class U>
        Vector3 operator+( const Vector3<U> &n) const {
            Vector3 c;
            c.x = x+n.x;
            c.y = y+n.y;
            c.z = z+n.z;
            return c;
        }

        template<class U>
        Vector3 operator-( const Vector3<U> &n) const {
            Vector3 c;
            c.x = x-n.x;
            c.y = y-n.y;
            c.z = z-n.z;
            return c;
        }

        template<class U>
        Vector3 operator*( const Vector3<U> &n) const {
            Vector3 c;
            c.x = x*n.x;
            c.y = y*n.y;
            c.z = z*n.z;
            return c;
        }

        template<class U>
        Vector3 operator/( const Vector3<U> &n) const {
            Vector3 c;
            c.x = x/n.x;
            c.y = y/n.y;
            c.z = z/n.z;
            return c;
        }

        // Operations

        template<class U>
        Tn DotProduct( const Vector3<U> &n ) const {
            Tn result;
            result = x*n.x + y*n.y + z*n.z;
            return result;
        }

        template<class U>
        Vector3 CrossProduct( const Vector3<U> &n) const {
            return Vector3(y*n.z-z*n.y,z*n.x-x*n.z,x*n.y-y*n.x);
        }

        Tn Trace() const {
            return x + y +z;
        }

        Tn Magnitude_Squared() const {
            return x*x + y*y + z*z;
        }

        Tn Magnitude() const {
            return (x*x + y*y + z*z).sqrt();
        }
    };
}
