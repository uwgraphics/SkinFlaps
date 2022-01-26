//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#ifdef __INTEL_COMPILER
#pragma warning( disable : 592 )
#endif

// #define USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION
// #define PERFORM_STRICT_QUATERNION_RENORMALIZATION

  { // Begin block : Scope of qV (if not maintained)

#ifndef COMPUTE_V_AS_QUATERNION
    Tn Vqvs;
    Tn Vqvvx;
    Tn Vqvvy;
    Tn Vqvvz;
#endif

    { // Begin block : Symmetric eigenanalysis

      Tn Vs11;
      Tn Vs21;
      Tn Vs31;
      Tn Vs22;
      Tn Vs32;
      Tn Vs33;

      Vqvs=Vone;
    
      //###########################################################
      // Compute normal equations matrix
      //###########################################################

      Vs11=(Va11 * Va11) + (Va21 * Va21) + (Va31 * Va31);
      Vs21=(Va12 * Va11) + (Va22 * Va21) + (Va32 * Va31);
      Vs31=(Va13 * Va11) + (Va23 * Va21) + (Va33 * Va31);
      Vs22=(Va12 * Va12) + (Va22 * Va22) + (Va32 * Va32);
      Vs32=(Va13 * Va12) + (Va23 * Va22) + (Va33 * Va32);
      Vs33=(Va13 * Va13) + (Va23 * Va23) + (Va33 * Va33);

      //###########################################################
      // Solve symmetric eigenproblem using Jacobi iteration
      //###########################################################

      for(int sweep=1;sweep<=4;sweep++){

        // First Jacobi conjugation

#define SS11 Ss11
#define SS21 Ss21
#define SS31 Ss31
#define SS22 Ss22
#define SS32 Ss32
#define SS33 Ss33
#define SQVVX Sqvvx
#define SQVVY Sqvvy
#define SQVVZ Sqvvz
#define STMP1 Stmp1
#define STMP2 Stmp2
#define STMP3 Stmp3

#define VS11 Vs11
#define VS21 Vs21
#define VS31 Vs31
#define VS22 Vs22
#define VS32 Vs32
#define VS33 Vs33
#define VQVVX Vqvvx
#define VQVVY Vqvvy
#define VQVVZ Vqvvz
#define VTMP1 Vtmp1
#define VTMP2 Vtmp2
#define VTMP3 Vtmp3

#include "Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.hpp"

#undef SS11
#undef SS21
#undef SS31
#undef SS22
#undef SS32
#undef SS33
#undef SQVVX
#undef SQVVY
#undef SQVVZ
#undef STMP1
#undef STMP2
#undef STMP3

#undef VS11
#undef VS21
#undef VS31
#undef VS22
#undef VS32
#undef VS33
#undef VQVVX
#undef VQVVY
#undef VQVVZ
#undef VTMP1
#undef VTMP2
#undef VTMP3

        // Second Jacobi conjugation

#define SS11 Ss22
#define SS21 Ss32
#define SS31 Ss21
#define SS22 Ss33
#define SS32 Ss31
#define SS33 Ss11
#define SQVVX Sqvvy
#define SQVVY Sqvvz
#define SQVVZ Sqvvx
#define STMP1 Stmp2
#define STMP2 Stmp3
#define STMP3 Stmp1

#define VS11 Vs22
#define VS21 Vs32
#define VS31 Vs21
#define VS22 Vs33
#define VS32 Vs31
#define VS33 Vs11
#define VQVVX Vqvvy
#define VQVVY Vqvvz
#define VQVVZ Vqvvx
#define VTMP1 Vtmp2
#define VTMP2 Vtmp3
#define VTMP3 Vtmp1

#include "Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.hpp"

#undef SS11
#undef SS21
#undef SS31
#undef SS22
#undef SS32
#undef SS33
#undef SQVVX
#undef SQVVY
#undef SQVVZ
#undef STMP1
#undef STMP2
#undef STMP3

#undef VS11
#undef VS21
#undef VS31
#undef VS22
#undef VS32
#undef VS33
#undef VQVVX
#undef VQVVY
#undef VQVVZ
#undef VTMP1
#undef VTMP2
#undef VTMP3

        // Third Jacobi conjugation

#define SS11 Ss33
#define SS21 Ss31
#define SS31 Ss32
#define SS22 Ss11
#define SS32 Ss21
#define SS33 Ss22
#define SQVVX Sqvvz
#define SQVVY Sqvvx
#define SQVVZ Sqvvy
#define STMP1 Stmp3
#define STMP2 Stmp1
#define STMP3 Stmp2

#define VS11 Vs33
#define VS21 Vs31
#define VS31 Vs32
#define VS22 Vs11
#define VS32 Vs21
#define VS33 Vs22
#define VQVVX Vqvvz
#define VQVVY Vqvvx
#define VQVVZ Vqvvy
#define VTMP1 Vtmp3
#define VTMP2 Vtmp1
#define VTMP3 Vtmp2

#include "Singular_Value_Decomposition_Jacobi_Conjugation_Kernel.hpp"

#undef SS11
#undef SS21
#undef SS31
#undef SS22
#undef SS32
#undef SS33
#undef SQVVX
#undef SQVVY
#undef SQVVZ
#undef STMP1
#undef STMP2
#undef STMP3

#undef VS11
#undef VS21
#undef VS31
#undef VS22
#undef VS32
#undef VS33
#undef VQVVX
#undef VQVVY
#undef VQVVZ
#undef VTMP1
#undef VTMP2
#undef VTMP3
      }



    } // End block : Symmetric eigenanalysis

    //###########################################################
    // Normalize quaternion for matrix V
    //###########################################################

#if !defined(USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION) || defined(PERFORM_STRICT_QUATERNION_RENORMALIZATION)

    Vtmp2=Vqvs * Vqvs;
    Vtmp1=Vqvvx * Vqvvx;
    Vtmp2=Vtmp1 + Vtmp2;
    Vtmp1=Vqvvy * Vqvvy;
    Vtmp2=Vtmp1 + Vtmp2;
    Vtmp1=Vqvvz * Vqvvz;
    Vtmp2=Vtmp1 + Vtmp2;

    Vtmp1=Vtmp2.rsqrt();
    Vtmp4=Vtmp1 * Vone_half;
    Vtmp3=Vtmp1 * Vtmp4;
    Vtmp3=Vtmp1 * Vtmp3;
    Vtmp3=Vtmp2 * Vtmp3;
    Vtmp1=Vtmp1 + Vtmp4;
    Vtmp1=Vtmp1 - Vtmp3;

    Vqvs=Vqvs * Vtmp1;
    Vqvvx=Vqvvx * Vtmp1;
    Vqvvy=Vqvvy * Vtmp1;
    Vqvvz=Vqvvz * Vtmp1;

#endif

    { // Begin block : Conjugation with V

      //###########################################################
      // Transform quaternion to matrix V
      //###########################################################

#ifndef COMPUTE_V_AS_MATRIX
      Tn Vv11;
      Tn Vv21;
      Tn Vv31;
      Tn Vv12;
      Tn Vv22;
      Tn Vv32;
      Tn Vv13;
      Tn Vv23;
      Tn Vv33;
#endif

      Vtmp1=Vqvvx * Vqvvx;
      Vtmp2=Vqvvy * Vqvvy;
      Vtmp3=Vqvvz * Vqvvz;
      Vv11=Vqvs * Vqvs;
      Vv22=Vv11 - Vtmp1;
      Vv33=Vv22 - Vtmp2;
      Vv33=Vv33 + Vtmp3;
      Vv22=Vv22 + Vtmp2;
      Vv22=Vv22 - Vtmp3;
      Vv11=Vv11 + Vtmp1;
      Vv11=Vv11 - Vtmp2;
      Vv11=Vv11 - Vtmp3;
      Vtmp1=Vqvvx + Vqvvx;
      Vtmp2=Vqvvy + Vqvvy;
      Vtmp3=Vqvvz + Vqvvz;
      Vv32=Vqvs * Vtmp1;
      Vv13=Vqvs * Vtmp2;
      Vv21=Vqvs * Vtmp3;
      Vtmp1=Vqvvy * Vtmp1;
      Vtmp2=Vqvvz * Vtmp2;
      Vtmp3=Vqvvx * Vtmp3;
      Vv12=Vtmp1 - Vv21;
      Vv23=Vtmp2 - Vv32;
      Vv31=Vtmp3 - Vv13;
      Vv21=Vtmp1 + Vv21;
      Vv32=Vtmp2 + Vv32;
      Vv13=Vtmp3 + Vv13;



      //###########################################################
      // Multiply (from the right) with V
      //###########################################################

      Vtmp2=Va12;
      Vtmp3=Va13;
      Va12=Vv12 * Va11;
      Va13=Vv13 * Va11;
      Va11=Vv11 * Va11;
      Vtmp1=Vv21 * Vtmp2;
      Va11=Va11 + Vtmp1;
      Vtmp1=Vv31 * Vtmp3;
      Va11=Va11 + Vtmp1;
      Vtmp1=Vv22 * Vtmp2;
      Va12=Va12 + Vtmp1;
      Vtmp1=Vv32 * Vtmp3;
      Va12=Va12 + Vtmp1;
      Vtmp1=Vv23 * Vtmp2;
      Va13=Va13 + Vtmp1;
      Vtmp1=Vv33 * Vtmp3;
      Va13=Va13 + Vtmp1;

      Vtmp2=Va22;
      Vtmp3=Va23;
      Va22=Vv12 * Va21;
      Va23=Vv13 * Va21;
      Va21=Vv11 * Va21;
      Vtmp1=Vv21 * Vtmp2;
      Va21=Va21 + Vtmp1;
      Vtmp1=Vv31 * Vtmp3;
      Va21=Va21 + Vtmp1;
      Vtmp1=Vv22 * Vtmp2;
      Va22=Va22 + Vtmp1;
      Vtmp1=Vv32 * Vtmp3;
      Va22=Va22 + Vtmp1;
      Vtmp1=Vv23 * Vtmp2;
      Va23=Va23 + Vtmp1;
      Vtmp1=Vv33 * Vtmp3;
      Va23=Va23 + Vtmp1;

      Vtmp2=Va32;
      Vtmp3=Va33;
      Va32=Vv12 * Va31;
      Va33=Vv13 * Va31;
      Va31=Vv11 * Va31;
      Vtmp1=Vv21 * Vtmp2;
      Va31=Va31 + Vtmp1;
      Vtmp1=Vv31 * Vtmp3;
      Va31=Va31 + Vtmp1;
      Vtmp1=Vv22 * Vtmp2;
      Va32=Va32 + Vtmp1;
      Vtmp1=Vv32 * Vtmp3;
      Va32=Va32 + Vtmp1;
      Vtmp1=Vv23 * Vtmp2;
      Va33=Va33 + Vtmp1;
      Vtmp1=Vv33 * Vtmp3;
      Va33=Va33 + Vtmp1;



    } // End block : Conjugation with V

  } // End block : Scope of qV (if not maintained)

    //###########################################################
    // Permute columns such that the singular values are sorted
    //###########################################################

Vtmp1=Va11 * Va11;
Vtmp4=Va21 * Va21;
Vtmp1=Vtmp1 + Vtmp4;
Vtmp4=Va31 * Va31;
Vtmp1=Vtmp1 + Vtmp4;

Vtmp2=Va12 * Va12;
Vtmp4=Va22 * Va22;
Vtmp2=Vtmp2 + Vtmp4;
Vtmp4=Va32 * Va32;
Vtmp2=Vtmp2 + Vtmp4;

Vtmp3=Va13 * Va13;
Vtmp4=Va23 * Va23;
Vtmp3=Vtmp3 + Vtmp4;
Vtmp4=Va33 * Va33;
Vtmp3=Vtmp3 + Vtmp4;

// Swap columns 1-2 if necessary

Mtmp4=Vtmp1 < Vtmp2;
Vtmp5=Va11 ^ Va12;
Vtmp5=Vtmp5.mask(Mtmp4);
Va11=Va11 ^ Vtmp5;
Va12=Va12 ^ Vtmp5;

Vtmp5=Va21 ^ Va22;
Vtmp5=Vtmp5.mask(Mtmp4);
Va21=Va21 ^ Vtmp5;
Va22=Va22 ^ Vtmp5;

Vtmp5=Va31 ^ Va32;
Vtmp5=Vtmp5.mask(Mtmp4);
Va31=Va31 ^ Vtmp5;
Va32=Va32 ^ Vtmp5;

#ifdef COMPUTE_V_AS_MATRIX
Vtmp5=Vv11 ^ Vv12;
Vtmp5=Vtmp5.mask(Mtmp4);
Vv11=Vv11 ^ Vtmp5;
Vv12=Vv12 ^ Vtmp5;

Vtmp5=Vv21 ^ Vv22;
Vtmp5=Vtmp5.mask(Mtmp4);
Vv21=Vv21 ^ Vtmp5;
Vv22=Vv22 ^ Vtmp5;

Vtmp5=Vv31 ^ Vv32;
Vtmp5=Vtmp5.mask(Mtmp4);
Vv31=Vv31 ^ Vtmp5;
Vv32=Vv32 ^ Vtmp5;
#endif

Vtmp5=Vtmp1 ^ Vtmp2;
Vtmp5=Vtmp5.mask(Mtmp4);
Vtmp1=Vtmp1 ^ Vtmp5;
Vtmp2=Vtmp2 ^ Vtmp5;

// If columns 1-2 have been swapped, negate 2nd column of A and
// V so that V is still a rotation

Vtmp5.Load_Aligned(_VNegTwo);
Vtmp5=Vtmp5.mask(Mtmp4);
Vtmp4=Vone;
Vtmp4=Vtmp4 + Vtmp5;

Va12=Va12 * Vtmp4;
Va22=Va22 * Vtmp4;
Va32=Va32 * Vtmp4;

#ifdef COMPUTE_V_AS_MATRIX
Vv12=Vv12 * Vtmp4;
Vv22=Vv22 * Vtmp4;
Vv32=Vv32 * Vtmp4;
#endif

// If columns 1-2 have been swapped, also update quaternion
//representation of V (the quaternion may become un-normalized after this)

#ifdef COMPUTE_V_AS_QUATERNION
Vtmp4=Vtmp4 * Vone_half;
Vtmp4=Vtmp4 - Vone_half;

Vtmp5=Vtmp4 * Vqvvz;
Vtmp5=Vtmp5 + Vqvs;
Vqvs=Vqvs * Vtmp4;
Vqvvz=Vqvvz - Vqvs;
Vqvs=Vtmp5;

Vtmp5=Vtmp4 * Vqvvx;
Vtmp5=Vtmp5 + Vqvvy;
Vqvvy=Vqvvy * Vtmp4;
Vqvvx=Vqvvx - Vqvvy;
Vqvvy=Vtmp5;
#endif

// Swap columns 1-3 if necessary

Mtmp4=Vtmp1 < Vtmp3;
Vtmp5=Va11 ^ Va13;
Vtmp5=Vtmp5.mask(Mtmp4);
Va11=Va11 ^ Vtmp5;
Va13=Va13 ^ Vtmp5;

Vtmp5=Va21 ^ Va23;
Vtmp5=Vtmp5.mask(Mtmp4);
Va21=Va21 ^ Vtmp5;
Va23=Va23 ^ Vtmp5;

Vtmp5=Va31 ^ Va33;
Vtmp5=Vtmp5.mask(Mtmp4);
Va31=Va31 ^ Vtmp5;
Va33=Va33 ^ Vtmp5;

#ifdef COMPUTE_V_AS_MATRIX
Vtmp5=Vv11 ^ Vv13;
Vtmp5=Vtmp5.mask(Mtmp4);
Vv11=Vv11 ^ Vtmp5;
Vv13=Vv13 ^ Vtmp5;

Vtmp5=Vv21 ^ Vv23;
Vtmp5=Vtmp5.mask(Mtmp4);
Vv21=Vv21 ^ Vtmp5;
Vv23=Vv23 ^ Vtmp5;

Vtmp5=Vv31 ^ Vv33;
Vtmp5=Vtmp5.mask(Mtmp4);
Vv31=Vv31 ^ Vtmp5;
Vv33=Vv33 ^ Vtmp5;
#endif

Vtmp5=Vtmp1 ^ Vtmp3;
Vtmp5=Vtmp5.mask(Mtmp4);
Vtmp1=Vtmp1 ^ Vtmp5;
Vtmp3=Vtmp3 ^ Vtmp5;

// If columns 1-3 have been swapped, negate 1st column of
// A and V so that V is still a rotation

Vtmp5.Load_Aligned(_VNegTwo);
Vtmp5=Vtmp5.mask(Mtmp4);
Vtmp4=Vone;
Vtmp4=Vtmp4 + Vtmp5;

Va11=Va11 * Vtmp4;
Va21=Va21 * Vtmp4;
Va31=Va31 * Vtmp4;

#ifdef COMPUTE_V_AS_MATRIX
Vv11=Vv11 * Vtmp4;
Vv21=Vv21 * Vtmp4;
Vv31=Vv31 * Vtmp4;
#endif

// If columns 1-3 have been swapped, also update quaternion representation
// of V (the quaternion may become un-normalized after this)

#ifdef COMPUTE_V_AS_QUATERNION
Vtmp4=Vtmp4 * Vone_half;
Vtmp4=Vtmp4 - Vone_half;

Vtmp5=Vtmp4 * Vqvvy;
Vtmp5=Vtmp5 + Vqvs;
Vqvs=Vqvs * Vtmp4;
Vqvvy=Vqvvy - Vqvs;
Vqvs=Vtmp5;

Vtmp5=Vtmp4 * Vqvvz;
Vtmp5=Vtmp5 + Vqvvx;
Vqvvx=Vqvvx * Vtmp4;
Vqvvz=Vqvvz - Vqvvx;
Vqvvx=Vtmp5;
#endif

// Swap columns 2-3 if necessary

Mtmp4=Vtmp2 < Vtmp3;
Vtmp5=Va12 ^ Va13;
Vtmp5=Vtmp5.mask(Mtmp4);
Va12=Va12 ^ Vtmp5;
Va13=Va13 ^ Vtmp5;

Vtmp5=Va22 ^ Va23;
Vtmp5=Vtmp5.mask(Mtmp4);
Va22=Va22 ^ Vtmp5;
Va23=Va23 ^ Vtmp5;

Vtmp5=Va32 ^ Va33;
Vtmp5=Vtmp5.mask(Mtmp4);
Va32=Va32 ^ Vtmp5;
Va33=Va33 ^ Vtmp5;

#ifdef COMPUTE_V_AS_MATRIX
Vtmp5=Vv12 ^ Vv13;
Vtmp5=Vtmp5.mask(Mtmp4);
Vv12=Vv12 ^ Vtmp5;
Vv13=Vv13 ^ Vtmp5;

Vtmp5=Vv22 ^ Vv23;
Vtmp5=Vtmp5.mask(Mtmp4);
Vv22=Vv22 ^ Vtmp5;
Vv23=Vv23 ^ Vtmp5;

Vtmp5=Vv32 ^ Vv33;
Vtmp5=Vtmp5.mask(Mtmp4);
Vv32=Vv32 ^ Vtmp5;
Vv33=Vv33 ^ Vtmp5;
#endif

Vtmp5=Vtmp2 ^ Vtmp3;
Vtmp5=Vtmp5.mask(Mtmp4);
Vtmp2=Vtmp2 ^ Vtmp5;
Vtmp3=Vtmp3 ^ Vtmp5;

// If columns 2-3 have been swapped, negate 3rd column of
// A and V so that V is still a rotation

Vtmp5.Load_Aligned(_VNegTwo);
Vtmp5=Vtmp5.mask(Mtmp4);
Vtmp4=Vone;
Vtmp4=Vtmp4 + Vtmp5;

Va13=Va13 * Vtmp4;
Va23=Va23 * Vtmp4;
Va33=Va33 * Vtmp4;

#ifdef COMPUTE_V_AS_MATRIX
Vv13=Vv13 * Vtmp4;
Vv23=Vv23 * Vtmp4;
Vv33=Vv33 * Vtmp4;
#endif

// If columns 2-3 have been swapped, also update quaternion
// representation of V (the quaternion may become un-normalized after this)

#ifdef COMPUTE_V_AS_QUATERNION
Vtmp4=Vtmp4 * Vone_half;
Vtmp4=Vtmp4 - Vone_half;

Vtmp5=Vtmp4 * Vqvvx;
Vtmp5=Vtmp5 + Vqvs;
Vqvs=Vqvs * Vtmp4;
Vqvvx=Vqvvx - Vqvs;
Vqvs=Vtmp5;

Vtmp5=Vtmp4 * Vqvvy;
Vtmp5=Vtmp5 + Vqvvz;
Vqvvz=Vqvvz * Vtmp4;
Vqvvy=Vqvvy - Vqvvz;
Vqvvz=Vtmp5;
#endif


//###########################################################
// Re-normalize quaternion for matrix V
//###########################################################

#ifdef COMPUTE_V_AS_QUATERNION
Vtmp2=Vqvs * Vqvs;
Vtmp1=Vqvvx * Vqvvx;
Vtmp2=Vtmp1 + Vtmp2;
Vtmp1=Vqvvy * Vqvvy;
Vtmp2=Vtmp1 + Vtmp2;
Vtmp1=Vqvvz * Vqvvz;
Vtmp2=Vtmp1 + Vtmp2;
Vtmp1=Vtmp2.rsqrt();

#ifdef PERFORM_STRICT_QUATERNION_RENORMALIZATION
Vtmp4=Vtmp1 * Vone_half;
Vtmp3=Vtmp1 * Vtmp4;
Vtmp3=Vtmp1 * Vtmp3;
Vtmp3=Vtmp2 * Vtmp3;
Vtmp1=Vtmp1 + Vtmp4;
Vtmp1=Vtmp1 - Vtmp3;
#endif

Vqvs=Vqvs * Vtmp1;
Vqvvx=Vqvvx * Vtmp1;
Vqvvy=Vqvvy * Vtmp1;
Vqvvz=Vqvvz * Vtmp1;

#endif


//###########################################################
// Construct QR factorization of A*V (=U*D) using Givens rotations
//###########################################################

#ifdef COMPUTE_U_AS_MATRIX
Vu11=Vone;
Vu21=Vu21 ^ Vu21;
Vu31=Vu31 ^ Vu31;
Vu12=Vu12 ^ Vu12;
Vu22=Vone;
Vu32=Vu32 ^ Vu32;
Vu13=Vu13 ^ Vu13;
Vu23=Vu23 ^ Vu23;
Vu33=Vone;
#endif

#ifdef COMPUTE_U_AS_QUATERNION
Vqus=Vone;
Vquvx=Vquvx ^ Vquvx;
Vquvy=Vquvy ^ Vquvy;
Vquvz=Vquvz ^ Vquvz;
#endif

// First Givens rotation

#define SAPIVOT Sa11
#define SANPIVOT Sa21
#define SA11 Sa11
#define SA21 Sa21
#define SA12 Sa12
#define SA22 Sa22
#define SA13 Sa13
#define SA23 Sa23
#define SU11 Su11
#define SU12 Su12
#define SU21 Su21
#define SU22 Su22
#define SU31 Su31
#define SU32 Su32

#define VAPIVOT Va11
#define VANPIVOT Va21
#define VA11 Va11
#define VA21 Va21
#define VA12 Va12
#define VA22 Va22
#define VA13 Va13
#define VA23 Va23
#define VU11 Vu11
#define VU12 Vu12
#define VU21 Vu21
#define VU22 Vu22
#define VU31 Vu31
#define VU32 Vu32

#include "Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp"

#undef SAPIVOT
#undef SANPIVOT
#undef SA11
#undef SA21
#undef SA12
#undef SA22
#undef SA13
#undef SA23
#undef SU11
#undef SU12
#undef SU21
#undef SU22
#undef SU31
#undef SU32

#undef VAPIVOT
#undef VANPIVOT
#undef VA11
#undef VA21
#undef VA12
#undef VA22
#undef VA13
#undef VA23
#undef VU11
#undef VU12
#undef VU21
#undef VU22
#undef VU31
#undef VU32

// Update quaternion representation of U

#ifdef COMPUTE_U_AS_QUATERNION
Vqus=Vch;
Vquvx=Vquvx ^ Vquvx;
Vquvy=Vquvy ^ Vquvy;
Vquvz=Vsh;
#endif

// Second Givens rotation

#define SAPIVOT Sa11
#define SANPIVOT Sa31
#define SA11 Sa11
#define SA21 Sa31
#define SA12 Sa12
#define SA22 Sa32
#define SA13 Sa13
#define SA23 Sa33
#define SU11 Su11
#define SU12 Su13
#define SU21 Su21
#define SU22 Su23
#define SU31 Su31
#define SU32 Su33

#define VAPIVOT Va11
#define VANPIVOT Va31
#define VA11 Va11
#define VA21 Va31
#define VA12 Va12
#define VA22 Va32
#define VA13 Va13
#define VA23 Va33
#define VU11 Vu11
#define VU12 Vu13
#define VU21 Vu21
#define VU22 Vu23
#define VU31 Vu31
#define VU32 Vu33

#include "Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp"

#undef SAPIVOT
#undef SANPIVOT
#undef SA11
#undef SA21
#undef SA12
#undef SA22
#undef SA13
#undef SA23
#undef SU11
#undef SU12
#undef SU21
#undef SU22
#undef SU31
#undef SU32

#undef VAPIVOT
#undef VANPIVOT
#undef VA11
#undef VA21
#undef VA12
#undef VA22
#undef VA13
#undef VA23
#undef VU11
#undef VU12
#undef VU21
#undef VU22
#undef VU31
#undef VU32

// Update quaternion representation of U

#ifdef COMPUTE_U_AS_QUATERNION
Vquvx=Vsh * Vquvz;
Vsh=Vsh * Vqus;
Vquvy=Vquvy - Vsh;
Vqus=Vch * Vqus;
Vquvz=Vch * Vquvz;
#endif

// Third Givens rotation

#define SAPIVOT Sa22
#define SANPIVOT Sa32
#define SA11 Sa21
#define SA21 Sa31
#define SA12 Sa22
#define SA22 Sa32
#define SA13 Sa23
#define SA23 Sa33
#define SU11 Su12
#define SU12 Su13
#define SU21 Su22
#define SU22 Su23
#define SU31 Su32
#define SU32 Su33

#define VAPIVOT Va22
#define VANPIVOT Va32
#define VA11 Va21
#define VA21 Va31
#define VA12 Va22
#define VA22 Va32
#define VA13 Va23
#define VA23 Va33
#define VU11 Vu12
#define VU12 Vu13
#define VU21 Vu22
#define VU22 Vu23
#define VU31 Vu32
#define VU32 Vu33

#include "Singular_Value_Decomposition_Givens_QR_Factorization_Kernel.hpp"

#undef SAPIVOT
#undef SANPIVOT
#undef SA11
#undef SA21
#undef SA12
#undef SA22
#undef SA13
#undef SA23
#undef SU11
#undef SU12
#undef SU21
#undef SU22
#undef SU31
#undef SU32

#undef VAPIVOT
#undef VANPIVOT
#undef VA11
#undef VA21
#undef VA12
#undef VA22
#undef VA13
#undef VA23
#undef VU11
#undef VU12
#undef VU21
#undef VU22
#undef VU31
#undef VU32

// Update quaternion representation of U

#ifdef COMPUTE_U_AS_QUATERNION
Vtmp1=Vsh * Vquvx;
Vtmp2=Vsh * Vquvy;
Vtmp3=Vsh * Vquvz;
Vsh=Vsh * Vqus;
Vqus=Vch * Vqus;
Vquvx=Vch * Vquvx;
Vquvy=Vch * Vquvy;
Vquvz=Vch * Vquvz;
Vquvx=Vquvx + Vsh;
Vqus=Vqus - Vtmp1;
Vquvy=Vquvy + Vtmp3;
Vquvz=Vquvz - Vtmp2;
#endif


#ifdef __INTEL_COMPILER
#pragma warning( default : 592 )
#endif
