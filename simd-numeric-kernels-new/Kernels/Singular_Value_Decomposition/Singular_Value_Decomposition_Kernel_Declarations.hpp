//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################


//###########################################################
// Local variable declarations
//###########################################################

const typename Tarch::Scalar Four_Gamma_Squared=(typename Tarch::Scalar) (sqrt(8.)+3.);
const typename Tarch::Scalar Sine_Pi_Over_Eight= (typename Tarch::Scalar) (.5*sqrt(2.-sqrt(2.)));
const typename Tarch::Scalar Cosine_Pi_Over_Eight= (typename Tarch::Scalar) (.5*sqrt(2.+sqrt(2.)));

#if 0
BUILD_CONSTANT(_VFour_Gamma_Squared,Four_Gamma_Squared);
BUILD_CONSTANT(_VSine_Pi_Over_Eight,Sine_Pi_Over_Eight);
BUILD_CONSTANT(_VCosine_Pi_Over_Eight,Cosine_Pi_Over_Eight);
BUILD_CONSTANT(_VOne_Half,.5f);
BUILD_CONSTANT(_VOne,1.0f);
BUILD_CONSTANT(_VTiny_Number,1.e-20f);
BUILD_CONSTANT(_VSmall_Number,1.e-12f);
BUILD_CONSTANT(_VNegTwo,-2.0f);
#endif

#define DECLARE_CONSTANT(name) typename Tarch::Scalar alignas(sizeof(Tarch::ScalarRegister)) name[Tarch::Width];
DECLARE_CONSTANT(_VFour_Gamma_Squared);
DECLARE_CONSTANT(_VSine_Pi_Over_Eight);
DECLARE_CONSTANT(_VCosine_Pi_Over_Eight);
DECLARE_CONSTANT(_VOne_Half);
DECLARE_CONSTANT(_VOne);
DECLARE_CONSTANT(_VTiny_Number);
DECLARE_CONSTANT(_VSmall_Number);
DECLARE_CONSTANT(_VNegTwo);
#undef DECLARE_CONSTANT

#define INITIALIZE_CONSTANT(name,value) \
        for(int i = 0; i < Tarch::Width; i++) name[i] = value;
INITIALIZE_CONSTANT(_VFour_Gamma_Squared,Four_Gamma_Squared);
INITIALIZE_CONSTANT(_VSine_Pi_Over_Eight,Sine_Pi_Over_Eight);
INITIALIZE_CONSTANT(_VCosine_Pi_Over_Eight,Cosine_Pi_Over_Eight);
INITIALIZE_CONSTANT(_VOne_Half,.5f);
INITIALIZE_CONSTANT(_VOne,1.0f);
INITIALIZE_CONSTANT(_VTiny_Number,1.e-20f);
INITIALIZE_CONSTANT(_VSmall_Number,1.e-12f);
INITIALIZE_CONSTANT(_VNegTwo,-2.0f);
#undef INITIALIZE_CONSTANT

typedef Number<Tarch> Tn;
using Tm = Mask<Tarch>;

Tn Vfour_gamma_squared;
Tn Vsine_pi_over_eight;
Tn Vcosine_pi_over_eight;
Tn Vone_half;
Tn Vone;
Tn Vtiny_number;
Tn Vsmall_number;

Vfour_gamma_squared.Load_Aligned(_VFour_Gamma_Squared);
Vsine_pi_over_eight.Load_Aligned(_VSine_Pi_Over_Eight);
Vcosine_pi_over_eight.Load_Aligned(_VCosine_Pi_Over_Eight);
Vone_half.Load_Aligned(_VOne_Half);
Vone.Load_Aligned(_VOne);
Vtiny_number.Load_Aligned(_VTiny_Number);
Vsmall_number.Load_Aligned(_VSmall_Number);

Tn Va11;
Tn Va21;
Tn Va31;
Tn Va12;
Tn Va22;
Tn Va32;
Tn Va13;
Tn Va23;
Tn Va33;

#ifdef COMPUTE_V_AS_MATRIX
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

#ifdef COMPUTE_V_AS_QUATERNION
Tn Vqvs;
Tn Vqvvx;
Tn Vqvvy;
Tn Vqvvz;
#endif

#ifdef COMPUTE_U_AS_MATRIX
Tn Vu11;
Tn Vu21;
Tn Vu31;
Tn Vu12;
Tn Vu22;
Tn Vu32;
Tn Vu13;
Tn Vu23;
Tn Vu33;
#endif

#ifdef COMPUTE_U_AS_QUATERNION
Tn Vqus;
Tn Vquvx;
Tn Vquvy;
Tn Vquvz;
#endif

Tn Vc;
Tn Vs;
Tn Vch;
Tn Vsh;
Tn Vtmp1;
Tn Vtmp2;
Tn Vtmp3;
Tn Vtmp4;
Tn Vtmp5;
Tm Mtmp1;
Tm Mtmp2;
Tm Mtmp3;
Tm Mtmp4;
Tm Mtmp5;
