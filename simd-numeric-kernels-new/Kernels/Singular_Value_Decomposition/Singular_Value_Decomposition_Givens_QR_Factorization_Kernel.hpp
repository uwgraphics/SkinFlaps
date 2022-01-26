//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


//###########################################################
// Compute the Givens half-angle, construct the Givens quaternion and the rotation sine/cosine (for the full angle)
//###########################################################

    Vsh=VANPIVOT * VANPIVOT;
    Mtmp1=Vsh >= Vsmall_number;
    Vsh=VANPIVOT.mask(Mtmp1);

    Vtmp5=Vtmp5 ^ Vtmp5;
    Vch=Vtmp5 - VAPIVOT;
    Vch=max(Vch,VAPIVOT);
    Vch=max(Vch,Vsmall_number);
    Mtmp5=VAPIVOT >= Vtmp5;

    Vtmp1=Vch * Vch;
    Vtmp2=Vsh * Vsh;
    Vtmp2=Vtmp1 + Vtmp2;
    Vtmp1=Vtmp2.rsqrt();

    Vtmp4=Vtmp1 * Vone_half;
    Vtmp3=Vtmp1 * Vtmp4;
    Vtmp3=Vtmp1 * Vtmp3;
    Vtmp3=Vtmp2 * Vtmp3;
    Vtmp1=Vtmp1 + Vtmp4;
    Vtmp1=Vtmp1 - Vtmp3;
    Vtmp1=Vtmp1 * Vtmp2;

    Vch=Vch + Vtmp1;

    Vtmp1 = Vsh;
    Vtmp5 = Vch;
    Vch = blend(Mtmp5,Vtmp1,Vtmp5);
    Vsh = blend(Mtmp5,Vtmp5,Vtmp1);

    Vtmp1=Vch * Vch;
    Vtmp2=Vsh * Vsh;
    Vtmp2=Vtmp1 + Vtmp2;
    Vtmp1=Vtmp2.rsqrt();

    Vtmp4=Vtmp1 * Vone_half;
    Vtmp3=Vtmp1 * Vtmp4;
    Vtmp3=Vtmp1 * Vtmp3;
    Vtmp3=Vtmp2 * Vtmp3;
    Vtmp1=Vtmp1 + Vtmp4;
    Vtmp1=Vtmp1 - Vtmp3;

    Vch=Vch * Vtmp1;
    Vsh=Vsh * Vtmp1;

    Vc=Vch * Vch;
    Vs=Vsh * Vsh;
    Vc=Vc - Vs;
    Vs=Vsh * Vch;
    Vs=Vs + Vs;

//###########################################################
// Rotate matrix A
//###########################################################

    Vtmp1=Vs * VA11;
    Vtmp2=Vs * VA21;
    VA11=Vc * VA11;
    VA21=Vc * VA21;
    VA11=VA11 + Vtmp2;
    VA21=VA21 - Vtmp1;

    Vtmp1=Vs * VA12;
    Vtmp2=Vs * VA22;
    VA12=Vc * VA12;
    VA22=Vc * VA22;
    VA12=VA12 + Vtmp2;
    VA22=VA22 - Vtmp1;

    Vtmp1=Vs * VA13;
    Vtmp2=Vs * VA23;
    VA13=Vc * VA13;
    VA23=Vc * VA23;
    VA13=VA13 + Vtmp2;
    VA23=VA23 - Vtmp1;

//###########################################################
// Update matrix U
//###########################################################

#ifdef COMPUTE_U_AS_MATRIX
    Vtmp1=Vs * VU11;
    Vtmp2=Vs * VU12;
    VU11=Vc * VU11;
    VU12=Vc * VU12;
    VU11=VU11 + Vtmp2;
    VU12=VU12 - Vtmp1;

    Vtmp1=Vs * VU21;
    Vtmp2=Vs * VU22;
    VU21=Vc * VU21;
    VU22=Vc * VU22;
    VU21=VU21 + Vtmp2;
    VU22=VU22 - Vtmp1;

    Vtmp1=Vs * VU31;
    Vtmp2=Vs * VU32;
    VU31=Vc * VU31;
    VU32=Vc * VU32;
    VU31=VU31 + Vtmp2;
    VU32=VU32 - Vtmp1;
#endif
