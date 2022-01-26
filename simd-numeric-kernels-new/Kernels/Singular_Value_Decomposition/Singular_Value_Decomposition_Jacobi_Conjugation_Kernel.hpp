//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


//###########################################################
// Compute the Givens angle (and half-angle)
//###########################################################

Vsh=VS21 * Vone_half;
Vtmp5=VS11 - VS22;

Vtmp2=Vsh * Vsh;
Mtmp1=Vtmp2 >= Vtiny_number;
Vsh = Vsh.mask(Mtmp1);
Vch = blend(Mtmp1,Vone,Vtmp5);

Vtmp1=Vsh * Vsh;
Vtmp2=Vch * Vch;
Vtmp3=Vtmp1 + Vtmp2;
Vtmp4=Vtmp3.rsqrt();

#ifdef USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION
Vs=Vtmp4 * Vone_half;
Vc=Vtmp4 * Vs;
Vc=Vtmp4 * Vc;
Vc=Vtmp3 * Vc;
Vtmp4=Vtmp4 + Vs;
Vtmp4=Vtmp4 - Vc;
#endif

Vsh=Vtmp4 * Vsh;
Vch=Vtmp4 * Vch;

Vtmp1=Vfour_gamma_squared * Vtmp1;
Mtmp1=Vtmp2 <= Vtmp1;
                      
Vsh = blend(Mtmp1,Vsh,Vsine_pi_over_eight);
Vch = blend(Mtmp1, Vch, Vcosine_pi_over_eight);

Vtmp1=Vsh * Vsh;
Vtmp2=Vch * Vch;
Vc=Vtmp2 - Vtmp1;
Vs=Vch * Vsh;
Vs=Vs + Vs;

//###########################################################
// Perform the actual Givens conjugation
//###########################################################

#ifndef USE_ACCURATE_RSQRT_IN_JACOBI_CONJUGATION
Vtmp3=Vtmp1 + Vtmp2;
VS33=VS33 * Vtmp3;
VS31=VS31 * Vtmp3;
VS32=VS32 * Vtmp3;
VS33=VS33 * Vtmp3;
#endif

Vtmp1=Vs * VS31;
Vtmp2=Vs * VS32;
VS31=Vc * VS31;
VS32=Vc * VS32;
VS31=Vtmp2 + VS31;
VS32=VS32 - Vtmp1;

Vtmp2=Vs * Vs;
Vtmp1=VS22 * Vtmp2;
Vtmp3=VS11 * Vtmp2;
Vtmp4=Vc * Vc;
VS11=VS11 * Vtmp4;
VS22=VS22 * Vtmp4;
VS11=VS11 + Vtmp1;
VS22=VS22 + Vtmp3;
Vtmp4=Vtmp4 - Vtmp2;
Vtmp2=VS21 + VS21;
VS21=VS21 * Vtmp4;
Vtmp4=Vc * Vs;
Vtmp2=Vtmp2 * Vtmp4;
Vtmp5=Vtmp5 * Vtmp4;
VS11=VS11 + Vtmp2;
VS21=VS21 - Vtmp5;
VS22=VS22 - Vtmp2;

//###########################################################
// Compute the cumulative rotation, in quaternion form
//###########################################################

Vtmp1=Vsh * Vqvvx;
Vtmp2=Vsh * Vqvvy;
Vtmp3=Vsh * Vqvvz;
Vsh=Vsh * Vqvs;

Vqvs=Vch * Vqvs;
Vqvvx=Vch * Vqvvx;
Vqvvy=Vch * Vqvvy;
Vqvvz=Vch * Vqvvz;

VQVVZ=VQVVZ + Vsh;
Vqvs=Vqvs - VTMP3;
VQVVX=VQVVX + VTMP2;
VQVVY=VQVVY - VTMP1;
