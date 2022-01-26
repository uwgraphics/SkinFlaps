//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


template<class T>
void Singular_Value_Decomposition_Reference(const T A[9], T U[9], T Sigma[3], T V[9]);

template<class T>
bool Singular_Value_Decomposition_Compare(const T U[9], const T Sigma[3], const T V[9],
                                          const T U_reference[9], const T Sigma_reference[3], const T V_reference[9] );
