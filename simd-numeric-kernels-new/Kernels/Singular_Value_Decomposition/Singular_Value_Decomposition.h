//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################


template<class Tarch,class T_DATA=void>
void Singular_Value_Decomposition(const T_DATA (&A)[9], T_DATA (&U)[9], T_DATA (&S)[3], T_DATA (&V)[9]);
