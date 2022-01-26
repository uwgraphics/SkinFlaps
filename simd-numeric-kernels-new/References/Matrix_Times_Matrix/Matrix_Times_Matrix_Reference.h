//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao, Qisi Wang.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################


template<class T>
void Matrix_Times_Matrix_Reference(const T A[9], const T B[9], T C[9]);

template<class T>
bool Matrix_Times_Matrix_Compare(const T C[9], const T C_reference[9]);
