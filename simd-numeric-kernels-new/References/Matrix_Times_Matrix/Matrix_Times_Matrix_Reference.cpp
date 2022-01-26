//#####################################################################
//  Copyright (c) 2011-2019 Nathan Mitchell, Eftychios Sifakis, Yutian Tao.
//  This file is covered by the FreeBSD license. Please refer to the
//  license.txt file for more information.
//#####################################################################


//#include <PhysBAM_Tools/Matrices/MATRIX_3X3.h>
#include <iostream>
#include <iomanip>
#include <bitset> // for bitset
#include <type_traits> // for is_same
#include <stdexcept> // for logic error
#include <Eigen/Dense>

using namespace Eigen;

namespace{

    namespace {
        typedef union {
            int i;
            float f;
        } floatConverter;

        typedef union {
            long long int i;
            double f;
        } doubleConverter;
    }

    template<class T>
        void Print_Bdiff(const T C[9], const T C_reference[9], std::ostream& output);

    template<>
        void Print_Bdiff(const float C[9], const float C_reference[9], std::ostream& output)
    {
        floatConverter cvt;
        floatConverter cvt_ref;
        for (int i=0; i<9; i++) {
            cvt.f = C[i];
            cvt_ref.f = C_reference[i];
                output << std::bitset<32>(cvt.i^cvt_ref.i)<<std::endl;
        }
    }

    template<>
        void Print_Bdiff(const double C[9], const double C_reference[9], std::ostream& output)
    {
        doubleConverter cvt;
        doubleConverter cvt_ref;
        for (int i=0; i<9; i++) {
            cvt.f = C[i];
            cvt_ref.f = C_reference[i];
            output << std::bitset<64>(cvt.i^cvt_ref.i)<<std::endl;
        }
    }

template<class T_MATRIX>
    void Print_Formatted(const T_MATRIX& A,std::ostream& output)
{
    for(int i=0;i<A.rows();i++){
        for(int j=0;j<A.cols();j++){
            output<<std::setw(12)<<A(i,j);
            if(j<A.cols()-1) output<<" ";}
        output<<std::endl;}
}
}

template<class T>
void Matrix_Times_Matrix_Reference(const T A[9], const T B[9], T C[9])
{
    Map<const Matrix<T,3,3>> mA=Map<const Matrix<T,3,3>>(A);
    Map<const Matrix<T,3,3>> mB=Map<const Matrix<T,3,3>>(B);
    Map<Matrix<T,3,3>> mC=Map<Matrix<T,3,3>>(C);

    mC=mA*mB;
}

template<class T>
bool Matrix_Times_Matrix_Compare(const T C[9], const T C_reference[9])
{
    Map<const Matrix<T,3,3>> mC=Map<const Matrix<T,3,3>>(C);
    Map<const Matrix<T,3,3>> mC_reference=Map<const Matrix<T,3,3>>(C_reference);

    std::cout<<"Computed matrix C :"<<std::endl;Print_Formatted(mC,std::cout);
    std::cout<<"Reference matrix C :"<<std::endl;Print_Formatted(mC_reference,std::cout);
    std::cout<<"Difference = "<<(mC-mC_reference).norm()<<std::endl;
    Print_Bdiff(C, C_reference, std::cout);
    if( (mC-mC_reference).norm() < 0.00001 )
        return true;
    else
        return false;
}

template void Matrix_Times_Matrix_Reference(const float A[9], const float B[9], float C[9]);
template bool Matrix_Times_Matrix_Compare(const float C[9], const float C_reference[9]);

template void Matrix_Times_Matrix_Reference(const double A[9], const double B[9], double C[9]);
template bool Matrix_Times_Matrix_Compare(const double C[9], const double C_reference[9]);
