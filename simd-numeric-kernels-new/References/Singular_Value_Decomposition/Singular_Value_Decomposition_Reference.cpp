//#####################################################################
//  Copyright (c) 2011-2013 Nathan Mitchell, Eftychios Sifakis.
//  This file is covered by the FreeBSD license. Please refer to the 
//  license.txt file for more information.
//#####################################################################


#include <iomanip>
#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

namespace{
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
void Singular_Value_Decomposition_Reference(const T A[9], T U[9], T Sigma[3], T V[9])
{
    Map<const Matrix<T,3,3>> mA=Map<const Matrix<T,3,3>>(A);

    Map<Matrix<T,3,3>> mU=Map<Matrix<T,3,3>>(U);
    Map<Matrix<T,3,1>> mSigma=Map<Matrix<T,3,1>>(Sigma);
    Map<Matrix<T,3,3>> mV=Map<Matrix<T,3,3>>(V);

    JacobiSVD<Matrix<T,3,3>> svd(mA, ComputeFullU|ComputeFullV);
    mU=svd.matrixU();
    mSigma=svd.singularValues();
    mV=svd.matrixV();

    if(mU.determinant() < 0.) {
        mU.col(2) *= -1.;
        mSigma(2) *= -1.;
    }

    if(mV.determinant() < 0.) {
        mV.col(2) *= -1.;
        mSigma(2) *= -1.;
    }
}

template<class T>
bool Singular_Value_Decomposition_Compare(const T U[9], const T Sigma[3], const T V[9],
                                          const T U_reference[9], const T Sigma_reference[3], const T V_reference[9])
{

    Map<const Matrix<T,3,3>> mU=Map<const Matrix<T,3,3>>(U);
    Map<const Matrix<T,3,1>> mSigma=Map<const Matrix<T,3,1>>(Sigma);
    Map<const Matrix<T,3,3>> mV=Map<const Matrix<T,3,3>>(V);

    Map<const Matrix<T,3,3>> mU_reference=Map<const Matrix<T,3,3>>(U_reference);
    Map<const Matrix<T,3,1>> mSigma_reference=Map<const Matrix<T,3,1>>(Sigma_reference);
    Map<const Matrix<T,3,3>> mV_reference=Map<const Matrix<T,3,3>>(V_reference);

    Matrix<T,3,3> mU_adjusted = mU_reference;
    Matrix<T,3,1> mSigma_adjusted = mSigma_reference;
    Matrix<T,3,3> mV_adjusted = mV_reference;

    for (int i = 0; i < 2; i++) {

      if (mU.col(i).dot(mU_adjusted.col(i)) < 0.) {
          mU_adjusted.col(i) *= -1.;
          mSigma_adjusted(i) *= -1.;
          mU_adjusted.col(2) *= -1.;
          mSigma_adjusted(2) *= -1.;
      }
	
      if (mV.col(i).dot(mV_adjusted.col(i)) < 0.) {
          mV_adjusted.col(i) *= -1.;
          mSigma_adjusted(i) *= -1.;
          mV_adjusted.col(2) *= -1.;
          mSigma_adjusted(2) *= -1.;
      }
	
    }
    
    std::cout<<"Computed matrix U :"<<std::endl;Print_Formatted(mU,std::cout);
    std::cout<<"Reference matrix U :"<<std::endl;Print_Formatted(mU_reference,std::cout);
    std::cout<<"Adjusted matrix U :"<<std::endl;Print_Formatted(mU_adjusted,std::cout);
    std::cout<<"Difference = "<< (mU-mU_adjusted).norm()  <<std::endl;

    std::cout<<std::endl;
    std::cout<<"Computed matrix Sigma :"<<std::endl;Print_Formatted(mSigma,std::cout);
    std::cout<<"Reference matrix Sigma :"<<std::endl;Print_Formatted(mSigma_reference,std::cout);
    std::cout<<"Adjusted matrix Sigma :"<<std::endl;Print_Formatted(mSigma_adjusted,std::cout);
    std::cout<<"Difference = "<<  (mSigma-mSigma_adjusted).norm()  <<std::endl;
    
    std::cout<<std::endl;
    std::cout<<"Computed matrix V :"<<std::endl;Print_Formatted(mV,std::cout); 
    std::cout<<"Reference matrix V :"<<std::endl;Print_Formatted(mV_reference,std::cout);
    std::cout<<"Adjusted matrix V :"<<std::endl;Print_Formatted(mV_adjusted,std::cout);
    std::cout<<"Difference = "<< (mV-mV_adjusted).norm()  <<std::endl;   

    Matrix<T,3,3> mA=mU*mSigma.asDiagonal()*mV.transpose();
    Matrix<T,3,3> mA_reference=mU_reference*mSigma_reference.asDiagonal()*mV_reference.transpose();

    std::cout<<std::endl;
    std::cout<<"Computed matrix A=U*Sigma*V^T :"<<std::endl;Print_Formatted(mA,std::cout);
    std::cout<<"Reference matrix A=U*Sigma*V^T :"<<std::endl;Print_Formatted(mA_reference,std::cout);
    std::cout<<"Difference = "<< (mA-mA_reference).norm()  <<std::endl;


    if(!((mU-mU_adjusted).norm() < 0.00001) )
        return false;
    

    if(!((mSigma-mSigma_adjusted).norm() < 0.00001) )
        return false;
    

    if(!((mV-mV_adjusted).norm() < 0.00001) )
        return false;
    

    






    return true;
}

template void Singular_Value_Decomposition_Reference(const float A[9], float U[9], float Sigma[3], float V[9]);
template bool Singular_Value_Decomposition_Compare(const float U[9], const float Sigma[3], const float V[9],
                                                   const float U_reference[9], const float Sigma_reference[3], const float V_reference[9]);
 
