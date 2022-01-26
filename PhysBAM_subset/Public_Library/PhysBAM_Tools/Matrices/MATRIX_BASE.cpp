//#####################################################################
// Copyright 2002-2013, Ronald Fedkiw, Geoffrey Irving, Igor Neverov, Craig Schroeder, Rahul Sheth, Tamar Shinar, Eftychios Sifakis, Huamin Wang, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_BASE
//#####################################################################
#include <PhysBAM_Tools/Matrices/MATRIX.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Matrices/MATRIX_MXN.h>
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Tools/Math_Tools/robust_givens_rotation.h>
#include <limits>
using namespace PhysBAM;
//#####################################################################
// Function Gram_Schmidt_QR_Factorization
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_Gram_Schmidt_QR_Factorization(MATRIX_BASE<T,T_MATRIX2>& R)
{
    R.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=1;j<=Columns();j++){ // for each column
        for(int i=1;i<=Rows();i++) R(j,j)+=sqr((*this)(i,j));R(j,j)=sqrt(R(j,j)); // compute the L2 norm
        T one_over_Rjj=1/R(j,j);
        for(int i=1;i<=Rows();i++) (*this)(i,j)*=one_over_Rjj; // orthogonalize the column
        for(int k=j+1;k<=Columns();k++){ // subtract this columns contributution from the rest of the columns
            for(int i=1;i<=Rows();i++) R(j,k)+=(*this)(i,j)*(*this)(i,k);
            for(int i=1;i<=Rows();i++) (*this)(i,k)-=R(j,k)*(*this)(i,j);}}
}
//#####################################################################
// Function Householder_QR_Factorization
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2,class T_MATRIX3> void MATRIX_BASE<T,T_MATRIX>::
Householder_QR_Factorization(MATRIX_BASE<T,T_MATRIX2>& V,MATRIX_BASE<T,T_MATRIX3>& R)
{
    V.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());R.Derived()=T_MATRIX3(Columns(),Columns());T_MATRIX temp(*this);LEFT_VECTOR a(Rows()),v,new_a;
    for(int j=1;j<=Columns();j++){ // for each column
        for(int i=1;i<=Rows();i++) a(i)=temp(i,j);
        v=a.Householder_Vector(j);v.Normalize();for(int i=1;i<=Rows();i++) V(i,j)=v(i); // store the v's in V
        for(int k=j;k<=Columns();k++){ // Householder transform each column
            for(int i=1;i<=Rows();i++) a(i)=temp(i,k);
            new_a=a.Householder_Transform(v);for(int i=1;i<=Rows();i++) temp(i,k)=new_a(i);}}
    for(int i=1;i<=Columns();i++) for(int j=1;j<=Columns();j++) R(i,j)=temp(i,j); // store R
}
//#####################################################################
// Function In_Place_Robust_Householder_QR_Solve
//#####################################################################
template<class T,class T_MATRIX> template<class T_VECTOR1,class T_VECTOR2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_Robust_Householder_QR_Solve(VECTOR_BASE<T,T_VECTOR1>& b,VECTOR_BASE<int,T_VECTOR2>& p)
{
    assert(Rows()==b.Size() && Columns()==p.Size());
    VECTOR_ND<T> a((INITIAL_SIZE)Rows());for(int i=1;i<=Columns();i++) p(i)=i; // TODO: This should not assume VECTOR_ND.
    VECTOR_ND<T> column_norm(Columns());for(int j=1;j<=Columns();j++) for(int i=1;i<=Rows();i++) column_norm(j)+=sqr((*this)(i,j));
    for(int j=1;j<=Columns();j++){
        int max_column=0;T max_column_norm=0;for(int k=j;k<=Columns();k++) if(column_norm(k)>max_column_norm){max_column_norm=column_norm(k);max_column=k;}
        if(max_column_norm<FLT_MIN) return;
        if(max_column!=j){exchange(column_norm(j),column_norm(max_column));exchange(p(j),p(max_column));for(int i=1;i<=Rows();i++) exchange((*this)(i,j),(*this)(i,max_column));}
        if(j==Rows()) return;
        Get_Column(j,a);VECTOR_ND<T> v=a.Householder_Vector(j);T two_over_v_dot_v=(T)2/v.Magnitude_Squared();
        if((*this)(j,j)>=0)(*this)(j,j)=-sqrt(max_column_norm);else (*this)(j,j)=sqrt(max_column_norm);for(int i=j+1;i<=Rows();i++)(*this)(i,j)=(T)0;
        for(int k=j+1;k<=Columns();k++){
            T v_dot_a=0;for(int i=j;i<=Rows();i++) v_dot_a+=v(i)*(*this)(i,k);T coefficient=v_dot_a*two_over_v_dot_v;for(int i=j;i<=Rows();i++) (*this)(i,k)-=coefficient*v(i);}
        T v_dot_b=0;for(int i=j;i<=Rows();i++) v_dot_b+=v(i)*b(i);T coefficient=v_dot_b*two_over_v_dot_v;for(int i=j;i<=Rows();i++) b(i)-=coefficient*v(i);
        for(int k=j+1;k<=Columns();k++) column_norm(k)-=sqr((*this)(j,k));}
}
//#####################################################################
// Function In_Place_Householder_QR_With_Column_Pivoting
// Input: Matrix R (this=input A, any size), Q (empty, correct size), Matrix P (empty, correct size)
// Output: Matrix R (this, upper triangular, destroys A), Q (rotation), P (permutation) such that A*P=Q*R, or A=Q*R*P^T
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX1,class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_Householder_QR_Factorization_With_Column_Pivoting(MATRIX_BASE<T,T_MATRIX1>& Q,MATRIX_BASE<T,T_MATRIX2>& P,const T tolerance) //defaults: tolerance=(T)1e-4
{
    MATRIX_BASE<T,T_MATRIX>& R=(*this);
    //TODO: Make more efficient by producing the economy size decomposition
    assert(R.Rows()==Q.Columns() && R.Rows()==Q.Rows());
    Q.Set_Identity_Matrix();P.Set_Zero_Matrix();
    VECTOR_ND<T> norms,a,v;VECTOR_ND<int> permutation;MATRIX_MXN<T> outer_product,householder_matrix;
    norms.Resize(R.Columns());permutation.Resize(R.Columns()); //also sets to zero
    T max;int max_index;
    householder_matrix.Resize(R.Rows(),R.Rows());
    for(int j=1;j<=R.Columns();j++){permutation(j)=j;for(int i=1;i<=R.Rows();i++){norms(j)+=R(i,j)*R(i,j);}}
    for(int k=1;k<=R.Columns();k++){
        if(k>R.Rows()){break;}
        if(k!=1){for(int i=k;i<=R.Columns();i++){norms(i)-=R(k-1,i)*R(k-1,i);}}
        max=norms(k);max_index=k;for(int i=k+1;i<=norms.Size();i++){if(norms(i)>max){max=norms(i);max_index=i;}}
        if(max<tolerance){break;}
        if(max_index!=k){for(int i=1;i<=R.Rows();i++){exchange(R(i,max_index),R(i,k));}
            exchange(norms(k),norms(max_index));exchange(permutation(k),permutation(max_index));}
        a.Resize(R.Rows()-k+1);v.Resize(R.Rows()-k+1);
        for(int i=k;i<=R.Rows();i++){v(i-k+1)=a(i-k+1)=R(i,k);}
        v(1)+=sign_nonzero(v(1))*sqrt(max);
        for(int j=k;j<=R.Columns();j++){
            if(j!=k){for(int i=k;i<=R.Rows();i++){a(i-k+1)=R(i,j);}}
            a=a.Householder_Transform(v);
            for(int i=k;i<=R.Rows();i++){R(i,j)=a(i-k+1);}}
        outer_product.Resize(v.Size(),v.Size());
        outer_product.Set_Identity_Matrix();
        outer_product-=(2/v.Magnitude_Squared())*MATRIX_MXN<T>::Outer_Product(v,v);
        householder_matrix.Set_Identity_Matrix();
        householder_matrix.Set_Submatrix(k,k,outer_product);
        Q.Derived()=Q.Derived()*householder_matrix;}
    for(int i=1;i<=R.Columns();i++){P(permutation(i),i)=(T)1;}
}
//#####################################################################
// Function In_Place_Householder_QR_With_Full_Pivoting
// Input: Matrix R (this=input A, any size), Q (empty, correct size), Matrix P (empty, correct size)
// Output: Matrix R (this, upper triangular, destroys A), Q (rotation), P (permutation) such that A*P=Q*R, or A=Q*R*P^T
// Slower than with column pivoting, but more robust
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX1,class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_Householder_QR_Factorization_With_Full_Pivoting(MATRIX_BASE<T,T_MATRIX1>& Q,MATRIX_BASE<T,T_MATRIX2>& P,const T tolerance) //defaults: tolerance=(T)1e-4
{
    MATRIX_BASE<T,T_MATRIX>& R=(*this);
    //TODO: Make more efficient by producing the economy size decomposition
    assert(R.Rows()==Q.Columns() && R.Rows()==Q.Rows());
    Q.Set_Identity_Matrix();P.Set_Zero_Matrix();
    VECTOR_ND<T> a,v;VECTOR_ND<int> row_permutation,col_permutation;MATRIX_MXN<T> outer_product,householder_matrix;
    row_permutation.Resize(R.Rows());col_permutation.Resize(R.Columns()); //also sets to zero
    T max,col_norm;int max_row_index,max_col_index;
    householder_matrix.Resize(R.Rows(),R.Rows());
    for(int i=1;i<=R.Rows();i++){row_permutation(i)=i;}for(int j=1;j<=R.Columns();j++){col_permutation(j)=j;}
    for(int k=1;k<=R.Columns();k++){
        if(k>R.Rows()){break;}
        max=abs(R(k,k));max_row_index=max_col_index=k;
        for(int i=k+1;i<=R.Rows();i++){for(int j=k+1;j<=R.Columns();j++){
            if(abs(R(i,j))>max){max=abs(R(i,j));max_row_index=i;max_col_index=j;}}}
        if(max<tolerance){break;}
        if(max_col_index!=k){for(int i=1;i<=R.Rows();i++){exchange(R(i,max_col_index),R(i,k));}
            exchange(col_permutation(k),col_permutation(max_col_index));}
        if(max_row_index!=k){for(int j=1;j<=R.Columns();j++){exchange(R(max_row_index,j),R(k,j));}
            for(int i=1;i<=Q.Rows();i++){exchange(Q(i,max_row_index),Q(i,k));}}
        a.Resize(R.Rows()-k+1);v.Resize(R.Rows()-k+1);
        col_norm=(T)0;
        for(int i=k;i<=R.Rows();i++){v(i-k+1)=a(i-k+1)=R(i,k);col_norm+=R(i,k)*R(i,k);}
        v(1)+=sign_nonzero(v(1))*sqrt(col_norm);
        for(int j=k;j<=R.Columns();j++){
            if(j!=k){for(int i=k;i<=R.Rows();i++){a(i-k+1)=R(i,j);}}
            a=a.Householder_Transform(v);
            for(int i=k;i<=R.Rows();i++){R(i,j)=a(i-k+1);}}
        outer_product.Resize(v.Size(),v.Size());
        outer_product.Set_Identity_Matrix();
        outer_product-=(2/v.Magnitude_Squared())*MATRIX_MXN<T>::Outer_Product(v,v);
        householder_matrix.Set_Identity_Matrix();
        householder_matrix.Set_Submatrix(k,k,outer_product);
        Q.Derived()=Q.Derived()*householder_matrix;}
    for(int j=1;j<=R.Columns();j++){P(col_permutation(j),j)=(T)1;}
}
//#####################################################################
// Function Robust_Householder_QR_Solve
// Input: Matrix R (this=input A, any size), Vector b (length of A.Rows()), Vector x (length of A.Columns())
// Output: Vector x, such that Ax=b
// Handles regular and rank deficient (minimum-norm) least-squares problems
//#####################################################################
template<class T,class T_MATRIX> template<class T_VECTOR1,class T_VECTOR2> void MATRIX_BASE<T,T_MATRIX>::
Robust_Householder_QR_Solve(VECTOR_BASE<T,T_VECTOR1>& x, const VECTOR_BASE<T,T_VECTOR2>& b,const bool use_full_pivoting,const T tolerance) //defaults: tolerance=(T)1e-4
{
    const MATRIX_BASE<T,T_MATRIX>& A=(*this);
    assert(A.Rows()==b.Size() && A.Columns()==x.Size());
    MATRIX_MXN<T> Q,R,P;R=A;
    Q.Resize(A.Rows(),A.Rows());P.Resize(A.Columns(),A.Columns());
    if(use_full_pivoting){R.In_Place_Householder_QR_Factorization_With_Full_Pivoting(Q,P);}
    else{R.In_Place_Householder_QR_Factorization_With_Column_Pivoting(Q,P);}
    int row_count=PhysBAM::min(R.Rows(),R.Columns());
    for(int i=row_count;i>=1;i--){if(abs(R(i,i))>=tolerance){break;}row_count--;}
    //full rank
    if(row_count==PhysBAM::min(R.Rows(),R.Columns())){x.Derived()=P*(R.Upper_Triangular_Solve(Q.Transpose_Times(b)));return;}
    //rank deficient
    else{//TODO: Replace Copies with MATRIX_VIEW for speed
        MATRIX_MXN<T> Q1,R1,P2,Q2,R2;
        R1.Resize(row_count,R.Columns());
        for(int j=1;j<=row_count;j++){for(int i=1;i<=R.Columns();i++){R1(j,i)=R(j,i);}}
        R1.Transpose();R2=R1;Q2.Resize(R1.Rows(),R1.Rows());P2.Resize(R1.Columns(),R1.Columns());
        //TODO: Probably don't need pivoting here
        R2.In_Place_Householder_QR_Factorization_With_Column_Pivoting(Q2,P2);
        VECTOR_ND<T> b2=Q.Transpose_Times(b);
        R1.Resize(row_count,row_count);Q1.Resize(row_count,Q2.Rows());
        VECTOR_ND<T> b_short;b_short.Resize(row_count);
        for(int i=1;i<=row_count;i++){b_short(i)=b2(i);
            for(int j=1;j<=row_count;j++){R1(i,j)=R2(j,i);}
            for(int j=1;j<=Q2.Rows();j++){Q1(i,j)=Q2(j,i);}}
        x.Derived()=P*Q1.Transposed()*(R1.Lower_Triangular_Solve(P2.Transposed()*b_short));}
}
//#####################################################################
// Function In_Place_PLU_Factorization
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_PLU_Factorization(MATRIX_BASE<T,T_MATRIX2>& L,COLUMN_PERMUTATION& p)
{
    assert((INITIAL_SIZE)Rows()==(INITIAL_SIZE)Columns());
    L.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    p=COLUMN_PERMUTATION(INITIAL_SIZE(Columns()));for(int i=1;i<=Columns();i++) p(i)=i; // initialize p
    for(int j=1;j<=Columns();j++){ // for each column
        // find the largest element and switch rows
        int row=j;T value=abs((*this)(j,j));
        for(int i=j+1;i<=Columns();i++) if(abs((*this)(i,j))>value){row=i;value=abs((*this)(i,j));}
        if(row!=j){ // need to switch rows
            exchange(p(j),p(row)); // update permutation matrix
            for(int k=1;k<=j-1;k++) exchange(L(j,k),L(row,k)); // update L
            for(int k=j;k<=Columns();k++) exchange((*this)(j,k),(*this)(row,k));} // update U
        // standard LU factorization steps
        T diagonal_inverse=1/(*this)(j,j);for(int i=j;i<=Columns();i++) L(i,j)=(*this)(i,j)*diagonal_inverse; // fill in the column for L
        for(int i=j+1;i<=Columns();i++) for(int k=j;k<=Columns();k++) (*this)(i,k)-=L(i,j)*(*this)(j,k);} // sweep across each row below row j  TODO: can order be changed?
}
//#####################################################################
// Function In_Place_Cholesky_Factorization
//#####################################################################
template<class T,class T_MATRIX> void MATRIX_BASE<T,T_MATRIX>::
In_Place_Cholesky_Factorization()
{
    assert(Rows()==Columns());
    for(int j=1;j<=Columns();j++){ // for each column
        for(int k=1;k<=j-1;k++) for(int i=j;i<=Rows();i++) (*this)(i,j)-=(*this)(j,k)*(*this)(i,k); // subtract off the known stuff in previous columns
        (*this)(j,j)=sqrt((*this)(j,j));T diagonal_inverse=1/(*this)(j,j);for(int i=j+1;i<=Columns();i++) (*this)(i,j)*=diagonal_inverse;} // update L
    for(int i=1;i<=Rows();i++) for(int j=i+1;j<=Columns();j++) (*this)(i,j)=0; // zero out upper triangular part  TODO: Loop the other way around
}
//#####################################################################
// Function In_Place_Cholesky_Factorization
//#####################################################################
template<class T,class T_MATRIX> template<class T_MATRIX2> void MATRIX_BASE<T,T_MATRIX>::
In_Place_LU_Factorization(MATRIX_BASE<T,T_MATRIX2>& L)
{
    assert(Rows()==Columns());
    L.Derived()=T_MATRIX2((INITIAL_SIZE)Rows(),(INITIAL_SIZE)Columns());
    for(int j=1;j<=Columns();j++){ // for each column
        T diagonal_inverse=1/(*this)(j,j);for(int i=j;i<=Columns();i++) L(i,j)=(*this)(i,j)*diagonal_inverse; // fill in the column for L
        for(int i=j+1;i<=Columns();i++) for(int k=j;k<=Columns();k++) (*this)(i,k)-=L(i,j)*(*this)(j,k);} // sweep across each row below row j  TODO: can the order of these loops be swapped?
}
//####################################################################################
// Function Number_Of_Nonzero_Rows
//####################################################################################
template<class T,class T_MATRIX> int MATRIX_BASE<T,T_MATRIX>::
Number_Of_Nonzero_Rows(const T threshold) const
{
    T threshold_squared=sqr(threshold);int nonzero_rows=0;
    for(int i=1;i<=Rows();i++){
        T row_norm_squared=0;for(int j=1;j<=Columns();j++) row_norm_squared+=sqr((*this)(i,j));
        if(row_norm_squared>threshold_squared) nonzero_rows++;}
    return nonzero_rows;
}
//####################################################################################
// Function Get_Singular_Values
//####################################################################################
template<class T, class T_MATRIX> template<class T_VECTOR1> void MATRIX_BASE<T,T_MATRIX>::
Get_Singular_Values(VECTOR_BASE<T,T_VECTOR1>& v,const T tolerance,const unsigned max_iterations) //default tolerance and max_iterations in header/definition
{
    int m,n;
    m=Rows();n=Columns();
    //sanity checks on incoming vector
    if(IS_SAME<T_VECTOR1,VECTOR_ND<T> >::value){v.Derived().Resize(min(Rows(),Columns()));}
    else{assert(v.Size()==min(Rows(),Columns()));}
    if(m<n) PHYSBAM_FATAL_ERROR("Get_Singular_Values: No support for wide (n>m) matrices!\n");
    //bidiagonalize first
    MATRIX_MXN<T> U,E,V;
    U.Resize(Rows(),Rows());
    V.Resize(Columns(),Columns());
    E.Resize(Rows(),Columns());
    Upper_Bidiagonalization(U,E,V);
    VECTOR_ND<T> d,e,lambda,mu;
    E.Get_Diagonal(d);E.Get_Diagonal(e,1);
    unsigned d_len = d.Size();
    int i,j,iUp,iLow;
    unsigned iter;
    lambda.Resize(d_len);lambda.Set_Zero();
    mu.Resize(d_len);mu.Set_Zero();
    T sigma_low = (T)0, threshold = (T)0;
    
    //estimate the first singular value
    lambda(d_len)=abs(d(d_len));
    for(j=d_len-1;j>=1;j--){lambda(j)=abs(d(j))*lambda(j+1)/(lambda(j+1)+abs(e(j)));}
    mu(1)=abs(d(1));
    for(j=1;j<=n-1;j++){mu(j+1)=abs(d(j+1))*mu(j)/(mu(j)+abs(e(j)));}
    sigma_low=min(lambda.Min(),mu.Min());
    threshold=max(tolerance*sigma_low,max_iterations*std::numeric_limits<T>::min());
    
    //use vector-type sweep - better than matrix-type if only singular values are needed
    T c_old,s_old,c,s,r,h;
    iUp=d_len-1;iLow=1;s_old=0;
    for(iter=1;iter<=max_iterations;iter++){
        for(i=iUp;i>=1;i--){iUp=i;
            if(abs(e(i))>threshold) break;}
        j=iUp;
        for(i=iLow;i<=iUp;i++){
            if(abs(e(i))>threshold){j=i;break;}}
        iLow=j;
        if((iUp==iLow && abs(e(iUp))<=threshold) || (iUp<iLow)){
            for(i=1;i<=d.Size();i++){v(i)=abs(d(i));}
            v.Derived().Sort_Descending();return;}
        c_old=(T)1;c=(T)1;
        d_len=iUp+1-iLow+1;
        for(i=iLow;i<=iUp;i++){
            robust_givens_rotation(c,s,r,c*d(i),e(i));
            if(i!=1) e(i-1)=r*s_old;
            robust_givens_rotation(c_old,s_old,d(i),c_old*r,d(i+1)*s);}
        h=c*d(d_len);
        e(d_len-1)=h*s_old;
        d(d_len)=h*c_old;
    }
    PHYSBAM_FATAL_ERROR("Get_Singular_Values: Failed, too many iterations! Exiting.\n");
}

//####################################################################################
// Function Singular_Value_Decomposition
//####################################################################################
template<class T, class T_MATRIX> template<class T_MATRIX2, class T_MATRIX3, class T_MATRIX4> void MATRIX_BASE<T,T_MATRIX>::
Singular_Value_Decomposition(MATRIX_BASE<T,T_MATRIX2>& U, MATRIX_BASE<T,T_MATRIX3>& E, MATRIX_BASE<T,T_MATRIX4>& V,
const T tolerance, const unsigned max_iterations) //default tolerance and max_iterations in header/definition
{
    int m,n;
    m=Rows();n=Columns();
    //sanity checks on incoming matrices
    assert(U.Rows()==U.Columns());assert(V.Rows()==V.Columns());
    assert(U.Columns()==m);assert(V.Rows()==n);
    assert(E.Rows()==m);assert(E.Columns()==n);
    
    if(m<n){//need to process the transpose instead
        E.Derived().Transpose();
        Derived().Transposed().Singular_Value_Decomposition(V,E,U);
        MATRIX_MXN<T> V_old;V_old.Resize(V.Rows(),V.Columns());
        V_old=V;int i,j;
        for(i=1;i<=m;i++){for(j=1;j<=n;j++){V.Derived()(j,i)=V_old(j,m-i+1);}}
        for(i=m+1;i<=n;i++){for(j=1;j<=n;j++){V.Derived()(j,i)=V_old(j,i);}}
        MATRIX_MXN<T> E_sub;E_sub.Resize(m,m);
        for(i=1;i<=m;i++){for(j=1;j<=m;j++){E_sub(j,i)=E(j,i);}}
        E_sub.Transpose();E_sub.Flip_Vertical();E_sub.Flip_Horizontal();
        E.Set_Zero_Matrix();E.Derived().Transpose();
        for(i=1;i<=m;i++){for(j=1;j<=m;j++){E(j,i)=E_sub(j,i);}}
        U.Flip_Horizontal();
        //reverse the order of singular values for this case
        VECTOR_ND<T> su,sv;su.Resize(m);sv.Resize(n);
        T s;int dim=min(E.Rows(),E.Columns());
        for(i=1;i<=dim/2;i++){
            s=E(i,i);E(i,i)=E(dim-i+1,dim-i+1);E(dim-i+1,dim-i+1)=s;
            for(j=1;j<=m;j++){su(j)=U(j,i);U(j,i)=U(j,dim-i+1);U(j,dim-i+1)=su(j);}
            for(j=1;j<=n;j++){sv(j)=V(j,i);V(j,i)=V(j,dim-i+1);V(j,dim-i+1)=sv(j);}}
        return;}

    //bidiagonalize first
    Upper_Bidiagonalization(U,E,V);
    V.Derived().Transpose();
    T_MATRIX2 U2;T_MATRIX4 V2;
    if(IS_SAME<T_MATRIX2,MATRIX_MXN<T> >::value) U2.Resize(U.Rows(),U.Columns());
    if(IS_SAME<T_MATRIX4,MATRIX_MXN<T> >::value) V2.Resize(V.Rows(),V.Columns());
    VECTOR_ND<T> d,e,lambda,mu;
    E.Get_Diagonal(d);E.Get_Diagonal(e,1);
    int d_len = d.Size();
    int i,j,iUp,iLow;
    unsigned iter;
    lambda.Resize(d_len);lambda.Set_Zero();
    mu.Resize(d_len);mu.Set_Zero();
    T sigma_low = (T)0, threshold = (T)0;
    
    //estimate the first singular value
    lambda(d_len)=abs(d(d_len));
    for(j=d_len-1;j>=1;j--){lambda(j)=abs(d(j))*lambda(j+1)/(lambda(j+1)+abs(e(j)));}
    mu(1)=abs(d(1));
    for(j=1;j<=n-1;j++){mu(j+1)=abs(d(j+1))*mu(j)/(mu(j)+abs(e(j)));}
    sigma_low=min(lambda.Min(),mu.Min());
    threshold=max(tolerance*sigma_low,max_iterations*std::numeric_limits<T>::min());
    
    T c,s,r;
    iUp=d_len-1;iLow=1;
    for(iter=1;iter<=max_iterations;iter++){
        for(i=iUp;i>=1;i--){iUp=i;
            if(abs(E(i,i+1))>threshold) break;}
        j=iUp;
        for(i=iLow;i<=iUp;i++){
            if(abs(E(i,i+1))>threshold){j=i;break;}}
        iLow=j;
        if((iUp==iLow && abs(E(iUp,iUp+1))<=threshold) || (iUp<iLow)){
            //sort the singular values, move the columns/rows accordingly
            //using shell sort
            int a=0,b=0,c=0,s=0,inc=1;T sw=0;VECTOR_ND<T> su,sv;
            su.Resize(m);sv.Resize(n);
            do {inc*=3;inc++;} while(inc<=n);
            do {inc/=3;for(a=inc;a<n;a++){
                sw=E(a+1,a+1);
                for(c=0;c<m;c++) su(c+1)=U(c+1,a+1);
                for(c=0;c<n;c++) sv(c+1)=V(c+1,a+1);
                b=a;
                while(E(b-inc+1,b-inc+1)<sw){
                    U(b+1,b+1)=U(b-inc+1,b-inc+1);
                    for(c=0;c<m;c++) U(c+1,b+1)=U(c+1,b-inc+1);
                    for(c=0;c<n;c++) V(c+1,b+1)=V(c+1,b-inc+1);
                    b-=inc;if(b<inc) break;}
                U(b+1,b+1)=sw;
                for(c=0;c<m;c++) U(c+1,b+1)=su(c+1);
                for(c=0;c<n;c++) V(c+1,b+1)=sv(c+1);}
            }while(inc>1);
            //remove all the negatives on the singular values
            for(c=0;c<n;c++){
                s=0;
                for(a=0;a<m;a++) if(U(a+1,c+1)<0) s++;
                for(c=0;c<n;c++) if(V(b+1,c+1)<0) s++;
                if(s>(m+n)/2){
                    for(a=0;a<m;a++) U(a+1,c+1)=-U(a+1,c+1);
                    for(c=0;c<n;c++) V(b+1,c+1)=-V(b+1,c+1);}}
            V.Derived().Transpose();return;}
        
        //start zero-shift QR downward sweep - need to get U and V matrices
        for(i=iLow;i<=iUp;i++){
            robust_givens_rotation(c,s,r,E(i,i),E(i,i+1));
            V2.Set_Identity_Matrix();
            V2(i,i)=c;V2(i,i+1)=s;
            V2(i+1,i)=-s;V2(i+1,i+1)=c;
            E.Derived()=E.Derived()*V2.Transposed();
            V.Derived()=V2*V.Derived();

            robust_givens_rotation(c,s,r,E(i,i),E(i+1,i));
            U2.Set_Identity_Matrix();
            U2(i,i)=c;U2(i,i+1)=s;
            U2(i+1,i)=-s;U2(i+1,i+1)=c;
            E.Derived()=U2*E.Derived();
            U.Derived()=U.Derived()*U2.Transposed();
        }
    }
    PHYSBAM_FATAL_ERROR("Singular_Value_Decomposition: Failed, too many iterations! Exiting.\n");
}
//####################################################################################
// Function Upper_Bidiagonalization
//####################################################################################
template<class T, class T_MATRIX> template<class T_MATRIX2, class T_MATRIX3, class T_MATRIX4> void MATRIX_BASE<T,T_MATRIX>::
Upper_Bidiagonalization(MATRIX_BASE<T,T_MATRIX2>& U, MATRIX_BASE<T,T_MATRIX3>& B, MATRIX_BASE<T,T_MATRIX4>& V)
{
    int m=Rows();
    int n=Columns();
    //sanity checks on incoming matrices
    assert(U.Rows()==U.Columns());assert(V.Rows()==V.Columns());
    assert(U.Columns()==m);assert(V.Rows()==n);
    assert(B.Rows()==m);assert(B.Columns()==n);
    
    if(m<n){//need to process the transpose instead
        B.Derived().Transpose();
        Derived().Transposed().Upper_Bidiagonalization(V,B,U);
        MATRIX_MXN<T> V_old;V_old.Resize(V.Rows(),V.Columns());
        V_old=V;int i,j;
        for(i=1;i<=m;i++){for(j=1;j<=n;j++){V.Derived()(j,i)=V_old(j,m-i+1);}}
        for(i=m+1;i<=n;i++){for(j=1;j<=n;j++){V.Derived()(j,i)=V_old(j,i);}}
        MATRIX_MXN<T> B_sub;B_sub.Resize(m,m);
        for(i=1;i<=m;i++){for(j=1;j<=m;j++){B_sub(j,i)=B(j,i);}}
        B_sub.Transpose();B_sub.Flip_Vertical();B_sub.Flip_Horizontal();
        B.Set_Zero_Matrix();B.Derived().Transpose();
        for(i=1;i<=m;i++){for(j=1;j<=m;j++){B(j,i)=B_sub(j,i);}}
        U.Flip_Horizontal();
        return;}

    B.Set_Zero_Matrix();
    U.Set_Identity_Matrix();
    V.Set_Identity_Matrix();
    VECTOR_ND<T> v1, v2, betaU, betaV;
    betaU.Resize(m); betaV.Resize(n);
    betaU.Set_Zero();betaV.Set_Zero();
    
    //uses Householder reflections to bidiagonalize matrix
    MATRIX_MXN<T> submatrix;
    int j,k;
    for(k=1;k<=n;k++){
        v1.Resize(m-k+1);v2.Resize(m-k+1);
        for(j=1;j<=m-k+1;j++){v1(j)=(*this)(k-1+j,k);}
        v1.Generalized_Householder_Transform(B(k,k),betaU(k),v2);
        for(j=1;j<=m-k+1;j++){(*this)(k-1+j,k)=v2(j);}
        submatrix.Resize(m-k+1,n-k);this->Get_Submatrix(k,k+1,submatrix);
        submatrix=submatrix-betaU(k)*(MATRIX_MXN<T>::Outer_Product(v2,submatrix.Left_Multiply_By_Vector(v2)));
        this->Set_Submatrix(k,k+1,submatrix);
        if(k<n-1){v1.Resize(n-k);v2.Resize(n-k);
            v2.Set_Zero();for(j=1;j<=n-k;j++){v1(j)=(*this)(k,k+j);}
            v1.Generalized_Householder_Transform(B(k,k+1),betaV(k),v2);
            for(j=1;j<=n-k;j++){(*this)(k,k+j)=v2(j);}
            submatrix.Resize(m-k,n-k);this->Get_Submatrix(k+1,k+1,submatrix);
            submatrix=(submatrix.Transposed()-betaV(k)*(MATRIX_MXN<T>::Outer_Product(v2,(submatrix.Transposed().Left_Multiply_By_Vector(v2))))).Transposed();
            this->Set_Submatrix(k+1,k+1,submatrix);}
        else if(k==n-1){
            B(n-1,n)=(*this)(n-1,n);}
    }
    //set U and V matrices from saved transforms
    for(k=n;k>=1;k--){
        v1.Resize(m-k+1);v1.Set_Zero();for(j=1;j<=m-k+1;j++){v1(j)=(*this)(k-1+j,k);}
        submatrix.Resize(m-k+1,m-k+1);U.Get_Submatrix(k,k,submatrix);
        submatrix=submatrix-betaU(k)*(MATRIX_MXN<T>::Outer_Product(v1,submatrix.Left_Multiply_By_Vector(v1)));
        U.Set_Submatrix(k,k,submatrix);
    }
    for(k=n-2;k>=1;k--){
        v1.Resize(n-k);v1.Set_Zero();for(j=1;j<=n-k;j++){v1(j)=(*this)(k,k+j);}
        submatrix.Resize(n-k,n-k+1);V.Get_Submatrix(k+1,k,submatrix);
        submatrix=submatrix-betaV(k)*(MATRIX_MXN<T>::Outer_Product(v1,submatrix.Left_Multiply_By_Vector(v1)));
        V.Set_Submatrix(k+1,k,submatrix);
    }
    return;
}
//####################################################################################
// Function Flip_Horizontal
//####################################################################################
template<class T, class T_MATRIX> void MATRIX_BASE<T,T_MATRIX>::Flip_Horizontal()
{
    VECTOR_ND<T> v;
    v.Resize(Rows());
    for(int i=1;i<=Columns()/2;i++){for(int j=1;j<=Rows();j++){
            v(j)=(*this)(j,i);(*this)(j,i)=(*this)(j,Columns()-i+1);(*this)(j,Columns()-i+1)=v(j);}}
}
//####################################################################################
// Function Flip_Vertical
//####################################################################################
template<class T, class T_MATRIX> void MATRIX_BASE<T,T_MATRIX>::Flip_Vertical()
{
    VECTOR_ND<T> v;
    v.Resize(Columns());
    for(int j=1;j<=Rows()/2;j++){for(int i=1;i<=Columns();i++){
            v(i)=(*this)(j,i);(*this)(j,i)=(*this)(Rows()-j+1,i);(*this)(Rows()-j+1,i)=v(i);}}
}
//####################################################################################
#define INSTANTIATION_HELPER(T) \
    template void MATRIX_BASE<T,MATRIX_MXN<T> >::In_Place_Robust_Householder_QR_Solve<VECTOR_ND<T>,VECTOR_ND<int> >(VECTOR_BASE<T,VECTOR_ND<T> >&,VECTOR_BASE<int,VECTOR_ND<int> >&); \
    template void MATRIX_BASE<T,MATRIX<T,6,6> >::In_Place_Gram_Schmidt_QR_Factorization<MATRIX<T,6,6> >(MATRIX_BASE<T,MATRIX<T,6,6> >&); \
    template void MATRIX_BASE<T,MATRIX_MXN<T> >::In_Place_Gram_Schmidt_QR_Factorization<MATRIX_MXN<T> >(MATRIX_BASE<T,MATRIX_MXN<T> >&); \
    template void MATRIX_BASE<T,MATRIX<T,1,1> >::In_Place_Cholesky_Factorization(); \
    template void MATRIX_BASE<T,MATRIX<T,3,3> >::In_Place_Cholesky_Factorization(); \
    template void MATRIX_BASE<T,MATRIX<T,4,4> >::In_Place_Cholesky_Factorization(); \
template void MATRIX_BASE<T,MATRIX<T,6,6> >::In_Place_Cholesky_Factorization(); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::In_Place_Cholesky_Factorization(); \
template void MATRIX_BASE<T,MATRIX<T,6,6> >::In_Place_PLU_Factorization<MATRIX<T,6,6> >(MATRIX_BASE<T,MATRIX<T,6,6> >&,VECTOR<int,6>&); \
template void MATRIX_BASE<T,MATRIX<T,12,12> >::In_Place_PLU_Factorization<MATRIX<T,12,12> >(MATRIX_BASE<T,MATRIX<T,12,12> >&,VECTOR<int,12>&); \
template void MATRIX_BASE<T,MATRIX<T,4,4> >::In_Place_PLU_Factorization<MATRIX<T,4,4> >(MATRIX_BASE<T,MATRIX<T,4,4> >&,VECTOR<int,4>&); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::In_Place_PLU_Factorization<MATRIX_MXN<T> >(MATRIX_BASE<T,MATRIX_MXN<T> >&,VECTOR_ND<int>&); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::In_Place_LU_Factorization<MATRIX_MXN<T> >(MATRIX_BASE<T,MATRIX_MXN<T> >&); \
template int MATRIX_BASE<T,MATRIX_MXN<T> >::Number_Of_Nonzero_Rows(const T threshold) const; \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::Householder_QR_Factorization<MATRIX_MXN<T>,MATRIX_MXN<T> >(MATRIX_BASE<T,MATRIX_MXN<T> >&,MATRIX_BASE<T,MATRIX_MXN<T> >&); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::In_Place_Householder_QR_Factorization_With_Column_Pivoting(MATRIX_BASE<T,MATRIX_MXN<T> >&,MATRIX_BASE<T,MATRIX_MXN<T> >& P,const T); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::In_Place_Householder_QR_Factorization_With_Full_Pivoting(MATRIX_BASE<T,MATRIX_MXN<T> >&,MATRIX_BASE<T,MATRIX_MXN<T> >& P,const T); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::Robust_Householder_QR_Solve(VECTOR_BASE<T,VECTOR_ND<T> >&,const VECTOR_BASE<T,VECTOR_ND<T> >&,const bool,const T); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::Upper_Bidiagonalization<MATRIX_MXN<T>,MATRIX_MXN<T>,MATRIX_MXN<T> >(MATRIX_BASE<T,MATRIX_MXN<T> >&,MATRIX_BASE<T,MATRIX_MXN<T> >&,MATRIX_BASE<T,MATRIX_MXN<T> >&); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::Singular_Value_Decomposition<MATRIX_MXN<T>,MATRIX_MXN<T>,MATRIX_MXN<T> >(MATRIX_BASE<T,MATRIX_MXN<T> >&,MATRIX_BASE<T,MATRIX_MXN<T> >&,MATRIX_BASE<T,MATRIX_MXN<T> >&,const T, const unsigned); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::Get_Singular_Values<VECTOR_ND<T> >(VECTOR_BASE<T,VECTOR_ND<T> >&,const T,const unsigned); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::Flip_Horizontal(); \
template void MATRIX_BASE<T,MATRIX_MXN<T> >::Flip_Vertical();

INSTANTIATION_HELPER(float);
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
INSTANTIATION_HELPER(double);
#endif
