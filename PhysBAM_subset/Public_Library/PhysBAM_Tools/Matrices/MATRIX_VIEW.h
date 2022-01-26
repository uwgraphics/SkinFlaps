//#####################################################################
// Copyright 2013, Rahul Sheth, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MATRIX_VIEW
//#####################################################################
// This class allows submatrices to be accessed and used without copying, like ARRAY_VIEW
// From a logical view, this matrix's row and column indices behave in the same order as the blocks were specified, only with irrelevant rows/columns collapsed.
// Internally, the base pointers, and ranges under each base pointer are traversed in the optimal sequential memory access pattern.
//#####################################################################
#ifndef __MATRIX_VIEW__
#define __MATRIX_VIEW__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Matrices/MATRIX_BASE.h>
#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Utilities/TYPE_UTILITIES.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Vectors/VECTOR_FORWARD.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Tools/Vectors/VECTOR_VIEW.h>
namespace PhysBAM{

//TODO: Change to FIELD_COMPARE?
struct SECOND_COMPARE
{
    template<class T> bool operator()(const T& a1,const T& a2) const
    {return a1.y<a2.y;}
    template<class T,class T_MATRIX,class T_VECTOR> bool operator()(PAIR<int,GENERAL_MATRIX<T,T_MATRIX,T_VECTOR>*>const& a1,PAIR<int,GENERAL_MATRIX<T,T_MATRIX,T_VECTOR>*>const& a2) const
    {return (*a1.y)<(*a2.y);}
};

struct SOURCE_COLUMN_FIRST_COMPARE
{
typedef VECTOR<int,2> INT2;
typedef PAIR<RANGE<INT2>,RANGE<INT2> > RANGE_MAPPING;
    bool operator()(const PAIR<int,RANGE_MAPPING>& a1,const PAIR<int,RANGE_MAPPING>& a2) const
    {const VECTOR<int,2>& mc1=a1.y.x.min_corner,mc2=a2.y.x.min_corner;
        return mc1(2)<mc2(2)||(mc1(2)==mc2(2)&&mc1(1)<mc2(1));}
};

template<class T,class T_MATRIX,class T_VECTOR>
class GENERAL_MATRIX
{
public:
    T_MATRIX* matrix;
    T_VECTOR* vector;
    T_VECTOR* negative_vector;
    int mode; // 0=matrix, 1=omega vector
    T my_zero;

    GENERAL_MATRIX(T_MATRIX& matrix_input)
    : matrix(&matrix_input),vector(NULL),negative_vector(NULL),mode(0),my_zero(0)
    {}

    GENERAL_MATRIX(T_VECTOR& vector_input)
    : matrix(NULL),vector(&vector_input),negative_vector(new T_VECTOR(-*vector)),mode(1),my_zero(0)
    {assert(vector_input.Size()==3);}

    ~GENERAL_MATRIX()
    {delete negative_vector;}

    const T& operator()(int i,int j) const{
        if(mode){
            assert(i>=1&&i<=3&&j>=1&&j<=3);
            if(i==j) return my_zero;
            else if(i==j%3+1) return vector->operator()(i%3+1);
            else /*if(j==i%3+1)*/ return negative_vector->operator()(j%3+1);
        }
        else return matrix->operator()(i,j);
    }

    T& operator()(int i,int j){
        if(mode){
            assert(i>=1&&i<=3&&j>=1&&j<=3);
            if(i==j) return my_zero;
            else if(i==j%3+1) return vector->operator()(i%3+1);
            else /*if(j==i%3+1)*/ return negative_vector->operator()(j%3+1); // warning: positive and negative components can go out of sync...
        }
        else return matrix->operator()(i,j);
    }

    int Rows() const
    {return mode?3:matrix->Rows();}

    int Columns() const
    {return mode?3:matrix->Columns();}

    bool operator<(const GENERAL_MATRIX<T,T_MATRIX,T_VECTOR>& other) const
    {return (mode?(void*)vector:matrix)<(other.mode?(void*)(other.vector):other.matrix);}

};

template<class T,class T_MATRIX>
class MATRIX_VIEW
{

typedef VECTOR<int,2> INT2;
typedef PAIR<RANGE<INT2>,RANGE<INT2> > RANGE_MAPPING;
typedef typename MATRIX_INFO<T_MATRIX>::RIGHT_VECTOR RIGHT_VECTOR;typedef typename MATRIX_INFO<T_MATRIX>::LEFT_VECTOR LEFT_VECTOR;

public:
    ARRAY<PAIR<T_MATRIX*,ARRAY<RANGE_MAPPING> > > blocks; // list of (&matrix, list of (src range,dst range))

    ARRAY<int> matrix_indirect;           // matrices and blocks in sequential memory order
    ARRAY<ARRAY<int> > block_indirect;

    INT2 size;
    bool auto_expand;
    T my_zero;

public:
    MATRIX_VIEW(ARRAY<PAIR<T_MATRIX*,ARRAY<RANGE_MAPPING> > > blocks_input,INT2 size_input,bool auto_expand_input=false)
        :blocks(blocks_input),size(size_input),auto_expand(auto_expand_input),my_zero(0)
    {
        ARRAY<RANGE<INT2> > dst_ranges;
        for(int i=1;i<=blocks.m;++i){ // check for out of bounds
            assert(blocks(i).x);
            RANGE<INT2> source_matrix_size(INT2(1,1),INT2(blocks(i).x->Rows(),blocks(i).x->Columns()));
            const ARRAY<RANGE_MAPPING>& ranges=blocks(i).y;
            for(int j=1;j<=ranges.m;++j){
                assert(source_matrix_size.Contains(ranges(j).x)); // check src bound
                assert(ranges(j).x.Edge_Lengths()==ranges(j).y.Edge_Lengths()); // src and dst must equal in size
                for(int k=1;k<=dst_ranges.Size();++k){assert(!ranges(j).y.Intersection(dst_ranges(k)));}// prevent overlap
                dst_ranges.Append(ranges(j).y);
                assert(ranges(j).y.min_corner(1)>=1 && ranges(j).y.min_corner(2)>=1);   // check dst bound
                if(auto_expand){
                    size(1)=max(size(1),ranges(j).y.max_corner(1));
                    size(2)=max(size(2),ranges(j).y.max_corner(2));}
                else{assert(ranges(j).y.max_corner(1)<=size(1) && ranges(j).y.max_corner(2)<=size(2));}
            }}
        Build_Indirection();
    }
    void Build_Indirection()
    {
        // sort to optimize memory access 
        matrix_indirect.Remove_All();
        ARRAY<PAIR<int,T_MATRIX*> > sorted_matrices;
        for(int i=1;i<=blocks.m;++i) sorted_matrices.Append(PAIR<int,T_MATRIX*>(i,blocks(i).x));
        Stable_Sort(sorted_matrices,SECOND_COMPARE());
        for(int i=1;i<=blocks.m;++i) matrix_indirect.Append(sorted_matrices(i).x);

        block_indirect.Remove_All();
        for(int i=1;i<=blocks.m;++i){
            const ARRAY<RANGE_MAPPING>& ranges=blocks(matrix_indirect(i)).y;
            ARRAY<PAIR<int,RANGE_MAPPING> > sorted_ranges;
            for(int j=1;j<=ranges.m;++j) sorted_ranges.Append(PAIR<int,RANGE_MAPPING>(j,ranges(j)));
            Stable_Sort(sorted_ranges,SOURCE_COLUMN_FIRST_COMPARE());
            ARRAY<int> block_indirect_element;
            for(int j=1;j<=ranges.m;++j) block_indirect_element.Append(sorted_ranges(j).x);
            block_indirect.Append(block_indirect_element);
        }
        /*LOG::cout<<"Matrix indirect = "<<matrix_indirect<<std::endl;
          LOG::cout<<"Block indirect = "<<block_indirect<<std::endl;*/
    }

    int Binary_Search(const ARRAY<TRIPLE<int,int,int> >& array,int query) const
    {assert(array.m);
        int begin=1,end=array.m;
        while(begin<end){
            int test=(begin+end)/2;
            if(test!=1&&query<array(test).x){end=test-1;}
            else if(test!=array.m&&query>=array(test+1).x){begin=test+1;}
            else return test;}
        return begin;
    }

    int Rows() const
    {return size(1);}

    int Columns() const
    {return size(2);}

    T& operator()(const int ii,const int jj) //TODO this is O(n) random access. A kd-tree will make this O(log n)
    {assert(RANGE<INT2>(INT2(1,1),size).Lazy_Inside(INT2(ii,jj)));
        for(int i=1;i<=blocks.m;++i){
            const ARRAY<RANGE_MAPPING>& ranges=blocks(i).y;
            for(int j=1;j<=ranges.m;++j)
                if(ranges(j).y.Contains(INT2(ii,jj))){
                    INT2 offset=ranges(j).x.min_corner-ranges(j).y.min_corner;
                    return blocks(i).x->operator()(ii+offset(1),jj+offset(2));
                }
        }
        return my_zero;
    }

    const T& operator()(const int ii,const int jj) const //TODO this is O(n) random access. A kd-tree will make this O(log n)
    {assert(RANGE<INT2>(INT2(1,1),size).Lazy_Inside(INT2(ii,jj)));
        for(int i=1;i<=blocks.m;++i){
            const ARRAY<RANGE_MAPPING>& ranges=blocks(i).y;
            for(int j=1;j<=ranges.m;++j)
                if(ranges(j).y.Contains(INT2(ii,jj))){
                    INT2 offset=ranges(j).x.min_corner-ranges(j).y.min_corner;
                    return blocks(i).x->operator()(ii+offset(1),jj+offset(2));
                }
        }
        return my_zero;
    }

    inline INTERVAL<int> Interval_From_Range_Axis(const RANGE<INT2>& range,int axis) const
    {
        return INTERVAL<int>(range.min_corner(axis),range.max_corner(axis));
    }

    template<class IN_VECTOR,class OUT_VECTOR>
    void Right_Multiply_By_Vector(IN_VECTOR& v,OUT_VECTOR& result) const
    {assert(v.Size()==Columns());
     assert(result.Size()==Rows());
        for(int i=1;i<=blocks.m;++i){
            const T_MATRIX* matrix=blocks(matrix_indirect(i)).x;
            const ARRAY<RANGE_MAPPING >& ranges=blocks(matrix_indirect(i)).y;
            for(int j=1;j<=ranges.m;++j){
                const RANGE_MAPPING& range_mapping=ranges(block_indirect(i)(j));
                const RANGE<INT2>& src_range=range_mapping.x;
                const RANGE<INT2>& range=range_mapping.y;
                INT2 offset=src_range.min_corner-range.min_corner;
                Add_Multiply_Helper(v,result,range,matrix,offset);
            }
        }
    }

    void Add_Multiply_Helper(const RIGHT_VECTOR& v,RIGHT_VECTOR& result,const RANGE<INT2>& range,const T_MATRIX* matrix,const INT2 offset) const{
        for(int jj=range.min_corner(2);jj<=range.max_corner(2);++jj) for(int ii=range.min_corner(1);ii<=range.max_corner(1);++ii) result(ii)+=(*matrix)(ii+offset(1),jj+offset(2))*v(jj);
    }//TODO deal with constness
    void Add_Multiply_Helper(/*const*/ VECTOR_VIEW<T,RIGHT_VECTOR>& v,RIGHT_VECTOR& result,const RANGE<INT2>& range,const T_MATRIX* matrix,const INT2 offset) const{
        for(VECTOR_VIEW_ITERATOR<T,RIGHT_VECTOR> jtr(v,Interval_From_Range_Axis(range,2));jtr.Valid();jtr.Next()) for(int ii=range.min_corner(1);ii<=range.max_corner(1);++ii) result(ii)+=(*matrix)(ii+offset(1),jtr.current_index+offset(2))*(*jtr);
    }
    void Add_Multiply_Helper(const RIGHT_VECTOR& v,VECTOR_VIEW<T,RIGHT_VECTOR>& result,const RANGE<INT2>& range,const T_MATRIX* matrix,const INT2 offset) const{
        VECTOR_VIEW_ITERATOR<T,RIGHT_VECTOR> itr_init(result,Interval_From_Range_Axis(range,1));
        for(int jj=range.min_corner(2);jj<=range.max_corner(2);++jj) for(VECTOR_VIEW_ITERATOR<T,RIGHT_VECTOR> itr(itr_init);itr.Valid();itr.Next()) (*itr)+=(*matrix)(itr.current_index+offset(1),jj+offset(2))*v(jj);
    }
    void Add_Multiply_Helper(/*const*/ VECTOR_VIEW<T,RIGHT_VECTOR>& v,VECTOR_VIEW<T,RIGHT_VECTOR>& result,const RANGE<INT2>& range,const T_MATRIX* matrix,const INT2 offset) const{
        VECTOR_VIEW_ITERATOR<T,RIGHT_VECTOR> itr_init(result,Interval_From_Range_Axis(range,1));
        for(VECTOR_VIEW_ITERATOR<T,RIGHT_VECTOR> jtr(v,Interval_From_Range_Axis(range,2));jtr.Valid();jtr.Next()) for(VECTOR_VIEW_ITERATOR<T,RIGHT_VECTOR> itr(itr_init);itr.Valid();itr.Next()) (*itr)+=(*matrix)(itr.current_index+offset(1),jtr.current_index+offset(2))*(*jtr);
    }
    //void Left_Multiply_By_Vector(const LEFT_VECTOR& v,LEFT_VECTOR& result) const
    //{assert(v.Size()==Rows());
    // assert(result.Size()==Columns());
    //    for(int i=1;i<=blocks.m;++i){
    //        const T_MATRIX* matrix=blocks(matrix_indirect(i)).x;
    //        const ARRAY<RANGE_MAPPING >& ranges=blocks(matrix_indirect(i)).y;
    //        for(int j=1;j<=ranges.m;++j){
    //            const RANGE_MAPPING& range_mapping=ranges(block_indirect(i)(j));
    //            const RANGE<INT2>& src_range=range_mapping.x;
    //            const RANGE<INT2>& range=range_mapping.y;
    //            INT2 offset=src_range.min_corner-range.min_corner; //TODO combine with right multiply
    //            for(int jj=range.min_corner(2);jj<=range.max_corner(2);++jj) for(int ii=range.min_corner(1);ii<=range.max_corner(1);++ii) result(jj)+=(*matrix)(ii+offset(1),jj+offset(2))*v(ii);
    //        }
    //    }
    //}

//#####################################################################
};

template<class T, class T_MATRIX> inline std::ostream& operator<<(std::ostream& output_stream,const MATRIX_VIEW<T,T_MATRIX>& A)
{output_stream<<std::endl;for(int i=1;i<=A.Rows();i++){output_stream<<"|";for(int j=1;j<=A.Columns();j++){output_stream<<std::setw(16)<<A(i,j);}output_stream<<"|"<<std::endl;}return output_stream;}

}
#endif
