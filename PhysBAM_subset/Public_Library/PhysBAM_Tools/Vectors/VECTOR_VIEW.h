//#####################################################################
// Copyright 2013, Rahul Sheth, Yue Yu.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VECTOR_VIEW
//#####################################################################
#ifndef __VECTOR_VIEW__
#define __VECTOR_VIEW__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Data_Structures/TRIPLE.h>
#include <PhysBAM_Tools/Math_Tools/INTERVAL.h>
#include <PhysBAM_Tools/Vectors/VECTOR_ND.h>
#include <PhysBAM_Tools/Arrays_Computations/SORT.h>
#include <PhysBAM_Tools/Log/LOG.h>
namespace PhysBAM{

//TODO: Change to FIELD_COMPARE?
struct FIRST_COMPARE
{
    template<class T> bool operator()(const T& a1,const T& a2) const
    {return a1.x<a2.x;}
};

template<class T,class T_VECTOR>
class VECTOR_VIEW
{
typedef PAIR<INTERVAL<int>,INTERVAL<int> > INTERVAL_MAPPING;
public:
    ARRAY<PAIR<T_VECTOR*,ARRAY<INTERVAL_MAPPING> > > intervals; // list of (&vector, list of (src interval,dst interval))
    ARRAY<PAIR<int,TRIPLE<int,int,int> > > binary_search_list; // logical begin index, (logical end index, physical array index, array's physical begin index)
    int size;
    bool auto_expand;
    T my_zero;
    VECTOR_VIEW(ARRAY<PAIR<T_VECTOR*,ARRAY<INTERVAL_MAPPING> > > intervals_input,int size_input,bool auto_expand_input=false)
        :intervals(intervals_input),size(size_input),auto_expand(auto_expand_input),my_zero(0)
    {
        ARRAY<INTERVAL<int> >dst_intervals;
        for(int i=1;i<=intervals.m;++i){ // check for out of bounds
            assert(intervals(i).x);
            INTERVAL<int> source_vector_size(1,intervals(i).x->Size());
            const ARRAY<INTERVAL_MAPPING>& ranges=intervals(i).y;
            for(int j=1;j<=ranges.m;++j){
                assert(source_vector_size.Contains(ranges(j).x)); // check src bound
                assert(ranges(j).x.Size()==ranges(j).y.Size()); // src and dst must equal in size
                for(int k=1;k<=dst_intervals.Size();++k){assert(!ranges(j).y.Intersection(dst_intervals(k)));}// prevent overlap
                dst_intervals.Append(ranges(j).y);
                assert(ranges(j).y.min_corner>=1);   // check dst bound
                if(auto_expand){
                    size=max(size,ranges(j).y.max_corner);}
                else{assert(ranges(j).y.max_corner<=size);}
                binary_search_list.Append(Tuple(ranges(j).y.min_corner,Tuple(ranges(j).y.max_corner,i,ranges(j).x.min_corner)));
            }}
        Stable_Sort(binary_search_list,FIRST_COMPARE());
    }
    int Size() const
    {return size;}
    int Binary_Search(const ARRAY<PAIR<int,TRIPLE<int,int,int> > >& array,int query) const
    {assert(array.m);
        int begin=1,end=array.m;
        while(begin<end){
            int test=(begin+end)/2;
            if(test!=1&&query<array(test).x){end=test-1;}
            else if(test!=array.m&&query>=array(test+1).x){begin=test+1;}
            else return test;}
        return begin;
    }
    T& operator()(const int ii) //TODO this is O(log n). Need O(1) sequential access...
    {
        if(ii<binary_search_list(1).x) return my_zero;
        int search_index=Binary_Search(binary_search_list,ii);
        if(ii>binary_search_list(search_index).y.x) return my_zero;
        return intervals(binary_search_list(search_index).y.y).x->operator()(ii+binary_search_list(search_index).y.z-binary_search_list(search_index).x);
    }
    void Write_To_Screen()
    {for(int i=1;i<=size;i++)printf("%5g ",(*this)(i));LOG::cout << std::endl;}
    void Write_To_Screen_Long()
    {for(int i=1;i<=size;i++)printf("%7g ",(*this)(i));LOG::cout << std::endl;}
};

template<class T,class T_VECTOR>
class VECTOR_VIEW_ITERATOR
{
public:
    VECTOR_VIEW<T,T_VECTOR>& vector;
    INTERVAL<int> query_interval;
    bool valid;
    int current_index;
    int current_search_index;
    int min_search_index;
    int max_search_index;

    VECTOR_VIEW_ITERATOR(VECTOR_VIEW<T,T_VECTOR>& vector_input,INTERVAL<int> query_interval_input)
    :vector(vector_input),query_interval(query_interval_input),valid(true)
    {
        assert(query_interval.Size()>=0 && INTERVAL<int>(1,vector.Size()).Contains(query_interval));
        current_index=query_interval.min_corner;
        if(query_interval.min_corner<vector.binary_search_list(1).x) current_index=query_interval.min_corner=vector.binary_search_list(1).x;
        if(query_interval.max_corner<vector.binary_search_list(1).x) valid=false;
        min_search_index=vector.Binary_Search(vector.binary_search_list,query_interval.min_corner);
        max_search_index=vector.Binary_Search(vector.binary_search_list,query_interval.max_corner);
        current_search_index=min_search_index;
        if(!valid) return;
        Update_Validity();
    }
    VECTOR_VIEW_ITERATOR(const VECTOR_VIEW_ITERATOR& other)
    :vector(other.vector),query_interval(other.query_interval),valid(other.valid),current_index(other.current_index),current_search_index(other.current_search_index),min_search_index(other.min_search_index),max_search_index(other.max_search_index)
    {
    }
    void Next()
    {assert(Valid());
        ++current_index;
        Update_Validity();}

    bool Valid() const
    {return valid;}

    T& operator*()
    {assert(Valid());
        return vector(current_index);}

    inline void Update_Validity()
    {
        if(current_search_index==max_search_index && current_index>query_interval.max_corner){valid=false;return;}
        if(current_index>vector.binary_search_list(current_search_index).y.x){
            if(current_search_index==vector.binary_search_list.Size()||current_search_index==max_search_index) valid=false;
            else current_index=vector.binary_search_list(++current_search_index).x;}
    }
//#####################################################################
};
}
#endif
