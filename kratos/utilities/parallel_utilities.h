//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_PARALLEL_UTILITIES_H_INCLUDED )
#define  KRATOS_PARALLEL_UTILITIES_H_INCLUDED


// System includes
#include<iostream>
#include<array>
#include<vector>
#include<cmath>
#include<limits>
#include<omp.h>

// External includes


// Project includes
#include "includes/define.h"


namespace Kratos
{
///@addtogroup KratosCore

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<class TContainerType, class TIteratorType=typename TContainerType::iterator>
class BlockPartition
{
public:

    static constexpr int max_allowed_threads = 128;

    BlockPartition(TIteratorType it_begin,
                   TIteratorType it_end,
                   int Nchunks = omp_get_max_threads())
        : mit_begin(it_begin), mit_end(it_end), mNchunks(Nchunks)
    {
        int NumTerms = it_end-it_begin;
        int mBlockPartitionize = NumTerms / mNchunks;
        mBlockPartition[0] = 0;
        mBlockPartition[mNchunks] = NumTerms;
        for(int i = 1; i < mNchunks; i++)
            mBlockPartition[i] = mBlockPartition[i-1] + mBlockPartitionize ;

    }

    BlockPartition(TContainerType& rData,
                   int Nchunks = omp_get_max_threads())
        : BlockPartition(rData.begin(), rData.end(), Nchunks)
    {}

    virtual ~BlockPartition() {}

    //version enabling reduction
    template <class TFunctionWithReduction, class TReduceHelper>
    inline void for_each(
        TReduceHelper& reduce_helper,
        TFunctionWithReduction&& f
    )
    {
        #pragma omp parallel for
        for(int i=0; i<mNchunks; ++i)
        {
            TReduceHelper local_helper;
            for (auto p = mit_begin + mBlockPartition[i], e = mit_begin + mBlockPartition[i+1]; p != e; ++p)
            {
                f(*p, local_helper); //note that we pass a reference to the value, not the iterator
            }
            reduce_helper.ThreadSafeMerge(local_helper);
        }
    }

    //simple version
    template <class TUnaryFunction>
    inline void for_each(TUnaryFunction&& f)
    {
        #pragma omp parallel for
        for(int i=0; i<mNchunks; ++i)
        {
            for (int k = mBlockPartition[i]; k < mBlockPartition[i+1]; ++k)
            {
                f(*(mit_begin + k)); //note that we pass a reference to the value, not the iterator
            }
        }
    }

private:

    TIteratorType mit_begin, mit_end;
    int mNchunks;
    std::array<int, max_allowed_threads> mBlockPartition;

};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
///This class is useful for iteration over containers
template<class TIndexType>
class IndexPartition
{
public:

    static constexpr int max_allowed_threads = 128;

    IndexPartition(TIndexType Size,
                   int Nchunks = omp_get_max_threads())
        : mSize(Size), mNchunks(Nchunks)
    {
        int mBlockPartitionize = mSize / mNchunks;
        mBlockPartition[0] = 0;
        mBlockPartition[mNchunks] = mSize;
        for(int i = 1; i < mNchunks; i++)
            mBlockPartition[i] = mBlockPartition[i-1] + mBlockPartitionize ;

    }

    virtual ~IndexPartition() {}

    //version enabling reduction
    template <class TFunctionWithReduction, class TReduceHelper>
    inline void for_each(
        TReduceHelper& reduce_helper,
        TFunctionWithReduction&& f
    )
    {
        #pragma omp parallel for
        for(int i=0; i<mNchunks; ++i)
        {
            TReduceHelper local_helper;
            for (TIndexType k = mBlockPartition[i]; k < mBlockPartition[i+1]; ++k)
            {
                f( k , local_helper); //note that we pass a reference to the value, not the iterator
            }
            reduce_helper.ThreadSafeMerge(local_helper);
        }
    }

    //simple version
    template <class TUnaryFunction>
    inline void for_each(TUnaryFunction &&f)
    {
        #pragma omp parallel for
        for(int i=0; i<mNchunks; ++i)
        {
            for (auto k = mBlockPartition[i]; k < mBlockPartition[i+1]; ++k)
            {
                f(k); //note that we pass a reference to the value, not the iterator
            }
        }
    }

private:

    TIndexType mSize;
    int mNchunks;

    std::array<TIndexType, max_allowed_threads> mBlockPartition;

};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<class TDataType>
class SumReduction
{
public:
    TDataType mvalue = TDataType(); //i am deliberately making the member value public, to allow one to change it as needed

    //a user could define a ThreadSafeMerge for his specific reduction case
    void ThreadSafeMerge(SumReduction<TDataType>& rOther)
    {
        #pragma omp atomic
        mvalue += rOther.mvalue;
    }
};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<class TDataType>
class SubReduction
{
public:
    TDataType mvalue = TDataType(); //i am deliberately making the member value public, to allow one to change it as needed

    //a user could define a ThreadSafeMerge for his specific reduction case
    void ThreadSafeMerge(SumReduction<TDataType>& rOther)
    {
        #pragma omp atomic
        mvalue -= rOther.mvalue;
    }
};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<class TDataType>
class MaxReduction
{
public:
    TDataType mvalue = -std::numeric_limits<TDataType>::max(); //i am deliberately making the member value public, to allow one to change it as needed

    //a user could define a ThreadSafeMerge for his specific reduction case
    void ThreadSafeMerge(SumReduction<TDataType>& rOther)
    {
        #pragma omp critical
        mvalue = std::max(mvalue,rOther.mvalue);
    }
};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<class TDataType>
class MinReduction
{
public:
    TDataType mvalue = std::numeric_limits<TDataType>::max(); //i am deliberately making the member value public, to allow one to change it as needed

    //a user could define a ThreadSafeMerge for his specific reduction case
    void ThreadSafeMerge(SumReduction<TDataType>& rOther)
    {
        #pragma omp critical
        mvalue = std::min(mvalue,rOther.mvalue);
    }
};


}  // namespace Kratos.

#endif // KRATOS_PARALLEL_UTILITIES_H_INCLUDED  defined

