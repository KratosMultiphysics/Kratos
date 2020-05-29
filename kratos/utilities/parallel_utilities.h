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
//                   Denis Demidov
//

#if !defined(KRATOS_PARALLEL_UTILITIES_H_INCLUDED )
#define  KRATOS_PARALLEL_UTILITIES_H_INCLUDED


// System includes
#include<iostream>
#include<array>
#include<vector>
#include<tuple>
#include<cmath>
#include<limits>
#include<omp.h>

#include <future>
#include <thread>
#include "utilities/reduction_utilities.h"

// External includes


// Project includes
#include "includes/define.h"
#include "includes/global_variables.h"

namespace Kratos
{
///@addtogroup KratosCore

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<
        class TContainerType,
        class TIteratorType=typename TContainerType::iterator,
        int TMaxThreads=Globals::MaxAllowedThreads
        >
class BlockPartition
{
public:


    BlockPartition(TIteratorType it_begin,
                   TIteratorType it_end,
                   int Nchunks = omp_get_max_threads())
        : mit_begin(it_begin), mit_end(it_end), mNchunks(Nchunks)
    {
        ptrdiff_t mBlockPartitionSize = (it_end-it_begin) / mNchunks;
        mBlockPartition[0] = it_begin;
        mBlockPartition[mNchunks] = it_end;
        for(int i = 1; i < mNchunks; i++)
            mBlockPartition[i] = mBlockPartition[i-1] + mBlockPartitionSize ;
    }

    BlockPartition(TContainerType& rData,
                   int Nchunks = omp_get_max_threads())
        : BlockPartition(rData.begin(), rData.end(), Nchunks)
    {}

    virtual ~BlockPartition() {}

    //version enabling reduction
    //simple version
    template <class TUnaryFunction>
    inline void for_each(TUnaryFunction&& f)
    {
        #pragma omp parallel for
        for(int i=0; i<mNchunks; ++i)
        {
            for (auto it = mBlockPartition[i]; it != mBlockPartition[i+1]; ++it)
            {
                f(*it); //note that we pass a reference to the value, not the iterator
            }
        }
    }

    //version with reduction
    template <class TReducer, class TUnaryFunction>
    inline typename TReducer::value_type for_reduce(TUnaryFunction &&f)
    {
        TReducer global_reducer;
        #pragma omp parallel for
        for(int i=0; i<mNchunks; ++i)
        {
            TReducer local_reducer;
            for (auto it = mBlockPartition[i]; it != mBlockPartition[i+1]; ++it)
            {
                local_reducer.LocalReduce(f(*it));
            }
            global_reducer.ThreadSafeReduce(local_reducer);
        }
        return global_reducer.GetValue();
    }

private:

    TIteratorType mit_begin, mit_end;
    int mNchunks;
    std::array<TIteratorType, TMaxThreads> mBlockPartition;

};

template <class TContainerType, class TFunctionType>
void block_for_each(TContainerType &&v, TFunctionType &&func) {
    BlockPartition<typename std::decay<TContainerType>::type>(std::forward<TContainerType>(v)).for_each(std::forward<TFunctionType>(func));
}

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
///This class is useful for iteration over containers
template<class TIndexType, int TMaxThreads=Globals::MaxAllowedThreads>
class IndexPartition
{
public:

    IndexPartition(TIndexType Size,
                   int Nchunks = omp_get_max_threads())
        : mSize(Size), mNchunks(Nchunks)
    {
        int mBlockPartitionSize = mSize / mNchunks;
        mBlockPartition[0] = 0;
        mBlockPartition[mNchunks] = mSize;
        for(int i = 1; i < mNchunks; i++)
            mBlockPartition[i] = mBlockPartition[i-1] + mBlockPartitionSize ;

    }

    virtual ~IndexPartition() {}

    //pure c++11 version (can handle exceptions)
    template <class TUnaryFunction>
    inline void for_pure_c11(TUnaryFunction &&f)
    {

        std::vector< std::future<void> > runners(mNchunks);
        const auto& partition = mBlockPartition;
        for(int i=0; i<mNchunks; ++i)
        {
            runners[i] = std::async(std::launch::async, [&partition, i,  &f]()
                {
                    for (auto k = partition[i]; k < partition[i+1]; ++k)
                        f(k);
                });
        }

        //here we impose a syncronization and we check the exceptions
        for(int i=0; i<mNchunks; ++i)
        {
            try {
                    runners[i].get();
            }
            catch(Exception& e) {
                KRATOS_ERROR << std::endl << "THREAD number: " << i << " caught exception " << e.what() << std::endl;
            } catch(std::exception& e) {
                KRATOS_ERROR << std::endl << "THREAD number: " << i << " caught exception " << e.what() << std::endl;
            } catch(...) {
                KRATOS_ERROR << std::endl << "unknown error" << std::endl;
            }
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

    //version with reduction
    template <class TReducer, class TUnaryFunction>
    inline typename TReducer::value_type for_reduce(TUnaryFunction &&f)
    {
        TReducer global_reducer;
        #pragma omp parallel for
        for(int i=0; i<mNchunks; ++i)
        {
            TReducer local_reducer;
            for (auto k = mBlockPartition[i]; k < mBlockPartition[i+1]; ++k)
            {
                local_reducer.LocalReduce(f(k));
            }
            global_reducer.ThreadSafeReduce(local_reducer);
        }
        return global_reducer.GetValue();
    }

private:

    TIndexType mSize;
    int mNchunks;

    std::array<TIndexType, TMaxThreads> mBlockPartition;

};

template <class TContainerType, class TFunctionType>
void indexbased_for_each(TContainerType &&v, TFunctionType &&func) {
    IndexPartition<typename std::decay<TContainerType>::type>(std::forward<TContainerType>(v)).for_each(std::forward<TFunctionType>(func));
}


}  // namespace Kratos.

#endif // KRATOS_PARALLEL_UTILITIES_H_INCLUDED  defined

