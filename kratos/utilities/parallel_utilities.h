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
                local_reducer.LocalMerge(f(*(mit_begin + k)));
            }
            global_reducer.ThreadSafeMerge(local_reducer);
        }
        return global_reducer.GetValue();
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

    //pure c++11 version (can handle exceptions)
    template <class TUnaryFunction>
    inline void for_pure_c11(TUnaryFunction &&f)
    {
        KRATOS_TRY
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
            } catch(const std::exception& e) {
                std::cout << std::endl << "THREAD number: " << i << " caught exception " << e.what() << std::endl;
        }

        }

        KRATOS_CATCH("")
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
                local_reducer.LocalMerge(f(k));
            }
            global_reducer.ThreadSafeMerge(local_reducer);
        }
        return global_reducer.GetValue();
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
    typedef TDataType value_type;
    TDataType mvalue = TDataType(); //i am deliberately making the member value public, to allow one to change it as needed

    TDataType GetValue() const
    {
        return mvalue;
    }

    void LocalMerge(const TDataType value){
        mvalue += value;
    }

    //a user could define a ThreadSafeMerge for his specific reduction case
    void ThreadSafeMerge(const SumReduction<TDataType>& rOther)
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
    typedef TDataType value_type;
    TDataType mvalue = TDataType(); //i am deliberately making the member value public, to allow one to change it as needed

    TDataType GetValue() const
    {
        return mvalue;
    }
    void LocalMerge(const TDataType value){
        mvalue -= value;
    }

    void ThreadSafeMerge(const SubReduction<TDataType>& rOther)
    {
        #pragma omp atomic
        mvalue += rOther.mvalue;
    }
};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<class TDataType>
class MaxReduction
{
public:
    typedef TDataType value_type;
    TDataType mvalue = std::numeric_limits<TDataType>::lowest(); //i am deliberately making the member value public, to allow one to change it as needed
    TDataType GetValue() const
    {
        return mvalue;
    }
    void LocalMerge(const TDataType value){
        mvalue = std::max(mvalue,value);
    }
    void ThreadSafeMerge(const MaxReduction<TDataType>& rOther)
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
    typedef TDataType value_type;
    TDataType mvalue = std::numeric_limits<TDataType>::max(); //i am deliberately making the member value public, to allow one to change it as needed

    TDataType GetValue() const
    {
        return mvalue;
    }
    void LocalMerge(const TDataType value){
        mvalue = std::min(mvalue,value);
    }

    void ThreadSafeMerge(const MinReduction<TDataType>& rOther)
    {
        #pragma omp critical
        mvalue = std::min(mvalue,rOther.mvalue);
    }
};

template <class... Reducer>
struct CombinedReduction {
    typedef std::tuple<typename Reducer::value_type...> value_type;

    std::tuple<Reducer...> mChild;

    CombinedReduction() {}

    value_type GetValue(){
        value_type return_value;
        fill_value<0>(return_value);
        return return_value;
    }

    template <int I, class T>
    typename std::enable_if<(I < sizeof...(Reducer)), void>::type
    fill_value(T& v) {
        std::get<I>(v) = std::get<I>(mChild).GetValue();
        fill_value<I+1>(v);
        };

    template <int I, class T>
    typename std::enable_if<(I == sizeof...(Reducer)), void>::type
    fill_value(T& v) {}

    // template <int I>
    // auto GetValue() const -> decltype(std::get<I>(mChild).GetValue()) {
    //     return std::get<I>(mChild).GetValue();
    // }

    template <class... T>
    void LocalMerge(const std::tuple<T...> &&v) {
        // Static recursive loop over tuple elements
        reduce_local<0>(v);
    }

    void ThreadSafeMerge(const CombinedReduction &other) {
        reduce_global<0>(other);
    }

    private:

        template <int I, class T>
        typename std::enable_if<(I < sizeof...(Reducer)), void>::type
        reduce_local(T &&v) {
            std::get<I>(mChild).LocalMerge(std::get<I>(v));
            reduce_local<I+1>(std::forward<T>(v));
        };

        template <int I, class T>
        typename std::enable_if<(I == sizeof...(Reducer)), void>::type
        reduce_local(T &&v) {
            // Exit static recursion
        }

        template <int I>
        typename std::enable_if<(I < sizeof...(Reducer)), void>::type
        reduce_global(const CombinedReduction &other) {
            std::get<I>(mChild).ThreadSafeMerge(std::get<I>(other.mChild));
            reduce_global<I+1>(other);
        }

        template <int I>
        typename std::enable_if<(I == sizeof...(Reducer)), void>::type
        reduce_global(const CombinedReduction &other) {
            // Exit static recursion
        }
};

}  // namespace Kratos.

#endif // KRATOS_PARALLEL_UTILITIES_H_INCLUDED  defined

