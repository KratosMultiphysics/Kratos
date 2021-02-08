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

#if !defined(KRATOS_REDUCTION_UTILITIES_H_INCLUDED )
#define  KRATOS_REDUCTION_UTILITIES_H_INCLUDED


// System includes
#include <iostream>
#include <array>
#include <vector>
#include <tuple>
#include <cmath>
#include <limits>
#include <future>
#include <thread>

// External includes
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

// Project includes
#include "includes/define.h"
#include "includes/global_variables.h"

namespace Kratos
{
///@addtogroup KratosCore

/** @brief utility function to do a sum reduction
 */
template<class TDataType>
class SumReduction
{
public:
    typedef TDataType value_type;
    TDataType mvalue = TDataType(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TDataType GetValue() const
    {
        return mvalue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mvalue += value;
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const SumReduction<TDataType>& rOther)
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
    TDataType mvalue = TDataType(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TDataType GetValue() const
    {
        return mvalue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mvalue -= value;
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const SubReduction<TDataType>& rOther)
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
    TDataType mvalue = std::numeric_limits<TDataType>::lowest(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TDataType GetValue() const
    {
        return mvalue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mvalue = std::max(mvalue,value);
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const MaxReduction<TDataType>& rOther)
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
    TDataType mvalue = std::numeric_limits<TDataType>::max(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TDataType GetValue() const
    {
        return mvalue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mvalue = std::min(mvalue,value);
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const MinReduction<TDataType>& rOther)
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

    /// access to reduced value
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

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    template <class... T>
    void LocalReduce(const std::tuple<T...> &&v) {
        // Static recursive loop over tuple elements
        reduce_local<0>(v);
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const CombinedReduction &other) {
        reduce_global<0>(other);
    }

    private:

        template <int I, class T>
        typename std::enable_if<(I < sizeof...(Reducer)), void>::type
        reduce_local(T &&v) {
            std::get<I>(mChild).LocalReduce(std::get<I>(v));
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
            std::get<I>(mChild).ThreadSafeReduce(std::get<I>(other.mChild));
            reduce_global<I+1>(other);
        }

        template <int I>
        typename std::enable_if<(I == sizeof...(Reducer)), void>::type
        reduce_global(const CombinedReduction &other) {
            // Exit static recursion
        }
};

}  // namespace Kratos.

#endif // KRATOS_REDUCTION_UTILITIES_H_INCLUDED  defined
