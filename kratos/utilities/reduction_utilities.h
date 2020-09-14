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
#include "includes/global_variables.h"

namespace Kratos
{
///@addtogroup KratosCore

template <class TDataType>
class Reducer
{
public:
    typedef TDataType value_type;
    TDataType mReducedValue; // deliberately making the member value public, to allow one to change it as needed

    /// Custom constructor
    Reducer(std::function<void(TDataType&, const TDataType&)>&& f)
        : mF(std::forward<std::function<void(TDataType&, const TDataType&)>>(f))
    {
    }

    /// Custom constructor with initialization
    Reducer(const TDataType& rInitialValue,
            std::function<void(TDataType&, const TDataType&)>&& f)
        : mReducedValue(rInitialValue),
          mF(std::forward<std::function<void(TDataType&, const TDataType&)>>(f))

    {
    }

    /// access to reduced value
    virtual TDataType GetValue() const
    {
        // this method is made virtual so, if required they can perform some thread independent
        // work here such as sorting, so we can optimize global sorting when we have locally thread sorted lists.
        return mReducedValue;
    }

    /// Destructor. Do nothing!!!
    virtual ~Reducer() {}

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType& value)
    {
        mF(mReducedValue, value);
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const Reducer<TDataType>& rOther)
    {
        const TDataType& r_other_value = rOther.GetValue();
#pragma omp critical
        {
            mF(mReducedValue, r_other_value);
        }
    }

private:
    const std::function<void(TDataType&, const TDataType&)> mF;
};

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
template <class TDataType>
class MaxReduction : public Reducer<TDataType>
{
public:
    typedef TDataType value_type;
    MaxReduction()
        : Reducer<TDataType>(std::numeric_limits<TDataType>::lowest(),
                             [](TDataType& rReducedValue, const TDataType& rValue) {
                                 rReducedValue = std::max(rReducedValue, rValue);
                             })
    {
    }
};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template <class TDataType>
class MinReduction : public Reducer<TDataType>
{
public:
    typedef TDataType value_type;
    MinReduction()
        : Reducer<TDataType>(std::numeric_limits<TDataType>::max(),
                             [](TDataType& rReducedValue, const TDataType& rValue) {
                                 rReducedValue = std::min(rReducedValue, rValue);
                             })
    {
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
