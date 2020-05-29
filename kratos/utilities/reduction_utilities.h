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

    void LocalReduce(const TDataType value){
        mvalue += value;
    }

    //a user could define a ThreadSafeReduce for his specific reduction case
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
    TDataType mvalue = TDataType(); //i am deliberately making the member value public, to allow one to change it as needed

    TDataType GetValue() const
    {
        return mvalue;
    }
    void LocalReduce(const TDataType value){
        mvalue -= value;
    }

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
    TDataType mvalue = std::numeric_limits<TDataType>::lowest(); //i am deliberately making the member value public, to allow one to change it as needed
    TDataType GetValue() const
    {
        return mvalue;
    }
    void LocalReduce(const TDataType value){
        mvalue = std::max(mvalue,value);
    }
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
    TDataType mvalue = std::numeric_limits<TDataType>::max(); //i am deliberately making the member value public, to allow one to change it as needed

    TDataType GetValue() const
    {
        return mvalue;
    }
    void LocalReduce(const TDataType value){
        mvalue = std::min(mvalue,value);
    }

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
    void LocalReduce(const std::tuple<T...> &&v) {
        // Static recursive loop over tuple elements
        reduce_local<0>(v);
    }

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

