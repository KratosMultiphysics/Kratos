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
//                   Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes
#include <tuple>
#include <limits>
#include <algorithm>
#include <mutex>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/atomic_utilities.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

namespace Internals
{
/** @brief Helper class for null-initializiation
 */
template <class TObjectType>
struct NullInitialized
{
    static TObjectType Get()
    {
        return TObjectType();
    }
};

template <class TValueType, std::size_t ArraySize>
struct NullInitialized<array_1d<TValueType,ArraySize>>
{
    static array_1d<TValueType,ArraySize> Get()
    {
        array_1d<TValueType,ArraySize> array;
        std::fill_n(array.begin(), ArraySize, NullInitialized<TValueType>::Get());
        return array;
    }
};
} // namespace Internals

///@addtogroup KratosCore

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

/** @brief utility function to do a sum reduction
 */
template<class TDataType, class TReturnType = TDataType>
class SumReduction
{
public:
    typedef TDataType   value_type;
    typedef TReturnType return_type;

    TReturnType mValue = Internals::NullInitialized<TReturnType>::Get(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TReturnType GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mValue += value;
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const SumReduction<TDataType, TReturnType>& rOther)
    {
        AtomicAdd(mValue, rOther.mValue);
    }
};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<class TDataType, class TReturnType = TDataType>
class SubReduction
{
public:
    typedef TDataType   value_type;
    typedef TReturnType return_type;

    TReturnType mValue = Internals::NullInitialized<TReturnType>::Get(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TReturnType GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mValue -= value;
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const SubReduction<TDataType, TReturnType>& rOther)
    {
        AtomicAdd(mValue, rOther.mValue);
    }
};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<class TDataType, class TReturnType = TDataType>
class MaxReduction
{
public:
    typedef TDataType   value_type;
    typedef TReturnType return_type;

    TReturnType mValue = std::numeric_limits<TReturnType>::lowest(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TReturnType GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mValue = std::max(mValue,value);
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const MaxReduction<TDataType, TReturnType>& rOther)
    {
        KRATOS_CRITICAL_SECTION
        LocalReduce(rOther.mValue);
    }
};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<class TDataType, class TReturnType = TDataType>
class AbsMaxReduction
{
public:
    typedef TDataType   value_type;
    typedef TReturnType return_type;

    TReturnType mValue = std::numeric_limits<TReturnType>::lowest(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TReturnType GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mValue = (std::abs(mValue) < std::abs(value)) ? value : mValue;
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const AbsMaxReduction<TDataType, TReturnType>& rOther)
    {
        KRATOS_CRITICAL_SECTION
        LocalReduce(rOther.mValue);
    }
};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************
template<class TDataType, class TReturnType = TDataType>
class MinReduction
{
public:
    typedef TDataType   value_type;
    typedef TReturnType return_type;

    TReturnType mValue = std::numeric_limits<TReturnType>::max(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TReturnType GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mValue = std::min(mValue,value);
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const MinReduction<TDataType, TReturnType>& rOther)
    {
        KRATOS_CRITICAL_SECTION
        LocalReduce(rOther.mValue);
    }
};


//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

template<class TDataType, class TReturnType = TDataType>
class AbsMinReduction
{
public:
    typedef TDataType   value_type;
    typedef TReturnType return_type;

    TReturnType mValue = std::numeric_limits<TReturnType>::max(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TReturnType GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mValue = (std::abs(mValue) < std::abs(value)) ? mValue : value;
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const AbsMinReduction<TDataType, TReturnType>& rOther)
    {
        KRATOS_CRITICAL_SECTION
        LocalReduce(rOther.mValue);
    }
};

//***********************************************************************************
//***********************************************************************************
//***********************************************************************************

template<class TDataType, class TReturnType = std::vector<TDataType>>
class AccumReduction
{
public:
    typedef TDataType   value_type;
    typedef TReturnType return_type;

    TReturnType mValue = TReturnType(); // deliberately making the member value public, to allow one to change it as needed

    /// access to reduced value
    TReturnType GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const TDataType value){
        mValue.insert(mValue.end(), value);
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const AccumReduction<TDataType, TReturnType>& rOther)
    {
        KRATOS_CRITICAL_SECTION
        std::copy(rOther.mValue.begin(), rOther.mValue.end(), std::inserter(mValue, mValue.end()));
    }
};

template<class MapType>
class MapReduction
{
public:
    using value_type = typename MapType::value_type;
    using return_type = MapType;

    return_type mValue;

    /// access to reduced value
    return_type GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const value_type rValue){
        mValue.emplace(rValue);
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(MapReduction<MapType>& rOther)
    {
        KRATOS_CRITICAL_SECTION
        mValue.merge(rOther.mValue);
    }
};

template <class... Reducer>
struct CombinedReduction {
    typedef std::tuple<typename Reducer::value_type...> value_type;
    typedef std::tuple<typename Reducer::return_type...> return_type;

    std::tuple<Reducer...> mChild;

    CombinedReduction() {}

    /// access to reduced value
    return_type GetValue(){
        return_type return_value;
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
