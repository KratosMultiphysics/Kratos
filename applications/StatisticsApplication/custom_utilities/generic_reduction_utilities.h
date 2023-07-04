//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#pragma once

// System includes
#include <tuple>

// Project includes
#include "utilities/atomic_utilities.h"
#include "utilities/parallel_utilities.h"

// Application incldues
#include "custom_utilities/data_type_traits.h"

namespace Kratos
{

namespace GenericReductionHelperUtilities
{

template<class TDataType>
void ResizeAndInitialize(
    TDataType& rLocalValue,
    const TDataType& rTargetValue,
    const double InitializationValue)
{
    if constexpr(std::is_same_v<TDataType, Vector>) {
        if (rLocalValue.size() != rTargetValue.size()) {
            rLocalValue = Vector(rTargetValue.size(), InitializationValue);
        }
    } else if constexpr(std::is_same_v<TDataType, Matrix>) {
        if (rLocalValue.size1() != rTargetValue.size1() || rLocalValue.size2() != rTargetValue.size2()) {
            rLocalValue = Matrix(rTargetValue.size1(), rTargetValue.size2(), InitializationValue);
        }
    }
}

template<class TDataType>
void ThreadSafeResizeAndInitialize(
    TDataType& rLocalValue,
    const TDataType& rTargetValue,
    const double InitializationValue)
{
    if constexpr(std::is_same_v<TDataType, Vector>) {
        KRATOS_CRITICAL_SECTION
        if (rLocalValue.size() != rTargetValue.size()) {
            rLocalValue = Vector(rTargetValue.size(), InitializationValue);
        }
    } else if constexpr(std::is_same_v<TDataType, Matrix>) {
        KRATOS_CRITICAL_SECTION
        if (rLocalValue.size1() != rTargetValue.size1() ||
            rLocalValue.size2() != rTargetValue.size2()) {
            rLocalValue = Matrix(rTargetValue.size1(), rTargetValue.size2(), InitializationValue);
        }
    }
}

} // namespace GenericReductionHelperUtilities

template<class TDataType>
class GenericSumReduction
{
public:
    using return_type = TDataType;

    static constexpr double mInitializationValue = 0.0;

    return_type mValue;

    /// Constructor
    GenericSumReduction()
    {
        DataTypeTraits<TDataType>::Initialize(mValue, mInitializationValue);
    }

    /// access to reduced value
    return_type GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const return_type& value){
        if (DataTypeTraits<TDataType>::Resize(mValue, value)) {
            DataTypeTraits<TDataType>::Initialize(mValue, mInitializationValue);
        }
        mValue += value;
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const GenericSumReduction<TDataType>& rOther)
    {
        {
            KRATOS_CRITICAL_SECTION
            if (DataTypeTraits<TDataType>::Resize(mValue, rOther.mValue)) {
                DataTypeTraits<TDataType>::Initialize(mValue, mInitializationValue);
            }
        }

        if constexpr(std::is_same_v<return_type, Vector>) {
            AtomicAddVector(mValue, rOther.mValue);
        } else if constexpr(std::is_same_v<return_type, Matrix>) {
            AtomicAddMatrix(mValue, rOther.mValue);
        } else {
            AtomicAdd(mValue, rOther.mValue);
        }

    }
};

template<class TDataType>
class GenericMinReduction
{
public:
    using data_type = TDataType;

    using indices_type = std::vector<unsigned int>;

    using return_type = std::tuple<TDataType, indices_type>;

    static constexpr double mInitializationValue = std::numeric_limits<double>::max();

    return_type mValue;

    /// Constructor
    GenericMinReduction()
    {
        auto& min_values = std::get<0>(mValue);
        auto& min_indices = std::get<1>(mValue);

        DataTypeTraits<TDataType>::Initialize(min_values, mInitializationValue);
        DataTypeTraits<indices_type>::Resize(min_indices, indices_type(DataTypeTraits<TDataType>::Size(min_values)));
        DataTypeTraits<indices_type>::Initialize(min_indices, std::numeric_limits<unsigned int>::max());
    }

    /// access to reduced value
    return_type GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const return_type& value){
        auto& min_values = std::get<0>(mValue);
        auto& min_indices = std::get<1>(mValue);

        const auto& values = std::get<0>(value);
        const auto& indices = std::get<1>(value);

        if (DataTypeTraits<TDataType>::Resize(min_values, values)) {
            // if resized, then initialize values
            DataTypeTraits<TDataType>::Initialize(min_values, mInitializationValue);
            // now initialize the indices vector size
            DataTypeTraits<indices_type>::Resize(min_indices, indices_type(DataTypeTraits<TDataType>::Size(min_values)));
            DataTypeTraits<indices_type>::Initialize(min_indices, std::numeric_limits<unsigned int>::max());
        }

        if constexpr(std::is_arithmetic_v<data_type>) {
            if (min_values > values) {
                min_values = values;
                min_indices = indices;
            }
        } else if constexpr(std::is_same_v<data_type, Matrix>) {
            for (IndexType i = 0; i < values.size1() * values.size2(); ++i) {
                if (min_values.data()[i] > values.data()[i]) {
                    min_values.data()[i] = value.data()[i];
                    min_indices[i] = indices[i];
                }
            }
        } else {
            for (IndexType i = 0; i < value.size(); ++i) {
                if (min_values[i] > values[i]) {
                    min_values[i] = value[i];
                    min_indices[i] = indices[i];
                }
            }
        }
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const GenericMinReduction<TDataType>& rOther)
    {
        KRATOS_CRITICAL_SECTION

        auto& min_values = std::get<0>(mValue);
        auto& min_indices = std::get<1>(mValue);

        const auto& values = std::get<0>(rOther.mValue);
        const auto& indices = std::get<1>(rOther.mValue);

        if (DataTypeTraits<TDataType>::Resize(min_values, values)) {
            // if resized, then initialize values
            DataTypeTraits<TDataType>::Initialize(min_values, mInitializationValue);
            // now initialize the indices vector size
            DataTypeTraits<indices_type>::Resize(min_indices, indices_type(DataTypeTraits<TDataType>::Size(min_values)));
            DataTypeTraits<indices_type>::Initialize(min_indices, std::numeric_limits<unsigned int>::max());
        }

        if constexpr(std::is_arithmetic_v<data_type>) {
            if (min_values > values) {
                min_values = values;
                min_indices = indices;
            }
        } else if constexpr(std::is_same_v<data_type, Matrix>) {
            for (IndexType i = 0; i < values.size1() * values.size2(); ++i) {
                if (min_values.data()[i] > values.data()[i]) {
                    min_values.data()[i] = values.data()[i];
                    min_indices[i] = indices[i];
                }
            }
        } else {
            for (IndexType i = 0; i < values.size(); ++i) {
                if (min_values[i] > values[i]) {
                    min_values[i] = values[i];
                    min_indices[i] = indices[i];
                }
            }
        }

    }
};

template<class TDataType>
class GenericMaxReduction
{
public:
    using data_type = TDataType;

    using indices_type = std::vector<unsigned int>;

    using return_type = std::tuple<TDataType, indices_type>;

    static constexpr double mInitializationValue = std::numeric_limits<double>::lowest();

    // initialize values and indices containers to min
    return_type mValue;

    /// Constructor
    GenericMaxReduction()
    {
        auto& max_values = std::get<0>(mValue);
        auto& max_indices = std::get<1>(mValue);

        DataTypeTraits<TDataType>::Initialize(max_values, mInitializationValue);
        DataTypeTraits<indices_type>::Resize(max_indices, indices_type(DataTypeTraits<TDataType>::Size(max_values)));
        DataTypeTraits<indices_type>::Initialize(max_indices, std::numeric_limits<unsigned int>::max());
    }

    /// access to reduced value
    return_type GetValue() const
    {
        return mValue;
    }

    /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
    void LocalReduce(const return_type& value){
        auto& max_values = std::get<0>(mValue);
        auto& max_indices = std::get<1>(mValue);

        const auto& values = std::get<0>(value);
        const auto& indices = std::get<1>(value);

        if (DataTypeTraits<TDataType>::Resize(max_values, values)) {
            // if resized, then initialize values
            DataTypeTraits<TDataType>::Initialize(max_values, mInitializationValue);
            // now initialize the indices vector size
            DataTypeTraits<indices_type>::Resize(max_indices, indices_type(DataTypeTraits<TDataType>::Size(max_values)));
            DataTypeTraits<indices_type>::Initialize(max_indices, std::numeric_limits<unsigned int>::max());
        }

        if constexpr(std::is_arithmetic_v<data_type>) {
            if (max_values < values) {
                max_values = values;
                max_indices = indices;
            }
        } else if constexpr(std::is_same_v<data_type, Matrix>) {
            for (IndexType i = 0; i < values.size1() * values.size2(); ++i) {
                if (max_values.data()[i] < values.data()[i]) {
                    max_values.data()[i] = values.data()[i];
                    max_indices[i] = indices[i];
                }
            }
        } else {
            for (IndexType i = 0; i < values.size(); ++i) {
                if (max_values[i] < values[i]) {
                    max_values[i] = values[i];
                    max_indices[i] = indices[i];
                }
            }
        }
    }

    /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
    void ThreadSafeReduce(const GenericMinReduction<TDataType>& rOther)
    {
        KRATOS_CRITICAL_SECTION

        auto& max_values = std::get<0>(mValue);
        auto& max_indices = std::get<1>(mValue);

        const auto& values = std::get<0>(rOther.mValue);
        const auto& indices = std::get<1>(rOther.mValue);

        if (DataTypeTraits<TDataType>::Resize(max_values, values)) {
            // if resized, then initialize values
            DataTypeTraits<TDataType>::Initialize(max_values, mInitializationValue);
            // now initialize the indices vector size
            DataTypeTraits<indices_type>::Resize(max_indices, indices_type(DataTypeTraits<TDataType>::Size(max_values)));
            DataTypeTraits<indices_type>::Initialize(max_indices, std::numeric_limits<unsigned int>::max());
        }

        if constexpr(std::is_arithmetic_v<data_type>) {
            if (max_values < values) {
                max_values = values;
                max_indices = indices;
            }
        } else if constexpr(std::is_same_v<data_type, Matrix>) {
            for (IndexType i = 0; i < values.size1() * values.size2(); ++i) {
                if (max_values.data()[i] < values.data()[i]) {
                    max_values.data()[i] = values.data()[i];
                    max_indices[i] = indices[i];
                }
            }
        } else {
            for (IndexType i = 0; i < values.size(); ++i) {
                if (max_values[i] < values[i]) {
                    max_values[i] = values[i];
                    max_indices[i] = indices[i];
                }
            }
        }

    }
};

} // namespace Kratos