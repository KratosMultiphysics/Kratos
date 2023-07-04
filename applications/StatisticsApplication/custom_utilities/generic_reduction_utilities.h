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
#include "includes/data_communicator.h"
#include "utilities/atomic_utilities.h"
#include "utilities/parallel_utilities.h"

// Application incldues
#include "custom_utilities/data_type_traits.h"

namespace Kratos
{

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

    static void Synchronize(
        return_type& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        if (DataTypeTraits<TDataType>::SynchronizeSize(rValue, rDataCommunicator)) {
            DataTypeTraits<TDataType>::Initialize(rValue, mInitializationValue);
        }

        typename DataTypeTraits<TDataType>::VectorType local_values, global_values;
        DataTypeTraits<TDataType>::FillToVector(local_values, rValue);
        rDataCommunicator.SumAll(local_values, global_values);
        DataTypeTraits<TDataType>::FillFromVector(rValue, global_values);
    }
};

template<class TDataType>
class GenericMinReduction
{
public:
    using data_type = TDataType;

    using indices_type = std::vector<int>;

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
            DataTypeTraits<indices_type>::Initialize(min_indices, std::numeric_limits<int>::max());
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
            DataTypeTraits<indices_type>::Initialize(min_indices, std::numeric_limits<int>::max());
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

    static void Synchronize(
        return_type& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        auto& min_values = std::get<0>(rValue);
        auto& min_indices = std::get<1>(rValue);

        if (DataTypeTraits<TDataType>::SynchronizeSize(min_values, rDataCommunicator)) {
            DataTypeTraits<TDataType>::Initialize(min_values, mInitializationValue);
            DataTypeTraits<indices_type>::Resize(min_indices, indices_type(DataTypeTraits<TDataType>::Size(min_values)));
            DataTypeTraits<indices_type>::Initialize(min_indices, std::numeric_limits<int>::max());
        }

        const int world_size = rDataCommunicator.Size();

        // first get the values from all ranks
        typename DataTypeTraits<TDataType>::VectorType local_values;
        DataTypeTraits<TDataType>::FillToVector(local_values, min_values);
        std::vector<typename DataTypeTraits<TDataType>::VectorType> global_values(world_size);
        global_values = rDataCommunicator.AllGatherv(local_values);

        // now get indices corresponding to above values from all ranks
        std::vector<indices_type> global_indices(world_size);
        global_indices = rDataCommunicator.AllGatherv(min_indices);

        // now get the min
        TDataType current_values;
        DataTypeTraits<TDataType>::Resize(current_values, min_values);
        for (int rank = 0; rank < world_size; ++rank) {
            const auto& rank_values = global_values[rank];
            const auto& rank_indices = global_indices[rank];
            DataTypeTraits<TDataType>::FillFromVector(current_values, rank_values);

            for (IndexType i = 0; i < DataTypeTraits<TDataType>::Size(rank_values); ++i) {
                if (DataTypeTraits<TDataType>::GetComponent(current_values, i) < DataTypeTraits<TDataType>::GetComponent(min_values, i)) {
                    DataTypeTraits<TDataType>::GetComponent(min_values, i) = DataTypeTraits<TDataType>::GetComponent(current_values, i);
                    DataTypeTraits<TDataType>::GetComponent(min_indices, i) = rank_indices[i];
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

    using indices_type = std::vector<int>;

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

    static void Synchronize(
        return_type& rValue,
        const DataCommunicator& rDataCommunicator)
    {
        auto& max_values = std::get<0>(rValue);
        auto& max_indices = std::get<1>(rValue);

        if (DataTypeTraits<TDataType>::SynchronizeSize(max_values, rDataCommunicator)) {
            DataTypeTraits<TDataType>::Initialize(max_values, mInitializationValue);
            DataTypeTraits<indices_type>::Resize(max_indices, indices_type(DataTypeTraits<TDataType>::Size(max_values)));
            DataTypeTraits<indices_type>::Initialize(max_indices, std::numeric_limits<unsigned int>::max());
        }

        const int world_size = rDataCommunicator.Size();

        // first get the values from all ranks
        typename DataTypeTraits<TDataType>::VectorType local_values;
        DataTypeTraits<TDataType>::FillToVector(local_values, max_values);
        std::vector<typename DataTypeTraits<TDataType>::VectorType> global_values(world_size);
        global_values = rDataCommunicator.AllGatherv(local_values);

        // now get indices corresponding to above values from all ranks
        std::vector<indices_type> global_indices(world_size);
        global_indices = rDataCommunicator.AllGatherv(max_indices);

        // now get the min
        TDataType current_values;
        DataTypeTraits<TDataType>::Resize(current_values, max_values);
        for (int rank = 0; rank < world_size; ++rank) {
            const auto& rank_values = global_values[rank];
            const auto& rank_indices = global_indices[rank];
            DataTypeTraits<TDataType>::FillFromVector(current_values, rank_values);

            for (IndexType i = 0; i < DataTypeTraits<TDataType>::Size(rank_values); ++i) {
                if (DataTypeTraits<TDataType>::GetComponent(current_values, i) > DataTypeTraits<TDataType>::GetComponent(max_values, i)) {
                    DataTypeTraits<TDataType>::GetComponent(max_values, i) = DataTypeTraits<TDataType>::GetComponent(current_values, i);
                    DataTypeTraits<TDataType>::GetComponent(max_indices, i) = rank_indices[i];
                }
            }
        }
    }
};

} // namespace Kratos