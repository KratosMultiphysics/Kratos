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

// System includes
#include <algorithm>
#include <cmath>
#include <numeric>
#include <tuple>
#include <type_traits>
#include <variant>
#include <vector>

// External includes

// Project includes
#include "includes/communicator.h"
#include "includes/define.h"
#include "includes/global_variables.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/data_containers.h"
#include "custom_utilities/generic_reduction_utilities.h"
#include "custom_utilities/method_utilities.h"
#include "custom_utilities/norms.h"

// Include base h
#include "spatial_methods.h"

namespace Kratos
{

namespace SpatialMethodHelperUtilities
{

using IndexType = unsigned int;

using IndicesType = std::vector<IndexType>;

using DataLocation = Globals::DataLocation;

IndexType GetDataLocationSize(
    const ModelPart& rModelPart,
    const DataLocation& rLocation)
{
    switch (rLocation) {
        case DataLocation::NodeHistorical:
        case DataLocation::NodeNonHistorical:
            return rModelPart.GetCommunicator().GlobalNumberOfNodes();
        case DataLocation::Condition:
            return rModelPart.GetCommunicator().GlobalNumberOfConditions();
        case DataLocation::Element:
            return rModelPart.GetCommunicator().GlobalNumberOfElements();
        default:
            KRATOS_ERROR << "Unsupported data location requested. Only supprts node/conditions/elements.";
            return 0;
    }
}

template<class TDataType>
class SumOperation
{
public:
    ///@name Type definitions
    ///@{

    using OperationTraits = DataTypeTraits<TDataType>;

    ///@}
    ///@name Life cycle
    ///@{

    SumOperation() { Initialize(); }

    SumOperation(const TDataType& rValue) { mValue = rValue; }

    ///@}
    ///@name Public operations
    ///@{

    TDataType GetValue() const { return mValue; }

    void Execute(const SumOperation<TDataType>& rOther)
    {
        if constexpr(OperationTraits::HasDynamicMemoryAllocation){
            if (!mIsInitialized) {
                mValue = rOther.mValue;
                Initialize();
                mIsInitialized = true;
            }
        }

        mValue += rOther.mValue;
    }

    void Synchronize(const DataCommunicator& rDataCommunicator)
    {
        if (rDataCommunicator.SynchronizeShape(mValue)) {
            Initialize();
        }

        mValue = rDataCommunicator.SumAll(mValue);
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    TDataType mValue;

    bool mIsInitialized = !OperationTraits::HasDynamicMemoryAllocation;

    ///@}
    ///@name Private operations
    ///@{

    void Initialize() {
        OperationTraits::Initialize(mValue, 0.0);
    }

    ///@}
};

template<class TDataType, class TSubOperation>
class MinMaxOperation
{
public:
    ///@name Type definitions
    ///@{

    static constexpr bool IsPrimitiveType = std::is_arithmetic_v<TDataType>;

    using OperationTraits = DataTypeTraits<TDataType>;

    using IndexType = unsigned int;

    using IndicesType = std::vector<IndexType>;

    using IndexReturnType = std::conditional_t<IsPrimitiveType, IndexType, IndicesType>;

    ///@}
    ///@name Life cycle
    ///@{

    MinMaxOperation() { Initialize(); }

    MinMaxOperation(
        const TDataType& rValue,
        const IndexType rId)
    {
        mValue = rValue;
        mIndices.resize(OperationTraits::Size(mValue), rId);
    }

    ///@}
    ///@name Public operations
    ///@{

    std::tuple<TDataType, IndexReturnType> GetValue() const
    {
        if constexpr(IsPrimitiveType) {
            return std::make_tuple(mValue, mIndices[0]);
        } else {
            return std::make_tuple(mValue, mIndices);
        }
    }

    void Execute(const MinMaxOperation<TDataType, TSubOperation>& rOther)
    {
        if constexpr(OperationTraits::HasDynamicMemoryAllocation){
            if (!mIsInitialized) {
                mValue = rOther.mValue;
                Initialize();
                mIsInitialized = true;
            }
        }

        for (IndexType i_comp = 0; i_comp < OperationTraits::Size(mValue); ++i_comp) {
            const auto other_value = OperationTraits::GetComponent(rOther.mValue, i_comp);
            const auto other_index = rOther.mIndices[i_comp];
            auto& current_value = OperationTraits::GetComponent(mValue, i_comp);
            auto& current_index = mIndices[i_comp];
            if (TSubOperation::IsCurrentValueNotValid(current_value, other_value)) {
                current_value = other_value;
                current_index = other_index;
            } else if (current_value == other_value) {
                current_index = std::min(current_index, other_index);
            }
        }
    }

    void Synchronize(const DataCommunicator& rDataCommunicator)
    {
        if (rDataCommunicator.SynchronizeShape(mValue)) {
            Initialize();
        }

        const auto& all_gather_values = rDataCommunicator.AllGatherv(std::vector<TDataType>{mValue});
        const auto& all_gather_indices = rDataCommunicator.AllGatherv(mIndices);

        for (IndexType i_comp = 0; i_comp < OperationTraits::Size(mValue); ++i_comp) {
            auto& current_value = OperationTraits::GetComponent(mValue, i_comp);
            auto& current_index = mIndices[i_comp];

            for (IndexType i_rank = 0; i_rank < all_gather_values.size(); ++i_rank) {
                const auto& other_value = OperationTraits::GetComponent(all_gather_values[i_rank][0], i_comp);
                const auto& other_index = all_gather_indices[i_rank][i_comp];

                if (TSubOperation::IsCurrentValueNotValid(current_value, other_value)) {
                    current_value = other_value;
                    current_index = other_index;
                } else if (current_value == other_value) {
                    current_index = std::min(current_index, other_index);
                }
            }
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    TDataType mValue;

    IndicesType mIndices;

    bool mIsInitialized = !OperationTraits::HasDynamicMemoryAllocation;

    ///@}
    ///@name Private operations
    ///@{

    void Initialize()
    {
        OperationTraits::Initialize(mValue, TSubOperation::InitialValue);
        mIndices.resize(OperationTraits::Size(mValue), std::numeric_limits<IndexType>::max());
    }

    ///@}
};

template<class TDataType>
struct MinSubOperation
{
    using RawDataType = typename DataTypeTraits<TDataType>::RawDataType;

    static constexpr RawDataType InitialValue = std::numeric_limits<RawDataType>::max();

    static bool IsCurrentValueNotValid(
        const RawDataType rCurrentValue,
        const RawDataType rOtherValue)
    {
        return rCurrentValue > rOtherValue;
    }
};

template<class TDataType>
struct MaxSubOperation
{
    using RawDataType = typename DataTypeTraits<TDataType>::RawDataType;

    static constexpr RawDataType InitialValue = std::numeric_limits<RawDataType>::lowest();

    static bool IsCurrentValueNotValid(
        const RawDataType rCurrentValue,
        const RawDataType rOtherValue)
    {
        return rCurrentValue < rOtherValue;
    }
};

template<class TDataType> class MinOperation: public MinMaxOperation<TDataType, MinSubOperation<TDataType>> { public: using MinMaxOperation<TDataType, MinSubOperation<TDataType>>::MinMaxOperation; };
template<class TDataType> class MaxOperation: public MinMaxOperation<TDataType, MaxSubOperation<TDataType>> { public: using MinMaxOperation<TDataType, MaxSubOperation<TDataType>>::MinMaxOperation; };

template<class TDataType>
class MedianOperation
{
public:
    ///@name Type definitions
    ///@{

    static constexpr bool IsPrimitiveType = std::is_arithmetic_v<TDataType>;

    using OperationTraits = DataTypeTraits<TDataType>;

    using IndexType = unsigned int;

    using IndicesType = std::vector<IndexType>;

    using IndexReturnType = std::conditional_t<IsPrimitiveType, IndexType, IndicesType>;

    using RawDataType = typename OperationTraits::RawDataType;

    using DataValuesType = std::vector<std::tuple<RawDataType, IndexType>>;

    ///@}
    ///@name Life cycle
    ///@{

    MedianOperation()
    {
        const IndexType data_size = OperationTraits::Size(TDataType{});
        mValues.resize(data_size);
        mResultantIndex.resize(data_size);
    }

    MedianOperation(
        const TDataType& rValue,
        const IndexType rId)
    {
        mResultantValue = rValue;

        const IndexType data_size = OperationTraits::Size(rValue);
        if (mValues.size() != data_size) {
            mValues.resize(data_size);
            mResultantIndex.resize(data_size);
        }

        for (IndexType i = 0; i < data_size; ++i) {
            mValues[i].push_back(std::make_tuple(OperationTraits::GetComponent(rValue, i), rId));
        }

    }

    ///@}
    ///@name Public operations
    ///@{

    std::tuple<TDataType, IndexReturnType> GetValue() const
    {
        if constexpr(IsPrimitiveType) {
            return std::make_tuple(mResultantValue, mResultantIndex[0]);
        } else {
            return std::make_tuple(mResultantValue, mResultantIndex);
        }
    }

    void Execute(const MedianOperation<TDataType>& rOther)
    {
        // first check the sizes. Requires resizing if dynamic types
        // such as Vector and Matrices are used.
        if (mValues.size() != rOther.mValues.size()) {
            mResultantValue = rOther.mResultantValue;
            mValues.resize(rOther.mValues.size());
            mResultantIndex.resize(rOther.mValues.size());
        }

        // iterate through each component
        for (IndexType i_comp = 0; i_comp < rOther.mValues.size(); ++i_comp) {
            auto& current_component_values = mValues[i_comp];
            const auto& r_other_component_values = rOther.mValues[i_comp];

            for (const auto& r_other_value_info : r_other_component_values) {
                current_component_values.push_back(r_other_value_info);
            }
        }
    }

    void Synchronize(const DataCommunicator& rDataCommunicator)
    {
        rDataCommunicator.SynchronizeShape(mResultantValue);

        const IndexType data_size = OperationTraits::Size(mResultantValue);

        if (mValues.size() != data_size) {
            mValues.resize(data_size);
            mResultantIndex.resize(data_size);
        }

        for (IndexType i_comp = 0; i_comp < data_size; ++i_comp) {
            auto& current_median_value = OperationTraits::GetComponent(mResultantValue, i_comp);
            auto& current_median_index = mResultantIndex[i_comp];
            const auto& current_values = mValues[i_comp];

            // get the values in rank 0
            std::vector<RawDataType> local_values(current_values.size());
            std::transform(current_values.begin(), current_values.end(), local_values.begin(), [](const auto& rV) { return std::get<0>(rV); });
            const auto& global_values = rDataCommunicator.Gatherv(local_values, 0);

            // get the indices in rank 0
            IndicesType local_indices(current_values.size());
            std::transform(current_values.begin(), current_values.end(), local_indices.begin(), [](const auto& rV) { return std::get<1>(rV); });
            const auto& global_indices = rDataCommunicator.Gatherv(local_indices, 0);

            if (rDataCommunicator.Rank() == 0) {
                const IndexType number_of_values = std::accumulate(global_values.begin(), global_values.end(), 0, [](const auto& rV1, const auto& rV2) { return rV1 + rV2.size();});
                std::vector<std::tuple<RawDataType, IndexType>> global_values_info;
                global_values_info.reserve(number_of_values);

                for (IndexType rank = 0; rank < global_values.size(); ++rank) {
                    const auto& rank_values = global_values[rank];
                    const auto& rank_indices = global_indices[rank];
                    for (IndexType i_value = 0; i_value < rank_values.size(); ++i_value) {
                        global_values_info.push_back(std::make_tuple(rank_values[i_value], rank_indices[i_value]));
                    }
                }

                // now sort everything
                std::sort(global_values_info.begin(), global_values_info.end(), [](const auto& rV1, const auto& rV2) { return std::get<0>(rV1) < std::get<0>(rV2); });

                if (number_of_values > 0) {
                    const IndexType mid_point = number_of_values / 2;
                    current_median_index = std::get<1>(global_values_info[mid_point]);
                    if (number_of_values % 2 != 0) {
                        current_median_value = std::get<0>(global_values_info[mid_point]);
                    } else {
                        const IndexType adjacent_mid_point = (number_of_values - 1) / 2;
                        current_median_value = (std::get<0>(global_values_info[adjacent_mid_point]) + std::get<0>(global_values_info[mid_point])) * 0.5;
                    }
                }
            }

            rDataCommunicator.Broadcast(current_median_value, 0);
            rDataCommunicator.Broadcast(current_median_index, 0);
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    // this is storing the sorted items
    std::vector<DataValuesType> mValues;

    TDataType mResultantValue;

    IndicesType mResultantIndex;

    ///@}
};

template<class TDataType, class TNormType, int TPower = 1>
typename TNormType::ResultantValueType<TDataType> GenericSumReduction(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    const TNormType& rNorm)
{
    KRATOS_TRY

    const auto data_container = DataContainers::GetDataContainer(rModelPart, rVariable, rLocation);

    return std::visit([&rModelPart, &rNorm](auto& rDataContainer) {
        using data_container_type = std::decay_t<decltype(rDataContainer)>;
        return GenericReductionUtilities::GenericReduction<data_container_type, TNormType, SumOperation, false, TPower>(
                rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, rNorm)
            .GetValue();
    }, data_container);

    KRATOS_CATCH("");
}

template<class TDataType, class TNormType, template <class T1> class OperationType>
std::tuple<typename TNormType::ResultantValueType<TDataType>, typename SpatialMethods::ItemPositionType<typename TNormType::ResultantValueType<TDataType>>> GenericReductionWithIndices(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    const TNormType& rNorm)
{
    KRATOS_TRY

    const auto data_container = DataContainers::GetDataContainer(rModelPart, rVariable, rLocation);

    return std::visit([&rModelPart, &rNorm](auto& rDataContainer) {
        using data_container_type = std::decay_t<decltype(rDataContainer)>;
        using norm_type = std::decay_t<decltype(rNorm)>;
        return GenericReductionUtilities::GenericReduction<data_container_type, norm_type, OperationType, true>(
                rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, rNorm)
            .GetValue();
    }, data_container);

    KRATOS_CATCH("");
}

template <class TDataType, class TNormType>
SpatialMethods::DistributionInfo<typename TNormType::ResultantValueType<TDataType>> GenericDistribution(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    Parameters Params,
    const TNormType& rNorm)
{
    KRATOS_TRY

    using norm_return_type = typename TNormType::ResultantValueType<TDataType>;

    using distribution_info_type = SpatialMethods::DistributionInfo<norm_return_type>;

    using data_type_traits = DataTypeTraits<norm_return_type>;

    using ReturnIndicesType = typename SpatialMethods::DistributionInfo<norm_return_type>::IndicesType;

    const auto data_container = DataContainers::GetDataContainer(rModelPart, rVariable, rLocation);

    return std::visit(
        [&](auto& rDataContainer) {
            using data_container_type = std::decay_t<decltype(rDataContainer)>;

            Parameters default_parameters = Parameters(R"(
            {
                "number_of_value_groups" : 10,
                "min_value"              : "min",
                "max_value"              : "max"
            })");

            // first fix the min value
            if (Params.Has("min_value") && Params["min_value"].IsDouble()) {
                default_parameters["min_value"].SetDouble(0.0);
            }
            if (Params.Has("max_value") && Params["max_value"].IsDouble()) {
                default_parameters["max_value"].SetDouble(0.0);
            }
            Params.RecursivelyValidateAndAssignDefaults(default_parameters);

            distribution_info_type distribution_info;

            auto& min_value = distribution_info.mMin;
            if (Params["min_value"].IsDouble()) {
                data_type_traits::Initialize(min_value, Params["min_value"].GetDouble());
            } else if (Params["min_value"].IsString() && Params["min_value"].GetString() == "min") {
                min_value = std::get<0>(GenericReductionUtilities::GenericReduction<data_container_type, TNormType, MinOperation, true>(rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, rNorm).GetValue());
            } else {
                KRATOS_ERROR
                    << "Unknown min_value. Allowed only double or \"min\" "
                       "string as a value. [ min_value = "
                    << Params["min_value"] << " ]\n.";
            }

            // first fix the max value
            auto& max_value = distribution_info.mMax;
            if (Params["max_value"].IsDouble()) {
                data_type_traits::Initialize(max_value, Params["max_value"].GetDouble());
            } else if (Params["max_value"].IsString() && Params["max_value"].GetString() == "max") {
                max_value = std::get<0>(GenericReductionUtilities::GenericReduction<data_container_type, TNormType, MaxOperation, true>(rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, rNorm).GetValue());
            } else {
                KRATOS_ERROR
                    << "Unknown max_value. Allowed only double or \"max\" "
                       "string as a value. [ max_value = "
                    << Params["max_value"] << " ]\n.";
            }

            auto& group_limits = distribution_info.mGroupUpperValues;
            const IndexType number_of_groups = Params["number_of_value_groups"].GetInt();

            // we need additional two groups to store values below the specified
            // minimum and values above the specified maximum.
            const IndexType number_of_all_groups = number_of_groups + 2;
            group_limits.resize(number_of_all_groups);
            for (IndexType i = 0; i < number_of_groups + 1; ++i) {
                group_limits[i] = min_value + (max_value - min_value) * static_cast<double>(i) / static_cast<double>(number_of_groups);
            }
            norm_return_type additional_max_value(max_value);
            data_type_traits::Initialize(additional_max_value, std::numeric_limits<typename DataTypeTraits<norm_return_type>::RawDataType>::max());
            group_limits.back() = additional_max_value;

            /// reduction class
            class DistributionReduction {
            public:
                using data_type = std::tuple<IndexType, IndicesType, norm_return_type>;

                using return_type = std::tuple<
                                        std::vector<ReturnIndicesType>, // group value counts
                                        std::vector<norm_return_type>, // group means
                                        std::vector<norm_return_type>  // group variances
                                    >;

                return_type mValue;

                /// access to reduced value
                return_type GetValue() const
                {
                    return mValue;
                }

                /// NON-THREADSAFE (fast) value of reduction, to be used within a single thread
                void LocalReduce(const data_type& rValue)
                {
                    ResizeAndInitialize(rValue);

                    auto& r_current_group_counts = std::get<0>(mValue);
                    auto& r_current_group_means = std::get<1>(mValue);
                    auto& r_current_group_variances = std::get<2>(mValue);

                    const auto& r_group_indices = std::get<1>(rValue);
                    const auto& r_value = std::get<2>(rValue);
                    const IndexType number_of_components = DataTypeTraits<norm_return_type>::Size(r_value);

                    for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                        const auto group_index = r_group_indices[i_comp];
                        const auto raw_value = DataTypeTraits<norm_return_type>::GetComponent(r_value, i_comp);
                        ++DataTypeTraits<ReturnIndicesType>::GetComponent(r_current_group_counts[group_index], i_comp);
                        DataTypeTraits<norm_return_type>::GetComponent(r_current_group_means[group_index], i_comp) += raw_value;
                        DataTypeTraits<norm_return_type>::GetComponent(r_current_group_variances[group_index], i_comp) += std::pow(raw_value, 2);
                    }
                }

                /// THREADSAFE (needs some sort of lock guard) reduction, to be used to sync threads
                void ThreadSafeReduce(const DistributionReduction& rOther)
                {
                    KRATOS_CRITICAL_SECTION

                    auto& r_group_counts = std::get<0>(mValue);
                    auto& r_group_means = std::get<1>(mValue);
                    auto& r_group_variances = std::get<2>(mValue);

                    const auto& r_other_group_counts = std::get<0>(rOther.mValue);
                    const auto& r_other_group_means = std::get<1>(rOther.mValue);
                    const auto& r_other_group_variances = std::get<2>(rOther.mValue);

                    const IndexType number_of_groups = r_other_group_counts.size();

                    // first check and resize everything
                    if (r_group_counts.size() != number_of_groups) {
                        r_group_counts.resize(number_of_groups);
                        r_group_means.resize(number_of_groups);
                        r_group_variances.resize(number_of_groups);
                    }

                    for (IndexType i_group = 0; i_group < number_of_groups; ++i_group) {
                        auto& r_group_count = r_group_counts[i_group];
                        auto& r_group_mean = r_group_means[i_group];
                        auto& r_group_variance = r_group_variances[i_group];

                        const auto& r_other_group_count = r_other_group_counts[i_group];
                        const auto& r_other_group_mean = r_other_group_means[i_group];
                        const auto& r_other_group_variance =
                            r_other_group_variances[i_group];

                        const IndexType number_of_components = DataTypeTraits<ReturnIndicesType>::Size(r_other_group_count);

                        if (DataTypeTraits<ReturnIndicesType>::Size(r_group_count) != number_of_components) {
                            r_group_count = r_other_group_count;
                            r_group_mean = r_other_group_mean;
                            r_group_variance = r_other_group_variance;
                            DataTypeTraits<ReturnIndicesType>::Initialize(r_group_count, 0);
                            data_type_traits::Initialize(r_group_mean, 0.0);
                            data_type_traits::Initialize(r_group_variance, 0.0);
                        }

                        for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                            DataTypeTraits<ReturnIndicesType>::GetComponent(r_group_count, i_comp) += DataTypeTraits<ReturnIndicesType>::GetComponent(r_other_group_count, i_comp);
                        }
                        r_group_mean += r_other_group_mean;
                        r_group_variance += r_other_group_variance;
                    }
                }

                void ResizeAndInitialize(const data_type& rValue)
                {
                    const auto number_of_groups = std::get<0>(rValue);
                    const auto& r_value = std::get<2>(rValue);
                    const auto number_of_components = DataTypeTraits<norm_return_type>::Size(r_value);

                    auto& r_current_group_counts = std::get<0>(mValue);
                    auto& r_current_group_means = std::get<1>(mValue);
                    auto& r_current_group_variances = std::get<2>(mValue);

                    if (r_current_group_counts.size() != number_of_groups) {
                        // resize
                        r_current_group_counts.resize(number_of_groups);
                        r_current_group_means.resize(number_of_groups);
                        r_current_group_variances.resize(number_of_groups);

                        // initialize
                        for (IndexType i_group = 0; i_group < number_of_groups; ++i_group) {
                            if constexpr(std::is_arithmetic_v<norm_return_type>) {
                                r_current_group_counts[i_group] = 0;
                            } else {
                                r_current_group_counts[i_group].resize(number_of_components);
                                std::fill(r_current_group_counts[i_group].begin(), r_current_group_counts[i_group].end(), 0);
                            }
                            r_current_group_means[i_group] = r_value;
                            r_current_group_variances[i_group] = r_value;
                            data_type_traits::Initialize(r_current_group_means[i_group], 0.0);
                            data_type_traits::Initialize(r_current_group_variances[i_group], 0.0);
                        }
                    }
                }
            };

            const IndexType number_of_components = DataTypeTraits<norm_return_type>::Size(max_value);

            const auto& reuduced_values =
                IndexPartition<IndexType>(rDataContainer.Size()).for_each<DistributionReduction>(TDataType{}, [&rDataContainer, &rNorm, &group_limits, number_of_components, number_of_groups](const IndexType Index, TDataType& rTLS) {
                    rDataContainer.GetValue(rTLS, Index);
                    auto norm_value = rNorm.Evaluate(rTLS);

                    IndicesType indices(number_of_components);
                    for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                        const auto& comp_value = DataTypeTraits<norm_return_type>::GetComponent(norm_value, i_comp);
                        IndexType i_group;
                        for (i_group = 0; i_group < group_limits.size() - 2; ++i_group) {
                            if (comp_value < DataTypeTraits<norm_return_type>::GetComponent(group_limits[i_group], i_comp)) {
                                indices[i_comp] = i_group;
                                i_group = std::numeric_limits<IndexType>::max();
                                break;
                            }
                        }

                        for (; i_group < group_limits.size(); ++i_group) {
                            if (comp_value <= DataTypeTraits<norm_return_type>::GetComponent(group_limits[i_group], i_comp)) {
                                indices[i_comp] = i_group;
                                break;
                            }
                        }
                    }

                    return std::make_tuple(group_limits.size(), indices, norm_value);
                });

            // now do the mpi communicaton
            const auto& r_data_communicator = rModelPart.GetCommunicator().GetDataCommunicator();
            distribution_info.mGroupNumberOfValues.resize(number_of_all_groups);
            double number_of_items = 0;
            for (IndexType i = 0; i < number_of_all_groups; ++i) {
                const auto& group_values = std::get<0>(reuduced_values)[i];
                distribution_info.mGroupNumberOfValues[i] = r_data_communicator.SumAll(group_values);
                number_of_items += DataTypeTraits<ReturnIndicesType>::GetComponent(group_values, 0);
            }
            number_of_items = std::max(r_data_communicator.SumAll(number_of_items), 1.0);
            distribution_info.mGroupMean     = r_data_communicator.SumAll(std::get<1>(reuduced_values));
            distribution_info.mGroupVariance = r_data_communicator.SumAll(std::get<2>(reuduced_values));

            // now revert back to the group values
            distribution_info.mGroupValueDistributionPercentage.resize(number_of_all_groups, distribution_info.mGroupMean[0]);
            for (IndexType i_group = 0; i_group < number_of_all_groups; ++i_group) {
                auto& current_number_of_values = distribution_info.mGroupNumberOfValues[i_group];
                auto& current_distribution_percentage = distribution_info.mGroupValueDistributionPercentage[i_group];
                auto& current_mean = distribution_info.mGroupMean[i_group];
                auto& current_variance = distribution_info.mGroupVariance[i_group];

                // post processing of values
                for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                    const auto n = DataTypeTraits<ReturnIndicesType>::GetComponent(current_number_of_values, i_comp);
                    if (n > 0) {
                        DataTypeTraits<norm_return_type>::GetComponent(current_distribution_percentage, i_comp) = n / number_of_items;
                        DataTypeTraits<norm_return_type>::GetComponent(current_mean, i_comp) /= n;
                        DataTypeTraits<norm_return_type>::GetComponent(current_variance, i_comp) /= n;
                    } else {
                        DataTypeTraits<norm_return_type>::GetComponent(current_distribution_percentage, i_comp) = 0.0;
                        DataTypeTraits<norm_return_type>::GetComponent(current_mean, i_comp) = 0.0;
                        DataTypeTraits<norm_return_type>::GetComponent(current_variance, i_comp) = 0.0;
                    }
                    DataTypeTraits<norm_return_type>::GetComponent(current_variance, i_comp) -= std::pow(DataTypeTraits<norm_return_type>::GetComponent(current_mean, i_comp), 2);
                }
            }

            // reversing group limit extention
            group_limits[group_limits.size() - 1] = max_value;

            return distribution_info;
        }, data_container);

    KRATOS_CATCH("");
}


} // namespace SpatialMethodHelperUtilities

SpatialMethods::IndexType SpatialMethods::Sum(
    const ModelPart& rModelPart,
    const Flags& rFlag,
    const DataLocation& rLocation)
{
    const auto data_container = DataContainers::GetDataContainer(rModelPart, rFlag, rLocation);

    return std::visit([&rModelPart](auto& rDataContainer) {
        using data_container_type = std::decay_t<decltype(rDataContainer)>;
        return GenericReductionUtilities::GenericReduction<data_container_type, Norms::Value, SpatialMethodHelperUtilities::SumOperation, false, 1>(
                rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, Norms::Value())
            .GetValue();
    }, data_container);
}

template<class TDataType>
TDataType SpatialMethods::Sum(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation)
{
    return SpatialMethodHelperUtilities::GenericSumReduction<TDataType, Norms::Value, 1>(
        rModelPart, rVariable, rLocation, Norms::Value());
}

template<class TDataType>
double SpatialMethods::Sum(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    const typename Norms::NormType<TDataType>::type& rNorm)
{
    return std::visit([&rModelPart, &rVariable, &rLocation](const auto& rNorm) {
        return SpatialMethodHelperUtilities::GenericSumReduction<TDataType, std::decay_t<decltype(rNorm)>, 1>(
            rModelPart, rVariable, rLocation, rNorm);
    }, rNorm);
}

template<class TDataType>
TDataType SpatialMethods::Mean(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation)
{
    const double number_of_items = SpatialMethodHelperUtilities::GetDataLocationSize(rModelPart, rLocation);
    return Sum(rModelPart, rVariable, rLocation) / std::max(number_of_items, 1.0);
}

template<class TDataType>
double SpatialMethods::Mean(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    const typename Norms::NormType<TDataType>::type& rNorm)
{
    const double number_of_items = SpatialMethodHelperUtilities::GetDataLocationSize(rModelPart, rLocation);
    return Sum(rModelPart, rVariable, rLocation, rNorm) / std::max(number_of_items, 1.0);
}

template<class TDataType>
TDataType SpatialMethods::RootMeanSquare(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation)
{
    const double number_of_items = SpatialMethodHelperUtilities::GetDataLocationSize(rModelPart, rLocation);
    auto sum_square = SpatialMethodHelperUtilities::GenericSumReduction<TDataType, Norms::Value, 2>(rModelPart, rVariable, rLocation, Norms::Value());

    const auto number_of_components = DataTypeTraits<TDataType>::Size(sum_square);
    for (IndexType i = 0; i < number_of_components; ++i) {
        auto& v = DataTypeTraits<TDataType>::GetComponent(sum_square, i);
        v = std::pow(v / std::max(number_of_items, 1.0), 0.5);
    }
    return sum_square;
}

template<class TDataType>
double SpatialMethods::RootMeanSquare(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    const typename Norms::NormType<TDataType>::type& rNorm)
{
    const double number_of_items = SpatialMethodHelperUtilities::GetDataLocationSize(rModelPart, rLocation);
    return std::visit([&rModelPart, &rVariable, &rLocation, number_of_items](const auto& rNorm){
        return std::pow(SpatialMethodHelperUtilities::GenericSumReduction<TDataType, std::decay_t<decltype(rNorm)>, 2>(rModelPart, rVariable, rLocation, rNorm) / std::max(number_of_items, 1.0), 0.5);
    }, rNorm);
}

template<class TDataType>
std::tuple<TDataType, TDataType> SpatialMethods::Variance(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation)
{
    const auto& mean = Mean(rModelPart, rVariable, rLocation);
    auto rms = RootMeanSquare(rModelPart, rVariable, rLocation);

    const auto number_of_components = DataTypeTraits<TDataType>::Size(rms);
    for (IndexType i = 0; i < number_of_components; ++i) {
        auto& v = DataTypeTraits<TDataType>::GetComponent(rms, i);
        v = std::pow(v, 2) - std::pow(DataTypeTraits<TDataType>::GetComponent(mean, i), 2);
    }
    return std::tuple<TDataType, TDataType>(mean, rms);
}

template<class TDataType>
std::tuple<double, double> SpatialMethods::Variance(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    const typename Norms::NormType<TDataType>::type& rNorm)
{
    const double mean = Mean(rModelPart, rVariable, rLocation, rNorm);
    const double rms = RootMeanSquare(rModelPart, rVariable, rLocation, rNorm);
    return std::make_tuple(mean, std::pow(rms, 2) - std::pow(mean, 2));
}

template<class TDataType>
std::tuple<TDataType, SpatialMethods::ItemPositionType<TDataType>> SpatialMethods::Min(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation)
{
    return SpatialMethodHelperUtilities::GenericReductionWithIndices<TDataType, Norms::Value, SpatialMethodHelperUtilities::MinOperation>(rModelPart, rVariable, rLocation, Norms::Value());
}

template<class TDataType>
std::tuple<double, SpatialMethods::IndexType> SpatialMethods::Min(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    const typename Norms::NormType<TDataType>::type& rNorm)
{
    return std::visit([&rModelPart, &rVariable, &rLocation](const auto& rNorm) {
        return SpatialMethodHelperUtilities::GenericReductionWithIndices<TDataType, std::decay_t<decltype(rNorm)>, SpatialMethodHelperUtilities::MinOperation>(rModelPart, rVariable, rLocation, rNorm);
    }, rNorm);
}

template<class TDataType>
std::tuple<TDataType, SpatialMethods::ItemPositionType<TDataType>> SpatialMethods::Max(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation)
{
    return SpatialMethodHelperUtilities::GenericReductionWithIndices<TDataType, Norms::Value, SpatialMethodHelperUtilities::MaxOperation>(rModelPart, rVariable, rLocation, Norms::Value());
}

template<class TDataType>
std::tuple<double, SpatialMethods::IndexType> SpatialMethods::Max(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    const typename Norms::NormType<TDataType>::type& rNorm)
{
    return std::visit([&rModelPart, &rVariable, &rLocation](const auto& rNorm) {
        return SpatialMethodHelperUtilities::GenericReductionWithIndices<TDataType, std::decay_t<decltype(rNorm)>, SpatialMethodHelperUtilities::MaxOperation>(rModelPart, rVariable, rLocation, rNorm);
    }, rNorm);
}

template<class TDataType>
std::tuple<TDataType, SpatialMethods::ItemPositionType<TDataType>> SpatialMethods::Median(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation)
{
    return SpatialMethodHelperUtilities::GenericReductionWithIndices<TDataType, Norms::Value, SpatialMethodHelperUtilities::MedianOperation>(rModelPart, rVariable, rLocation, Norms::Value());
}

template<class TDataType>
std::tuple<double, SpatialMethods::IndexType> SpatialMethods::Median(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    const typename Norms::NormType<TDataType>::type& rNorm)
{
    return std::visit([&rModelPart, &rVariable, &rLocation](const auto& rNorm) {
        return SpatialMethodHelperUtilities::GenericReductionWithIndices<TDataType, std::decay_t<decltype(rNorm)>, SpatialMethodHelperUtilities::MedianOperation>(rModelPart, rVariable, rLocation, rNorm);
    }, rNorm);
}

template<class TDataType>
SpatialMethods::DistributionInfo<TDataType> SpatialMethods::Distribution(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    Parameters Params)
{
    return SpatialMethodHelperUtilities::GenericDistribution<TDataType, Norms::Value>(rModelPart, rVariable, rLocation, Params, Norms::Value());
}

template<class TDataType>
SpatialMethods::DistributionInfo<double> SpatialMethods::Distribution(
    const ModelPart& rModelPart,
    const Variable<TDataType>& rVariable,
    const DataLocation& rLocation,
    Parameters Params,
    const typename Norms::NormType<TDataType>::type& rNorm)
{
    return std::visit([&rModelPart, &rVariable, &rLocation, Params](const auto& rNorm) {
        return SpatialMethodHelperUtilities::GenericDistribution<TDataType, std::decay_t<decltype(rNorm)>>(rModelPart, rVariable, rLocation, Params, rNorm);
    }, rNorm);
}

// template instantiations
#define KRATOS_TEMPLATE_VARIABLE_METHOD_INSTANTIATION(...)                                                                                                  \
    template __VA_ARGS__ SpatialMethods::Sum(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&);                 \
    template double SpatialMethods::Sum(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&, const typename Norms::NormType<__VA_ARGS__>::type&);                 \
    template __VA_ARGS__ SpatialMethods::Mean(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&);                 \
    template double SpatialMethods::Mean(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&, const typename Norms::NormType<__VA_ARGS__>::type&);                 \
    template __VA_ARGS__ SpatialMethods::RootMeanSquare(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&);                 \
    template double SpatialMethods::RootMeanSquare(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&, const typename Norms::NormType<__VA_ARGS__>::type&);                 \
    template std::tuple<__VA_ARGS__, __VA_ARGS__> SpatialMethods::Variance(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&);            \
    template std::tuple<double, double> SpatialMethods::Variance(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&, const typename Norms::NormType<__VA_ARGS__>::type&);            \
    template std::tuple<__VA_ARGS__, SpatialMethods::ItemPositionType<__VA_ARGS__>> SpatialMethods::Min(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&); \
    template std::tuple<double, SpatialMethods::IndexType> SpatialMethods::Min(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&, const typename Norms::NormType<__VA_ARGS__>::type&); \
    template std::tuple<__VA_ARGS__, SpatialMethods::ItemPositionType<__VA_ARGS__>> SpatialMethods::Max(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&); \
    template std::tuple<double, SpatialMethods::IndexType> SpatialMethods::Max(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&, const typename Norms::NormType<__VA_ARGS__>::type&); \
    template std::tuple<__VA_ARGS__, SpatialMethods::ItemPositionType<__VA_ARGS__>> SpatialMethods::Median(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&); \
    template std::tuple<double, SpatialMethods::IndexType> SpatialMethods::Median(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&, const typename Norms::NormType<__VA_ARGS__>::type&); \
    template SpatialMethods::DistributionInfo<__VA_ARGS__> SpatialMethods::Distribution(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&, Parameters);\
    template SpatialMethods::DistributionInfo<double> SpatialMethods::Distribution(const ModelPart&, const Variable<__VA_ARGS__>&, const DataLocation&, Parameters, const typename Norms::NormType<__VA_ARGS__>::type&);

KRATOS_TEMPLATE_VARIABLE_METHOD_INSTANTIATION(double)
KRATOS_TEMPLATE_VARIABLE_METHOD_INSTANTIATION(array_1d<double, 3>)
KRATOS_TEMPLATE_VARIABLE_METHOD_INSTANTIATION(array_1d<double, 4>)
KRATOS_TEMPLATE_VARIABLE_METHOD_INSTANTIATION(array_1d<double, 6>)
KRATOS_TEMPLATE_VARIABLE_METHOD_INSTANTIATION(array_1d<double, 9>)
KRATOS_TEMPLATE_VARIABLE_METHOD_INSTANTIATION(Vector)
KRATOS_TEMPLATE_VARIABLE_METHOD_INSTANTIATION(Matrix)

#undef KRATOS_TEMPLATE_VARIABLE_METHOD_INSTANTIATION

} // namespace Kratos