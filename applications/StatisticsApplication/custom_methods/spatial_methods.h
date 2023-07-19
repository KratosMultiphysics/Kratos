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
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <tuple>
#include <vector>
#include <type_traits>
#include <variant>

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
#include "custom_utilities/norms.h"
#include "custom_utilities/method_utilities.h"
#include "custom_utilities/generic_reduction_utilities.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

class SpatialMethods
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = unsigned int;

    using IndicesType = std::vector<IndexType>;

    using DataLocation = Globals::DataLocation;

    ///@}
    ///@name Class definitions
    ///@{

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
            if (OperationTraits::Resize(mValue, rOther.mValue)) {
                Initialize();
            }
            mValue += rOther.mValue;
        }

        void Synchronize(const DataCommunicator& rDataCommunicator)
        {
            if (OperationTraits::SynchronizeSize(mValue, rDataCommunicator)) {
                Initialize();
            }

            typename OperationTraits::VectorType local_values, global_values;
            OperationTraits::FillToVector(local_values, mValue);
            rDataCommunicator.SumAll(local_values, global_values);
            OperationTraits::FillFromVector(mValue, global_values);
        }

        ///@}

    private:
        ///@name Private member variables
        ///@{

        TDataType mValue;

        ///@}
        ///@name Private operations
        ///@{

        void Initialize() { OperationTraits::Initialize(mValue, 0.0); }

        ///@}
    };

    template<class TDataType>
    class MinOperation
    {
    public:
        ///@name Type definitions
        ///@{

        using IndicesTraits = DataTypeTraits<IndicesType>;

        using OperationTraits = DataTypeTraits<TDataType>;

        ///@}
        ///@name Life cycle
        ///@{

        MinOperation() { Initialize(); }

        MinOperation(
            const TDataType& rValue,
            const IndexType rId)
        {
            mValue = rValue;
            IndicesTraits::Resize(mIndices, IndicesType(OperationTraits::Size(mValue)));
            std::fill(mIndices.begin(), mIndices.end(), rId);
        }

        ///@}
        ///@name Public operations
        ///@{

        std::tuple<TDataType, IndicesType> GetValue() const
        {
            return std::make_tuple(mValue, mIndices);
        }

        void Execute(const MinOperation<TDataType>& rOther)
        {
            if (OperationTraits::Resize(mValue, rOther.mValue)) {
                Initialize();
            }

            for (IndexType i = 0; i < OperationTraits::Size(mValue); ++i) {
                if (OperationTraits::GetComponent(mValue, i) > OperationTraits::GetComponent(rOther.mValue, i)) {
                    OperationTraits::GetComponent(mValue, i) = OperationTraits::GetComponent(rOther.mValue, i);
                    mIndices[i] = rOther.mIndices[i];
                }
            }
        }

        void Synchronize(const DataCommunicator& rDataCommunicator)
        {
            if (OperationTraits::SynchronizeSize(mValue, rDataCommunicator)) {
                Initialize();
            }

            typename OperationTraits::VectorType local_values;
            OperationTraits::FillToVector(local_values, mValue);
            auto global_values = rDataCommunicator.AllGatherv(local_values);
            auto global_indices = rDataCommunicator.AllGatherv(mIndices);


            for (IndexType i = 0; i < OperationTraits::Size(mValue); ++i) {
                auto& current_value = local_values[i];
                auto& current_index = mIndices[i];
                for (IndexType rank = 0; rank < global_values.size(); ++rank) {
                    if (current_value > global_values[rank][i]) {
                        current_value = global_values[rank][i];
                        current_index = global_indices[rank][i];
                    }
                }
            }

            OperationTraits::FillFromVector(mValue, local_values);
        }

        ///@}

    private:
        ///@name Private member variables
        ///@{

        TDataType mValue;

        IndicesType mIndices;

        ///@}
        ///@name Private operations
        ///@{

        void Initialize()
        {
            OperationTraits::Initialize(mValue, std::numeric_limits<double>::max());
            IndicesTraits::Resize(mIndices, IndicesType(OperationTraits::Size(mValue)));
            IndicesTraits::Initialize(mIndices, std::numeric_limits<IndexType>::max());
        }

        ///@}
    };

    template<class TDataType>
    class MaxOperation
    {
    public:
        ///@name Type definitions
        ///@{

        using IndicesTraits = DataTypeTraits<IndicesType>;

        using OperationTraits = DataTypeTraits<TDataType>;

        ///@}
        ///@name Life cycle
        ///@{

        MaxOperation() { Initialize(); }

        MaxOperation(
            const TDataType& rValue,
            const IndexType rId)
        {
            mValue = rValue;
            IndicesTraits::Resize(mIndices, IndicesType(OperationTraits::Size(mValue)));
            std::fill(mIndices.begin(), mIndices.end(), rId);
        }

        ///@}
        ///@name Public operations
        ///@{

        std::tuple<TDataType, IndicesType> GetValue() const
        {
            return std::make_tuple(mValue, mIndices);
        }

        void Execute(const MaxOperation<TDataType>& rOther)
        {
            if (OperationTraits::Resize(mValue, rOther.mValue)) {
                Initialize();
            }

            for (IndexType i = 0; i < OperationTraits::Size(mValue); ++i) {
                if (OperationTraits::GetComponent(mValue, i) < OperationTraits::GetComponent(rOther.mValue, i)) {
                    OperationTraits::GetComponent(mValue, i) = OperationTraits::GetComponent(rOther.mValue, i);
                    mIndices[i] = rOther.mIndices[i];
                }
            }
        }

        void Synchronize(const DataCommunicator& rDataCommunicator)
        {
            if (OperationTraits::SynchronizeSize(mValue, rDataCommunicator)) {
                Initialize();
            }

            typename OperationTraits::VectorType local_values;
            OperationTraits::FillToVector(local_values, mValue);
            auto global_values = rDataCommunicator.AllGatherv(local_values);
            auto global_indices = rDataCommunicator.AllGatherv(mIndices);


            for (IndexType i = 0; i < OperationTraits::Size(mValue); ++i) {
                auto& current_value = local_values[i];
                auto& current_index = mIndices[i];
                for (IndexType rank = 0; rank < global_values.size(); ++rank) {
                    if (current_value < global_values[rank][i]) {
                        current_value = global_values[rank][i];
                        current_index = global_indices[rank][i];
                    }
                }
            }

            OperationTraits::FillFromVector(mValue, local_values);
        }

        ///@}

    private:
        ///@name Private member variables
        ///@{

        TDataType mValue;

        IndicesType mIndices;

        ///@}
        ///@name Private operations
        ///@{

        void Initialize()
        {
            OperationTraits::Initialize(mValue, std::numeric_limits<double>::lowest());
            IndicesTraits::Resize(mIndices, IndicesType(OperationTraits::Size(mValue)));
            IndicesTraits::Initialize(mIndices, std::numeric_limits<IndexType>::max());
        }

        ///@}
    };

    template<class TDataType>
    class MedianOperation
    {
    public:
        ///@name Type definitions
        ///@{

        using IndicesTraits = DataTypeTraits<IndicesType>;

        using OperationTraits = DataTypeTraits<TDataType>;

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

        std::tuple<TDataType, IndicesType> GetValue() const
        {
            return std::make_tuple(mResultantValue, mResultantIndex);
        }

        void Execute(const MedianOperation<TDataType>& rOther)
        {
            // first check the sizes. Requires resizing if dynamic types
            // such as Vector and Matrices are used.
            if (mValues.size() != rOther.mValues.size()) {
                mValues.resize(rOther.mValues.size());
                mResultantIndex.resize(rOther.mValues.size());
            }

            // iterate through each component
            for (IndexType i_comp = 0; i_comp < rOther.mValues.size(); ++i_comp) {
                auto& current_component_values = mValues[i_comp];
                const auto& r_other_component_values = rOther.mValues[i_comp];

                // assumes all the values in r_component_values is sorted
                // in the acending order of values.
                for (const auto& r_other_value_info : r_other_component_values) {
                    auto i_lower_bound = std::lower_bound(
                        current_component_values.begin(),
                        current_component_values.end(), r_other_value_info,
                        [](const auto& rV1, const auto& rV2) {
                            return std::get<0>(rV1) < std::get<0>(rV2);
                        });
                    current_component_values.insert(i_lower_bound, r_other_value_info);
                }
            }
        }

        void Synchronize(const DataCommunicator& rDataCommunicator)
        {
            std::vector<RawDataType> results(mValues.size());
            const IndexType data_size = rDataCommunicator.MaxAll(OperationTraits::Size(mResultantValue));
            if (mValues.size() != data_size) {
                mValues.resize(data_size);
                OperationTraits::SynchronizeSize(mResultantValue, rDataCommunicator);
                mResultantIndex.resize(data_size);
                results.resize(data_size);
            }

            for (IndexType i_comp = 0; i_comp < data_size; ++i_comp) {
                auto& current_median_value = results[i_comp];
                auto& current_median_index = mResultantIndex[i_comp];
                const auto& current_values = mValues[i_comp];

                // get the values in rank 0
                std::vector<RawDataType> local_values(current_values.size());
                std::transform(current_values.begin(), current_values.end(), local_values.begin(), [](const auto& rV) { return std::get<0>(rV); });
                const auto& global_values = rDataCommunicator.Gatherv(local_values, 0);

                // get the indices in rank 0
                std::vector<IndexType> local_indices(current_values.size());
                std::transform(current_values.begin(), current_values.end(), local_indices.begin(), [](const auto& rV) { return std::get<1>(rV); });
                const auto& global_indices = rDataCommunicator.Gatherv(local_indices, 0);

                if (rDataCommunicator.Rank() == 0) {
                    const IndexType number_of_values = std::accumulate(global_values.begin(), global_values.end(), 0, [](const auto& rV1, const auto& rV2) { return rV1 + rV2.size();});
                    std::vector<RawDataType> sorted_values;
                    sorted_values.reserve(number_of_values);
                    std::vector<IndexType> sorted_indices;
                    sorted_indices.reserve(number_of_values);

                    for (IndexType rank = 0; rank < global_values.size(); ++rank) {
                        const auto& rank_values = global_values[rank];
                        const auto& rank_indices = global_indices[rank];
                        for (IndexType i_value = 0; i_value < rank_values.size(); ++i_value) {
                            const auto lower_bound = std::lower_bound(sorted_values.begin(), sorted_values.end(), rank_values[i_value]);
                            sorted_indices.insert(sorted_indices.begin() + std::distance(sorted_values.begin(), lower_bound), rank_indices[i_value]);
                            sorted_values.insert(lower_bound, rank_values[i_value]);
                        }
                    }

                    if (number_of_values > 0) {
                        const IndexType mid_point = number_of_values / 2;
                        current_median_index = sorted_indices[mid_point];
                        if (number_of_values % 2 != 0) {
                            current_median_value = sorted_values[mid_point];
                        } else {
                            const IndexType adjacent_mid_point = (number_of_values - 1) / 2;
                            current_median_value = (sorted_values[adjacent_mid_point] + sorted_values[mid_point]) * 0.5;
                        }
                    }
                }

                rDataCommunicator.Broadcast(current_median_value, 0);
                rDataCommunicator.Broadcast(current_median_index, 0);
                OperationTraits::GetComponent(mResultantValue, i_comp) = current_median_value;
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

    template<class TDataType>
    class DistributionInfo
    {
    public:
        ///@name Type definitions
        ///@{

        using DataTypeTrait = DataTypeTraits<TDataType>;

        using RawDataType = typename DataTypeTrait::RawDataType;

        ///@}
        ///@name Life cycle
        ///@{

        DistributionInfo()
        {
            DataTypeTrait::Initialize(mMin, std::numeric_limits<RawDataType>::max());
            DataTypeTrait::Initialize(mMax, std::numeric_limits<RawDataType>::lowest());
        }

        ///@}
        ///@name Public operations
        ///@{

        TDataType GetMin() const { return mMin; }

        TDataType GetMax() const { return mMax; }

        std::vector<TDataType> GetGroupUpperValues() const { return mGroupUpperValues; }

        std::vector<std::vector<IndexType>> GetGroupNumberOfValues() const { return mGroupNumberOfValues; }

        std::vector<TDataType> GetGroupValueDistributionPercentage() const { return mGroupValueDistributionPercentage; }

        std::vector<TDataType> GetGroupMeans() const { return mGroupMean; }

        std::vector<TDataType> GetGroupVariances() const { return mGroupVariance; }

        ///@}

        ///@name Private member variables
        ///@{

        TDataType mMin;

        TDataType mMax;

        std::vector<TDataType> mGroupUpperValues;

        std::vector<std::vector<IndexType>> mGroupNumberOfValues;

        std::vector<TDataType> mGroupValueDistributionPercentage;

        std::vector<TDataType> mGroupMean;

        std::vector<TDataType> mGroupVariance;

        ///@}
        ///@name Friends
        ///@{


        ///@}
    };

    ///@}
    ///@name Dependent type definitions
    ///@{

    using DistributionInfoType = std::variant<
                                        DistributionInfo<int>,
                                        DistributionInfo<double>,
                                        DistributionInfo<array_1d<double, 3>>,
                                        DistributionInfo<array_1d<double, 4>>,
                                        DistributionInfo<array_1d<double, 6>>,
                                        DistributionInfo<array_1d<double, 9>>,
                                        DistributionInfo<Vector>,
                                        DistributionInfo<Matrix>
                                    >;

    ///@}
    ///@name Static operations
    ///@{

    template<class TDataType>
    static DistributionInfoType Distribution(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const std::string& rNormType,
        const DataLocation& rLocation,
        Parameters Params)
    {
        KRATOS_TRY

        const auto data_container = DataContainers::GetDataContainer(rModelPart, rVariable, rLocation);

        const auto r_norm_type = Norms::GetNorm<TDataType>(rNormType);

        return std::visit([&](auto& rDataContainer, auto& rNorm) -> DistributionInfoType {
            using data_container_type = std::decay_t<decltype(rDataContainer)>;
            using norm_type = std::decay_t<decltype(rNorm)>;
            using norm_return_type = typename norm_type::ResultantValueType<TDataType>;
            using data_type_traits = DataTypeTraits<norm_return_type>;

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

            DistributionInfo<norm_return_type> distribution_info;

            auto& min_value = distribution_info.mMin;
            if (Params["min_value"].IsDouble()) {
                data_type_traits::Initialize(min_value, Params["min_value"].GetDouble());
            } else if (Params["min_value"].IsString() && Params["min_value"].GetString() == "min") {
                min_value = std::get<0>(GenericReductionUtilities::GenericReduction<data_container_type, norm_type, MinOperation, true>(rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, rNorm).GetValue());
            } else {
                KRATOS_ERROR << "Unknown min_value. Allowed only double or \"min\" "
                                "string as a value. [ min_value = "
                                << Params["min_value"] << " ]\n.";
            }

            // first fix the max value
            auto& max_value = distribution_info.mMax;
            if (Params["max_value"].IsDouble()) {
                data_type_traits::Initialize(max_value, Params["max_value"].GetDouble());
            } else if (Params["max_value"].IsString() && Params["max_value"].GetString() == "max") {
                max_value = std::get<0>(GenericReductionUtilities::GenericReduction<data_container_type, norm_type, MaxOperation, true>(rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, rNorm).GetValue());
            } else {
                KRATOS_ERROR << "Unknown max_value. Allowed only double or \"max\" "
                                "string as a value. [ max_value = "
                                << Params["max_value"] << " ]\n.";
            }

            auto& group_limits = distribution_info.mGroupUpperValues;
            const IndexType number_of_groups = Params["number_of_value_groups"].GetInt();

            // we need additional two groups to store values below the specified minimum and values above the
            // specified maximum.
            const IndexType number_of_all_groups = number_of_groups + 2;
            group_limits.resize(number_of_all_groups);
            for (IndexType i = 0; i < number_of_groups + 1; ++i) {
                group_limits[i] = min_value + (max_value - min_value) * static_cast<double>(i) / static_cast<double>(number_of_groups);
            }

            const IndexType number_of_components = DataTypeTraits<norm_return_type>::Size(max_value);

            // final group limit is extended by a small amount. epsilon in numeric
            // limits cannot be used since testing also need to have the same
            // extending value in python. Therefore hard coded value is used
            auto& last_group_limit = group_limits[group_limits.size() - 2];
            for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                DataTypeTraits<norm_return_type>::GetComponent(last_group_limit, i_comp) += 1e-16;
            }
            norm_return_type additional_max_value;
            DataTypeTraits<norm_return_type>::Resize(additional_max_value, max_value);
            DataTypeTraits<norm_return_type>::Initialize(additional_max_value, std::numeric_limits<typename DataTypeTraits<norm_return_type>::RawDataType>::max());
            group_limits.back() = additional_max_value;

            /// reduction class
            class DistributionReduction
            {
            public:
                using data_type = std::tuple<
                                        IndexType,
                                        IndicesType,
                                        norm_return_type
                                    >;

                using return_type = std::tuple<
                                            std::vector<IndicesType>,       // group value counts
                                            std::vector<norm_return_type>,  // group means
                                            std::vector<norm_return_type>   // group variances
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
                        ++r_current_group_counts[group_index][i_comp];
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
                        const auto& r_other_group_variance = r_other_group_variances[i_group];

                        const IndexType number_of_components = r_other_group_count.size();

                        if (r_group_count.size() != number_of_components) {
                            r_group_count.resize(number_of_components);
                            std::fill(r_group_count.begin(), r_group_count.end(), 0);
                            DataTypeTraits<norm_return_type>::Resize(r_group_mean, r_other_group_mean);
                            DataTypeTraits<norm_return_type>::Initialize(r_group_mean, 0.0);
                            DataTypeTraits<norm_return_type>::Resize(r_group_variance, r_other_group_variance);
                            DataTypeTraits<norm_return_type>::Initialize(r_group_variance, 0.0);
                        }

                        for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                            r_group_count[i_comp] += r_other_group_count[i_comp];
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
                            r_current_group_counts[i_group].resize(number_of_components);
                            std::fill(r_current_group_counts[i_group].begin(), r_current_group_counts[i_group].end(), 0);

                            DataTypeTraits<norm_return_type>::Resize(r_current_group_means[i_group], r_value);
                            DataTypeTraits<norm_return_type>::Initialize(r_current_group_means[i_group], 0.0);
                            DataTypeTraits<norm_return_type>::Resize(r_current_group_variances[i_group], r_value);
                            DataTypeTraits<norm_return_type>::Initialize(r_current_group_variances[i_group], 0.0);
                        }
                    }
                }
            };

            const auto& reuduced_values = IndexPartition<IndexType>(rDataContainer.Size()).for_each<DistributionReduction>(TDataType{}, [&rDataContainer, &rNorm, &group_limits, number_of_components, number_of_groups](const IndexType Index, TDataType& rTLS) {
                rDataContainer.GetValue(rTLS, Index);
                auto norm_value = rNorm.Evaluate(rTLS);

                IndicesType indices(number_of_components);
                for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                    const auto& comp_value = DataTypeTraits<norm_return_type>::GetComponent(norm_value, i_comp);
                    for (IndexType i_group = 0; i_group < group_limits.size(); ++i_group) {
                        if (comp_value < DataTypeTraits<norm_return_type>::GetComponent(group_limits[i_group], i_comp)) {
                            indices[i_comp] = i_group;
                            break;
                        }
                    }
                }

                return std::make_tuple(group_limits.size(), indices, norm_value);
            });

            // now prepare data for mpi communication
            std::vector<typename DataTypeTraits<norm_return_type>::RawDataType> local_values;
            local_values.resize(number_of_components * number_of_all_groups * 2);
            auto values_begin = local_values.begin();
            IndicesType local_distribution;
            local_distribution.resize(number_of_components * number_of_all_groups);
            auto indices_begin = local_distribution.begin();
            for (IndexType i_group = 0; i_group < number_of_all_groups; ++i_group) {
                DataTypeTraits<IndicesType>::FillToVector(indices_begin, std::get<0>(reuduced_values)[i_group]); // get the group counts
                DataTypeTraits<norm_return_type>::FillToVector(values_begin, std::get<1>(reuduced_values)[i_group]); // put the means
                DataTypeTraits<norm_return_type>::FillToVector(values_begin + number_of_components, std::get<2>(reuduced_values)[i_group]); // put the variances
                values_begin += 2 * number_of_components;
                indices_begin += number_of_components;
            }

            // now do the mpi communicaton
            const auto& global_distribution = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_distribution);
            const auto& global_values = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_values);

            const double number_of_items = static_cast<double>(std::max(
                std::accumulate(global_distribution.begin(), global_distribution.end(), 0), 1)) / number_of_components;

            // now revert back to the group values
            distribution_info.mGroupNumberOfValues.resize(number_of_all_groups);
            distribution_info.mGroupValueDistributionPercentage.resize(number_of_all_groups);
            distribution_info.mGroupMean.resize(number_of_all_groups);
            distribution_info.mGroupVariance.resize(number_of_all_groups);
            auto global_indices_begin = global_distribution.begin();
            auto global_values_begin = global_values.begin();
            for (IndexType i_group = 0; i_group < number_of_all_groups; ++i_group) {
                auto& current_number_of_values = distribution_info.mGroupNumberOfValues[i_group];
                current_number_of_values.resize(number_of_components);
                DataTypeTraits<IndicesType>::FillFromVector(current_number_of_values, global_indices_begin, global_indices_begin + number_of_components);

                auto& current_distribution_percentage = distribution_info.mGroupValueDistributionPercentage[i_group];
                DataTypeTraits<norm_return_type>::Resize(current_distribution_percentage, max_value);
                DataTypeTraits<norm_return_type>::FillFromVector(current_distribution_percentage, global_indices_begin, global_indices_begin + number_of_components);

                auto& current_mean = distribution_info.mGroupMean[i_group];
                DataTypeTraits<norm_return_type>::Resize(current_mean, max_value);
                DataTypeTraits<norm_return_type>::FillFromVector(current_mean, global_values_begin, global_values_begin + number_of_components);

                auto& current_variance = distribution_info.mGroupVariance[i_group];
                DataTypeTraits<norm_return_type>::Resize(current_variance, max_value);
                DataTypeTraits<norm_return_type>::FillFromVector(current_variance, global_values_begin + number_of_components, global_values_begin + 2 * number_of_components);

                global_indices_begin += number_of_components;
                global_values_begin += number_of_components * 2;

                // post processing of values
                current_distribution_percentage /= number_of_items;
                for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                    const auto n = current_number_of_values[i_comp];
                    if (n > 0) {
                        DataTypeTraits<norm_return_type>::GetComponent(current_mean, i_comp) /= n;
                        DataTypeTraits<norm_return_type>::GetComponent(current_variance, i_comp) /= n;
                    }
                    DataTypeTraits<norm_return_type>::GetComponent(current_variance, i_comp) -= std::pow(DataTypeTraits<norm_return_type>::GetComponent(current_mean, i_comp), 2);
                }
            }

            // reversing group limit extention
            for (IndexType i_comp = 0; i_comp < number_of_components; ++i_comp) {
                DataTypeTraits<norm_return_type>::GetComponent(group_limits[group_limits.size() - 2], i_comp) += 1e-16;
            }
            group_limits[group_limits.size() - 1] = max_value;

            return distribution_info;

        }, data_container, r_norm_type);

        KRATOS_CATCH("");
    }

    template<class TDataType, int TPower = 1>
    static SupportedDataType GenericSumReduction(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const std::string& rNormType,
        const DataLocation& rLocation)
    {
        KRATOS_TRY

        const auto data_container = DataContainers::GetDataContainer(rModelPart, rVariable, rLocation);

        const auto r_norm_type = Norms::GetNorm<TDataType>(rNormType);

        return std::visit([&rModelPart](auto& rDataContainer, auto& rNorm) -> SupportedDataType {
            using data_container_type = std::decay_t<decltype(rDataContainer)>;
            using norm_type = std::decay_t<decltype(rNorm)>;
            return GenericReductionUtilities::GenericReduction<data_container_type, norm_type, SumOperation, false, TPower>(
                    rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, rNorm)
                .GetValue();
        }, data_container, r_norm_type);

        KRATOS_CATCH("");
    }

    template<class TDataType, template <class T1> class OperationType>
    static std::tuple<SupportedDataType, IndicesType> GenericReductionWithIndices(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const std::string& rNormType,
        const DataLocation& rLocation)
    {
        KRATOS_TRY

        const auto data_container = DataContainers::GetDataContainer(rModelPart, rVariable, rLocation);

        const auto r_norm_type = Norms::GetNorm<TDataType>(rNormType);

        return std::visit([&rModelPart](auto& rDataContainer, auto& rNorm) {
            using data_container_type = std::decay_t<decltype(rDataContainer)>;
            using norm_type = std::decay_t<decltype(rNorm)>;
            const auto& value_pair = GenericReductionUtilities::GenericReduction<data_container_type, norm_type, OperationType, true>(
                    rModelPart.GetCommunicator().GetDataCommunicator(), rDataContainer, rNorm)
                .GetValue();

            const SupportedDataType value = std::get<0>(value_pair);
            const IndicesType indices = std::get<1>(value_pair);

            return std::make_pair(value, indices);
        }, data_container, r_norm_type);

        KRATOS_CATCH("");
    }

    template<class T>
    static DataLocation GetDataLocation()
    {
        if constexpr(std::is_same_v<T, MethodUtilities::HistoricalDataValueRetrievalFunctor<ModelPart::NodeType>>) {
            return DataLocation::NodeHistorical;
        } else if constexpr(std::is_same_v<T, MethodUtilities::NonHistoricalDataValueRetrievalFunctor<ModelPart::NodeType>>) {
            return DataLocation::NodeNonHistorical;
        } else if constexpr(std::is_same_v<T, MethodUtilities::NonHistoricalDataValueRetrievalFunctor<ModelPart::ConditionType>>) {
            return DataLocation::Condition;
        } else if constexpr(std::is_same_v<T, MethodUtilities::NonHistoricalDataValueRetrievalFunctor<ModelPart::ElementType>>) {
            return DataLocation::Element;
        } else {
            KRATOS_ERROR << "Unsupported type";
            return DataLocation::NodeHistorical;
        }
    }

    template <class TContainerType, class TContainerItemType, template <class T> class TDataRetrievalFunctor>
    class ContainerSpatialMethods
    {
    public:
        // special overloaded method for flags
        int static CalculateSum(const ModelPart& rModelPart, const Flags& rVariable)
        {
            const TContainerType& r_container =
                MethodUtilities::GetLocalDataContainer<TContainerType>(rModelPart);
            int sum = 0;

    #pragma omp parallel for reduction(+ : sum)
            for (int i = 0; i < static_cast<int>(r_container.size()); ++i)
            {
                const TContainerItemType& r_item = *(r_container.begin() + i);
                sum += r_item.Is(rVariable);
            }

            sum = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(sum);
            return sum;
        }

        template <class TDataType>
        TDataType static CalculateSum(const ModelPart& rModelPart, const Variable<TDataType>& rVariable)
        {
            KRATOS_TRY

            const auto value = GenericSumReduction<TDataType, 1>(rModelPart, rVariable, "value", GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
            return std::get<TDataType>(value);

            KRATOS_CATCH("");
        }

        template <class TDataType>
        double static CalculateNormSum(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            KRATOS_TRY

            const auto value = GenericSumReduction<TDataType, 1>(rModelPart, rVariable, rNormType, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
            return std::get<double>(value);

            KRATOS_CATCH("");
        }

        template <class TDataType>
        TDataType static CalculateRootMeanSquare(const ModelPart& rModelPart, const Variable<TDataType>& rVariable)
        {
            KRATOS_TRY

            const auto value = GenericSumReduction<TDataType, 2>(rModelPart, rVariable, "value", GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
            const auto square_sum = std::get<TDataType>(value);

            const auto data_container = DataContainers::GetDataContainer(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());

            const int local_size = std::visit([](const auto& rDataContainer){ return rDataContainer.Size(); }, data_container);
            const int global_size = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_size);

            return MethodUtilities::RaiseToPower<TDataType>(square_sum * (1.0 / std::max(global_size, 1)), 0.5);

            KRATOS_CATCH("");
        }

        template <class TDataType>
        double static CalculateNormRootMeanSquare(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            KRATOS_TRY

            const auto value = GenericSumReduction<TDataType, 2>(rModelPart, rVariable, rNormType, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
            const double square_sum = std::get<double>(value);

            const auto data_container = DataContainers::GetDataContainer(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
            const int local_size = std::visit([](const auto& rDataContainer){ return rDataContainer.Size(); }, data_container);
            const int global_size = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_size);

            return MethodUtilities::RaiseToPower<double>(square_sum * (1.0 / std::max(global_size, 1)), 0.5);

            KRATOS_CATCH("");
        }

        template <class TDataType>
        TDataType static CalculateMean(const ModelPart& rModelPart, const Variable<TDataType>& rVariable)
        {
            const TDataType& sum = CalculateSum<TDataType>(rModelPart, rVariable);
            const TContainerType& r_container =
                MethodUtilities::GetLocalDataContainer<TContainerType>(rModelPart);

            const double number_of_items =
                rModelPart.GetCommunicator().GetDataCommunicator().SumAll(
                    static_cast<double>(r_container.size()));

            return sum * (1.0 / std::max(number_of_items, 1.0));
        }

        template <class TDataType>
        double static CalculateNormMean(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            const double sum =
                CalculateNormSum<TDataType>(rModelPart, rVariable, rNormType, Params);
            const TContainerType& r_container =
                MethodUtilities::GetLocalDataContainer<TContainerType>(rModelPart);

            const double number_of_items =
                rModelPart.GetCommunicator().GetDataCommunicator().SumAll(
                    static_cast<double>(r_container.size()));

            if (number_of_items > 0)
            {
                return sum * (1.0 / number_of_items);
            }

            return 0.0;
        }

        template <class TDataType>
        std::tuple<TDataType, TDataType> static CalculateVariance(
            const ModelPart& rModelPart, const Variable<TDataType>& rVariable)
        {
            TDataType mean = CalculateMean<TDataType>(rModelPart, rVariable);
            TDataType rms = CalculateRootMeanSquare<TDataType>(rModelPart, rVariable);
            TDataType global_variance = MethodUtilities::RaiseToPower<TDataType>(rms, 2) - MethodUtilities::RaiseToPower<TDataType>(mean, 2);

            return std::make_tuple<TDataType, TDataType>(
                std::forward<TDataType>(mean), std::forward<TDataType>(global_variance));
        }

        template <class TDataType>
        std::tuple<double, double> static CalculateNormVariance(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            double mean = CalculateNormMean<TDataType>(rModelPart, rVariable, rNormType, Params);
            double rms = CalculateNormRootMeanSquare<TDataType>(rModelPart, rVariable, rNormType, Params);
            double global_variance = MethodUtilities::RaiseToPower<double>(rms, 2) - MethodUtilities::RaiseToPower<double>(mean, 2);

            return std::make_tuple<double, double>(
                std::forward<double>(mean), std::forward<double>(global_variance));
        }

        template <class TDataType>
        std::tuple<double, IndexType> static GetNormMax(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            KRATOS_TRY

            const auto value_info = GenericReductionWithIndices<TDataType, MaxOperation>(rModelPart, rVariable, rNormType, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
            const double value = std::get<double>(std::get<0>(value_info));
            const IndexType index = std::get<1>(value_info)[0];
            return std::make_tuple(value, index);

            KRATOS_CATCH("");
        }

        template <class TDataType>
        std::tuple<double, std::size_t> static GetNormMin(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            KRATOS_TRY

            const auto value_info = GenericReductionWithIndices<TDataType, MinOperation>(rModelPart, rVariable, rNormType, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
            const double value = std::get<double>(std::get<0>(value_info));
            const IndexType index = std::get<1>(value_info)[0];
            return std::make_tuple(value, index);

            KRATOS_CATCH("");
        }

        template <class TDataType>
        double static GetNormMedian(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            KRATOS_TRY

            const auto value_info = GenericReductionWithIndices<TDataType, MedianOperation>(rModelPart, rVariable, rNormType, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
            const double value = std::get<double>(std::get<0>(value_info));
            // const IndexType index = std::get<1>(value_info)[0];
            return value;

            KRATOS_CATCH("");
        }

        template <class TDataType>
        std::tuple<double, double, std::vector<double>, std::vector<IndexType>, std::vector<double>, std::vector<double>, std::vector<double>> static GetNormDistribution(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            KRATOS_TRY

            const auto& distribution_info = Distribution(rModelPart, rVariable, rNormType, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), Params);
            const auto& values = std::get<DistributionInfo<double>>(distribution_info);

            const auto& r_group_number_of_values = values.GetGroupNumberOfValues();
            std::vector<IndexType> number_of_values(r_group_number_of_values.size());
            std::transform(r_group_number_of_values.begin(), r_group_number_of_values.end(), number_of_values.begin(), [](const auto& V1) { return V1[0]; });

            return std::make_tuple(values.GetMin(), values.GetMax(), values.GetGroupUpperValues(),
                                   number_of_values,
                                   values.GetGroupValueDistributionPercentage(),
                                   values.GetGroupMeans(), values.GetGroupVariances());

            KRATOS_CATCH("");
        }
    };

using NodeType = ModelPart::NodeType;
using ElementType = ModelPart::ElementType;
using ConditionType = ModelPart::ConditionType;

using NodesContainerType = ModelPart::NodesContainerType;
using ElementsContainerType = ModelPart::ElementsContainerType;
using ConditionsContainerType = ModelPart::ConditionsContainerType;

class HistoricalSpatialMethods
    : public SpatialMethods::ContainerSpatialMethods<NodesContainerType, NodeType, MethodUtilities::HistoricalDataValueRetrievalFunctor>
{
};

class NodalNonHistoricalSpatialMethods
    : public SpatialMethods::ContainerSpatialMethods<NodesContainerType, NodeType, MethodUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

class ConditionNonHistoricalSpatialMethods
    : public SpatialMethods::ContainerSpatialMethods<ConditionsContainerType, ConditionType, MethodUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

class ElementNonHistoricalSpatialMethods
    : public SpatialMethods::ContainerSpatialMethods<ElementsContainerType, ElementType, MethodUtilities::NonHistoricalDataValueRetrievalFunctor>
{
};

};

///@}

///@} addtogroup block

} // namespace Kratos.
