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

    template<class TDataType>
    using MethodReturnType = std::conditional_t<std::is_arithmetic_v<TDataType>, std::variant<TDataType>, std::variant<TDataType, double>>;

    ///@}
    ///@name Class definitions
    ///@{

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
    static TDataType Sum(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation);

    template<class TDataType>
    static double Sum(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation,
        const std::string& rNormType);

    template<class TDataType>
    static TDataType Mean(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation);

    template<class TDataType>
    static double Mean(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation,
        const std::string& rNormType);

    template<class TDataType>
    static TDataType RootMeanSquare(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation);

    template<class TDataType>
    static double RootMeanSquare(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation,
        const std::string& rNormType);

    template<class TDataType>
    static std::tuple<TDataType, TDataType> Variance(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation);

    template<class TDataType>
    static std::tuple<double, double> Variance(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation,
        const std::string& rNormType);

    // template<class TDataType>
    // static SupportedDataType Median(
    //     const ModelPart& rModelPart,
    //     const Variable<TDataType>& rVariable,
    //     const std::string& rNormType,
    //     const DataLocation& rLocation);

    template<class TDataType>
    static std::tuple<SpatialMethods::MethodReturnType<TDataType>, SpatialMethods::IndicesType> Min(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const std::string& rNormType,
        const DataLocation& rLocation);

    // template<class TDataType>
    // static SupportedDataType Max(
    //     const ModelPart& rModelPart,
    //     const Variable<TDataType>& rVariable,
    //     const std::string& rNormType,
    //     const DataLocation& rLocation);

    template<class TDataType>
    static DistributionInfoType Distribution(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const std::string& rNormType,
        const DataLocation& rLocation,
        Parameters Params);

    ///@}
    ///@name Class definitions
    ///@{

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

    ///@name Static operations
    ///@{

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
            return Sum(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
        }

        template <class TDataType>
        double static CalculateNormSum(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            return Sum(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), rNormType);
        }

        template <class TDataType>
        TDataType static CalculateRootMeanSquare(const ModelPart& rModelPart, const Variable<TDataType>& rVariable)
        {
            return RootMeanSquare(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
        }

        template <class TDataType>
        double static CalculateNormRootMeanSquare(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            return RootMeanSquare(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), rNormType);
        }

        template <class TDataType>
        TDataType static CalculateMean(const ModelPart& rModelPart, const Variable<TDataType>& rVariable)
        {
            return Mean(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
        }

        template <class TDataType>
        double static CalculateNormMean(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            return Mean(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), rNormType);
        }

        template <class TDataType>
        std::tuple<TDataType, TDataType> static CalculateVariance(
            const ModelPart& rModelPart, const Variable<TDataType>& rVariable)
        {
            return Variance(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
        }

        template <class TDataType>
        std::tuple<double, double> static CalculateNormVariance(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            return Variance(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), rNormType);
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
