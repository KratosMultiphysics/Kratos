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
    using ItemPositionType = std::conditional_t<std::is_arithmetic_v<TDataType>, IndexType, IndicesType>;

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
    };

    ///@}
    ///@name Static operations
    ///@{

    static IndexType Sum(
        const ModelPart& rModelPart,
        const Flags& rFlag,
        const DataLocation& rLocation);

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
        const typename Norms::NormType<TDataType>::type& rNorm);

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
        const typename Norms::NormType<TDataType>::type& rNorm);

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
        const typename Norms::NormType<TDataType>::type& rNorm);

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
        const typename Norms::NormType<TDataType>::type& rNorm);

    template<class TDataType>
    static std::tuple<TDataType, ItemPositionType<TDataType>> Min(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation);

    template<class TDataType>
    static std::tuple<double, IndexType> Min(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation,
        const typename Norms::NormType<TDataType>::type& rNorm);

    template<class TDataType>
    static std::tuple<TDataType, ItemPositionType<TDataType>> Max(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation);

    template<class TDataType>
    static std::tuple<double, IndexType> Max(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation,
        const typename Norms::NormType<TDataType>::type& rNorm);

    template<class TDataType>
    static std::tuple<TDataType, ItemPositionType<TDataType>> Median(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation);

    template<class TDataType>
    static std::tuple<double, IndexType> Median(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation,
        const typename Norms::NormType<TDataType>::type& rNorm);

    template<class TDataType>
    static DistributionInfo<TDataType> Distribution(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation,
        Parameters Params);

    template<class TDataType>
    static DistributionInfo<double> Distribution(
        const ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const DataLocation& rLocation,
        Parameters Params,
        const typename Norms::NormType<TDataType>::type& rNorm);

    ///@}
    ///@name Static operations
    ///@{

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
        int static CalculateSum(const ModelPart& rModelPart, const Flags& rFlag)
        {
            return Sum(rModelPart, rFlag, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>());
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
            return Sum(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), Norms::GetNorm<TDataType>(rNormType));
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
            return RootMeanSquare(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), Norms::GetNorm<TDataType>(rNormType));
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
            return Mean(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), Norms::GetNorm<TDataType>(rNormType));
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
            return Variance(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), Norms::GetNorm<TDataType>(rNormType));
        }

        template <class TDataType>
        std::tuple<double, IndexType> static GetNormMax(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            return Max(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), Norms::GetNorm<TDataType>(rNormType));
        }

        template <class TDataType>
        std::tuple<double, std::size_t> static GetNormMin(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            return Min(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), Norms::GetNorm<TDataType>(rNormType));
        }

        template <class TDataType>
        double static GetNormMedian(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            return std::get<0>(Median(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), Norms::GetNorm<TDataType>(rNormType)));
        }

        template <class TDataType>
        std::tuple<double, double, std::vector<double>, std::vector<IndexType>, std::vector<double>, std::vector<double>, std::vector<double>> static GetNormDistribution(
            const ModelPart& rModelPart,
            const Variable<TDataType>& rVariable,
            const std::string& rNormType,
            Parameters Params)
        {
            KRATOS_TRY

            const auto& values = Distribution(rModelPart, rVariable, GetDataLocation<TDataRetrievalFunctor<TContainerItemType>>(), Params, Norms::GetNorm<TDataType>(rNormType));

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
