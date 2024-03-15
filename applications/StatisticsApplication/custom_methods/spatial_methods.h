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
#include <vector>
#include <variant>

// External includes

// Project includes
#include "includes/communicator.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "expression/expression.h"
#include "expression/container_expression.h"

// Application includes
#include "custom_utilities/norms.h"
#include "custom_utilities/data_type_traits.h"

namespace Kratos
{
///@addtogroup RANSApplication
///@{

///@name Kratos Globals
///@{

class KRATOS_API(STATISTICS_APPLICATION) SpatialMethods
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

        using IndicesType = std::conditional_t<std::is_arithmetic_v<TDataType>, IndexType, std::vector<IndexType>>;

        using RawDataType = typename DataTypeTrait::RawDataType;

        KRATOS_CLASS_POINTER_DEFINITION(DistributionInfo);

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

        std::vector<IndicesType> GetGroupNumberOfValues() const { return mGroupNumberOfValues; }

        std::vector<TDataType> GetGroupValueDistributionPercentage() const { return mGroupValueDistributionPercentage; }

        std::vector<TDataType> GetGroupMeans() const { return mGroupMean; }

        std::vector<TDataType> GetGroupVariances() const { return mGroupVariance; }

        ///@}
        ///@name Private member variables
        ///@{

        TDataType mMin;

        TDataType mMax;

        std::vector<TDataType> mGroupUpperValues;

        std::vector<IndicesType> mGroupNumberOfValues;

        std::vector<TDataType> mGroupValueDistributionPercentage;

        std::vector<TDataType> mGroupMean;

        std::vector<TDataType> mGroupVariance;

        ///@}
    };

    ///@}
    ///@name Dependent types for expressions
    ///@{

    using ContainerExpressionType = std::variant<
                                            ContainerExpression<ModelPart::NodesContainerType>,
                                            ContainerExpression<ModelPart::ConditionsContainerType>,
                                            ContainerExpression<ModelPart::ElementsContainerType>
                                        >;

    using ExpressionReturnType = std::variant<
                                        double,
                                        array_1d<double, 3>,
                                        array_1d<double, 4>,
                                        array_1d<double, 6>,
                                        array_1d<double, 9>,
                                        Vector,
                                        Matrix
                                    >;

    using ExpressionReturnTypeWithIndices = std::variant<
                                                std::tuple<double, IndexType>,
                                                std::tuple<array_1d<double, 3>, std::vector<IndexType>>,
                                                std::tuple<array_1d<double, 4>, std::vector<IndexType>>,
                                                std::tuple<array_1d<double, 6>, std::vector<IndexType>>,
                                                std::tuple<array_1d<double, 9>, std::vector<IndexType>>,
                                                std::tuple<Vector, std::vector<IndexType>>,
                                                std::tuple<Matrix, std::vector<IndexType>>
                                            >;

    using ExpressionDistributionReturnType = std::variant<
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

    static ExpressionReturnType Sum(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator);

    static ExpressionReturnType Sum(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator,
        const Norms::AllNormTypes& rNorm);

    static ExpressionReturnType Sum(
        const ContainerExpressionType& rContainerExpression);

    static ExpressionReturnType Sum(
        const ContainerExpressionType& rContainerExpression,
        const Norms::AllNormTypes& rNorm);

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

    static ExpressionReturnType Mean(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator);

    static ExpressionReturnType Mean(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator,
        const Norms::AllNormTypes& rNorm);

    static ExpressionReturnType Mean(
        const ContainerExpressionType& rContainerExpression);

    static ExpressionReturnType Mean(
        const ContainerExpressionType& rContainerExpression,
        const Norms::AllNormTypes& rNorm);

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

    static ExpressionReturnType RootMeanSquare(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator);

    static ExpressionReturnType RootMeanSquare(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator,
        const Norms::AllNormTypes& rNorm);

    static ExpressionReturnType RootMeanSquare(
        const ContainerExpressionType& rContainerExpression);

    static ExpressionReturnType RootMeanSquare(
        const ContainerExpressionType& rContainerExpression,
        const Norms::AllNormTypes& rNorm);

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

    static std::tuple<ExpressionReturnType, ExpressionReturnType> Variance(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator);

    static std::tuple<ExpressionReturnType, ExpressionReturnType> Variance(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator,
        const Norms::AllNormTypes& rNorm);

    static std::tuple<ExpressionReturnType, ExpressionReturnType> Variance(
        const ContainerExpressionType& rContainerExpression);

    static std::tuple<ExpressionReturnType, ExpressionReturnType> Variance(
        const ContainerExpressionType& rContainerExpression,
        const Norms::AllNormTypes& rNorm);

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

    static ExpressionReturnTypeWithIndices Min(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator);

    static ExpressionReturnTypeWithIndices Min(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator,
        const Norms::AllNormTypes& rNorm);

    static ExpressionReturnTypeWithIndices Min(
        const ContainerExpressionType& rContainerExpression);

    static ExpressionReturnTypeWithIndices Min(
        const ContainerExpressionType& rContainerExpression,
        const Norms::AllNormTypes& rNorm);

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

    static ExpressionReturnTypeWithIndices Max(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator);

    static ExpressionReturnTypeWithIndices Max(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator,
        const Norms::AllNormTypes& rNorm);

    static ExpressionReturnTypeWithIndices Max(
        const ContainerExpressionType& rContainerExpression);

    static ExpressionReturnTypeWithIndices Max(
        const ContainerExpressionType& rContainerExpression,
        const Norms::AllNormTypes& rNorm);

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

    static ExpressionReturnTypeWithIndices Median(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator);

    static ExpressionReturnTypeWithIndices Median(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator,
        const Norms::AllNormTypes& rNorm);

    static ExpressionReturnTypeWithIndices Median(
        const ContainerExpressionType& rContainerExpression);

    static ExpressionReturnTypeWithIndices Median(
        const ContainerExpressionType& rContainerExpression,
        const Norms::AllNormTypes& rNorm);

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

    static ExpressionDistributionReturnType Distribution(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator,
        Parameters Params);

    static ExpressionDistributionReturnType Distribution(
        const Expression& rExpression,
        const DataCommunicator& rDataCommunicator,
        Parameters Params,
        const Norms::AllNormTypes& rNorm);

    static ExpressionDistributionReturnType Distribution(
        const ContainerExpressionType& rContainerExpression,
        Parameters Params);

    static ExpressionDistributionReturnType Distribution(
        const ContainerExpressionType& rContainerExpression,
        Parameters Params,
        const Norms::AllNormTypes& rNorm);

    ///@}
};

///@}

///@} addtogroup block

} // namespace Kratos.
