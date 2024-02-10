//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/data_communicator.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) OptimizationUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    template<class TContainerType>
    static GeometryData::KratosGeometryType GetContainerEntityGeometryType(
        const TContainerType& rContainer,
        const DataCommunicator& rDataCommunicator);

    template<class TContainerType, class TDataType>
    static bool IsVariableExistsInAllContainerProperties(
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const DataCommunicator& rDataCommunicator);

    template<class TContainerType, class TDataType>
    static bool IsVariableExistsInAtLeastOneContainerProperties(
        const TContainerType& rContainer,
        const Variable<TDataType>& rVariable,
        const DataCommunicator& rDataCommunicator);

    template<class TContainerType>
    static void CreateEntitySpecificPropertiesForContainer(
        ModelPart& rModelPart,
        TContainerType& rContainer);

    template<class TDataType>
    static IndexType GetVariableDimension(
        const Variable<TDataType>& rVariable,
        const IndexType DomainSize);

    static void CopySolutionStepVariablesList(
        ModelPart& rDestinationModelPart,
        const ModelPart& rOriginModelPart);

    static IndexType Factorial(const IndexType N);

    static IndexType NChooseK(
        const IndexType N,
        const IndexType K);

    template<class TContainerType>
    static ContainerExpression<TContainerType> SmoothClamp(
        const ContainerExpression<TContainerType>& rInputExpression,
        const double MinValue,
        const double MaxValue,
        const IndexType NumberOfContinuousDerivatives);

    template<class TContainerType>
    static ContainerExpression<TContainerType> SmoothClampGradient(
        const ContainerExpression<TContainerType>& rInputExpression,
        const double MinValue,
        const double MaxValue,
        const IndexType NumberOfContinuousDerivatives);

    ///@}
};

///@}
}