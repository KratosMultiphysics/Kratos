//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/model_part.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) ControlUtils
{
public:
    ///@name Public static operations
    ///@{

    template<class TContainerType>
    static void AssignEquivalentProperties(
        TContainerType& rSourceContainer,
        TContainerType& rDestinationContainer);

    template<class TContainerType>
    static void ClipContainerExpression(
        ContainerExpression<TContainerType>& rContainerExpression,
        const double Min,
        const double Max);

    template<class TContainerType>
    static void GetIntegrationPoints(
        std::vector<Point>& rOutput,
        const TContainerType& rContainer);

    template<class TContainerType>
    static void GetIntegrationPointAreas(
        std::vector<double>& rOutput,
        const TContainerType& rContainer);

    template<class EntityType, class TDataType>
    static void EvaluateAtPoints(
        std::vector<TDataType>& rOutput,
        const Variable<TDataType>& rVariable,
        ModelPart& rModelPart,
        const std::vector<Point>& rCoordinates);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/