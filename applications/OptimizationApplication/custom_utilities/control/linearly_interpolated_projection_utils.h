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
#include <vector>

// Project includes
#include "includes/define.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) LinearlyInterpolatedProjectionUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    template<class TContainerType>
    static ContainerExpression<TContainerType> ProjectForward(
        const ContainerExpression<TContainerType>& rInputExpression,
        const std::vector<double>& rXValues,
        const std::vector<double>& rYValues,
        const double Beta,
        const int PenaltyFactor);

    template<class TContainerType>
    static ContainerExpression<TContainerType> ProjectBackward(
        const ContainerExpression<TContainerType>& rInputExpression,
        const std::vector<double>& rXValues,
        const std::vector<double>& rYValues,
        const double Beta,
        const int PenaltyFactor);

    template<class TContainerType>
    static ContainerExpression<TContainerType> CalculateForwardProjectionGradient(
        const ContainerExpression<TContainerType>& rInputExpression,
        const std::vector<double>& rXValues,
        const std::vector<double>& rYValues,
        const double Beta,
        const int PenaltyFactor);

    ///@}
};

///@}
}