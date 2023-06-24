//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <variant>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) SigmoidalProjectionUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    using GeometryType = ModelPart::ElementType::GeometryType;

    ///@}
    ///@name Static operations
    ///@{

    template<class TContainerType>
    static ContainerExpression<TContainerType> ProjectForward(
        ContainerExpression<TContainerType>& rInputExpression,
        std::vector<double> xValues,
        std::vector<double> yValues,
        double beta,
        int penalFac);

    template<class TContainerType>
    static ContainerExpression<TContainerType> ProjectBackward(
        ContainerExpression<TContainerType>& rInputExpression,
        std::vector<double> xValues,
        std::vector<double> yValues,
        double beta,
        int penalFac);

    template<class TContainerType>
    static ContainerExpression<TContainerType> ComputeFirstDerivative(
        ContainerExpression<TContainerType>& rInputExpression,
        std::vector<double> xValues,
        std::vector<double> yValues,
        double beta,
        int penalFac);

    ///@}
};

///@}
}