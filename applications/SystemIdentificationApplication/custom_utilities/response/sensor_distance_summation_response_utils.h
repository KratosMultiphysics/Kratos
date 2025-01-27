//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: SystemIdentificationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "expression/container_expression.h"

// Application includes
#include "custom_utilities/distance_matrix.h"

namespace Kratos
{

///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorDistanceSummationResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Static operations
    ///@{

    static double CalculateValue(
        ModelPart& rModelPart,
        const DistanceMatrix& rDistanceMatrix);

    static ContainerExpression<ModelPart::NodesContainerType> CalculateGradient(
        ModelPart& rModelPart,
        const DistanceMatrix& rDistanceMatrix);

    ///@}
};

///@}
}