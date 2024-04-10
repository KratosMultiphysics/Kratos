//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "expression/container_expression.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class SensorLocalizationResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorLocalizationResponseUtils);

    ///@}
    ///@name Life cycle
    ///@{

    SensorLocalizationResponseUtils(
        ModelPart& rSensorModelPart,
        const double P);

    ///@}
    ///@name Public operations
    ///@{

    double CalculateValue(const std::vector<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>& rMasksList);

    ContainerExpression<ModelPart::NodesContainerType> CalculateGradient(const std::vector<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>& rMasksList);

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart* mpSensorModelPart;

    const double mP;

    ///@}
    ///@name Private operations
    ///@{

    void ComputeClusterDifference(
        Matrix& rOutput,
        const std::vector<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>& rMasksList);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/