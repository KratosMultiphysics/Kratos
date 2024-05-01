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

// External includes

// Project includes
#include "expression/container_expression.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/sensor_mask_status.h"

namespace Kratos {
///@name Kratos Classes
///@{

class SensorCoverageResponseUtils
{
public:
    ///@name Public static operations
    ///@{

    template<class TContainerType>
    static double CalculateValue(const SensorMaskStatus<TContainerType>& rSensorMaskStatus);

    template<class TContainerType>
    static ContainerExpression<ModelPart::NodesContainerType> CalculateGradient(const SensorMaskStatus<TContainerType>& rSensorMaskStatus);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/