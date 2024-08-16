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
#include "includes/define.h"
#include "includes/model_part.h"
#include "expression/container_expression.h"

// Application includes
#include "custom_utilities/sensor_mask_status.h"

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorCoverageResponseUtils
{
public:
    ///@name Public static operations
    ///@{

    static double CalculateValue(const SensorMaskStatus& rSensorMaskStatus);

    static ContainerExpression<ModelPart::NodesContainerType> CalculateGradient(const SensorMaskStatus& rSensorMaskStatus);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/