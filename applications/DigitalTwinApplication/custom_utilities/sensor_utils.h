//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
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
#include "custom_sensors/sensor.h"

namespace Kratos {
///@name Kratos Classes
///@{

class SensorUtils {
public:
    ///@name Public static operations
    ///@{

    static void SetSensor(ModelPart::NodeType& rNode, Sensor::Pointer pSensor);

    static Sensor::Pointer GetSensor(const ModelPart::NodeType& rNode);

    static void AddSensors(
        ModelPart& rModelPart,
        const std::vector<Sensor::Pointer>& rSensorsList);

    ///@}
};

} // namespace Kratos