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

// Application includes
#include "custom_sensors/sensor.h"

namespace Kratos {
///@name Kratos Classes
///@{

class SensorUtils {
public:
    ///@name Public static operations
    ///@{

    /**
     * @brief Set the Sensor to the given node.
     */
    static void SetSensor(ModelPart::NodeType& rNode, Sensor::Pointer pSensor);

    /**
     * @brief Get the Sensor from the node.
     */
    static Sensor::Pointer GetSensor(const ModelPart::NodeType& rNode);

    /**
     * @brief Create new node for the sensors and assign the respective sensor in the model part.
     */
    static void AddSensors(
        ModelPart& rModelPart,
        const std::vector<Sensor::Pointer>& rSensorsList);

    /**
     * @brief Assigns sensor ids to the given sensor list.
     */
    static void AssignSensorIds(std::vector<Sensor::Pointer>& rSensorsList);

    ///@}
};

} // namespace Kratos