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
#include <variant>

// External includes

// Project includes
#include "includes/node.h"
#include "geometries/geometry.h"

// Application includes
#include "custom_sensors/sensor_view.h"

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorUtils
{
public:
    ///@name Type definitions
    ///@{

    using SensorViewType = std::variant<
                                    SensorView<ModelPart::NodesContainerType>::Pointer,
                                    SensorView<ModelPart::ConditionsContainerType>::Pointer,
                                    SensorView<ModelPart::ElementsContainerType>::Pointer
                                >;


    ///@}
    ///@name Public static operations
    ///@{

    static bool IsPointInGeometry(
        const Point& rPoint,
        const Geometry<Node>& rGeometry);

    static SensorViewType CreateSensorView(
        Sensor::Pointer pSensor,
        const std::string& rExpressionName);

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/