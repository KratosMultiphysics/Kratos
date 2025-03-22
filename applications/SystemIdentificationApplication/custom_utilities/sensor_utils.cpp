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

// System includes

// External includes

// Project includes

// Application includes

// Include base h
#include "sensor_utils.h"

namespace Kratos {

bool SensorUtils::IsPointInGeometry(
    const Point& rPoint,
    const Geometry<Node>& rGeometry)
{
    Point result;
    return rGeometry.IsInside(rPoint, result);
}

SensorUtils::SensorViewType SensorUtils::CreateSensorView(
    Sensor::Pointer pSensor,
    const std::string& rExpressionName)
{
    KRATOS_TRY

    auto p_expression = pSensor->GetContainerExpression(rExpressionName);
    if (std::holds_alternative<ContainerExpression<ModelPart::NodesContainerType>::Pointer>(p_expression)) {
        return std::make_shared<SensorView<ModelPart::NodesContainerType>>(pSensor, rExpressionName);
    } else if (std::holds_alternative<ContainerExpression<ModelPart::ConditionsContainerType>::Pointer>(p_expression)) {
        return std::make_shared<SensorView<ModelPart::ConditionsContainerType>>(pSensor, rExpressionName);
    } else if (std::holds_alternative<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>(p_expression)) {
        return std::make_shared<SensorView<ModelPart::ElementsContainerType>>(pSensor, rExpressionName);
    } else {
        KRATOS_ERROR << "Unsupported expression type.";
        return std::make_shared<SensorView<ModelPart::ElementsContainerType>>(pSensor, rExpressionName);
    }

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/