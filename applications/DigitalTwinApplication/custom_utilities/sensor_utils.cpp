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

// System includes

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "digital_twin_application_variables.h"

// Include base h
#include "sensor_utils.h"


namespace Kratos {

void SensorUtils::SetSensor(ModelPart::NodeType& rNode, Sensor::Pointer pSensor)
{
    rNode.Coordinates() = pSensor->GetLocation();
    rNode.SetValue(SENSOR, pSensor);
}

Sensor::Pointer SensorUtils::GetSensor(const ModelPart::NodeType& rNode)
{
    return rNode.GetValue(SENSOR);
}

void SensorUtils::AddSensors(
    ModelPart& rModelPart,
    const std::vector<Sensor::Pointer>& rSensorsList)
{
    KRATOS_TRY

    const auto max_id = rModelPart.GetCommunicator().GetDataCommunicator().MaxAll(block_for_each<MaxReduction<int>>(rModelPart.Nodes(), [](const auto& rNode) -> int {
        return rNode.Id();
    }));

    for (IndexType i = 0; i < rSensorsList.size(); ++i) {
        SetSensor(*(rModelPart.CreateNewNode(i + max_id + 1, 0, 0, 0)), rSensorsList[i]);
    }

    KRATOS_CATCH("");
}

} // namespace Kratos