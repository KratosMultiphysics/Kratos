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
#include <algorithm>

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

void SensorUtils::AssignSensorIds(std::vector<Sensor::Pointer>& rSensorsList)
{
    KRATOS_TRY

    IndexPartition<int>(rSensorsList.size()).for_each([&](const auto Index) {
        (*(rSensorsList.begin() + Index))->SetValue(SENSOR_ID, Index + 1);
    });

    KRATOS_CATCH("");
}

IndexType SensorUtils::GetMostDistanced(
    const std::vector<Sensor::Pointer>& rOriginSensors,
    const std::vector<Sensor::Pointer>& rTestSensors)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(rOriginSensors.empty())
        << "No origin sensors provided.";

    std::vector<double> distances(rTestSensors.size());

    IndexPartition<IndexType>(rTestSensors.size()).for_each([&rOriginSensors, &rTestSensors, &distances](const auto Index) {
        auto& p_test_sensor = rTestSensors[Index];
        double min_distance = std::numeric_limits<double>::max();
        for (const auto& p_origin_sensor : rOriginSensors) {
            min_distance = std::min(min_distance, norm_2(p_origin_sensor->GetLocation() - p_test_sensor->GetLocation()));
        }
        distances[Index] = min_distance;
    });

    return std::distance(distances.begin(), std::max_element(distances.begin(), distances.end()));

    KRATOS_CATCH("");
}

} // namespace Kratos