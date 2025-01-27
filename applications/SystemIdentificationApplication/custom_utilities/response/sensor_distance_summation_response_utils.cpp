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

// System includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "expression/literal_flat_expression.h"

// Application includes
#include "system_identification_application_variables.h"

// Include base h
#include "sensor_distance_summation_response_utils.h"

namespace Kratos
{

double SensorDistanceSummationResponseUtils::CalculateValue(
    ModelPart &rModelPart,
    const DistanceMatrix &rDistanceMatrix)
{
    KRATOS_TRY

    return IndexPartition<IndexType>(rDistanceMatrix.GetEntriesSize()).for_each<SumReduction<double>>([&rModelPart, &rDistanceMatrix](const auto Index) {
        const auto& index_pair = rDistanceMatrix.GetIndexPair(Index);

        const auto i_index = std::get<0>(index_pair);
        const auto j_index = std::get<1>(index_pair);

        return (rModelPart.NodesBegin() + i_index)->GetValue(SENSOR_STATUS) * (rModelPart.NodesBegin() + j_index)->GetValue(SENSOR_STATUS) * rDistanceMatrix.GetDistance(Index);
    });

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorDistanceSummationResponseUtils::CalculateGradient(
    ModelPart &rModelPart,
    const DistanceMatrix &rDistanceMatrix)
{
    KRATOS_TRY

    auto p_expression = LiteralFlatExpression<double>::Create(rModelPart.NumberOfNodes(), {});

    IndexPartition<IndexType>(rModelPart.NumberOfNodes()).for_each([&p_expression, &rModelPart, &rDistanceMatrix](const auto k) {
        double& r_value = *(p_expression->begin() + k);
        r_value = 0.0;
        for (IndexType i = 0; i < rModelPart.NumberOfNodes(); ++i) {
            r_value += (rModelPart.NodesBegin() + i)->GetValue(SENSOR_STATUS) * rDistanceMatrix.GetDistance(i, k);
        }
    });

    ContainerExpression<ModelPart::NodesContainerType> result(rModelPart);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} // namespace Kratos