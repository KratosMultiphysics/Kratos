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
#include "sensor_inverse_distance_summation_response_utils.h"

namespace Kratos
{

double SensorInverseDistanceSummationResponseUtils::CalculateValue(
    ModelPart &rModelPart,
    const double P,
    const DistanceMatrix &rDistanceMatrix)
{
    KRATOS_TRY

    return IndexPartition<IndexType>(rDistanceMatrix.GetEntriesSize()).for_each<SumReduction<double>>([P, &rModelPart, &rDistanceMatrix](const auto Index) {
        const auto& index_pair = rDistanceMatrix.GetIndexPair(Index);

        const auto i_index = std::get<0>(index_pair);
        const auto j_index = std::get<1>(index_pair);

        if (i_index != j_index) {
            return std::pow((rModelPart.NodesBegin() + i_index)->GetValue(SENSOR_STATUS), P) * std::pow((rModelPart.NodesBegin() + j_index)->GetValue(SENSOR_STATUS), P) / rDistanceMatrix.GetDistance(Index);
        } else {
            return 0.0;
        }
    });

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorInverseDistanceSummationResponseUtils::CalculateGradient(
    ModelPart &rModelPart,
    const double P,
    const DistanceMatrix &rDistanceMatrix)
{
    KRATOS_TRY

    auto p_expression = LiteralFlatExpression<double>::Create(rModelPart.NumberOfNodes(), {});

    IndexPartition<IndexType>(rModelPart.NumberOfNodes()).for_each([P, &p_expression, &rModelPart, &rDistanceMatrix](const auto k) {
        double& r_value = *(p_expression->begin() + k);
        r_value = 0.0;
        const double coeff = std::pow((rModelPart.NodesBegin() + k)->GetValue(SENSOR_STATUS), P - 1) * P;
        for (IndexType i = 0; i < rModelPart.NumberOfNodes(); ++i) {
            if (i != k) {
                r_value += std::pow((rModelPart.NodesBegin() + i)->GetValue(SENSOR_STATUS), P) * coeff / rDistanceMatrix.GetDistance(i, k);
            }
        }
    });

    ContainerExpression<ModelPart::NodesContainerType> result(rModelPart);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} // namespace Kratos