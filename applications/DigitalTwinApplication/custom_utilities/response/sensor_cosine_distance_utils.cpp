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

// System includes
#include <cmath>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "expression/literal_flat_expression.h"
#include "expression/expression_utils.h"

// Application includes
#include "digital_twin_application_variables.h"

// Include base h
#include "sensor_cosine_distance_utils.h"

namespace Kratos {

SensorCosineDistanceResponseUtils::SensorCosineDistanceResponseUtils(
    ModelPart& rSensorModelPart,
    const double P)
    : mpSensorModelPart(&rSensorModelPart),
      mP(P)
{
}

void SensorCosineDistanceResponseUtils::Initialize(const std::vector<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>& rMasksList)
{
    KRATOS_TRY

    const IndexType n = rMasksList.size();
    mDistances.resize(n, n, false);
    mDistances.clear();

    for (IndexType i = 0; i < n; ++i) {
        for (IndexType j = i + 1; j < n; ++j) {
            mDistances(i, j) = ExpressionUtils::InnerProduct(*rMasksList[i], *rMasksList[j]);
            mDistances(j, i) = mDistances(i, j);
        }
    };

    KRATOS_CATCH("");
}

double SensorCosineDistanceResponseUtils::CalculateValue() const
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();

    // calculating the p-norm of the distances to reduce the maximum if required
    const auto summation = IndexPartition<IndexType>(n).for_each<SumReduction<double>>([&](const auto i) {
        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
        double value = 0.0;
        for (IndexType j = i + 1; j < n; ++j) {
            const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
            value += std::pow(sensor_status_i * sensor_status_j * mDistances(i ,j), mP);
        }
        return value;
    });

    return std::pow(summation, 1 / mP);

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorCosineDistanceResponseUtils::CalculateGradient() const
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();

    const auto summation = IndexPartition<IndexType>(n).for_each<SumReduction<double>>([&](const auto i) {
        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
        double value = 0.0;
        for (IndexType j = i + 1; j < n; ++j) {
            const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
            value += std::pow(sensor_status_i * sensor_status_j * mDistances(i ,j), mP);
        }
        return value;
    });

    const double coeff = std::pow(summation, 1 / mP - 1) / mP;

    auto p_expression = LiteralFlatExpression<double>::Create(n, {});

    IndexPartition<IndexType>(n).for_each([&](const auto k) {
        auto& value = *(p_expression->begin() + k);
        value = 0.0;
        const double sensor_status_k = (mpSensorModelPart->NodesBegin() + k)->GetValue(SENSOR_STATUS);
            for (IndexType i = 0; i < k; ++i) {
                const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
                value += mP * std::pow(sensor_status_i * sensor_status_k * mDistances(i, k), mP - 1) * sensor_status_i * mDistances(i, k);
            }
            for (IndexType i = k + 1; i < n; ++i) {
                const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
                value += mP * std::pow(sensor_status_i * sensor_status_k * mDistances(i, k), mP - 1) * sensor_status_i * mDistances(i, k);
            }
        value *= coeff;
    });

    ContainerExpression<ModelPart::NodesContainerType> result(*mpSensorModelPart);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/