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
#include <limits>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "expression/literal_flat_expression.h"

// Application includes
#include "digital_twin_application_variables.h"

// Include base h
#include "sensor_distance_p_norm_utils.h"

namespace Kratos {

SensorDistancePNormResponseUtils::SensorDistancePNormResponseUtils(
    ModelPart& rSensorModelPart,
    const double P)
    : mpSensorModelPart(&rSensorModelPart),
      mP(P)
{
}

void SensorDistancePNormResponseUtils::Initialize()
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();
    mDistances.resize(n, n, false);
    mDistances.clear();

    IndexPartition<IndexType>(n).for_each([&](const auto i) {
        for (IndexType j = i + 1; j < n; ++j) {
            const auto distance = norm_2((mpSensorModelPart->NodesBegin() + i)->Coordinates() - (mpSensorModelPart->NodesBegin() + j)->Coordinates());
            if (distance > std::numeric_limits<double>::epsilon()) {
                mDistances(i, j) = std::pow(distance, mP);
                mDistances(j, i) = mDistances(i, j);
            }
        }
    });

    KRATOS_CATCH("");
}

double SensorDistancePNormResponseUtils::CalculateValue() const
{
    KRATOS_TRY

    const IndexType n = mDistances.size1();

    const double summation = IndexPartition<IndexType>(n).for_each<SumReduction<double>>([&](const auto i) {
        double value = 0.0;
        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
            for (IndexType j = i + 1; j < n; ++j) {
                const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
                value += sensor_status_i * sensor_status_j * mDistances(i, j);
            }
        return value;
    });

    return std::pow(summation, 1 / mP);

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorDistancePNormResponseUtils::CalculateGradient() const
{
    KRATOS_TRY

    const IndexType n = mDistances.size1();

    const double summation = IndexPartition<IndexType>(n).for_each<SumReduction<double>>([&](const auto i) {
        double value = 0.0;
        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
            for (IndexType j = i + 1; j < n; ++j) {
                const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
                value += sensor_status_i * sensor_status_j * mDistances(i, j);
            }
        return value;
    });

    const double coeff = std::pow(summation, 1 / mP - 1) / mP;

    auto p_expression = LiteralFlatExpression<double>::Create(n, {});

    IndexPartition<IndexType>(n).for_each([&](const auto k) {
        auto& value = *(p_expression->begin() + k);
        value = 0.0;
        for (IndexType i = 0; i < k; ++i) {
            const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
            value += sensor_status_i * mDistances(i, k);
        }
        for (IndexType i = k + 1; i < n; ++i) {
            const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
            value += sensor_status_i *  mDistances(i, k);
        }
        value *= coeff;
    });

    ContainerExpression<ModelPart::NodesContainerType> result(*mpSensorModelPart);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/