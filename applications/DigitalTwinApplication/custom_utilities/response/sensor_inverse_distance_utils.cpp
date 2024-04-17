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

// Application includes
#include "digital_twin_application_variables.h"

// Include base h
#include "sensor_inverse_distance_utils.h"

namespace Kratos {

SensorInverseDistanceResponseUtils::SensorInverseDistanceResponseUtils(
    ModelPart& rSensorModelPart,
    const double P)
    : mpSensorModelPart(&rSensorModelPart),
      mP(P)
{
}

void SensorInverseDistanceResponseUtils::Initialize()
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();
    mDistances.resize(n, n, false);
    mDistances.clear();

    IndexPartition<IndexType>(n).for_each([&](const auto i) {
        for (IndexType j = i + 1; j < n; ++j) {
            mDistances(i, j) = norm_2((mpSensorModelPart->NodesBegin() + i)->Coordinates() - (mpSensorModelPart->NodesBegin() + j)->Coordinates());
            mDistances(j, i) = mDistances(i, j);
        }
    });

    KRATOS_CATCH("");
}

double SensorInverseDistanceResponseUtils::CalculateValue() const
{
    KRATOS_TRY

    const IndexType n = mDistances.size1();

    const double summation = IndexPartition<IndexType>(n).for_each<SumReduction<double>>([&](const auto i) {
        double value = 0.0;
        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
        if (sensor_status_i > 1e-8) {
            for (IndexType j = i + 1; j < n; ++j) {
                const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
                if (sensor_status_j > 1e-8) {
                    value += std::exp(-mP * sensor_status_i * sensor_status_j * mDistances(i, j));
                }
            }
        }
        return value;
    });

    return std::log(summation) / mP;

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorInverseDistanceResponseUtils::CalculateGradient() const
{
    KRATOS_TRY

    const IndexType n = mDistances.size1();

    const double summation = IndexPartition<IndexType>(n).for_each<SumReduction<double>>([&](const auto i) {
        double value = 0.0;
        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
        if (sensor_status_i > 1e-8) {
            for (IndexType j = i + 1; j < n; ++j) {
                const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
                if (sensor_status_j > 1e-8) {
                    value += std::exp(-mP * sensor_status_i * sensor_status_j * mDistances(i, j));
                }
            }
        }
        return value;
    });

    auto p_expression = LiteralFlatExpression<double>::Create(n, {});

    IndexPartition<IndexType>(n).for_each([&](const auto k) {
        auto& value = *(p_expression->begin() + k);
        value = 0.0;
        const double sensor_status_k = (mpSensorModelPart->NodesBegin() + k)->GetValue(SENSOR_STATUS);
        if (sensor_status_k > 1e-8) {
            for (IndexType i = 0; i < k; ++i) {
                const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
                if (sensor_status_i > 1e-8) {
                    value -= std::exp(-mP * sensor_status_i * sensor_status_k * mDistances(i, k)) * sensor_status_i * mDistances(i, k);
                }
            }
            for (IndexType i = k + 1; i < n; ++i) {
                const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
                if (sensor_status_i > 1e-8) {
                    value -= std::exp(-mP * sensor_status_i * sensor_status_k * mDistances(i, k)) * sensor_status_i * mDistances(i, k);
                }
            }
        }
        value /= summation;
    });

    ContainerExpression<ModelPart::NodesContainerType> result(*mpSensorModelPart);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/