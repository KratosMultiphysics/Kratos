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
#include "sensor_distance_variance_utils.h"

namespace Kratos {

SensorDistanceVarianceUtils::SensorDistanceVarianceUtils(
    ModelPart& rSensorModelPart)
    : mpSensorModelPart(&rSensorModelPart)
{
}

void SensorDistanceVarianceUtils::Initialize()
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();
    mDistances.resize(n, n, false);
    mDistances.clear();

    IndexPartition<IndexType>(n).for_each([&](const auto i) {
        for (IndexType j = i + 1; j < n; ++j) {
            const auto distance = norm_2((mpSensorModelPart->NodesBegin() + i)->Coordinates() - (mpSensorModelPart->NodesBegin() + j)->Coordinates());
            mDistances(i, j) = distance;
            mDistances(j, i) = mDistances(i, j);
        }
    });

    KRATOS_CATCH("");
}

double SensorDistanceVarianceUtils::CalculateValue()
{
    KRATOS_TRY

    const IndexType n = mDistances.size1();

    // first compute the mean
    const double sum = IndexPartition<IndexType>(n).for_each<SumReduction<double>>([&](const auto i) {
        double value = 0.0;
        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
            for (IndexType j = i + 1; j < n; ++j) {
                const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
                value += sensor_status_i * sensor_status_j * mDistances(i, j);
            }
        return value;
    });

    mMean = sum / (n * (n - 1) / 2);

    const double mean_square = std::pow(mMean, 2);

    const double variance = IndexPartition<IndexType>(n).for_each<SumReduction<double>>([&](const auto i) {
        double value = 0.0;
        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
            for (IndexType j = i + 1; j < n; ++j) {
                const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
                value += std::pow(sensor_status_i * sensor_status_j * mDistances(i, j), 2) - mean_square;
            }
        return value;
    });

    return variance / (n * (n - 1) / 2);

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorDistanceVarianceUtils::CalculateGradient() const
{
    KRATOS_TRY

    const IndexType n = mDistances.size1();

    const double coeff = 2.0 / (n * (n - 1) / 2);

    auto p_expression = LiteralFlatExpression<double>::Create(n, {});

    IndexPartition<IndexType>(n).for_each([&](const auto k) {
        const double sensor_status_k = (mpSensorModelPart->NodesBegin() + k)->GetValue(SENSOR_STATUS);

        auto& value = *(p_expression->begin() + k);
        value = 0.0;
        for (IndexType i = 0; i < k; ++i) {
            const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
            const double d_ik = mDistances(i, k);
            const double x_ik = sensor_status_i * sensor_status_k * d_ik;
            value += sensor_status_i * d_ik * (x_ik - mMean);
        }
        for (IndexType i = k + 1; i < n; ++i) {
            const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
            const double d_ik = mDistances(i, k);
            const double x_ik = sensor_status_i * sensor_status_k * d_ik;
            value += sensor_status_i * d_ik * (x_ik - mMean);
        }
        value *= coeff;
    });

    ContainerExpression<ModelPart::NodesContainerType> result(*mpSensorModelPart);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/