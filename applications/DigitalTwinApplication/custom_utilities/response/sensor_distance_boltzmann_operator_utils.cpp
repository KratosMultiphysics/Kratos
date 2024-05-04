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
#include "sensor_distance_boltzmann_operator_utils.h"

namespace Kratos {

SensorDistanceBoltzmannOperatorResponseUtils::SensorDistanceBoltzmannOperatorResponseUtils(
    ModelPart& rSensorModelPart,
    const double Beta)
    : mpSensorModelPart(&rSensorModelPart),
      mBeta(Beta)
{
}

void SensorDistanceBoltzmannOperatorResponseUtils::Initialize()
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();
    mDistances.resize(n, n, false);
    mDistances.clear();

    const double max_distance = IndexPartition<IndexType>(n).for_each<MaxReduction<double>>([&](const auto i) {
        double max_distance = 0.0;
        for (IndexType j = i + 1; j < n; ++j) {
            const auto distance = norm_2((mpSensorModelPart->NodesBegin() + i)->Coordinates() - (mpSensorModelPart->NodesBegin() + j)->Coordinates());
            mDistances(i, j) = distance;

            if (max_distance < distance) {
                max_distance = distance;
            }
        }
        return max_distance;
    });

    // For IEEE-compatible type double, overflow is guaranteed if 709.8 < num
    // and underflow is guaranteed if num < -708.4. hence we calculate the optimum scaling
    // to not to have overflow or under flow
    mScaling = max_distance;

    IndexPartition<IndexType>(n).for_each([&](const auto i) {
        for (IndexType j = i + 1; j < n; ++j) {
            const double distance_ratio = mDistances(i, j) / mScaling;
            mDistances(j, i) = std::exp(mBeta * distance_ratio);
            mDistances(i, j) = distance_ratio * mDistances(j, i);
        }
    });

    KRATOS_CATCH("");
}

double SensorDistanceBoltzmannOperatorResponseUtils::CalculateValue()
{
    KRATOS_TRY

    const IndexType n = mDistances.size1();

    std::tie(mNumerator, mDenominator) = IndexPartition<IndexType>(n).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([&](const auto i) {
        double numerator = 0.0;
        double denominator = 0.0;
        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
        for (IndexType j = i + 1; j < n; ++j) {
            const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
            numerator += sensor_status_i * sensor_status_j * mDistances(i, j);
            denominator += sensor_status_i * sensor_status_j * mDistances(j, i);
        }
        return std::make_tuple(numerator, denominator);
    });

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        return mNumerator * mScaling / mDenominator;
    } else {
        return 0.0;
    }

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorDistanceBoltzmannOperatorResponseUtils::CalculateGradient() const
{
    KRATOS_TRY

    const IndexType n = mDistances.size1();

    auto p_expression = LiteralFlatExpression<double>::Create(n, {});

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        const double coeff_1 = mScaling / mDenominator;
        const double coeff_2 = mNumerator * mScaling / (mDenominator * mDenominator);
        IndexPartition<IndexType>(n).for_each([&](const auto k) {
            double numerator_derivative = 0.0;
            double denominator_derivative = 0.0;

            for (IndexType i = 0; i < k; ++i) {
                const double sensor_status = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
                numerator_derivative += sensor_status * mDistances(i, k);
                denominator_derivative += sensor_status * mDistances(k, i);
            }

            for (IndexType i = k + 1; i < n; ++i) {
                const double sensor_status = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
                numerator_derivative += sensor_status * mDistances(k, i);
                denominator_derivative += sensor_status * mDistances(i, k);
            }

            *(p_expression->begin() + k) = coeff_1 * numerator_derivative - denominator_derivative * coeff_2;
        });
    } else {
        IndexPartition<IndexType>(n).for_each([&](const auto Index) {
            *(p_expression->begin() + Index) = 0.0;
        });
    }

    ContainerExpression<ModelPart::NodesContainerType> result(*mpSensorModelPart);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/