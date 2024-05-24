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
    const double P,
    const double Beta)
    : mpSensorModelPart(&rSensorModelPart),
      mP(P),
      mBeta(Beta)
{
}

void SensorDistanceBoltzmannOperatorResponseUtils::Initialize()
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();
    mDistances.resize(n, n, false);
    mSensorDistance.resize(n);
    mDistances.clear();

    IndexPartition<IndexType>(n).for_each([&](const auto i) {
        for (IndexType j = i + 1; j < n; ++j) {
            // 1.0 is added so that distances are never zero, and the mP of distance
            // is always greater than 1.0
            // here we multiply everything by n so that, at the point of response value calculation
            // the total summation will be divided by "n" which is good for scaling
            // in the boltzmann operator.
            mDistances(i, j) = n * std::pow(1 + norm_2((mpSensorModelPart->NodesBegin() + i)->Coordinates() - (mpSensorModelPart->NodesBegin() + j)->Coordinates()), mP);
        }
    });

    KRATOS_CATCH("");
}

double SensorDistanceBoltzmannOperatorResponseUtils::CalculateValue()
{
    KRATOS_TRY

    const IndexType n = mDistances.size1();

    std::tie(mNumerator, mDenominator) = IndexPartition<IndexType>(n).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([&](const auto i) {
        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);

        double& sensor_distance = mSensorDistance[i];

        sensor_distance = 0.0;
        for (IndexType j = 0; j < i; ++j) {
            const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
            sensor_distance += sensor_status_i * sensor_status_j / mDistances(j, i);
        }
        for (IndexType j = i + 1; j < n; ++j) {
            const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
            sensor_distance += sensor_status_i * sensor_status_j / mDistances(i, j);
        }

        return std::make_tuple(sensor_distance * std::exp(mBeta * sensor_distance), std::exp(mBeta * sensor_distance));
    });

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        return mNumerator / mDenominator;
    } else {
        return 0.0;
    }

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorDistanceBoltzmannOperatorResponseUtils::CalculateGradient() const
{
    KRATOS_TRY

    const IndexType n = mDistances.size1();

    auto p_expression = LiteralFlatExpression<double>::Create(mpSensorModelPart->NumberOfNodes(), {});
    IndexPartition<IndexType>(mpSensorModelPart->NumberOfNodes()).for_each([&p_expression](const auto Index) {
        *(p_expression->begin() + Index) = 0.0;
    });

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        IndexPartition<IndexType>(n).for_each([&](const auto k) {
            double& value = *(p_expression->begin() + k);
            const double coeff = (1.0 + mBeta * mSensorDistance[k] - mNumerator * mBeta / mDenominator) * std::exp(mBeta * mSensorDistance[k]) / mDenominator;

            for (IndexType i = 0; i < k; ++i) {
                const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);

                double temp = coeff;
                temp += (1.0 + mBeta * mSensorDistance[i] - mNumerator * mBeta / mDenominator) * std::exp(mBeta * mSensorDistance[i]) / mDenominator;

                value += sensor_status_i * temp / mDistances(i, k);
            }

            for (IndexType i = k + 1; i < n; ++i) {
                const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);

                double temp = coeff;
                temp += (1.0 + mBeta * mSensorDistance[i] - mNumerator * mBeta / mDenominator) * std::exp(mBeta * mSensorDistance[i]) / mDenominator;

                value += sensor_status_i * temp / mDistances(k, i);
            }
        });
    }

    ContainerExpression<ModelPart::NodesContainerType> result(*mpSensorModelPart);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/