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
#include "custom_utilities/control_utils.h"
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
      mBoltzmannOperator(Beta)
{
}

void SensorDistanceBoltzmannOperatorResponseUtils::Initialize()
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();
    const IndexType m = ControlUtils::GetDistVectorSize(n);

    mDistances.resize(m, false);
    mBoltzmannOperator.GetData().resize(m, false);

    mMaxDistance = IndexPartition<IndexType>(m).for_each<MaxReduction<double>>([&](const auto iDist) {
        IndexType i, j;
        std::tie(i, j) = ControlUtils::GetPairIndicesFromDistIndex(n, iDist);
        const double distance = norm_2((mpSensorModelPart->NodesBegin() + i)->Coordinates() - (mpSensorModelPart->NodesBegin() + j)->Coordinates());
        mDistances[iDist] = distance;
        return distance;
    });

    IndexPartition<IndexType>(m).for_each([&](const auto iDist) {
        mDistances[iDist] = n * std::pow(1 + mDistances[iDist] / mMaxDistance, mP);
    });

    KRATOS_CATCH("");
}

double SensorDistanceBoltzmannOperatorResponseUtils::CalculateValue()
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();
    const IndexType m = ControlUtils::GetDistVectorSize(n);

    Vector& r_data = mBoltzmannOperator.GetData();

    IndexPartition<IndexType>(m).for_each([&](const auto iDist) {
        IndexType i, j;
        std::tie(i, j) = ControlUtils::GetPairIndicesFromDistIndex(n, iDist);

        const double sensor_status_i = (mpSensorModelPart->NodesBegin() + i)->GetValue(SENSOR_STATUS);
        const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
        r_data[iDist] = sensor_status_i * sensor_status_j / mDistances[iDist];
    });

    mBoltzmannOperator.CalculateCoefficients();

    return mBoltzmannOperator.GetValue();

    KRATOS_CATCH("");
}

ContainerExpression<ModelPart::NodesContainerType> SensorDistanceBoltzmannOperatorResponseUtils::CalculateGradient() const
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();
    const IndexType m = ControlUtils::GetDistVectorSize(n);

    auto p_expression = LiteralFlatExpression<double>::Create(mpSensorModelPart->NumberOfNodes(), {});

    IndexPartition<IndexType>(n).for_each(Vector(m), [&](const auto k, auto& rGradients) {
        rGradients.clear();

        for (IndexType j = 0; j < k; ++j) {
            const IndexType i_dist = ControlUtils::GetDistIndexFromPairIndices(n, j, k);
            const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
            rGradients[i_dist] = sensor_status_j / mDistances[i_dist];
        }

        for (IndexType j = k + 1; j < n; ++j) {
            const IndexType i_dist = ControlUtils::GetDistIndexFromPairIndices(n, k, j);
            const double sensor_status_j = (mpSensorModelPart->NodesBegin() + j)->GetValue(SENSOR_STATUS);
            rGradients[i_dist] = sensor_status_j / mDistances[i_dist];
        }

        *(p_expression->begin() + k) = mBoltzmannOperator.GetGradient(rGradients);
    });

    ContainerExpression<ModelPart::NodesContainerType> result(*mpSensorModelPart);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

} /* namespace Kratos.*/