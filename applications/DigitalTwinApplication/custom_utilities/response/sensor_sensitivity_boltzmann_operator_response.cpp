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
#include "expression/literal_flat_expression.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "digital_twin_application_variables.h"

// Include base h
#include "sensor_sensitivity_boltzmann_operator_response.h"

namespace Kratos {

template<class TContainerType>
SensorSensitivityBoltzmannOperatorResponseUtils<TContainerType>::SensorSensitivityBoltzmannOperatorResponseUtils(
    ModelPart& rSensorModelPart,
    const std::vector<typename ContainerExpression<TContainerType>::Pointer>& rSensitivityDistributions,
    const double Beta)
    : mpSensorModelPart(&rSensorModelPart),
      mBeta(Beta),
      mSensorSensitivities(rSensitivityDistributions)

{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mpSensorModelPart->NumberOfNodes() == mSensorSensitivities.size())
        << "Number of sensors in the " << mpSensorModelPart->FullName()
        << " and number of sensor distributions mismatch [ number of sensors = "
        << mpSensorModelPart->NumberOfNodes() << ", number of sensor distributions = "
        << mSensorSensitivities.size() << " ].";

    KRATOS_ERROR_IF(mSensorSensitivities.empty())
        << "Empty sensor masks list provided.";

    KRATOS_CATCH("");
}

template<class TContainerType>
double SensorSensitivityBoltzmannOperatorResponseUtils<TContainerType>::CalculateValue()
{
    KRATOS_TRY

    const auto& r_container = mSensorSensitivities.front()->GetContainer();
    mElementRedundancy.resize(r_container.size());

    IndexPartition<IndexType>(r_container.size()).for_each([&](const auto iEntity) {
        double& value = mElementRedundancy[iEntity];
        value = 0.0;

        for (IndexType i_sensor = 0; i_sensor < mSensorSensitivities.size(); ++i_sensor) {
            const double sensor_sensitivity = mSensorSensitivities[i_sensor]->GetExpression().Evaluate(iEntity, iEntity, 0);
            const double sensor_status = (mpSensorModelPart->NodesBegin() + i_sensor)->GetValue(SENSOR_STATUS);
            value += sensor_status * sensor_sensitivity;
        }

        value /= mSensorSensitivities.size();
    });

    std::tie(mNumerator, mDenominator) = IndexPartition<IndexType>(r_container.size()).for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>([&](const auto iEntity) {
        const double value = std::exp(mBeta * mElementRedundancy[iEntity]);
        return std::make_tuple(mElementRedundancy[iEntity] * value, value);
    });

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        return mNumerator / mDenominator;
    } else {
        return 0.0;
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<ModelPart::NodesContainerType> SensorSensitivityBoltzmannOperatorResponseUtils<TContainerType>::CalculateGradient() const
{
    KRATOS_TRY

    const IndexType n = mpSensorModelPart->NumberOfNodes();

    const auto& r_container = mSensorSensitivities.front()->GetContainer();

    auto p_expression = LiteralFlatExpression<double>::Create(n, {});

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        const double coeff_1 = 1.0 / (mDenominator) ;
        const double coeff_2 = mNumerator / (mDenominator * mDenominator);
        IndexPartition<IndexType>(n).for_each([&](const auto k) {
            const auto& r_sensor_sensitivity_expression = mSensorSensitivities[k]->GetExpression();
            double numerator_derivative = 0.0;
            double denominator_derivative = 0.0;

            for (IndexType i_entity = 0; i_entity < r_container.size(); ++i_entity) {
                const double value = std::exp(mBeta * mElementRedundancy[i_entity]);
                const double m_kj = r_sensor_sensitivity_expression.Evaluate(i_entity, i_entity, 0);
                numerator_derivative += value * (mBeta * mElementRedundancy[i_entity] + 1) * m_kj / n;
                denominator_derivative += mBeta * value * m_kj / n;
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

// template instantiations
template class KRATOS_API(DIGITAL_TWIN_APPLICATION) SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::NodesContainerType>;
template class KRATOS_API(DIGITAL_TWIN_APPLICATION) SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::ConditionsContainerType>;
template class KRATOS_API(DIGITAL_TWIN_APPLICATION) SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::ElementsContainerType>;

} /* namespace Kratos.*/