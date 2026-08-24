//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// --- Kratos SysId Includes ---
#include "custom_adjoint/sensor_aggregate_response.hpp"
#include "system_identification_application_variables.h"


namespace Kratos {


SensorAggregateResponse::SensorAggregateResponse() noexcept
    : SensorAggregateResponse(1.0)
{}


SensorAggregateResponse::SensorAggregateResponse(double Exponent) noexcept
    : mExponent(Exponent)
{}


double SensorAggregateResponse::ComputeValue(
    const ModelPart& rModelPart,
    int iBuffer) const {
    KRATOS_TRY
        double sum = 0.0;

        for (auto& p_sensor : mSensors) {
            const double sensor_value = p_sensor->ComputeValue(rModelPart, iBuffer);
            const double raw_sensor_error = sensor_value - p_sensor->GetNode()->GetValue(SENSOR_MEASURED_VALUE);
            const double current_sensor_error = std::abs(raw_sensor_error) < p_sensor->GetErrorThreshold() ? 0.0 : raw_sensor_error;
            const double normalized_sensor_error = current_sensor_error / p_sensor->GetNode()->GetValue(SENSOR_NORMALIZATION_FACTOR);

            //p_sensor->SetSensorValue(sensor_value);
            p_sensor->GetNode()->SetValue(SENSOR_ERROR, current_sensor_error);
            p_sensor->GetNode()->SetValue(SENSOR_RELATIVE_ERROR, current_sensor_error / p_sensor->GetNode()->GetValue(SENSOR_MEASURED_VALUE));

            sum += (std::pow(0.5 * std::pow(normalized_sensor_error, 2) * p_sensor->GetWeight(), mExponent));
        }

        return std::pow(sum, 1 / mExponent);
    KRATOS_CATCH("");
}


void SensorAggregateResponse::GetStateVariables(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    const Element& rElement,
    const ProcessInfo& rProcessInfo) const {
        KRATOS_TRY
            rVariables.clear();
            std::vector<IAdjoint::DynamicVariable> sensor_variables;
            for (const auto& p_sensor : mSensors) {
                sensor_variables.clear();
                p_sensor->GetStateVariables(sensor_variables, rElement, rProcessInfo);
                rVariables.insert(rVariables.end(), sensor_variables.begin(), sensor_variables.end());
            }

            std::sort(
                rVariables.begin(),
                rVariables.end(),
                [] (const IAdjoint::DynamicVariable& rLeft, const IAdjoint::DynamicVariable& rRight) {
                    if (rLeft.Key() != rRight.Key()) {
                        return rLeft.Key() < rRight.Key();
                    }
                    return rLeft.GetDynamicIndex() < rRight.GetDynamicIndex();
                });
            rVariables.erase(
                std::unique(
                    rVariables.begin(),
                    rVariables.end(),
                    [] (const IAdjoint::DynamicVariable& rLeft, const IAdjoint::DynamicVariable& rRight) {
                        return rLeft == rRight;
                    }),
                rVariables.end());
        KRATOS_CATCH("")
}


void SensorAggregateResponse::GetStateVariables(
    std::vector<IAdjoint::DynamicVariable>& rVariables,
    const Condition& rCondition,
    const ProcessInfo& rProcessInfo) const {
        KRATOS_TRY
            rVariables.clear();
            std::vector<IAdjoint::DynamicVariable> sensor_variables;
            for (const auto& p_sensor : mSensors) {
                sensor_variables.clear();
                p_sensor->GetStateVariables(sensor_variables, rCondition, rProcessInfo);
                rVariables.insert(rVariables.end(), sensor_variables.begin(), sensor_variables.end());
            }

            std::sort(
                rVariables.begin(),
                rVariables.end(),
                [] (const IAdjoint::DynamicVariable& rLeft, const IAdjoint::DynamicVariable& rRight) {
                    if (rLeft.Key() != rRight.Key()) {
                        return rLeft.Key() < rRight.Key();
                    }
                    return rLeft.GetDynamicIndex() < rRight.GetDynamicIndex();
                });
            rVariables.erase(
                std::unique(
                    rVariables.begin(),
                    rVariables.end(),
                    [] (const IAdjoint::DynamicVariable& rLeft, const IAdjoint::DynamicVariable& rRight) {
                        return rLeft == rRight;
                    }),
                rVariables.end());
        KRATOS_CATCH("")
}


void SensorAggregateResponse::ComputeDerivative(
    Vector& rOutput,
    const Element& rElement,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_ERROR;
}


void SensorAggregateResponse::ComputeDerivative(
    Vector& rOutput,
    const Condition& rCondition,
    std::span<const IAdjoint::DynamicVariable> Variables,
    const ProcessInfo& rProcessInfo,
    int iBuffer) const {
        KRATOS_ERROR;
}


} // namespace Kratos
