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

// --- Kratos Core Includes ---
#include "adjoint/sensor_response.hpp"


namespace Kratos {


class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorAggregateResponse : public ResponseFunction {
public:

    KRATOS_CLASS_POINTER_DEFINITION(SensorAggregateResponse);

    SensorAggregateResponse() noexcept;

    SensorAggregateResponse(double Exponent) noexcept;

    void AddSensor(SensorResponse::Pointer pSensor);

    /// @copydoc ResponseFunction::ComputeValue
    double ComputeValue(
        const ModelPart& rModelPart,
        int iBuffer) const override;

    /// @copydoc ResponseFunction::GetStateVariables(std::vector<IAdjoint::DynamicVariable>&,const Element&,const ProcessInfo&) const
    void GetStateVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Element& rElement,
        const ProcessInfo& rProcessInfo) const override;

    /// @copydoc ResponseFunction::GetStateVariables(std::vector<IAdjoint::DynamicVariable>&,const Condition&,const ProcessInfo&) const
    void GetStateVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Condition& rCondition,
        const ProcessInfo& rProcessInfo) const override;

    /// @copydoc ResponseFunction::ComputeDerivative(Vector&,const Element&,std::span<const IAdjoint::DynamicVariable>,const ProcessInfo&,int)
    void ComputeDerivative(
        Vector& rOutput,
        const Element& rElement,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const override;

    /// @copydoc ResponseFunction::ComputeDerivative(Vector&,const Condition&,std::span<const IAdjoint::DynamicVariable>,const ProcessInfo&,int)
    void ComputeDerivative(
        Vector& rOutput,
        const Condition& rCondition,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const override;

private:
    double mExponent;

    std::vector<SensorResponse::Pointer> mSensors;
};


} // namespace Kratos
