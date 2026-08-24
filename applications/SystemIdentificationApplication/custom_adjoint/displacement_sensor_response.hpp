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


class KRATOS_API(KRATOS_CORE) DisplacementSensorResponse final : public SensorResponse {
public:
    using SensorResponse::SensorResponse;

    /// @copydoc SensorResponse::Create
    SensorResponse::Pointer Create(
        const ModelPart& rDomainModelPart,
        const ModelPart& rSensorModelPart,
        IndexType Id,
        Parameters SensorParameters) const override;

    Parameters GetSensorParameters() const override;

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

    Parameters GetDefaultParameters() const override;

private:
    array_1d<double,3> mDirection;
}; // class DisplacementSensorResponse


} // namespace Kratos
