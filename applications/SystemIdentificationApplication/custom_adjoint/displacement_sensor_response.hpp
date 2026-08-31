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


/// @brief Class representing @ref DISPLACEMENT "displacement" measurements at a physical point.
/// @ingroup adjoints
class KRATOS_API(KRATOS_CORE) DisplacementSensorResponse final : public SensorResponse {
public:
    KRATOS_CLASS_POINTER_DEFINITION(DisplacementSensorResponse);

    using SensorResponse::SensorResponse;

    /// @copydoc SensorResponse::Create
    SensorResponse::Pointer Create(
        const ModelPart& rDomainModelPart,
        const ModelPart& rSensorModelPart,
        IndexType Id,
        Parameters SensorParameters) const override;

    /// @copydoc ResponseFunction::ComputeValue(const Element&,const ProcessInfo&,int) const
    [[nodiscard]] double ComputeValue(
        const Element& rElement,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const override;

    using SensorResponse::ComputeValue;

    /// @copydoc ResponseFunction::GetStateVariables(std::vector<IAdjoint::DynamicVariable>&,const Element&,const ProcessInfo&) const
    void GetStateVariables(
        std::vector<IAdjoint::DynamicVariable>& rVariables,
        const Element& rElement,
        const ProcessInfo& rProcessInfo) const override;

    using SensorResponse::GetStateVariables;

    /// @copydoc ResponseFunction::ComputeDerivative(Vector&,const Element&,std::span<const IAdjoint::DynamicVariable>,const ProcessInfo&,int)
    void ComputeDerivative(
        Vector& rOutput,
        const Element& rElement,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const override;

    using SensorResponse::ComputeDerivative;

    Parameters GetDefaultParameters() const override;

private:
    array_1d<double,3> mDirection;
}; // class DisplacementSensorResponse


} // namespace Kratos
