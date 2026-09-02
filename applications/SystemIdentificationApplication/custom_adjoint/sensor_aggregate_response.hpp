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

// --- STL Includes ---
#include <memory>
#include <span>


namespace Kratos {


/// @brief Response defining an @f$L_p@f$ norm of measurement gaps over a list of sensors.
/// @ingroup adjoints
class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) SensorAggregateResponse : public ResponseFunction {
public:

    KRATOS_CLASS_POINTER_DEFINITION(SensorAggregateResponse);

    SensorAggregateResponse() noexcept;

    SensorAggregateResponse(std::size_t Exponent);

    SensorAggregateResponse(SensorAggregateResponse&& rRhs) noexcept;

    SensorAggregateResponse(const SensorAggregateResponse& rRhs);

    ~SensorAggregateResponse();

    SensorAggregateResponse& operator=(SensorAggregateResponse&& rRhs) noexcept;

    SensorAggregateResponse& operator=(const SensorAggregateResponse& rRhs);

    /// @copydoc ResponseFunction::ComputeCache
    void ComputeCache(const ModelPart& rModelPart) override;

    /// @copydoc ResponseFunction::ClearCache
    void ClearCache() override;

    /// @brief Add a sensor to consider output from.
    void AddSensors(std::span<const SensorResponse::Pointer> Sensors);

    /// @copydoc ResponseFunction::ComputeInnerValue(const Element&,const ProcessInfo&,int) const
    [[nodiscard]] double ComputeInnerValue(
        const Element& rElement,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const override;

    /// @copydoc ResponseFunction::ComputeInnerValue(const Condition&,const ProcessInfo&,int) const
    [[nodiscard]] double ComputeInnerValue(
        const Condition& Condition,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const final override;

    /// @copydoc ResponseFunction::ComputeValue
    [[nodiscard]] double ComputeValue(
        double InnerSum,
        const ProcessInfo& rProcessInfo,
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

    /// @copydoc ResponseFunction::ComputeInnerDerivative(Vector&,const Element&,std::span<const IAdjoint::DynamicVariable>,const ProcessInfo&,int)
    void ComputeInnerDerivative(
        Vector& rOutput,
        const Element& rElement,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const override;

    /// @copydoc ResponseFunction::ComputeInnerDerivative(Vector&,const Condition&,std::span<const IAdjoint::DynamicVariable>,const ProcessInfo&,int)
    void ComputeInnerDerivative(
        Vector& rOutput,
        const Condition& rCondition,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const override;

    /// @copydoc ResponseFunction::ComputeDerivative
    void ComputeDerivative(
        Vector& rDerivative,
        std::span<const IAdjoint::DynamicVariable> Variables,
        const ProcessInfo& rProcessInfo,
        int iBuffer) const override;

private:
    struct Impl;
    std::unique_ptr<Impl> mpImpl;
};


} // namespace Kratos
