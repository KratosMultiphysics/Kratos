// KRATOS Multiphysics - IgaApplication
// License: BSD License

#pragma once

#include "custom_conditions/gap_sbm_load_solid_condition.h"

namespace Kratos
{

/// Gap-SBM Neumann condition using a finite-strain constitutive response.
class KRATOS_API(IGA_APPLICATION) GapSbmLoadSolidGreenLagrangeCondition
    : public GapSbmLoadSolidCondition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GapSbmLoadSolidGreenLagrangeCondition);

    GapSbmLoadSolidGreenLagrangeCondition() = default;
    GapSbmLoadSolidGreenLagrangeCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : GapSbmLoadSolidCondition(NewId, pGeometry) {}
    GapSbmLoadSolidGreenLagrangeCondition(IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : GapSbmLoadSolidCondition(NewId, pGeometry, pProperties) {}

    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& rNodes,
        PropertiesType::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    std::string Info() const override
    {
        return "GapSbmLoadSolidGreenLagrangeCondition #" + std::to_string(Id());
    }

private:
    void CalculateGreenKinematics(Matrix& rGradients, Matrix& rF, Matrix& rB);
    void SetFiniteStrainParameters(const Matrix& rF, ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rVariables, bool ComputeTangent) const;
};

} // namespace Kratos
