// KRATOS Multiphysics - IgaApplication
// License: BSD License

#pragma once

#include "custom_conditions/sbm_load_solid_condition.h"

namespace Kratos
{

/// Shifted boundary traction condition in the reference configuration.
class KRATOS_API(IGA_APPLICATION) SbmLoadSolidGreenLagrangeCondition : public SbmLoadSolidCondition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SbmLoadSolidGreenLagrangeCondition);

    SbmLoadSolidGreenLagrangeCondition() = default;
    SbmLoadSolidGreenLagrangeCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : SbmLoadSolidCondition(NewId, pGeometry) {}
    SbmLoadSolidGreenLagrangeCondition(IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : SbmLoadSolidCondition(NewId, pGeometry, pProperties) {}

    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& rNodes,
        PropertiesType::Pointer pProperties) const override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    std::string Info() const override
    {
        return "SbmLoadSolidGreenLagrangeCondition #" + std::to_string(Id());
    }

private:
    void CalculateReferenceGradients(Matrix& rSurrogateGradients, Matrix& rProjectedGradients);
    void CalculateDeformationGradientAndB(const Matrix& rGradients, Matrix& rF, Matrix& rB) const;
    void SetConstitutiveParameters(const Matrix& rF, ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rVariables, bool ComputeTangent) const;
    void CalculateNominalTraction(const Matrix& rF, const Vector& rPK2Stress,
        const array_1d<double, 3>& rNormal, Vector& rTraction) const;
    void CalculateNominalTractionTangent(const Matrix& rGradients, const Matrix& rF,
        const Matrix& rB, const ConstitutiveVariables& rVariables,
        const array_1d<double, 3>& rNormal, Matrix& rTangent) const;
};

} // namespace Kratos
