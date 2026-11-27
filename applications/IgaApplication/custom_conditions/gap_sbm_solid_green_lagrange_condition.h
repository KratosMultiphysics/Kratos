// KRATOS Multiphysics - IgaApplication
// License: BSD License

#pragma once

#include "custom_conditions/gap_sbm_solid_condition.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) GapSbmSolidGreenLagrangeCondition : public GapSbmSolidCondition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GapSbmSolidGreenLagrangeCondition);

    GapSbmSolidGreenLagrangeCondition() = default;
    GapSbmSolidGreenLagrangeCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : GapSbmSolidCondition(NewId, pGeometry) {}
    GapSbmSolidGreenLagrangeCondition(IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : GapSbmSolidCondition(NewId, pGeometry, pProperties) {}

    Condition::Pointer Create(IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) const override;
    Condition::Pointer Create(IndexType NewId, NodesArrayType const& rNodes,
        PropertiesType::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    std::string Info() const override
    {
        return "GapSbmSolidGreenLagrangeCondition #" + std::to_string(Id());
    }

private:
    void CalculateGapKinematics(Vector& rTaylorShapeFunctions, Matrix& rGradients,
        Matrix& rF, Matrix& rB);
    void SetConstitutiveParameters(const Matrix& rF, ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rVariables, bool ComputeTangent) const;
    void CalculateNominalTractionTangent(const Matrix& rGradients, const Matrix& rF,
        const Matrix& rB, const ConstitutiveVariables& rVariables, Matrix& rTangent) const;

    Matrix mInitialTractionTangent;
};

} // namespace Kratos
