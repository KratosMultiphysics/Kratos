//  KRATOS Multiphysics - IgaApplication
//  License: BSD License

#pragma once

#include "custom_conditions/sbm_solid_condition.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) SbmSolidGreenLagrangeCondition : public SbmSolidCondition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SbmSolidGreenLagrangeCondition);

    SbmSolidGreenLagrangeCondition() = default;

    SbmSolidGreenLagrangeCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : SbmSolidCondition(NewId, pGeometry) {}

    SbmSolidGreenLagrangeCondition(
        IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : SbmSolidCondition(NewId, pGeometry, pProperties) {}

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
        return "SbmSolidGreenLagrangeCondition #" + std::to_string(Id());
    }

private:
    void CalculateKinematics(Matrix& rDN_DX, Matrix& rF, Matrix& rB) const;
    void SetConstitutiveParameters(const Matrix& rF, double DeterminantF,
        ConstitutiveLaw::Parameters& rValues, ConstitutiveVariables& rVariables,
        bool ComputeTangent) const;
    void CalculateTractionAndTangent(const Matrix& rDN_DX, const Matrix& rF,
        const Matrix& rB, const ConstitutiveVariables& rVariables,
        Vector& rTraction, Matrix& rTractionTangent) const;
    Matrix mInitialTractionTangent;
};

} // namespace Kratos
