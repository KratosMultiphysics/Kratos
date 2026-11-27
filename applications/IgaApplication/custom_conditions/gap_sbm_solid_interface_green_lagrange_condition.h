// KRATOS Multiphysics - IgaApplication
// License: BSD License

#pragma once

#include "custom_conditions/gap_sbm_solid_interface_condition.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) GapSbmSolidInterfaceGreenLagrangeCondition
    : public GapSbmSolidInterfaceCondition
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GapSbmSolidInterfaceGreenLagrangeCondition);

    GapSbmSolidInterfaceGreenLagrangeCondition() = default;
    GapSbmSolidInterfaceGreenLagrangeCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : GapSbmSolidInterfaceCondition(NewId, pGeometry) {}
    GapSbmSolidInterfaceGreenLagrangeCondition(IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : GapSbmSolidInterfaceCondition(NewId, pGeometry, pProperties) {}

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
        return "GapSbmSolidInterfaceGreenLagrangeCondition #" + std::to_string(Id());
    }

private:
    void CalculateSideKinematics(const GeometryType& rGeometry, const Vector& rDistance,
        const Vector& rDisplacement, Matrix& rGradients, Matrix& rF, Matrix& rB);
    void SetFiniteStrainParameters(const Matrix& rF, ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rVariables, bool ComputeTangent) const;
    void CalculateNominalTractionAndTangent(const Matrix& rGradients, const Matrix& rF,
        const Matrix& rB, const ConstitutiveVariables& rVariables,
        Vector& rTraction, Matrix& rTangent) const;
    void UpdateFiniteStrainResponses(const ProcessInfo& rCurrentProcessInfo, bool Finalize);

    Matrix mInitialTractionTangentPlus;
    Matrix mInitialTractionTangentMinus;
};

} // namespace Kratos
