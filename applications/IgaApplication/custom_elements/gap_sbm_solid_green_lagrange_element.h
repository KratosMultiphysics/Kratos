// KRATOS Multiphysics - IgaApplication
// License: BSD License

#pragma once

#include "custom_elements/gap_sbm_solid_element.h"

namespace Kratos
{

class KRATOS_API(IGA_APPLICATION) GapSbmSolidGreenLagrangeElement : public GapSbmSolidElement
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(GapSbmSolidGreenLagrangeElement);

    GapSbmSolidGreenLagrangeElement() = default;
    GapSbmSolidGreenLagrangeElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : GapSbmSolidElement(NewId, pGeometry) {}
    GapSbmSolidGreenLagrangeElement(IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties)
        : GapSbmSolidElement(NewId, pGeometry, pProperties) {}

    Element::Pointer Create(IndexType NewId, NodesArrayType const& rNodes,
        PropertiesType::Pointer pProperties) const override;
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties) const override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;
    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;
    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    std::string Info() const override
    {
        return "GapSbmSolidGreenLagrangeElement #" + std::to_string(Id());
    }

private:
    void CalculateGreenKinematics(Vector& rTaylorShapeFunctions, Matrix& rGradients,
        Matrix& rF, Matrix& rB);
    void SetFiniteStrainParameters(const Matrix& rF, ConstitutiveLaw::Parameters& rValues,
        ConstitutiveVariables& rVariables, bool ComputeTangent) const;
    void CalculateAndAddGeometricStiffness(MatrixType& rLHS, const Matrix& rGradients,
        const Vector& rPK2Stress, double IntegrationWeight) const;
};

} // namespace Kratos
