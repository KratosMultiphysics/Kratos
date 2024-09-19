// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#pragma once

#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "integration_scheme.h"
#include "stress_state_policy.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) LineInterfaceElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(LineInterfaceElement);

    using Element::GeometryType;
    using Element::PropertiesType;

    LineInterfaceElement() = default;

    LineInterfaceElement(IndexType                      NewId,
                         const GeometryType::Pointer&   rGeometry,
                         const PropertiesType::Pointer& rProperties);

    LineInterfaceElement(IndexType NewId, const GeometryType::Pointer& rGeometry);
    Element::Pointer Create(IndexType NewId, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const override;
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo&) override;
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo&) override;
    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;
    void CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                      std::vector<ConstitutiveLaw::Pointer>&    rOutput,
                                      const ProcessInfo&) override;
    using Element::CalculateOnIntegrationPoints;
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override;
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    int  Check(const ProcessInfo& rCurrentProcessInfo) const override;

private:
    Element::DofsVectorType GetDofs() const;

    std::vector<Matrix> CalculateLocalBMatricesAtIntegrationPoints() const;
    std::vector<double> CalculateIntegrationCoefficients() const;
    std::vector<Matrix> CalculateConstitutiveMatricesAtIntegrationPoints();
    std::vector<Vector> CalculateRelativeDisplacementsAtIntegrationPoints(const std::vector<Matrix>& rLocalBMatrices) const;
    std::vector<Vector> CalculateTractionsAtIntegrationPoints(const std::vector<Vector>& rRelativeDisplacements);

    std::unique_ptr<IntegrationScheme>    mIntegrationScheme;
    std::unique_ptr<StressStatePolicy>    mStressStatePolicy;
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLaws;
};

} // namespace Kratos
