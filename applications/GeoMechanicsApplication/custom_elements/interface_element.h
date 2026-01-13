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
//                   Gennady Markelov
//

#pragma once

#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "integration_coefficients_calculator.h"
#include "integration_scheme.h"
#include "stress_state_policy.h"
#include <optional>

namespace Kratos
{
class KRATOS_API(GEO_MECHANICS_APPLICATION) InterfaceElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(InterfaceElement);

    using Element::GeometryType;
    using Element::PropertiesType;

    // Following the UPwBaseElement example, we will follow the rule of 5
    // to avoid the noexcept code smell.
    ~InterfaceElement() override                             = default;
    InterfaceElement(const InterfaceElement&)                = delete;
    InterfaceElement& operator=(const InterfaceElement&)     = delete;
    InterfaceElement(InterfaceElement&&) noexcept            = default;
    InterfaceElement& operator=(InterfaceElement&&) noexcept = default;

    InterfaceElement(IndexType                          NewId,
                     const GeometryType::Pointer&       rGeometry,
                     const PropertiesType::Pointer&     rProperties,
                     std::unique_ptr<StressStatePolicy> pStressStatePolicy);

    InterfaceElement(IndexType                          NewId,
                     const GeometryType::Pointer&       rGeometry,
                     std::unique_ptr<StressStatePolicy> pStressStatePolicy);
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

    void Calculate(const Variable<Vector>& rVariable, Vector& rOutput, const ProcessInfo& rProcessInfo) override;
    using Element::Calculate;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override;
    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;
    int  Check(const ProcessInfo& rCurrentProcessInfo) const override;

    const IntegrationScheme& GetIntegrationScheme() const;

    const Geometry<Node>& GetMidGeometry() const;

private:
    InterfaceElement() = default;

    Element::DofsVectorType GetDofs() const;

    std::vector<Matrix> CalculateLocalBMatricesAtIntegrationPoints() const;
    std::vector<double> CalculateIntegrationCoefficients() const;
    std::vector<Matrix> CalculateConstitutiveMatricesAtIntegrationPoints(const std::vector<Vector>& rRelativeDisplacements,
                                                                         const ProcessInfo& rProcessInfo);
    std::vector<Vector> CalculateRelativeDisplacementsAtIntegrationPoints(const std::vector<Matrix>& rLocalBMatrices) const;
    std::vector<Vector> CalculateTractionsAtIntegrationPoints(const std::vector<Vector>& rRelativeDisplacements,
                                                              const ProcessInfo& rProcessInfo);
    void ApplyRotationToBMatrix(Matrix& rBMatrix, const Matrix& rRotationMatrix) const;
    void MakeIntegrationSchemeAndAssignFunction();
    void InterpolateNodalStressesToInitialTractions(const std::vector<std::optional<Vector>>& rInterfaceNodalCauchyStresses) const;
    Vector InterpolateNodalStressToIntegrationPoints(const Geo::IntegrationPointType& rIntegrationPoint,
                                                     const std::vector<Vector>& rNodalStresses) const;
    Matrix RotateStressToLocalCoordinates(const Geo::IntegrationPointType& rIntegrationPoint,
                                          const Vector& rGlobalStressVector) const;
    Vector ConvertLocalStressToTraction(const Matrix& rLocalStress) const;
    std::function<Matrix(const Geometry<Node>&, const array_1d<double, 3>&)> mfpCalculateRotationMatrix;

    std::unique_ptr<IntegrationScheme>    mpIntegrationScheme;
    std::unique_ptr<StressStatePolicy>    mpStressStatePolicy;
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLaws;
    IntegrationCoefficientsCalculator     mIntegrationCoefficientsCalculator;

    friend class Serializer;
};

} // namespace Kratos
