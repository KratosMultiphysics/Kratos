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

#include "contribution_calculators/calculation_contribution.h"
#include "contribution_calculators/stiffness_calculator.hpp"
#include "custom_elements/contribution_calculators/up_coupling_calculator.hpp"
#include "geo_aliases.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "integration_coefficients_calculator.hpp"
#include "integration_scheme.h"
#include "stress_state_policy.h"

#include <optional>

#include "custom_utilities/element_utilities.hpp"

namespace Kratos
{
class RetentionLaw;

class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwInterfaceElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwInterfaceElement);

    using Element::GeometryType;
    using Element::PropertiesType;

    // Following the UPwBaseElement example, we will follow the rule of 5
    // to avoid the noexcept code smell.
    ~UPwInterfaceElement() override                                = default;
    UPwInterfaceElement(const UPwInterfaceElement&)                = delete;
    UPwInterfaceElement& operator=(const UPwInterfaceElement&)     = delete;
    UPwInterfaceElement(UPwInterfaceElement&&) noexcept            = default;
    UPwInterfaceElement& operator=(UPwInterfaceElement&&) noexcept = default;

    UPwInterfaceElement(IndexType                                   NewId,
                        const GeometryType::Pointer&                rpGeometry,
                        const PropertiesType::Pointer&              rpProperties,
                        std::unique_ptr<StressStatePolicy>          pStressStatePolicy,
                        IsDiffOrderElement                          IsDiffOrder,
                        const std::vector<CalculationContribution>& rContributions);

    UPwInterfaceElement(IndexType                                   NewId,
                        const GeometryType::Pointer&                rpGeometry,
                        std::unique_ptr<StressStatePolicy>          pStressStatePolicy,
                        IsDiffOrderElement                          IsDiffOrder,
                        const std::vector<CalculationContribution>& rContributions);
    Element::Pointer Create(IndexType NewId, const NodesArrayType& rNodes, PropertiesType::Pointer pProperties) const override;
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;
    void CalculateAndAssignStifnessMatrix(Element::MatrixType& rLeftHandSideMatrix, const ProcessInfo& rProcessInfo);
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo&) override;
    void CalculateAndAssignStifnessForceVector(Element::VectorType& rRightHandSideVector,
                                               const ProcessInfo&   rProcessInfo);
    void CalculateAndAssignCouplingMatrix(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rProcessInfo) const;
    void CalculateAndAssembleCouplingForceVector(Element::VectorType& rRightHandSideVector,
                                                 const ProcessInfo&   rProcessInfo) const;
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

    const Geometry<Node>& GetDisplacementMidGeometry() const;
    const Geometry<Node>& GetWaterPressureMidGeometry() const;

private:
    UPwInterfaceElement() = default;

    Element::DofsVectorType GetDofs() const;

    const GeometryType&       GetDisplacementGeometry() const;
    const GeometryType&       GetWaterPressureGeometry() const;
    [[nodiscard]] std::size_t NumberOfUDofs() const;

    std::vector<Matrix> CalculateLocalBMatricesAtIntegrationPoints() const;
    std::vector<double> CalculateIntegrationCoefficients() const;
    std::vector<Vector> CalculateRelativeDisplacementsAtIntegrationPoints(const std::vector<Matrix>& rLocalBMatrices) const;
    void ApplyRotationToBMatrix(Matrix& rBMatrix, const Matrix& rRotationMatrix) const;
    void MakeIntegrationSchemeAndAssignFunction();
    void InterpolateNodalStressesToInitialTractions(const std::vector<std::optional<Vector>>& rInterfaceNodalCauchyStresses) const;
    Vector InterpolateNodalStressToIntegrationPoints(const Geo::IntegrationPointType& rIntegrationPoint,
                                                     const std::vector<Vector>& rNodalStresses) const;
    Matrix RotateStressToLocalCoordinates(const Geo::IntegrationPointType& rIntegrationPoint,
                                          const Vector& rGlobalStressVector) const;
    Vector ConvertLocalStressToTraction(const Matrix& rLocalStress) const;

    Geo::BMatricesGetter                              CreateBMatricesGetter() const;
    Geo::StrainVectorsGetter                          CreateRelativeDisplacementsGetter() const;
    Geo::IntegrationCoefficientsGetter                CreateIntegrationCoefficientsGetter() const;
    Geo::PropertiesGetter                             CreatePropertiesGetter() const;
    Geo::ConstitutiveLawsGetter                       CreateConstitutiveLawsGetter() const;
    std::function<const Matrix()>                     CreateNpContainerGetter() const;
    std::function<Vector()>                           CreateVoigtVectorGetter() const;
    std::function<std::vector<double>()>              CreateBiotCoefficientsGetter() const;
    std::function<std::vector<double>(const Vector&)> CreateBishopCoefficientsGetter() const;
    std::function<Vector()> CreateIntegrationPointFluidPressuresGetter() const;

    std::vector<double> CalculateBiotCoefficients() const;
    std::vector<double> CalculateBishopCoefficients(const Vector& rFluidPressure) const;
    Vector              CalculateIntegrationPointFluidPressures() const;
    Vector GetWaterPressureGeometryNodalVariable(const Variable<double>& rVariable) const;
    Matrix GetNpContainer() const;

    template <unsigned int MatrixSize>
    typename StiffnessCalculator<MatrixSize>::InputProvider CreateStiffnessInputProvider(const ProcessInfo& rProcessInfo);

    template <unsigned int MatrixSize>
    auto CreateStiffnessCalculator(const ProcessInfo& rProcessInfo);

    template <unsigned int MatrixSize>
    void CalculateAndAssignStiffnessMatrix(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rProcessInfo);

    template <unsigned int MatrixSize>
    void CalculateAndAssignStiffnesForceVector(VectorType& rRightHandSideVector, const ProcessInfo& rProcessInfo);

    template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
    typename UPCouplingCalculator<NumberOfRows, NumberOfColumns>::InputProvider CreateCouplingInputProvider(
        const ProcessInfo& rProcessInfo) const;

    template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
    auto CreateCouplingCalculator(const ProcessInfo& rProcessInfo) const;

    template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
    void CalculateAndAssignCouplingMatrix(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rProcessInfo) const;

    template <unsigned int NumberOfRows, unsigned int NumberOfColumns>
    void CalculateAndAssembleCouplingForceVector(VectorType&        rRightHandSideVector,
                                                 const ProcessInfo& rProcessInfo) const;

    std::function<Matrix(const Geometry<Node>&, const array_1d<double, 3>&)> mfpCalculateRotationMatrix;

    std::unique_ptr<IntegrationScheme>    mpIntegrationScheme;
    std::unique_ptr<StressStatePolicy>    mpStressStatePolicy;
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLaws;
    IntegrationCoefficientsCalculator     mIntegrationCoefficientsCalculator;
    Geo::GeometryUniquePtr                mpOptionalPressureGeometry;
    std::vector<CalculationContribution>  mContributions;
    std::vector<RetentionLaw::Pointer>    mRetentionLawVector;

    friend class Serializer;
};
} // namespace Kratos
