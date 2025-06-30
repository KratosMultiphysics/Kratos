// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

// Project includes
#include "containers/array_1d.h"
#include "custom_retention/retention_law.h"
#include "custom_retention/retention_law_factory.h"
#include "geometries/geometry.h"
#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_utilities/element_utilities.hpp"
#include "integration_coefficients_calculator.h"
#include "stress_state_policy.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) UPwBaseElement : public Element
{
public:
    using SizeType = std::size_t;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(UPwBaseElement);

    using Element::Element;

    /// Constructor using an array of nodes
    UPwBaseElement(IndexType                                       NewId,
                   const NodesArrayType&                           ThisNodes,
                   std::unique_ptr<StressStatePolicy>              pStressStatePolicy,
                   std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : Element(NewId, ThisNodes),
          mpStressStatePolicy{std::move(pStressStatePolicy)},
          mIntegrationCoefficientsCalculator{std::move(pCoefficientModifier)}
    {
    }

    /// Constructor using Geometry
    UPwBaseElement(IndexType                                       NewId,
                   GeometryType::Pointer                           pGeometry,
                   std::unique_ptr<StressStatePolicy>              pStressStatePolicy,
                   std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : Element(NewId, pGeometry),
          mpStressStatePolicy{std::move(pStressStatePolicy)},
          mIntegrationCoefficientsCalculator{std::move(pCoefficientModifier)}
    {
    }

    /// Constructor using Properties
    UPwBaseElement(IndexType                                       NewId,
                   GeometryType::Pointer                           pGeometry,
                   PropertiesType::Pointer                         pProperties,
                   std::unique_ptr<StressStatePolicy>              pStressStatePolicy,
                   std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier = nullptr)
        : Element(NewId, pGeometry, pProperties),
          mpStressStatePolicy{std::move(pStressStatePolicy)},
          mIntegrationCoefficientsCalculator{std::move(pCoefficientModifier)}
    {
        // this is needed for interface elements
        mThisIntegrationMethod = this->GetIntegrationMethod();
    }

    ~UPwBaseElement() override                           = default;
    UPwBaseElement(const UPwBaseElement&)                = delete;
    UPwBaseElement& operator=(const UPwBaseElement&)     = delete;
    UPwBaseElement(UPwBaseElement&&) noexcept            = delete;
    UPwBaseElement& operator=(UPwBaseElement&&) noexcept = delete;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override;

    void ResetConstitutiveLaw() override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    void GetFirstDerivativesVector(Vector& rValues, int Step = 0) const override;

    void GetSecondDerivativesVector(Vector& rValues, int Step = 0) const override;

    void CalculateOnIntegrationPoints(const Variable<ConstitutiveLaw::Pointer>& rVariable,
                                      std::vector<ConstitutiveLaw::Pointer>&    rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rValues,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rValues,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rValues,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rValues,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    using Element::CalculateOnIntegrationPoints;

    void SetValuesOnIntegrationPoints(const Variable<double>&    rVariable,
                                      const std::vector<double>& rValues,
                                      const ProcessInfo&         rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Vector>&    rVariable,
                                      const std::vector<Vector>& rValues,
                                      const ProcessInfo&         rCurrentProcessInfo) override;

    void SetValuesOnIntegrationPoints(const Variable<Matrix>&    rVariable,
                                      const std::vector<Matrix>& rValues,
                                      const ProcessInfo&         rCurrentProcessInfo) override;

    using Element::SetValuesOnIntegrationPoints;

    std::string Info() const override
    {
        const std::string constitutive_info =
            !mConstitutiveLawVector.empty() ? mConstitutiveLawVector[0]->Info() : "not defined";
        return "U-Pw Base class Element #" + std::to_string(Id()) + "\nConstitutive law: " + constitutive_info;
    }

    void PrintInfo(std::ostream& rOStream) const override { rOStream << Info(); }

protected:
    /// Member Variables
    GeometryData::IntegrationMethod mThisIntegrationMethod;

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    std::vector<RetentionLaw::Pointer>    mRetentionLawVector;

    std::vector<Vector> mStressVector;
    std::vector<Vector> mStateVariablesFinalized;

    virtual void CalculateMaterialStiffnessMatrix(MatrixType& rStiffnessMatrix, const ProcessInfo& CurrentProcessInfo);

    virtual void CalculateAll(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& CurrentProcessInfo,
                              bool               CalculateStiffnessMatrixFlag,
                              bool               CalculateResidualVectorFlag);

    std::vector<double> CalculateIntegrationCoefficients(const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
                                                         const Vector& rDetJs) const;

    void CalculateDerivativesOnInitialConfiguration(double&      rDetJ,
                                                    Matrix&      rJ0,
                                                    Matrix&      rInvJ0,
                                                    Matrix&      rDNu_DX0,
                                                    unsigned int IntegrationPointIndex) const;

    void CalculateJacobianOnCurrentConfiguration(double& detJ, Matrix& rJ, Matrix& rInvJ, unsigned int GPoint) const;

    virtual std::size_t GetNumberOfDOF() const;

    StressStatePolicy& GetStressStatePolicy() const;

    std::unique_ptr<IntegrationCoefficientModifier> CloneIntegrationCoefficientModifier() const;

private:
    [[nodiscard]] virtual DofsVectorType GetDofs() const;

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    }

    void load(Serializer& rSerializer) override{KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)}

    std::unique_ptr<StressStatePolicy> mpStressStatePolicy;
    IntegrationCoefficientsCalculator mIntegrationCoefficientsCalculator;
};

// Class UPwBaseElement

} // namespace Kratos
