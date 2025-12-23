// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//                   Richard Faasse
//                   Gennady Markelov
//

#pragma once

#include "calculation_contribution.h"
#include "compressibility_calculator.hpp"
#include "custom_retention/retention_law_factory.h"
#include "filter_compressibility_calculator.hpp"
#include "fluid_body_flow_calculator.hpp"
#include "includes/element.h"
#include "includes/serializer.h"
#include "integration_coefficients_calculator.hpp"
#include "permeability_calculator.hpp"

#include <optional>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) PwElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(PwElement);

    explicit PwElement(IndexType NewId = 0);

    PwElement(IndexType                                       NewId,
              const GeometryType::Pointer&                    pGeometry,
              const std::vector<CalculationContribution>&     rContributions,
              std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier);

    PwElement(IndexType                                       NewId,
              const GeometryType::Pointer&                    pGeometry,
              const PropertiesType::Pointer&                  pProperties,
              const std::vector<CalculationContribution>&     rContributions,
              std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier);

    ~PwElement() override                      = default;
    PwElement(const PwElement&)                = delete;
    PwElement& operator=(const PwElement&)     = delete;
    PwElement(PwElement&&) noexcept            = delete;
    PwElement& operator=(PwElement&&) noexcept = delete;

    Element::Pointer Create(IndexType               NewId,
                            const NodesArrayType&   rThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override;

    void Initialize(const ProcessInfo&) override;

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override;

    using Element::CalculateOnIntegrationPoints;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    GeometryData::IntegrationMethod GetIntegrationMethod() const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override;

private:
    std::vector<CalculationContribution> mContributions;
    IntegrationCoefficientsCalculator    mIntegrationCoefficientsCalculator;
    std::vector<RetentionLaw::Pointer>   mRetentionLawVector;
    Vector                               mIntegrationCoefficients;
    Matrix                               mNContainer;
    Vector                               mDetJCcontainer;
    std::vector<double>                  mFluidPressures;

    std::vector<double> CalculateIntegrationCoefficients(const Vector& rDetJs) const;

    std::vector<Vector> CalculateProjectedGravityAtIntegrationPoints(const Matrix& rNContainer) const;

    std::unique_ptr<IntegrationCoefficientModifier> CloneIntegrationCoefficientModifier() const;

    std::unique_ptr<ContributionCalculator<TNumNodes>> CreateCalculator(const CalculationContribution& rContribution,
                                                                        const ProcessInfo& rCurrentProcessInfo);

    void CachingDataForCalculator();

    typename CompressibilityCalculator<TNumNodes>::InputProvider CreateCompressibilityInputProvider(
        const ProcessInfo& rCurrentProcessInfo);

    typename FilterCompressibilityCalculator<TNumNodes>::InputProvider CreateFilterCompressibilityInputProvider(
        const ProcessInfo& rCurrentProcessInfo);

    typename PermeabilityCalculator<TNumNodes>::InputProvider CreatePermeabilityInputProvider();

    typename FluidBodyFlowCalculator<TNumNodes>::InputProvider CreateFluidBodyFlowInputProvider();

    auto MakePropertiesGetter()
    {
        return [this]() -> const Properties& { return GetProperties(); };
    }

    auto MakeRetentionLawsGetter()
    {
        return [this]() -> const std::vector<RetentionLaw::Pointer>& { return mRetentionLawVector; };
    }

    auto GetNContainer()
    {
        return [this]() -> const Matrix& { return mNContainer; };
    }

    Matrix CalculateNContainer();

    auto GetIntegrationCoefficients()
    {
        return [this]() -> const Vector& { return mIntegrationCoefficients; };
    }

    Vector CalculateIntegrationCoefficients();

    auto GetFluidPressures()
    {
        return [this]() -> const std::vector<double>& { return mFluidPressures; };
    }

    std::vector<double> CalculateFluidPressure();

    auto MakeProjectedGravityForIntegrationPointsGetter() const
    {
        return [this]() -> std::vector<Vector> {
            return CalculateProjectedGravityAtIntegrationPoints(mNContainer);
        };
    }

    static auto MakeMatrixScalarFactorGetter(const ProcessInfo& rCurrentProcessInfo)
    {
        return [&rCurrentProcessInfo]() { return rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT]; };
    }

    auto MakeNodalVariableGetter() const
    {
        return [this](const Variable<double>& rVariable) -> Vector {
            return VariablesUtilities::GetNodalValuesOf<TNumNodes>(rVariable, this->GetGeometry());
        };
    }

    auto MakeShapeFunctionLocalGradientsGetter()
    {
        return [this]() {
            GeometryType::ShapeFunctionsGradientsType dN_dX_container;
            if (GetGeometry().LocalSpaceDimension() == 1) {
                dN_dX_container = GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
                std::transform(dN_dX_container.begin(), dN_dX_container.end(),
                               mDetJCcontainer.begin(), dN_dX_container.begin(), std::divides<>());
            } else {
                GetGeometry().ShapeFunctionsIntegrationPointsGradients(
                    dN_dX_container, mDetJCcontainer, this->GetIntegrationMethod());
            }

            return dN_dX_container;
        };
    }

    auto MakeLocalSpaceDimensionGetter() const
    {
        return [this]() -> std::size_t { return this->GetGeometry().LocalSpaceDimension(); };
    }

    [[nodiscard]] DofsVectorType GetDofs() const;

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;
};

} // namespace Kratos
