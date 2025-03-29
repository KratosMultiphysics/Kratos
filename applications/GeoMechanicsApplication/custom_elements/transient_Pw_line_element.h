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
//

#pragma once

#include "calculation_contribution.h"
#include "compressibility_calculator.h"
#include "custom_retention/retention_law_factory.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "filter_compressibility_calculator.h"
#include "fluid_body_flow_calculator.h"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "permeability_calculator.h"
#include "pw_line_integration_coefficients.h"
#include <numeric>
#include <optional>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientPwLineElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientPwLineElement);

    explicit TransientPwLineElement(IndexType NewId = 0) : Element(NewId) {}

    TransientPwLineElement(IndexType                                   NewId,
                           const GeometryType::Pointer&                pGeometry,
                           const std::vector<CalculationContribution>& rContributions)
        : Element(NewId, pGeometry), mContributions(rContributions)
    {
    }

    TransientPwLineElement(IndexType                                   NewId,
                           const GeometryType::Pointer&                pGeometry,
                           const PropertiesType::Pointer&              pProperties,
                           const std::vector<CalculationContribution>& rContributions)
        : Element(NewId, pGeometry, pProperties), mContributions(rContributions)
    {
    }

    Element::Pointer Create(IndexType NewId, const NodesArrayType& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<TransientPwLineElement>(NewId, GetGeometry().Create(rThisNodes),
                                                      pProperties, mContributions);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<TransientPwLineElement>(NewId, pGeom, pProperties, mContributions);
    }

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo&) const override
    {
        rElementalDofList = GetDofs();
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const override
    {
        rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
    }

    void Initialize(const ProcessInfo&) override
    {
        mRetentionLawVector.resize(GetGeometry().IntegrationPointsNumber(GetIntegrationMethod()));

        for (auto& r_retention_law : mRetentionLawVector) {
            r_retention_law = RetentionLawFactory::Clone(GetProperties());
        }
    }

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        rRightHandSideVector = ZeroVector{TNumNodes};
        rLeftHandSideMatrix  = ZeroMatrix{TNumNodes, TNumNodes};

        for (const auto& rContribution : mContributions) {
            const auto calculator = CreateCalculator(rContribution, rCurrentProcessInfo);
            const auto [LHSContribution, RHSContribution] = calculator->LocalSystemContribution();
            if (LHSContribution) rLeftHandSideMatrix += *LHSContribution;
            rRightHandSideVector += RHSContribution;
        }

        KRATOS_CATCH("")
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override
    {
        rRightHandSideVector = ZeroVector{TNumNodes};
        for (const auto& rContribution : mContributions) {
            const auto calculator = CreateCalculator(rContribution, rCurrentProcessInfo);
            noalias(rRightHandSideVector) += calculator->RHSContribution();
        }
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        rLeftHandSideMatrix = ZeroMatrix{TNumNodes, TNumNodes};
        for (const auto& rContribution : mContributions) {
            const auto calculator = CreateCalculator(rContribution, rCurrentProcessInfo);
            if (const auto LHSContribution = calculator->LHSContribution())
                rLeftHandSideMatrix += *LHSContribution;
        }
    }

    GeometryData::IntegrationMethod GetIntegrationMethod() const override
    {
        switch (this->GetGeometry().GetGeometryOrderType()) {
        case GeometryData::Kratos_Cubic_Order:
            return GeometryData::IntegrationMethod::GI_GAUSS_3;
        case GeometryData::Kratos_Quartic_Order:
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        default:
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override
    {
        KRATOS_TRY

        CheckElementLength();
        CheckHasSolutionStepsDataFor(WATER_PRESSURE);
        CheckHasSolutionStepsDataFor(DT_WATER_PRESSURE);
        CheckHasSolutionStepsDataFor(VOLUME_ACCELERATION);
        CheckHasDofsFor(WATER_PRESSURE);
        CheckProperties();
        CheckForNonZeroZCoordinateIn2D();
        CheckRetentionLaw(rCurrentProcessInfo);

        KRATOS_CATCH("")

        return 0;
    }

private:
    std::vector<RetentionLaw::Pointer>                 mRetentionLawVector;
    std::vector<CalculationContribution>               mContributions;
    std::unique_ptr<IntegrationCoefficientsCalculator> mpIntegrationCoefficientsCalculator =
        std::make_unique<PwLineIntegrationCoefficients>();

    void CheckElementLength() const
    {
        constexpr auto min_domain_size = 1.0e-15;
        KRATOS_ERROR_IF(GetGeometry().DomainSize() < min_domain_size)
            << "Length smaller than " << min_domain_size << " for element " << Id() << std::endl;
    }

    void CheckHasSolutionStepsDataFor(const VariableData& rVariable) const
    {
        for (const auto& node : GetGeometry()) {
            KRATOS_ERROR_IF_NOT(node.SolutionStepsDataHas(rVariable))
                << "Missing variable " << rVariable.Name() << " on node " << node.Id() << std::endl;
        }
    }

    void CheckHasDofsFor(const Variable<double>& rVariable) const
    {
        for (const auto& node : GetGeometry()) {
            KRATOS_ERROR_IF_NOT(node.HasDofFor(rVariable))
                << "Missing degree of freedom for " << rVariable.Name() << " on node " << node.Id()
                << std::endl;
        }
    }

    void CheckProperties() const
    {
        CheckProperty(DENSITY_WATER);
        CheckProperty(DENSITY_SOLID);
        constexpr auto max_value = 1.0;
        CheckProperty(POROSITY, max_value);
        CheckProperty(BULK_MODULUS_SOLID);
        CheckProperty(BULK_MODULUS_FLUID);
        CheckProperty(DYNAMIC_VISCOSITY);
        CheckProperty(BIOT_COEFFICIENT);
        CheckProperty(PERMEABILITY_XX);
    }

    void CheckProperty(const Kratos::Variable<double>& rVariable, std::optional<double> MaxValue = std::nullopt) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(rVariable))
            << rVariable.Name()
            << " does not exist in the material properties (Id = " << GetProperties().Id()
            << ") at element " << Id() << std::endl;
        constexpr auto min_value = 0.0;
        if (MaxValue.has_value()) {
            KRATOS_ERROR_IF(GetProperties()[rVariable] < min_value ||
                            GetProperties()[rVariable] > MaxValue.value())
                << rVariable.Name() << " of material Id = " << GetProperties().Id() << " at element "
                << Id() << " has an invalid value " << GetProperties()[rVariable] << " which is outside of the range [ "
                << min_value << ", " << MaxValue.value() << "]" << std::endl;
        } else {
            KRATOS_ERROR_IF(GetProperties()[rVariable] < min_value)
                << rVariable.Name() << " of material Id = " << GetProperties().Id()
                << " at element " << Id() << " has an invalid value " << GetProperties()[rVariable]
                << " which is below the minimum allowed value of " << min_value << std::endl;
        }
    }

    void CheckProperty(const Kratos::Variable<std::string>& rVariable, const std::string& rName) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(rVariable))
            << rVariable.Name() << " does not exist in the pressure element's properties" << std::endl;
        KRATOS_ERROR_IF_NOT(GetProperties()[rVariable] == rName)
            << rVariable.Name() << " has a value of (" << GetProperties()[rVariable]
            << ") instead of (" << rName << ") at element " << Id() << std::endl;
    }

    void CheckForNonZeroZCoordinateIn2D() const
    {
        if constexpr (TDim == 2) {
            const auto& r_geometry = GetGeometry();
            auto        pos        = std::find_if(r_geometry.begin(), r_geometry.end(),
                                                  [](const auto& node) { return node.Z() != 0.0; });
            KRATOS_ERROR_IF_NOT(pos == r_geometry.end())
                << "Node with non-zero Z coordinate found. Id: " << pos->Id() << std::endl;
        }
    }

    void CheckRetentionLaw(const ProcessInfo& rCurrentProcessInfo) const
    {
        if (!mRetentionLawVector.empty()) {
            mRetentionLawVector[0]->Check(this->GetProperties(), rCurrentProcessInfo);
        }
    }

    array_1d<double, TNumNodes> GetNodalValuesOf(const Variable<double>& rNodalVariable) const
    {
        auto        result     = array_1d<double, TNumNodes>{};
        const auto& r_geometry = GetGeometry();
        std::transform(r_geometry.begin(), r_geometry.end(), result.begin(), [&rNodalVariable](const auto& node) {
            return node.FastGetSolutionStepValue(rNodalVariable);
        });
        return result;
    }

    std::vector<Vector> CalculateProjectedGravityAtIntegrationPoints(const Matrix& rNContainer) const
    {
        const auto number_integration_points = GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
        GeometryType::JacobiansType J_container{number_integration_points};
        for (auto& j : J_container) {
            j.resize(GetGeometry().WorkingSpaceDimension(), GetGeometry().LocalSpaceDimension(), false);
        }
        GetGeometry().Jacobian(J_container, this->GetIntegrationMethod());

        array_1d<double, TNumNodes * TDim> volume_acceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(
            volume_acceleration, GetGeometry(), VOLUME_ACCELERATION);
        array_1d<double, TDim> body_acceleration;

        std::vector<Vector> projected_gravity;
        projected_gravity.reserve(number_integration_points);

        for (unsigned int integration_point_index = 0;
             integration_point_index < number_integration_points; ++integration_point_index) {
            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                body_acceleration, rNContainer, volume_acceleration, integration_point_index);
            array_1d<double, TDim> tangent_vector = column(J_container[integration_point_index], 0);
            tangent_vector /= norm_2(tangent_vector);
            projected_gravity.emplace_back(
                ScalarVector(1, std::inner_product(tangent_vector.begin(), tangent_vector.end(),
                                                   body_acceleration.begin(), 0.0)));
        }
        return projected_gravity;
    }

    std::unique_ptr<ContributionCalculator> CreateCalculator(const CalculationContribution& rContribution,
                                                             const ProcessInfo& rCurrentProcessInfo)
    {
        switch (rContribution) {
        case CalculationContribution::Permeability:
            return std::make_unique<PermeabilityCalculator>(CreatePermeabilityInputProvider());
        case CalculationContribution::Compressibility:
            if (GetProperties()[RETENTION_LAW] == "PressureFilterLaw") {
                return std::make_unique<FilterCompressibilityCalculator>(
                    CreateFilterCompressibilityInputProvider(rCurrentProcessInfo));
            }
            return std::make_unique<CompressibilityCalculator>(CreateCompressibilityInputProvider(rCurrentProcessInfo));
        case CalculationContribution::FluidBodyFlow:
            return std::make_unique<FluidBodyFlowCalculator>(CreateFluidBodyFlowInputProvider());
        default:
            KRATOS_ERROR << "Unknown contribution" << std::endl;
        }
    }

    CompressibilityCalculator::InputProvider CreateCompressibilityInputProvider(const ProcessInfo& rCurrentProcessInfo)
    {
        return CompressibilityCalculator::InputProvider(
            MakePropertiesGetter(), MakeRetentionLawsGetter(), MakeNContainerGetter(),
            MakeIntegrationCoefficientsGetter(), MakeMatrixScalarFactorGetter(rCurrentProcessInfo),
            MakeNodalVariableGetter());
    }

    FilterCompressibilityCalculator::InputProvider CreateFilterCompressibilityInputProvider(const ProcessInfo& rCurrentProcessInfo)
    {
        return FilterCompressibilityCalculator::InputProvider(
            MakePropertiesGetter(), MakeNContainerGetter(), MakeIntegrationCoefficientsGetter(),
            MakeProjectedGravityForIntegrationPointsGetter(),
            MakeMatrixScalarFactorGetter(rCurrentProcessInfo), MakeNodalVariableGetter());
    }

    PermeabilityCalculator::InputProvider CreatePermeabilityInputProvider()
    {
        return PermeabilityCalculator::InputProvider(
            MakePropertiesGetter(), MakeRetentionLawsGetter(), MakeIntegrationCoefficientsGetter(),
            MakeNodalVariableGetter(), MakeShapeFunctionLocalGradientsGetter());
    }

    FluidBodyFlowCalculator::InputProvider CreateFluidBodyFlowInputProvider()
    {
        return FluidBodyFlowCalculator::InputProvider(
            MakePropertiesGetter(), MakeRetentionLawsGetter(), MakeIntegrationCoefficientsGetter(),
            MakeProjectedGravityForIntegrationPointsGetter(),
            MakeShapeFunctionLocalGradientsGetter(), MakeLocalSpaceDimensionGetter());
    }

    auto MakePropertiesGetter()
    {
        return [this]() -> const Properties& { return GetProperties(); };
    }

    auto MakeRetentionLawsGetter()
    {
        return [this]() -> const std::vector<RetentionLaw::Pointer>& { return mRetentionLawVector; };
    }

    auto MakeNContainerGetter()
    {
        return [this]() -> const Matrix& {
            return GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());
        };
    }

    auto MakeIntegrationCoefficientsGetter()
    {
        return [this]() -> Vector {
            Vector det_J_container;
            GetGeometry().DeterminantOfJacobian(det_J_container, this->GetIntegrationMethod());
            return mpIntegrationCoefficientsCalculator->CalculateIntegrationCoefficients(
                GetGeometry().IntegrationPoints(GetIntegrationMethod()), det_J_container,
                GetProperties()[CROSS_AREA]);
        };
    }

    auto MakeProjectedGravityForIntegrationPointsGetter()
    {
        return [this]() -> std::vector<Vector> {
            return CalculateProjectedGravityAtIntegrationPoints(
                GetGeometry().ShapeFunctionsValues(GetIntegrationMethod()));
        };
    }

    static auto MakeMatrixScalarFactorGetter(const ProcessInfo& rCurrentProcessInfo)
    {
        return [&rCurrentProcessInfo]() { return rCurrentProcessInfo[DT_PRESSURE_COEFFICIENT]; };
    }

    auto MakeNodalVariableGetter() const
    {
        return
            [this](const Variable<double>& variable) -> Vector { return GetNodalValuesOf(variable); };
    }

    auto MakeShapeFunctionLocalGradientsGetter()
    {
        return [this]() {
            Vector det_J_container;
            GetGeometry().DeterminantOfJacobian(det_J_container, this->GetIntegrationMethod());
            GeometryType::ShapeFunctionsGradientsType dN_dX_container =
                GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
            std::transform(dN_dX_container.begin(), dN_dX_container.end(), det_J_container.begin(),
                           dN_dX_container.begin(), std::divides<>());

            return dN_dX_container;
        };
    }

    auto MakeLocalSpaceDimensionGetter() const
    {
        return [this]() -> std::size_t { return this->GetGeometry().LocalSpaceDimension(); };
    }

    [[nodiscard]] DofsVectorType GetDofs() const
    {
        return Geo::DofUtilities::ExtractDofsFromNodes(GetGeometry(), WATER_PRESSURE);
    }

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
        rSerializer.save("RetentionlawVector", mRetentionLawVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }
};

} // namespace Kratos
