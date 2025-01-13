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
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "permeability_calculator.h"

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

        rRightHandSideVector = CalculateFluidBodyVector();
        rLeftHandSideMatrix  = ZeroMatrix{TNumNodes, TNumNodes};

        for (const auto& rContribution : mContributions) {
            const auto calculator = CreateCalculator(rContribution, rCurrentProcessInfo);
            const auto [LHSContribution, RHSContribution] = calculator->LocalSystemContribution();
            rLeftHandSideMatrix += LHSContribution;
            rRightHandSideVector += RHSContribution;
        }

        KRATOS_CATCH("")
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override
    {
        rRightHandSideVector = CalculateFluidBodyVector();
        for (const auto& rContribution : mContributions) {
            const auto calculator = CreateCalculator(rContribution, rCurrentProcessInfo);
            rRightHandSideVector += calculator->RHSContribution();
        }
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        rLeftHandSideMatrix = ZeroMatrix{TNumNodes, TNumNodes};
        for (const auto& rContribution : mContributions) {
            const auto calculator = CreateCalculator(rContribution, rCurrentProcessInfo);
            rLeftHandSideMatrix += calculator->LHSContribution();
        }
    }

    GeometryData::IntegrationMethod GetIntegrationMethod() const override
    {
        switch (TNumNodes) {
        case 2:
        case 3:
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        case 4:
            return GeometryData::IntegrationMethod::GI_GAUSS_3;
        case 5:
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        default:
            KRATOS_ERROR << "Can't return integration method: unexpected number of nodes: " << TNumNodes
                         << std::endl;
        }
    }

    int Check(const ProcessInfo&) const override
    {
        KRATOS_TRY

        CheckDomainSize();
        CheckHasSolutionStepsDataFor(WATER_PRESSURE);
        CheckHasSolutionStepsDataFor(DT_WATER_PRESSURE);
        CheckHasDofsFor(WATER_PRESSURE);
        CheckProperties();
        CheckForNonZeroZCoordinateIn2D();

        KRATOS_CATCH("")

        return 0;
    }

private:
    std::vector<RetentionLaw::Pointer>   mRetentionLawVector;
    std::vector<CalculationContribution> mContributions;

    void CheckDomainSize() const
    {
        constexpr auto min_domain_size = 1.0e-15;
        KRATOS_ERROR_IF(GetGeometry().DomainSize() < min_domain_size)
            << "DomainSize smaller than " << min_domain_size << " for element " << Id() << std::endl;
    }

    void CheckHasSolutionStepsDataFor(const Variable<double>& rVariable) const
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
        CheckProperty(POROSITY);
        CheckProperty(BULK_MODULUS_SOLID);
        CheckProperty(BULK_MODULUS_FLUID);
        CheckProperty(DYNAMIC_VISCOSITY);
        CheckProperty(BIOT_COEFFICIENT);
        CheckProperty(PERMEABILITY_XX);
    }

    void CheckProperty(const Kratos::Variable<double>& rVariable) const
    {
        KRATOS_ERROR_IF_NOT(GetProperties().Has(rVariable))
            << rVariable.Name() << " does not exist in the pressure element's properties" << std::endl;
        KRATOS_ERROR_IF(GetProperties()[rVariable] < 0.0)
            << rVariable.Name() << " has an invalid value at element " << Id() << std::endl;
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
                << " Node with non-zero Z coordinate found. Id: " << pos->Id() << std::endl;
        }
    }

    Vector CalculateIntegrationCoefficients(const Vector& rDetJContainer) const
    {
        const auto& r_properties         = GetProperties();
        const auto& r_integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());

        auto result = Vector{r_integration_points.size()};
        std::transform(r_integration_points.begin(), r_integration_points.end(), rDetJContainer.begin(),
                       result.begin(), [&r_properties](const auto& rIntegrationPoint, const auto& rDetJ) {
            return rIntegrationPoint.Weight() * rDetJ * r_properties[CROSS_AREA];
        });
        return result;
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

    array_1d<double, TNumNodes> CalculateFluidBodyVector() const
    {
        Vector det_J_container;
        GetGeometry().DeterminantOfJacobian(det_J_container, this->GetIntegrationMethod());
        GeometryType::ShapeFunctionsGradientsType shape_function_gradients =
            GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
        std::transform(shape_function_gradients.begin(), shape_function_gradients.end(),
                       det_J_container.begin(), shape_function_gradients.begin(), std::divides<>());
        const Matrix& N_container = GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());

        const auto integration_coefficients = CalculateIntegrationCoefficients(det_J_container);

        const auto&                 r_properties = GetProperties();
        BoundedMatrix<double, 1, 1> constitutive_matrix;
        GeoElementUtilities::FillPermeabilityMatrix(constitutive_matrix, r_properties);

        RetentionLaw::Parameters RetentionParameters(GetProperties());

        const auto projected_gravity = CalculateProjectedGravityAtIntegrationPoints(N_container);

        array_1d<double, TNumNodes> fluid_body_vector = ZeroVector(TNumNodes);
        for (unsigned int integration_point_index = 0;
             integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
             ++integration_point_index) {
            const auto N = Vector{row(N_container, integration_point_index)};
            double     RelativePermeability =
                mRetentionLawVector[integration_point_index]->CalculateRelativePermeability(RetentionParameters);
            fluid_body_vector +=
                r_properties[DENSITY_WATER] * RelativePermeability *
                prod(prod(shape_function_gradients[integration_point_index], constitutive_matrix),
                     ScalarVector(1, projected_gravity[integration_point_index])) *
                integration_coefficients[integration_point_index] / r_properties[DYNAMIC_VISCOSITY];
        }
        return fluid_body_vector;
    }

    Vector CalculateProjectedGravityAtIntegrationPoints(const Matrix& rNContainer) const
    {
        const auto number_integration_points = GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
        GeometryType::JacobiansType J_container;
        J_container.resize(number_integration_points, false);
        for (std::size_t i = 0; i < number_integration_points; ++i) {
            J_container[i].resize(GetGeometry().WorkingSpaceDimension(),
                                  GetGeometry().LocalSpaceDimension(), false);
        }
        GetGeometry().Jacobian(J_container, this->GetIntegrationMethod());

        array_1d<double, TNumNodes * TDim> volume_acceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(
            volume_acceleration, GetGeometry(), VOLUME_ACCELERATION);
        array_1d<double, TDim> body_acceleration;

        Vector projected_gravity = ZeroVector(number_integration_points);

        for (unsigned int integration_point_index = 0;
             integration_point_index < number_integration_points; ++integration_point_index) {
            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                body_acceleration, rNContainer, volume_acceleration, integration_point_index);
            array_1d<double, TDim> tangent_vector = column(J_container[integration_point_index], 0);
            tangent_vector /= norm_2(tangent_vector);
            projected_gravity(integration_point_index) =
                MathUtils<double>::Dot(tangent_vector, body_acceleration);
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
            return CalculateIntegrationCoefficients(det_J_container);
        };
    }

    auto MakeProjectedGravityForIntegrationPointsGetter()
    {
        return [this]() -> Vector {
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
