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
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/hydraulic_discharge.h"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "custom_utilities/variables_utilities.hpp"
#include "filter_compressibility_calculator.hpp"
#include "fluid_body_flow_calculator.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "integration_coefficients_calculator.h"
#include "permeability_calculator.hpp"

#include <optional>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) PwElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(PwElement);

    explicit PwElement(IndexType NewId = 0) : Element(NewId) {}

    PwElement(IndexType                                       NewId,
              const GeometryType::Pointer&                    pGeometry,
              const std::vector<CalculationContribution>&     rContributions,
              std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
        : Element(NewId, pGeometry),
          mContributions(rContributions),
          mIntegrationCoefficientsCalculator{std::move(pCoefficientModifier)}
    {
    }

    PwElement(IndexType                                       NewId,
              const GeometryType::Pointer&                    pGeometry,
              const PropertiesType::Pointer&                  pProperties,
              const std::vector<CalculationContribution>&     rContributions,
              std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
        : Element(NewId, pGeometry, pProperties),
          mContributions(rContributions),
          mIntegrationCoefficientsCalculator{std::move(pCoefficientModifier)}
    {
    }

    ~PwElement() override                      = default;
    PwElement(const PwElement&)                = delete;
    PwElement& operator=(const PwElement&)     = delete;
    PwElement(PwElement&&) noexcept            = delete;
    PwElement& operator=(PwElement&&) noexcept = delete;

    Element::Pointer Create(IndexType NewId, const NodesArrayType& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<PwElement>(NewId, GetGeometry().Create(rThisNodes), pProperties,
                                         mContributions, this->CloneIntegrationCoefficientModifier());
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<PwElement>(NewId, pGeom, pProperties, mContributions,
                                         this->CloneIntegrationCoefficientModifier());
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

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        if (this->GetGeometry().LocalSpaceDimension() != 1) {
            for (auto& r_node : this->GetGeometry()) {
                r_node.FastGetSolutionStepValue(HYDRAULIC_DISCHARGE) = 0.0;
            }
        }
        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        if (this->GetGeometry().LocalSpaceDimension() != 1) {
            GeometryType& r_geometry = this->GetGeometry();
            const auto    number_of_integration_points =
                r_geometry.IntegrationPointsNumber(this->GetIntegrationMethod());
            Vector                                    det_J_container(number_of_integration_points);
            GeometryType::ShapeFunctionsGradientsType dN_dx_container;
            r_geometry.ShapeFunctionsIntegrationPointsGradients(dN_dx_container, det_J_container,
                                                                GetIntegrationMethod());
            const auto integration_coefficients = this->CalculateIntegrationCoefficients(det_J_container);
            std::vector<array_1d<double, 3>> fluid_flux;
            this->CalculateOnIntegrationPoints(FLUID_FLUX_VECTOR, fluid_flux, rCurrentProcessInfo);
            HydraulicDischarge::CalculateHydraulicDischarge(fluid_flux, integration_coefficients, dN_dx_container,
                                                            this->GetIntegrationMethod(), r_geometry);
        }

        KRATOS_CATCH("")
    }

    using Element::CalculateOnIntegrationPoints;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                                      std::vector<double>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const GeometryType& r_geometry = this->GetGeometry();
        const auto          number_of_integration_points =
            r_geometry.IntegrationPointsNumber(this->GetIntegrationMethod());

        auto& r_properties = this->GetProperties();
        rOutput.resize(number_of_integration_points);

        if (rVariable == DEGREE_OF_SATURATION || rVariable == EFFECTIVE_SATURATION || rVariable == BISHOP_COEFFICIENT ||
            rVariable == DERIVATIVE_OF_SATURATION || rVariable == RELATIVE_PERMEABILITY) {
            Matrix N_container(number_of_integration_points, TNumNodes);
            N_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
            RetentionLaw::Parameters    RetentionParameters(r_properties);
            Vector                      Np(TNumNodes);
            array_1d<double, TNumNodes> pressure_vector;
            VariablesUtilities::GetNodalValues(r_geometry, WATER_PRESSURE, pressure_vector.begin());

            for (unsigned int integration_point = 0;
                 integration_point < number_of_integration_points; ++integration_point) {
                noalias(Np) = row(N_container, integration_point);

                RetentionParameters.SetFluidPressure(
                    GeoTransportEquationUtilities::CalculateFluidPressure(Np, pressure_vector));
                rOutput[integration_point] = mRetentionLawVector[integration_point]->CalculateValue(
                    RetentionParameters, rVariable, rOutput[integration_point]);
            }
        } else if (rVariable == HYDRAULIC_HEAD) {
            // Defining the shape functions, the Jacobian and the shape functions local gradients containers
            const Matrix& N_container = r_geometry.ShapeFunctionsValues(GetIntegrationMethod());

            const auto nodal_hydraulic_head =
                GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(r_geometry, r_properties);

            for (unsigned int integration_point = 0;
                 integration_point < number_of_integration_points; ++integration_point) {
                const auto& shape_function_values = row(N_container, integration_point);
                rOutput[integration_point] =
                    std::inner_product(shape_function_values.begin(), shape_function_values.end(),
                                       nodal_hydraulic_head.begin(), 0.0);
            }
        } else if (r_properties.Has(rVariable)) {
            // Map initial material property to gauss points, as required for the output
            std::fill_n(rOutput.begin(), number_of_integration_points, r_properties.GetValue(rVariable));
        } else {
            std::fill(rOutput.begin(), rOutput.end(), 0.0);
        }

        KRATOS_CATCH("")
    }

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                      std::vector<array_1d<double, 3>>&    rOutput,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const GeometryType& r_geom = this->GetGeometry();
        const IndexType     number_of_integration_points =
            r_geom.IntegrationPointsNumber(this->GetIntegrationMethod());
        if (rOutput.size() != number_of_integration_points)
            rOutput.resize(number_of_integration_points);

        if (rVariable == FLUID_FLUX_VECTOR) {
            std::vector<double> permeability_update_factors(number_of_integration_points, 1.0);
            const auto fluid_fluxes = GeoTransportEquationUtilities::CalculateFluidFluxes<TDim, TNumNodes>(
                this->GetGeometry(), this->GetIntegrationMethod(), this->GetProperties(),
                mRetentionLawVector, permeability_update_factors);

            for (unsigned int integration_point = 0;
                 integration_point < number_of_integration_points; ++integration_point) {
                GeoElementUtilities::FillArray1dOutput(rOutput[integration_point],
                                                       fluid_fluxes[integration_point]);
            }
        } else {
            for (unsigned int integration_point = 0;
                 integration_point < number_of_integration_points; ++integration_point) {
                noalias(rOutput[integration_point]) = ZeroVector(3);
            }
        }

        KRATOS_CATCH("")
    }

    void CalculateLocalSystem(MatrixType&        rLeftHandSideMatrix,
                              VectorType&        rRightHandSideVector,
                              const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        CachingDataForCalculator();

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
        CachingDataForCalculator();
        for (const auto& rContribution : mContributions) {
            const auto calculator = CreateCalculator(rContribution, rCurrentProcessInfo);
            noalias(rRightHandSideVector) += calculator->RHSContribution();
        }
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        rLeftHandSideMatrix = ZeroMatrix{TNumNodes, TNumNodes};
        CachingDataForCalculator();
        for (const auto& rContribution : mContributions) {
            const auto calculator = CreateCalculator(rContribution, rCurrentProcessInfo);
            if (const auto LHSContribution = calculator->LHSContribution())
                rLeftHandSideMatrix += *LHSContribution;
        }
    }

    GeometryData::IntegrationMethod GetIntegrationMethod() const override
    {
        switch (this->GetGeometry().GetGeometryOrderType()) {
            using enum GeometryData::KratosGeometryOrderType;
            using enum GeometryData::IntegrationMethod;
        case Kratos_Cubic_Order:
            return GetGeometry().LocalSpaceDimension() == 1 ? GI_GAUSS_3 : GI_GAUSS_4;
        case Kratos_Quartic_Order:
            return GI_GAUSS_5;
        default:
            return GI_GAUSS_2;
        }
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override
    {
        KRATOS_TRY

        CheckUtilities::CheckDomainSize(
            GetGeometry().DomainSize(), Id(),
            GetGeometry().LocalSpaceDimension() == 1 ? "Length" : std::optional<std::string>{});

        const auto r_geometry = this->GetGeometry();
        CheckUtilities::CheckHasNodalSolutionStepData(
            r_geometry,
            {std::cref(WATER_PRESSURE), std::cref(DT_WATER_PRESSURE), std::cref(VOLUME_ACCELERATION)});
        CheckUtilities::CheckHasDofs(r_geometry, {std::cref(WATER_PRESSURE)});

        const auto            r_properties = this->GetProperties();
        const CheckProperties check_properties(r_properties, "material properties at element",
                                               this->Id(), CheckProperties::Bounds::AllInclusive);
        check_properties.Check(DENSITY_WATER);
        check_properties.Check(DENSITY_SOLID);
        constexpr auto max_value_porosity = 1.0;
        check_properties.Check(POROSITY, max_value_porosity);
        check_properties.Check(BULK_MODULUS_SOLID);
        check_properties.Check(BULK_MODULUS_FLUID);
        check_properties.Check(DYNAMIC_VISCOSITY);
        check_properties.Check(BIOT_COEFFICIENT);
        check_properties.Check(BULK_MODULUS_SOLID);

        check_properties.Check(PERMEABILITY_XX);
        if (GetGeometry().LocalSpaceDimension() > 1) {
            check_properties.Check(PERMEABILITY_YY);
            check_properties.Check(PERMEABILITY_XY);
            if constexpr (TDim > 2) {
                check_properties.Check(PERMEABILITY_ZZ);
                check_properties.Check(PERMEABILITY_YZ);
                check_properties.Check(PERMEABILITY_ZX);
            }
        }
        if (this->GetGeometry().WorkingSpaceDimension() == 2)
            CheckUtilities::CheckForNonZeroZCoordinateIn2D(this->GetGeometry());

        return RetentionLaw::Check(mRetentionLawVector, r_properties, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable,
                                      std::vector<Matrix>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const auto& r_geometry = this->GetGeometry();
        const auto number_of_integration_points = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());

        rOutput.resize(number_of_integration_points);

        if (rVariable == PERMEABILITY_MATRIX) {
            // If the permeability of the element is a given property
            BoundedMatrix<double, TDim, TDim> permeability_matrix;
            GeoElementUtilities::FillPermeabilityMatrix(permeability_matrix, this->GetProperties());
            std::fill_n(rOutput.begin(), number_of_integration_points, permeability_matrix);
        } else {
            for (unsigned int i = 0; i < number_of_integration_points; ++i) {
                rOutput[i].resize(TDim, TDim, false);
                noalias(rOutput[i]) = ZeroMatrix(TDim, TDim);
            }
        }

        KRATOS_CATCH("")
    }

private:
    std::vector<CalculationContribution> mContributions;
    IntegrationCoefficientsCalculator    mIntegrationCoefficientsCalculator;
    std::vector<RetentionLaw::Pointer>   mRetentionLawVector;
    Vector                               mIntegrationCoefficients;
    Matrix                               mNContainer;
    Vector                               mDetJCcontainer;
    std::vector<double>                  mFluidPressures;

    array_1d<double, TNumNodes> GetNodalValuesOf(const Variable<double>& rNodalVariable) const
    {
        auto        result     = array_1d<double, TNumNodes>{};
        const auto& r_geometry = GetGeometry();
        std::ranges::transform(r_geometry, result.begin(), [&rNodalVariable](const auto& node) {
            return node.FastGetSolutionStepValue(rNodalVariable);
        });
        return result;
    }

    std::vector<double> CalculateIntegrationCoefficients(const Vector& rDetJs) const
    {
        const GeometryType::IntegrationPointsArrayType& integration_points =
            this->GetGeometry().IntegrationPoints(GetIntegrationMethod());
        return mIntegrationCoefficientsCalculator.Run<>(integration_points, rDetJs, this);
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
        if (GetGeometry().LocalSpaceDimension() == 1) {
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
        } else {
            for (unsigned int integration_point_index = 0;
                 integration_point_index < number_integration_points; ++integration_point_index) {
                GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                    body_acceleration, rNContainer, volume_acceleration, integration_point_index);
                projected_gravity.emplace_back(body_acceleration);
            }
        }

        return projected_gravity;
    }

    std::unique_ptr<IntegrationCoefficientModifier> CloneIntegrationCoefficientModifier() const
    {
        return mIntegrationCoefficientsCalculator.CloneModifier();
    }

    std::unique_ptr<ContributionCalculator<TNumNodes>> CreateCalculator(const CalculationContribution& rContribution,
                                                                        const ProcessInfo& rCurrentProcessInfo)
    {
        switch (rContribution) {
            using enum CalculationContribution;
        case Permeability:
            return std::make_unique<PermeabilityCalculator<TNumNodes>>(CreatePermeabilityInputProvider());
        case Compressibility:
            if (GetProperties()[RETENTION_LAW] == "PressureFilterLaw") {
                return std::make_unique<FilterCompressibilityCalculator<TNumNodes>>(
                    CreateFilterCompressibilityInputProvider(rCurrentProcessInfo));
            }
            return std::make_unique<CompressibilityCalculator<TNumNodes>>(
                CreateCompressibilityInputProvider(rCurrentProcessInfo));
        case FluidBodyFlow:
            return std::make_unique<FluidBodyFlowCalculator<TNumNodes>>(CreateFluidBodyFlowInputProvider());
        default:
            KRATOS_ERROR << "Unknown contribution" << std::endl;
        }
    }

    void CachingDataForCalculator()
    {
        mIntegrationCoefficients = CalculateIntegrationCoefficients();
        mNContainer              = CalculateNContainer();
        mFluidPressures          = CalculateFluidPressure();
    }

    typename CompressibilityCalculator<TNumNodes>::InputProvider CreateCompressibilityInputProvider(const ProcessInfo& rCurrentProcessInfo)
    {
        return typename CompressibilityCalculator<TNumNodes>::InputProvider(
            MakePropertiesGetter(), MakeRetentionLawsGetter(), GetNContainer(),
            GetIntegrationCoefficients(), MakeMatrixScalarFactorGetter(rCurrentProcessInfo),
            MakeNodalVariableGetter(), GetFluidPressures());
    }

    typename FilterCompressibilityCalculator<TNumNodes>::InputProvider CreateFilterCompressibilityInputProvider(
        const ProcessInfo& rCurrentProcessInfo)
    {
        return typename FilterCompressibilityCalculator<TNumNodes>::InputProvider(
            MakePropertiesGetter(), GetNContainer(), GetIntegrationCoefficients(),
            MakeProjectedGravityForIntegrationPointsGetter(),
            MakeMatrixScalarFactorGetter(rCurrentProcessInfo), MakeNodalVariableGetter());
    }

    typename PermeabilityCalculator<TNumNodes>::InputProvider CreatePermeabilityInputProvider()
    {
        return typename PermeabilityCalculator<TNumNodes>::InputProvider(
            MakePropertiesGetter(), MakeRetentionLawsGetter(), GetIntegrationCoefficients(),
            MakeNodalVariableGetter(), MakeShapeFunctionLocalGradientsGetter(), GetFluidPressures());
    }

    typename FluidBodyFlowCalculator<TNumNodes>::InputProvider CreateFluidBodyFlowInputProvider()
    {
        return typename FluidBodyFlowCalculator<TNumNodes>::InputProvider(
            MakePropertiesGetter(), MakeRetentionLawsGetter(), GetIntegrationCoefficients(),
            MakeProjectedGravityForIntegrationPointsGetter(), MakeShapeFunctionLocalGradientsGetter(),
            MakeLocalSpaceDimensionGetter(), GetFluidPressures());
    }

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

    Matrix CalculateNContainer()
    {
        return GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());
    }

    auto GetIntegrationCoefficients()
    {
        return [this]() -> const Vector& { return mIntegrationCoefficients; };
    }

    Vector CalculateIntegrationCoefficients()
    {
        GetGeometry().DeterminantOfJacobian(mDetJCcontainer, this->GetIntegrationMethod());
        return mIntegrationCoefficientsCalculator.Run<Vector>(
            GetGeometry().IntegrationPoints(GetIntegrationMethod()), mDetJCcontainer, this);
    }

    auto GetFluidPressures()
    {
        return [this]() -> const std::vector<double>& { return mFluidPressures; };
    }

    std::vector<double> CalculateFluidPressure()
    {
        array_1d<double, TNumNodes> pressure_vector;
        VariablesUtilities::GetNodalValues(this->GetGeometry(), WATER_PRESSURE, pressure_vector.begin());
        return GeoTransportEquationUtilities::CalculateFluidPressures(mNContainer, pressure_vector);
    }

    auto MakeProjectedGravityForIntegrationPointsGetter()
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
        return
            [this](const Variable<double>& variable) -> Vector { return GetNodalValuesOf(variable); };
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

    [[nodiscard]] DofsVectorType GetDofs() const
    {
        return Geo::DofUtilities::ExtractDofsFromNodes(GetGeometry(), WATER_PRESSURE);
    }

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
        std::vector<size_t> contributions(mContributions.size());
        std::transform(mContributions.begin(), mContributions.end(), contributions.begin(),
                       [](CalculationContribution c) { return static_cast<size_t>(c); });
        rSerializer.save("Contribution", contributions);
        rSerializer.save("RetentionlawVector", mRetentionLawVector);
        rSerializer.save("IntegrationCoefficientsCalculator", mIntegrationCoefficientsCalculator);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
        std::vector<size_t> contributions;
        rSerializer.load("Contribution", contributions);
        mContributions.clear();
        std::transform(contributions.begin(), contributions.end(), std::back_inserter(mContributions),
                       [](size_t c) { return static_cast<CalculationContribution>(c); });
        rSerializer.load("RetentionlawVector", mRetentionLawVector);
        rSerializer.save("IntegrationCoefficientsCalculator", mIntegrationCoefficientsCalculator);
    }
};

} // namespace Kratos
