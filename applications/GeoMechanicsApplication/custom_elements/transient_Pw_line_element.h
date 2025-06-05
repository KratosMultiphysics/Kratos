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
#include "custom_utilities/check_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/dof_utilities.h"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "custom_utilities/variables_utilities.hpp"
#include "filter_compressibility_calculator.h"
#include "fluid_body_flow_calculator.h"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/constitutive_law.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "integration_coefficients_calculator.h"
#include "permeability_calculator.h"
#include "utilities/geometry_utilities.h"

#include <optional>

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
class KRATOS_API(GEO_MECHANICS_APPLICATION) TransientPwLineElement : public Element
{
public:
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(TransientPwLineElement);

    explicit TransientPwLineElement(IndexType NewId = 0) : Element(NewId) {}

    TransientPwLineElement(IndexType                                       NewId,
                           const GeometryType::Pointer&                    pGeometry,
                           const std::vector<CalculationContribution>&     rContributions,
                           std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
        : Element(NewId, pGeometry),
          mContributions(rContributions),
          mIntegrationCoefficientsCalculator{std::move(pCoefficientModifier)}
    {
    }

    TransientPwLineElement(IndexType                                       NewId,
                           const GeometryType::Pointer&                    pGeometry,
                           const PropertiesType::Pointer&                  pProperties,
                           const std::vector<CalculationContribution>&     rContributions,
                           std::unique_ptr<IntegrationCoefficientModifier> pCoefficientModifier)
        : Element(NewId, pGeometry, pProperties),
          mContributions(rContributions),
          mIntegrationCoefficientsCalculator{std::move(pCoefficientModifier)}
    {
    }

    Element::Pointer Create(IndexType NewId, const NodesArrayType& rThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<TransientPwLineElement>(NewId, GetGeometry().Create(rThisNodes),
                                                      pProperties, mContributions,
                                                      this->CloneIntegrationCoefficientModifier());
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return make_intrusive<TransientPwLineElement>(NewId, pGeom, pProperties, mContributions,
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
        if (GetGeometry().LocalSpaceDimension() != 1) {
            mConstitutiveLawVector.resize(GetGeometry().IntegrationPointsNumber(GetIntegrationMethod()));
            for (auto& constitutive_law : mConstitutiveLawVector) {
                constitutive_law = nullptr;
            }
        }
        mRetentionLawVector.resize(GetGeometry().IntegrationPointsNumber(GetIntegrationMethod()));

        for (auto& r_retention_law : mRetentionLawVector) {
            r_retention_law = RetentionLawFactory::Clone(GetProperties());
        }
    }

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY
        // Reset hydraulic discharge
        for (auto& r_node : this->GetGeometry()) {
            r_node.FastGetSolutionStepValue(HYDRAULIC_DISCHARGE) = 0.0;
        }
        KRATOS_CATCH("")
    }

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_TRY

        this->CalculateHydraulicDischarge(rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateHydraulicDischarge(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if (GetGeometry().LocalSpaceDimension() == 1) {
            // there is no need to calculate the discharge for a line element
            return;
        }
        std::vector<array_1d<double, 3>> fluid_flux;
        this->CalculateOnIntegrationPoints(FLUID_FLUX_VECTOR, fluid_flux, rCurrentProcessInfo);

        const GeometryType& r_geometry = this->GetGeometry();
        const IndexType     number_of_integration_points =
            r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
            r_geometry.IntegrationPoints(GetIntegrationMethod());

        // Gradient of shape functions and determinant of Jacobian
        Matrix                                    grad_Np_T(TNumNodes, TDim);
        Vector                                    det_J_container(number_of_integration_points);
        GeometryType::ShapeFunctionsGradientsType dN_dx_container;
        r_geometry.ShapeFunctionsIntegrationPointsGradients(dN_dx_container, det_J_container,
                                                            GetIntegrationMethod());
        const auto integration_coefficients =
            this->CalculateIntegrationCoefficients(IntegrationPoints, det_J_container);

        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            noalias(grad_Np_T) = dN_dx_container[integration_point];

            auto integration_coefficient = integration_coefficients[integration_point];

            for (unsigned int node = 0; node < TNumNodes; ++node) {
                double hydraulic_discharge = 0;
                for (unsigned int direction = 0; direction < TDim; ++direction) {
                    hydraulic_discharge +=
                        grad_Np_T(node, direction) * fluid_flux[integration_point][direction];
                }

                hydraulic_discharge *= integration_coefficient;
                hydraulic_discharge += r_geometry[node].FastGetSolutionStepValue(HYDRAULIC_DISCHARGE);
                ThreadSafeNodeWrite(this->GetGeometry()[node], HYDRAULIC_DISCHARGE, hydraulic_discharge);
            }
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
            RetentionLaw::Parameters RetentionParameters(r_properties);
            Vector                   Np(TNumNodes);
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

    void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                      std::vector<Vector>&    rOutput,
                                      const ProcessInfo&      rCurrentProcessInfo) override
    {
        KRATOS_TRY

        const auto& r_geometry = this->GetGeometry();
        const auto number_of_integration_points = r_geometry.IntegrationPointsNumber(GetIntegrationMethod());
        const auto& r_properties = this->GetProperties();

        rOutput.resize(number_of_integration_points);

        if (r_properties.Has(rVariable)) {
            // Map initial material property to integration points, as required for the output
            std::fill_n(rOutput.begin(), mConstitutiveLawVector.size(), r_properties.GetValue(rVariable));
        } else {
            for (unsigned int integration_point = 0;
                 integration_point < mConstitutiveLawVector.size(); ++integration_point)
                rOutput[integration_point] = mConstitutiveLawVector[integration_point]->GetValue(
                    rVariable, rOutput[integration_point]);
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
            const auto          fluid_fluxes =
                this->CalculateFluidFluxes(permeability_update_factors, rCurrentProcessInfo);

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

    std::vector<array_1d<double, TDim>> CalculateFluidFluxes(const std::vector<double>& rPermeabilityUpdateFactors,
                                                             const ProcessInfo& rCurrentProcessInfo)
    {
        const GeometryType& r_geometry = this->GetGeometry();
        const IndexType     number_of_integration_points =
            r_geometry.IntegrationPointsNumber(this->GetIntegrationMethod());

        std::vector<array_1d<double, TDim>> fluid_fluxes;
        fluid_fluxes.reserve(number_of_integration_points);
        array_1d<double, TNumNodes> pressure_vector;
        VariablesUtilities::GetNodalValues(r_geometry, WATER_PRESSURE, pressure_vector.begin());
        Matrix N_container(number_of_integration_points, TNumNodes);
        N_container = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
        const PropertiesType&             r_properties = this->GetProperties();
        BoundedMatrix<double, TDim, TDim> permeability_matrix;
        GeoElementUtilities::FillPermeabilityMatrix(permeability_matrix, r_properties);

        auto relative_permeability_values =
            this->CalculateRelativePermeabilityValues(GeoTransportEquationUtilities::CalculateFluidPressures(
                N_container, pressure_vector));
        std::transform(relative_permeability_values.cbegin(), relative_permeability_values.cend(),
                       rPermeabilityUpdateFactors.cbegin(), relative_permeability_values.begin(),
                       std::multiplies<>{});
        const auto dynamic_viscosity_inverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
        array_1d<double, TNumNodes * TDim> volume_acceleration;
        GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(
            volume_acceleration, r_geometry, VOLUME_ACCELERATION);
        array_1d<double, TDim> body_acceleration;
        Matrix grad_Np_T(TNumNodes, TDim);
        Vector                                    det_J_Container(number_of_integration_points);
        GeometryType::ShapeFunctionsGradientsType dN_dx_container;
        r_geometry.ShapeFunctionsIntegrationPointsGradients(
            dN_dx_container, det_J_Container, this->GetIntegrationMethod());
        for (unsigned int integration_point = 0; integration_point < number_of_integration_points;
             ++integration_point) {
            noalias(grad_Np_T) = dN_dx_container[integration_point];

            GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
                body_acceleration, N_container, volume_acceleration, integration_point);

            array_1d<double, TDim> GradPressureTerm = prod(trans(grad_Np_T), pressure_vector);
            GradPressureTerm += PORE_PRESSURE_SIGN_FACTOR * r_properties[DENSITY_WATER] * body_acceleration;

            fluid_fluxes.push_back(PORE_PRESSURE_SIGN_FACTOR * dynamic_viscosity_inverse * relative_permeability_values[integration_point] *
                                   prod(permeability_matrix, GradPressureTerm));
        }

        return fluid_fluxes;
    }

    std::vector<double> CalculateRelativePermeabilityValues(const std::vector<double>& rFluidPressures) const
    {
        KRATOS_ERROR_IF_NOT(rFluidPressures.size() == mRetentionLawVector.size());

        auto retention_law_params = RetentionLaw::Parameters{this->GetProperties()};

        auto result = std::vector<double>{};
        result.reserve(rFluidPressures.size());
        std::transform(mRetentionLawVector.begin(), mRetentionLawVector.end(),
                       rFluidPressures.begin(), std::back_inserter(result),
                       [&retention_law_params](const auto& pRetentionLaw, auto FluidPressure) {
            retention_law_params.SetFluidPressure(FluidPressure);
            return pRetentionLaw->CalculateRelativePermeability(retention_law_params);
        });
        return result;
    }

    void CalculateDerivativesOnInitialConfiguration(
        double& rDetJ, Matrix& rJ0, Matrix& rInvJ0, Matrix& rDNu_DX0, unsigned int IntegrationPointIndex) const
    {
        KRATOS_TRY

        const auto& r_geometry           = this->GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints(GetIntegrationMethod());

        GeometryUtils::JacobianOnInitialConfiguration(
            r_geometry, r_integration_points[IntegrationPointIndex], rJ0);
        const auto& r_dn_de =
            r_geometry.ShapeFunctionsLocalGradients(GetIntegrationMethod())[IntegrationPointIndex];
        MathUtils<>::InvertMatrix(rJ0, rInvJ0, rDetJ);
        GeometryUtils::ShapeFunctionsGradients(r_dn_de, rInvJ0, rDNu_DX0);

        KRATOS_CATCH("")
    }

    std::vector<double> CalculateIntegrationCoefficients(const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
                                                         const Vector& rDetJs) const
    {
        return mIntegrationCoefficientsCalculator.Run<>(rIntegrationPoints, rDetJs, this);
    }

    template <class TValueType>
    inline void ThreadSafeNodeWrite(NodeType& rNode, const Variable<TValueType>& Var, const TValueType Value)
    {
        rNode.SetLock();
        rNode.FastGetSolutionStepValue(Var) = Value;
        rNode.UnSetLock();
    }

    //////////////////////////////////////////////////

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
            return GetGeometry().LocalSpaceDimension() == 1 ? GeometryData::IntegrationMethod::GI_GAUSS_3
                                                            : IntegrationMethod::GI_GAUSS_4;
        case GeometryData::Kratos_Quartic_Order:
            return GeometryData::IntegrationMethod::GI_GAUSS_5;
        default:
            return GeometryData::IntegrationMethod::GI_GAUSS_2;
        }
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override
    {
        KRATOS_TRY

        CheckUtilities::CheckDomainSize(
            GetGeometry().DomainSize(), Id(),
            GetGeometry().LocalSpaceDimension() == 1 ? "Length" : std::optional<std::string>{});
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

    std::vector<ConstitutiveLaw::Pointer> GetConstitutiveLawVector() const
    {
        return mConstitutiveLawVector;
    }

    std::vector<RetentionLaw::Pointer> GetRetentionLawVector() const { return mRetentionLawVector; }

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
            for (unsigned int i = 0; i < mRetentionLawVector.size(); ++i) {
                rOutput[i].resize(TDim, TDim, false);
                noalias(rOutput[i]) = ZeroMatrix(TDim, TDim);
            }
        }

        KRATOS_CATCH("")
    }

private:
    std::vector<CalculationContribution>  mContributions;
    IntegrationCoefficientsCalculator     mIntegrationCoefficientsCalculator;
    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLawVector;
    std::vector<RetentionLaw::Pointer>    mRetentionLawVector;

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
        if (GetGeometry().LocalSpaceDimension() > 1) {
            CheckProperty(PERMEABILITY_YY);
            CheckProperty(PERMEABILITY_XY);
            if constexpr (TDim > 2) {
                CheckProperty(PERMEABILITY_ZZ);
                CheckProperty(PERMEABILITY_YZ);
                CheckProperty(PERMEABILITY_ZX);
            }
        }
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
            MakeNodalVariableGetter(), MakeShapeFunctionLocalGradientsGetter(), MakeNContainerGetter());
    }

    FluidBodyFlowCalculator::InputProvider CreateFluidBodyFlowInputProvider()
    {
        return FluidBodyFlowCalculator::InputProvider(
            MakePropertiesGetter(), MakeRetentionLawsGetter(), MakeIntegrationCoefficientsGetter(),
            MakeProjectedGravityForIntegrationPointsGetter(), MakeNContainerGetter(),
            MakeShapeFunctionLocalGradientsGetter(), MakeLocalSpaceDimensionGetter(),
            MakeNodalVariableGetter());
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
            return mIntegrationCoefficientsCalculator.Run<Vector>(
                GetGeometry().IntegrationPoints(GetIntegrationMethod()), det_J_container, this);
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
            GeometryType::ShapeFunctionsGradientsType dN_dX_container;
            if (GetGeometry().LocalSpaceDimension() == 1) {
                dN_dX_container = GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
                std::transform(dN_dX_container.begin(), dN_dX_container.end(),
                               det_J_container.begin(), dN_dX_container.begin(), std::divides<>());
            } else {
                GetGeometry().ShapeFunctionsIntegrationPointsGradients(
                    dN_dX_container, det_J_container, this->GetIntegrationMethod());
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
        rSerializer.save("RetentionlawVector", mRetentionLawVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    }
};

} // namespace Kratos
