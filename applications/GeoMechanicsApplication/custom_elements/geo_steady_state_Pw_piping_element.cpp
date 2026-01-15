// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#include "custom_elements/geo_steady_state_Pw_piping_element.h"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/dof_utilities.hpp"
#include "custom_utilities/element_utilities.hpp"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/transport_equation_utilities.hpp"
#include "custom_utilities/variables_utilities.hpp"
#include "geo_mechanics_application_variables.h"
#include "includes/cfd_variables.h"

#include <numeric>

namespace Kratos
{
using namespace std::string_literals;

template <unsigned int TDim, unsigned int TNumNodes>
GeoSteadyStatePwPipingElement<TDim, TNumNodes>::GeoSteadyStatePwPipingElement(IndexType NewId)
    : Element(NewId)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoSteadyStatePwPipingElement<TDim, TNumNodes>::GeoSteadyStatePwPipingElement(IndexType NewId,
                                                                              const GeometryType::Pointer& rpGeometry)
    : Element(NewId, rpGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoSteadyStatePwPipingElement<TDim, TNumNodes>::GeoSteadyStatePwPipingElement(
    IndexType NewId, const GeometryType::Pointer& rpGeometry, const PropertiesType::Pointer& rpProperties)
    : Element(NewId, rpGeometry, rpProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoSteadyStatePwPipingElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                        const NodesArrayType& rThisNodes,
                                                                        PropertiesType::Pointer pProperties) const
{
    return this->Create(NewId, GetGeometry().Create(rThisNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer GeoSteadyStatePwPipingElement<TDim, TNumNodes>::Create(IndexType NewId,
                                                                        GeometryType::Pointer pGeometry,
                                                                        PropertiesType::Pointer pProperties) const
{
    return make_intrusive<GeoSteadyStatePwPipingElement>(NewId, pGeometry, pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                const ProcessInfo&) const
{
    rElementalDofList = GetDofs();
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                                                      const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                                          VectorType& rRightHandSideVector,
                                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    Vector det_J_container;
    GetGeometry().DeterminantOfJacobian(det_J_container, this->GetIntegrationMethod());
    GeometryType::ShapeFunctionsGradientsType dN_dX_container =
        GetGeometry().ShapeFunctionsLocalGradients(this->GetIntegrationMethod());
    std::transform(dN_dX_container.begin(), dN_dX_container.end(), det_J_container.begin(),
                   dN_dX_container.begin(), std::divides<>());
    const Matrix& r_N_container = GetGeometry().ShapeFunctionsValues(GetIntegrationMethod());

    const auto integration_coefficients = CalculateIntegrationCoefficients(det_J_container);
    const auto permeability_matrix = CalculatePermeabilityMatrix(dN_dX_container, integration_coefficients);
    AddContributionsToLhsMatrix(rLeftHandSideMatrix, permeability_matrix);

    const auto fluid_body_vector =
        CalculateFluidBodyVector(r_N_container, dN_dX_container, integration_coefficients);
    AddContributionsToRhsVector(rRightHandSideVector, permeability_matrix, fluid_body_vector);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod GeoSteadyStatePwPipingElement<TDim, TNumNodes>::GetIntegrationMethod() const
{
    KRATOS_ERROR_IF(TNumNodes != 2)
        << "Can't return integration method: unexpected number of nodes: " << TNumNodes << std::endl;
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

template <unsigned int TDim, unsigned int TNumNodes>
int GeoSteadyStatePwPipingElement<TDim, TNumNodes>::Check(const ProcessInfo&) const
{
    KRATOS_TRY
    CheckUtilities::CheckDomainSize(GetGeometry().DomainSize(), Id());
    CheckUtilities::CheckHasNodalSolutionStepData(this->GetGeometry(), {std::cref(WATER_PRESSURE)});
    CheckUtilities::CheckHasDofs(this->GetGeometry(), {std::cref(WATER_PRESSURE)});

    const CheckProperties check_properties(this->GetProperties(), "material properties at element",
                                           this->Id(), CheckProperties::Bounds::AllInclusive);
    check_properties.Check(DENSITY_WATER);
    check_properties.Check(DYNAMIC_VISCOSITY);
    check_properties.Check(PIPE_HEIGHT);
    if constexpr (TDim == 3) {
        check_properties.Check(PIPE_WIDTH_FACTOR);
    }

    if constexpr (TDim == 2) {
        CheckUtilities::CheckForNonZeroZCoordinateIn2D(this->GetGeometry());
    }

    KRATOS_CATCH("")

    return 0;
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    Element::Initialize(rCurrentProcessInfo);
    // all these except the PIPE_ELEMENT_LENGTH seem to be in the erosion_process_strategy only. Why do this: it is used in output for dGeoFlow
    this->SetValue(PIPE_ELEMENT_LENGTH, CalculateLength(this->GetGeometry()));

    if (!mIsInitialized) {
        this->SetValue(PIPE_EROSION, false);
        constexpr double small_pipe_height = 1e-10;
        this->SetValue(PIPE_HEIGHT, small_pipe_height);
        this->SetValue(PREV_PIPE_HEIGHT, small_pipe_height);
        this->SetValue(DIFF_PIPE_HEIGHT, 0.);
        this->SetValue(PIPE_ACTIVE, false);
        mIsInitialized = true;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculateEquilibriumPipeHeight(
    const PropertiesType& rProp, const GeometryType& rGeom, double)
{
    // calculate head gradient over element ( now without abs in CalculateHeadGradient )
    const auto head_gradient = CalculateHeadGradient(rProp, rGeom);
    // return infinite when the head gradient dh/dx is 0
    if (std::abs(head_gradient) < std::numeric_limits<double>::epsilon()) return 1e10;

    const auto particle_d = GeoTransportEquationUtilities::CalculateParticleDiameter(rProp);
    // for a more generic element calculate slope of pipe (in degrees! see formula), currently pipe is assumed to be horizontal
    constexpr double pipe_slope = 0.0;
    const auto       theta      = rProp[PIPE_THETA];

    return rProp[PIPE_MODEL_FACTOR] * (Globals::Pi / 3.0) * particle_d *
           (rProp[DENSITY_SOLID] / rProp[DENSITY_WATER] - 1.0) * rProp[PIPE_ETA] *
           (std::sin(MathUtils<>::DegreesToRadians(theta + pipe_slope)) /
            std::cos(MathUtils<>::DegreesToRadians(theta))) /
           std::abs(head_gradient);
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable, std::vector<bool>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == PIPE_ACTIVE) {
        const auto number_of_integration_points =
            this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        rValues.resize(number_of_integration_points);
        std::fill_n(rValues.begin(), number_of_integration_points, this->GetValue(rVariable));
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == PIPE_HEIGHT) {
        const auto number_of_integration_points =
            this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        rValues.resize(number_of_integration_points);
        std::fill_n(rValues.begin(), number_of_integration_points, this->GetValue(rVariable));
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable, std::vector<array_1d<double, 3>>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == FLUID_FLUX_VECTOR || rVariable == LOCAL_FLUID_FLUX_VECTOR) {
        const auto number_of_integration_points =
            this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        rValues.resize(number_of_integration_points);
        const auto dynamic_viscosity_inverse = 1.0 / this->GetProperties()[DYNAMIC_VISCOSITY];
        const auto permeability_matrix       = FillPermeabilityMatrix(this->GetValue(PIPE_HEIGHT));
        auto       local_fluid_flux_vector   = array_1d<double, 3>(3, 0.0);
        local_fluid_flux_vector[0] = -dynamic_viscosity_inverse * permeability_matrix(0, 0) *
                                     CalculateHeadGradient(this->GetProperties(), this->GetGeometry());
        std::fill_n(rValues.begin(), number_of_integration_points, local_fluid_flux_vector);

        if (rVariable == LOCAL_FLUID_FLUX_VECTOR) return;

        // For the global fluid flux vector the local fluid flux vector should be rotated to the direction of the element
        const auto& r_integration_points =
            this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
        for (std::size_t i = 0; i < number_of_integration_points; ++i) {
            Matrix jacobian;
            this->GetGeometry().Jacobian(jacobian, r_integration_points[i]);
            const auto tangential_vector =
                GeoMechanicsMathUtilities::Normalized(Vector{column(jacobian, 0)});
            std::transform(tangential_vector.begin(), tangential_vector.end(), rValues[i].begin(),
                           [&local_fluid_flux_vector](auto component) {
                return component * local_fluid_flux_vector[0];
            });
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == PERMEABILITY_MATRIX || rVariable == LOCAL_PERMEABILITY_MATRIX) {
        // permeability matrix
        const auto number_of_integration_points =
            this->GetGeometry().IntegrationPointsNumber(this->GetIntegrationMethod());
        rValues.resize(number_of_integration_points);
        std::fill_n(rValues.begin(), number_of_integration_points,
                    FillPermeabilityMatrix(this->GetValue(PIPE_HEIGHT)));
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string GeoSteadyStatePwPipingElement<TDim, TNumNodes>::Info() const
{
    return "GeoSteadyStatePwPipingElement"s;
}

template <unsigned int TDim, unsigned int TNumNodes>
double GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculateLength(const GeometryType& Geom)
{
    // currently length is only calculated in x direction, to be changed for inclined pipes
    return std::abs(Geom.GetPoint(1)[0] - Geom.GetPoint(0)[0]);
}

template <unsigned int TDim, unsigned int TNumNodes>
double GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculateHeadGradient(const PropertiesType& rProp,
                                                                             const GeometryType& rGeom)
{
    const auto nodal_heads = GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(rGeom, rProp);
    Matrix shapefunctions_local_gradient;
    // iso-parametric derivative
    GetGeometry().ShapeFunctionsLocalGradients(shapefunctions_local_gradient, GetGeometry().Center());
    // local derivative
    shapefunctions_local_gradient /= GetGeometry().DeterminantOfJacobian(GetGeometry().Center());
    return Vector{prod(trans(shapefunctions_local_gradient), nodal_heads)}[0];
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::AddContributionsToLhsMatrix(
    MatrixType& rLeftHandSideMatrix, const BoundedMatrix<double, TNumNodes, TNumNodes>& rPermeabilityMatrix)
{
    rLeftHandSideMatrix = rPermeabilityMatrix;
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::AddContributionsToRhsVector(
    VectorType&                                        rRightHandSideVector,
    const BoundedMatrix<double, TNumNodes, TNumNodes>& rPermeabilityMatrix,
    const array_1d<double, TNumNodes>&                 rFluidBodyVector) const
{
    const auto permeability_vector = array_1d<double, TNumNodes>{-prod(
        rPermeabilityMatrix, VariablesUtilities::GetNodalValuesOf<TNumNodes>(WATER_PRESSURE, GetGeometry()))};
    rRightHandSideVector           = permeability_vector + rFluidBodyVector;
}

template <unsigned int TDim, unsigned int TNumNodes>
Vector GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculateIntegrationCoefficients(const Vector& rDetJContainer) const
{
    const auto& r_integration_points = GetGeometry().IntegrationPoints(GetIntegrationMethod());

    auto result = Vector{r_integration_points.size()};
    // all governed by PIPE_HEIGHT and element length so without CROSS_AREA
    std::ranges::transform(r_integration_points, rDetJContainer, result.begin(),
                           [](const auto& rIntegrationPoint, const auto& rDetJ) {
        return rIntegrationPoint.Weight() * rDetJ;
    });
    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
Matrix GeoSteadyStatePwPipingElement<TDim, TNumNodes>::FillPermeabilityMatrix(double pipe_height) const
{
    if constexpr (TDim == 2) {
        return ScalarMatrix{1, 1, std::pow(pipe_height, 3) / 12.0};
    }
    return ScalarMatrix{1, 1, this->GetProperties()[PIPE_WIDTH_FACTOR] * std::pow(pipe_height, 4) / 12.0};
}

template <unsigned int TDim, unsigned int TNumNodes>
BoundedMatrix<double, TNumNodes, TNumNodes> GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculatePermeabilityMatrix(
    const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients, const Vector& rIntegrationCoefficients) const
{
    const auto&  r_properties              = GetProperties();
    const double dynamic_viscosity_inverse = 1.0 / r_properties[DYNAMIC_VISCOSITY];
    auto         constitutive_matrix       = FillPermeabilityMatrix(this->GetValue(PIPE_HEIGHT));

    auto result = BoundedMatrix<double, TNumNodes, TNumNodes>{ZeroMatrix{TNumNodes, TNumNodes}};
    for (unsigned int integration_point_index = 0;
         integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
         ++integration_point_index) {
        result += GeoTransportEquationUtilities::CalculatePermeabilityMatrix<TDim, TNumNodes>(
            rShapeFunctionGradients[integration_point_index], dynamic_viscosity_inverse,
            constitutive_matrix, 1.0, rIntegrationCoefficients[integration_point_index]);
    }

    return result;
}

template <unsigned int TDim, unsigned int TNumNodes>
array_1d<double, TNumNodes> GeoSteadyStatePwPipingElement<TDim, TNumNodes>::CalculateFluidBodyVector(
    const Matrix&                                    rNContainer,
    const GeometryType::ShapeFunctionsGradientsType& rShapeFunctionGradients,
    const Vector&                                    rIntegrationCoefficients) const
{
    const std::size_t number_integration_points =
        GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
    GeometryType::JacobiansType J_container;
    J_container.resize(number_integration_points, false);
    for (std::size_t i = 0; i < number_integration_points; ++i) {
        J_container[i].resize(GetGeometry().WorkingSpaceDimension(), GetGeometry().LocalSpaceDimension(), false);
    }
    GetGeometry().Jacobian(J_container, this->GetIntegrationMethod());

    const auto& r_properties        = GetProperties();
    auto        constitutive_matrix = FillPermeabilityMatrix(r_properties[PIPE_HEIGHT]);

    array_1d<double, TNumNodes * TDim> volume_acceleration;
    GeoElementUtilities::GetNodalVariableVector<TDim, TNumNodes>(volume_acceleration, GetGeometry(),
                                                                 VOLUME_ACCELERATION);

    array_1d<double, TNumNodes> fluid_body_vector = ZeroVector(TNumNodes);
    for (unsigned int integration_point_index = 0;
         integration_point_index < GetGeometry().IntegrationPointsNumber(GetIntegrationMethod());
         ++integration_point_index) {
        array_1d<double, TDim> body_acceleration;
        GeoElementUtilities::InterpolateVariableWithComponents<TDim, TNumNodes>(
            body_acceleration, rNContainer, volume_acceleration, integration_point_index);

        array_1d<double, TDim> tangent_vector = column(J_container[integration_point_index], 0);
        tangent_vector /= norm_2(tangent_vector);

        const auto projected_gravity =
            array_1d<double, 1>(1, std::inner_product(tangent_vector.begin(), tangent_vector.end(),
                                                      body_acceleration.begin(), 0.0));
        const auto N = Vector{row(rNContainer, integration_point_index)};
        fluid_body_vector +=
            r_properties[DENSITY_WATER] *
            prod(prod(rShapeFunctionGradients[integration_point_index], constitutive_matrix), projected_gravity) *
            rIntegrationCoefficients[integration_point_index] / r_properties[DYNAMIC_VISCOSITY];
    }
    return fluid_body_vector;
}

template <unsigned int TDim, unsigned int TNumNodes>
Element::DofsVectorType GeoSteadyStatePwPipingElement<TDim, TNumNodes>::GetDofs() const
{
    return Geo::DofUtilities::ExtractDofsFromNodes(GetGeometry(), WATER_PRESSURE);
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element)
    rSerializer.save("mIsInitialized", mIsInitialized);
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoSteadyStatePwPipingElement<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element)
    rSerializer.load("mIsInitialized", mIsInitialized);
}
template class GeoSteadyStatePwPipingElement<2, 2>;
template class GeoSteadyStatePwPipingElement<3, 2>;

} // Namespace Kratos