//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "geometries/geometry_data.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "custom_elements/data_containers/k_omega_sst/element_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "omega_k_based_wall_condition_data.h"

namespace Kratos
{
namespace KOmegaSSTWallConditionData
{
template<unsigned int TDim>
const Variable<double>& OmegaKBasedWallConditionData<TDim>::GetScalarVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

template<unsigned int TDim>
void OmegaKBasedWallConditionData<TDim>::Check(
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = rCondition.GetGeometry();
    const auto& r_properties = rCondition.GetProperties();
    const int number_of_nodes = r_geometry.PointsNumber();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(VON_KARMAN))
        << "VON_KARMAN is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(r_properties.Has(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT))
        << "RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT is not found in condition properties [ Condition.Id() = "
        << rCondition.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}

template<unsigned int TDim>
GeometryData::IntegrationMethod OmegaKBasedWallConditionData<TDim>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_1;
}

template<unsigned int TDim>
void OmegaKBasedWallConditionData<TDim>::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const double beta_star = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    const double sigma_omega_1 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1];
    const double sigma_omega_2 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2];
    mKappa = rCurrentProcessInfo[VON_KARMAN];
    mCmu25 = std::pow(beta_star, 0.25);

    KRATOS_ERROR_IF(!(this->GetGeometry().Has(RANS_Y_PLUS)))
        << "RANS_Y_PLUS value is not set at " << this->GetGeometry() << "\n";

    mDensity = this->GetElementProperties()[DENSITY];

    const double y_plus_limit = this->GetConditionProperties()[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];

    mYPlus = std::max(this->GetGeometry().GetValue(RANS_Y_PLUS), y_plus_limit);

    // Blended sigma is computed on the center of the condition, therefore
    // it will be a constant for all gauss point evaluations
    const auto& parent_element_geometry = this->GetGeometry().GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry();
    Vector W;
    Matrix N;
    ShapeFunctionDerivativesArrayType dNdX;
    RansCalculationUtilities::CalculateGeometryData(
        parent_element_geometry, GeometryData::IntegrationMethod::GI_GAUSS_1, W, N, dNdX);

    const Vector& r_N = row(N, 0);

    auto& cl_parameters = this->GetConstitutiveLawParameters();
    cl_parameters.SetShapeFunctionsValues(r_N);

    double tke, omega, nu, wall_distance;
    this->GetConstitutiveLaw().CalculateValue(cl_parameters, EFFECTIVE_VISCOSITY, nu);
    nu /= mDensity;

    FluidCalculationUtilities::EvaluateInPoint(
        parent_element_geometry, r_N,
        std::tie(tke, TURBULENT_KINETIC_ENERGY),
        std::tie(omega, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(wall_distance, DISTANCE));

    array_1d<double, TDim> k_gradient, epsilon_gradient;
    FluidCalculationUtilities::EvaluateGradientInPoint(
        parent_element_geometry, dNdX[0],
        std::tie(k_gradient, TURBULENT_KINETIC_ENERGY),
        std::tie(epsilon_gradient, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));

    const double cross_diffusion =
        KOmegaSSTElementData::CalculateCrossDiffusionTerm<TDim>(
            sigma_omega_2, omega, k_gradient, epsilon_gradient);

    const double f1 = KOmegaSSTElementData::CalculateF1(
        tke, omega, nu, wall_distance, beta_star, cross_diffusion, sigma_omega_2);

    mBlendedSigmaOmega =
        KOmegaSSTElementData::CalculateBlendedPhi(sigma_omega_1, sigma_omega_2, f1);

    KRATOS_CATCH("");
}

template<unsigned int TDim>
bool OmegaKBasedWallConditionData<TDim>::IsWallFluxComputable() const
{
    return true;
}

template<unsigned int TDim>
double OmegaKBasedWallConditionData<TDim>::CalculateWallFlux(
    const Vector& rShapeFunctions)
{
    using namespace RansCalculationUtilities;

    auto& cl_parameters = this->GetConstitutiveLawParameters();
    cl_parameters.SetShapeFunctionsValues(rShapeFunctions);

    double nu, nu_t, tke;

    this->GetConstitutiveLaw().CalculateValue(cl_parameters, EFFECTIVE_VISCOSITY, nu);
    nu /= mDensity;

    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rShapeFunctions,
        std::tie(nu_t, TURBULENT_VISCOSITY),
        std::tie(tke, TURBULENT_KINETIC_ENERGY));

    const double u_tau = mCmu25 * std::sqrt(std::max(tke, 0.0));

    return (nu + mBlendedSigmaOmega * nu_t) * std::pow(u_tau, 3) /
           (mKappa * std::pow(mCmu25 * mYPlus * nu, 2));
}

// template instantiations
template class OmegaKBasedWallConditionData<2>;
template class OmegaKBasedWallConditionData<3>;

} // namespace KOmegaSSTWallConditionData

} // namespace Kratos