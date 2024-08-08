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
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"

// Application includes
#include "custom_conditions/data_containers/k_omega/condition_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "omega_vis_log_law_based_wall_condition_data.h"

namespace Kratos
{
namespace KOmegaWallConditionData
{
const Variable<double>& OmegaVisLogBasedWallConditionData::GetScalarVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

void OmegaVisLogBasedWallConditionData::Check(
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const auto& r_geometry = rCondition.GetGeometry();
    const int number_of_nodes = r_geometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}

void OmegaVisLogBasedWallConditionData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mKappa = rCurrentProcessInfo[VON_KARMAN];
    mCmu25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);

    const auto& r_condition_geometry = this->GetGeometry();

    // get parent element
    const auto& r_parent_element = r_condition_geometry.GetValue(NEIGHBOUR_ELEMENTS)[0];
    const auto& r_parent_element_geometry = r_parent_element.GetGeometry();

    // get fluid properties from parent element
    const auto& r_elem_properties = this->GetElementProperties();
    const double rho = r_elem_properties[DENSITY];
    const double nu = r_elem_properties[DYNAMIC_VISCOSITY] / rho;
    const array_1d<double, 3>& normal = r_condition_geometry.GetValue(NORMAL);

    // get surface properties from condition
    const auto& r_cond_properties = this->GetConditionProperties();
    const double beta = r_cond_properties.GetValue(WALL_SMOOTHNESS_BETA);

    // Get Shape function data
    RansCalculationUtilities::CalculateGeometryData(
        r_parent_element_geometry,
        GeometryData::IntegrationMethod::GI_GAUSS_1, mGaussPointWeights, mShapeFunctions, mShapeFunctionDerivatives);

    mWallDistance = r_condition_geometry.GetValue(DISTANCE);

    array_1d<double, 3> wall_velocity, fluid_velocity, mesh_velocity;
    FluidCalculationUtilities::EvaluateInPoint(
        r_parent_element_geometry, row(mShapeFunctions, 0),
        std::tie(fluid_velocity, VELOCITY),
        std::tie(mesh_velocity, MESH_VELOCITY));

    noalias(wall_velocity) = fluid_velocity - mesh_velocity;
    const double wall_velocity_magnitude = std::sqrt(std::pow(norm_2(wall_velocity), 2) - std::pow(inner_prod(wall_velocity, normal), 2));

    double mUTau, y_plus;
    RansCalculationUtilities::CalculateYPlusAndUtau(
            y_plus, mUTau, wall_velocity_magnitude, mWallDistance, nu, mKappa, beta);

    const double u_tau_vis = wall_velocity_magnitude / y_plus;
    const double u_tau_log = wall_velocity_magnitude / (std::log(y_plus) / mKappa + beta);
    mUTau = std::pow(std::pow(u_tau_vis, 4) + std::pow(u_tau_log, 4), 0.25);

    const double mOmegaViscous = 6.0 * nu / (0.0075 * mWallDistance * mWallDistance);
    const double mOmegaLog = mUTau / (mCmu25 * mKappa * mWallDistance);
    mOmegaBlended = std::sqrt(std::pow(mOmegaViscous, 2) + std::pow(mOmegaLog, 2));

    KRATOS_CATCH("");
}

double OmegaVisLogBasedWallConditionData::CalculateWallFlux(
    const Vector& rShapeFunctions,
    const ScalarWallFluxConditionData::Parameters& rParameters)
{
    const double omega_vis_contribution = 12.0 * rParameters.mKinematicViscosity / (0.0075 * std::pow(mWallDistance, 3));
    const double omega_log_contribution = mUTau / (mCmu25 * mKappa * std::pow(mWallDistance, 2));

    return (rParameters.mKinematicViscosity) *
           (mOmegaViscous * omega_vis_contribution + mOmegaLog * omega_log_contribution) /
           mOmegaBlended;
}

} // namespace KOmegaWallConditionData

} // namespace Kratos