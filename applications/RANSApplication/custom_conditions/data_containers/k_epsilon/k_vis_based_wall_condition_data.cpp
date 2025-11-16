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
#include "k_vis_based_wall_condition_data.h"

namespace Kratos
{
namespace KEpsilonWallConditionData
{
const Variable<double>& KVisBasedWallConditionData::GetScalarVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

void KVisBasedWallConditionData::Check(
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const auto& r_geometry = rCondition.GetGeometry();
    const int number_of_nodes = r_geometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_KINETIC_ENERGY, r_node);
    }

    KRATOS_CATCH("");
}

void KVisBasedWallConditionData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mWallDistance = this->GetGeometry().GetValue(DISTANCE);

    const auto& r_element_properties = this->GetElementProperties();
    mKinematicViscosity = r_element_properties[DYNAMIC_VISCOSITY] / r_element_properties[DENSITY];

    Vector ws;
    Matrix ns;
    Geometry<Node>::ShapeFunctionsGradientsType dn_dxs;

    // get parent element
    auto& r_parent_element = this->GetGeometry().GetValue(NEIGHBOUR_ELEMENTS)[0];

    // Get Shape function data
    RansCalculationUtilities::CalculateGeometryData(
        r_parent_element.GetGeometry(),
        GeometryData::IntegrationMethod::GI_GAUSS_1, ws, ns, dn_dxs);

    array_1d<double, 3> fluid_velocity, mesh_velocity;
    FluidCalculationUtilities::EvaluateInPoint(
        r_parent_element.GetGeometry(), row(ns, 0),
        std::tie(fluid_velocity, VELOCITY),
        std::tie(mesh_velocity, MESH_VELOCITY));

    mWallVelocityMagnitude = norm_2(fluid_velocity - mesh_velocity);

    KRATOS_CATCH("");
}

double KVisBasedWallConditionData::CalculateWallFlux(
    const Vector& rShapeFunctions,
    const ScalarWallFluxConditionData::Parameters& rParameters)
{
    KRATOS_TRY

    const double y_plus = std::sqrt(mWallVelocityMagnitude * mWallDistance / mKinematicViscosity);

    if (y_plus > 1e-12) {
        const double u_tau = y_plus * mKinematicViscosity / mWallDistance;
        const double u_tau_dy = -0.5 * std::sqrt(mWallVelocityMagnitude * mKinematicViscosity) * std::pow(mWallDistance, -1.5);
        const double y_plus_dy = 0.5 * std::sqrt(mWallVelocityMagnitude / (mWallDistance * mKinematicViscosity));

        const double C = 11.0;
        const double Ceps2_sqr = 1.9 * 1.9;

        const double Cf = (1 / std::pow(y_plus + C, 2) + 2 * y_plus / std::pow(C, 3) - 1 / std::pow(C, 2));
        const double Cf_dy = (-2 / std::pow(y_plus + C, 3) + 2/std::pow(C, 3)) * y_plus_dy;

        const double k_plus = 2400.0 / Ceps2_sqr * Cf;
        const double k_plus_dy = 2400.0 / Ceps2_sqr * Cf_dy;

        const double g_n = -(k_plus_dy * std::pow(u_tau, 2) + k_plus * 2 * u_tau * u_tau_dy);

        return (rParameters.mKinematicViscosity) * g_n;
    } else {
        return 0.0;
    }

    KRATOS_CATCH("");
}

} // namespace KEpsilonWallConditionData

} // namespace Kratos