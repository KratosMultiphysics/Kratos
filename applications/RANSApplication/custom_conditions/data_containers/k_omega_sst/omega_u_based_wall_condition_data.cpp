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

// Application includes
#include "custom_elements/data_containers/k_omega_sst/element_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "omega_u_based_wall_condition_data.h"

namespace Kratos
{
namespace KOmegaSSTWallConditionData
{
const Variable<double>& OmegaUBasedWallConditionData::GetScalarVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

void OmegaUBasedWallConditionData::Check(
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
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(WALL_CORRECTION_FACTOR))
        << "WALL_CORRECTION_FACTOR is not found in process info.\n";

    KRATOS_ERROR_IF_NOT(r_properties.Has(WALL_SMOOTHNESS_BETA))
        << "WALL_SMOOTHNESS_BETA is not found in condition properties [ Condition.Id() = "
        << rCondition.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    KRATOS_ERROR_IF_NOT(r_geometry.Has(DISTANCE))
        << "DISTANCE is not found in condition with id " << rCondition.Id() << ".\n";

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}

void OmegaUBasedWallConditionData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mBetaStar = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    mSigmaOmega1 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1];
    mSigmaOmega2 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2];
    mWallCorrectionFactor = rCurrentProcessInfo[WALL_CORRECTION_FACTOR];
    mCmu25 = std::pow(mBetaStar, 0.25);
    mBeta = this->GetConditionProperties()[WALL_SMOOTHNESS_BETA];

    mWallHeight = this->GetGeometry().GetValue(DISTANCE);

    KRATOS_CATCH("");
}

double OmegaUBasedWallConditionData::CalculateWallFlux(
    const Vector& rShapeFunctions,
    const ScalarWallFluxConditionData::Parameters& rParameters)
{
    double tke, omega;
    array_1d<double, 3> velocity, fluid_velocity, mesh_velocity;

    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rShapeFunctions,
        std::tie(velocity, VELOCITY),
        std::tie(tke, TURBULENT_KINETIC_ENERGY),
        std::tie(omega, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(fluid_velocity, VELOCITY),
        std::tie(mesh_velocity, MESH_VELOCITY));

    noalias(velocity) = fluid_velocity - mesh_velocity;

    const double u_tau = norm_2(velocity) / ((1.0 / rParameters.mKappa) * std::log(rParameters.mYPlus) + mBeta);

    const double f1 = KOmegaSSTElementData::CalculateF1(
        tke, omega, rParameters.mKinematicViscosity, mWallHeight, mBetaStar, 1e-10, mSigmaOmega2);

    const double blended_sigma_omega =
        KOmegaSSTElementData::CalculateBlendedPhi(mSigmaOmega1, mSigmaOmega2, f1);

    return (rParameters.mKinematicViscosity + mWallCorrectionFactor * blended_sigma_omega * rParameters.mWallTurbulentViscosity) * std::pow(u_tau, 3) /
           (rParameters.mKappa * std::pow(mCmu25 * rParameters.mYPlus * rParameters.mKinematicViscosity, 2));
}


} // namespace KOmegaSSTWallConditionData

} // namespace Kratos