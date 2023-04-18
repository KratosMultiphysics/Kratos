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
#include "custom_utilities/fluid_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "epsilon_u_based_wall_condition_data.h"

namespace Kratos
{
namespace KEpsilonWallConditionData
{
const Variable<double>& EpsilonUBasedWallConditionData::GetScalarVariable()
{
    return TURBULENT_ENERGY_DISSIPATION_RATE;
}

void EpsilonUBasedWallConditionData::Check(
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = rCondition.GetGeometry();
    const auto& r_properties = rCondition.GetProperties();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA))
        << "TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA is not found in process info.\n";

    KRATOS_ERROR_IF_NOT(r_properties.Has(WALL_SMOOTHNESS_BETA))
        << "WALL_SMOOTHNESS_BETA is not found in condition properties [ Condition.Id() = "
        << rCondition.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    for (const auto& r_node: r_geometry) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}

void EpsilonUBasedWallConditionData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mEpsilonSigma = rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    mBeta = this->GetConditionProperties()[WALL_SMOOTHNESS_BETA];

    KRATOS_CATCH("");
}

double EpsilonUBasedWallConditionData::CalculateWallFlux(
    const Vector& rShapeFunctions,
    const ScalarWallFluxConditionData::Parameters& rParameters)
{
    array_1d<double, 3> velocity;
    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rShapeFunctions,
        std::tie(velocity, VELOCITY));

    const double u_tau = norm_2(velocity) / ((1.0 / rParameters.mKappa) * std::log(rParameters.mYPlus) + mBeta);

    return (rParameters.mKinematicViscosity + rParameters.mWallTurbulentViscosity / mEpsilonSigma) * std::pow(u_tau, 5) /
           (rParameters.mKappa * std::pow(rParameters.mYPlus * rParameters.mKinematicViscosity, 2));
}

} // namespace KEpsilonWallConditionData

} // namespace Kratos