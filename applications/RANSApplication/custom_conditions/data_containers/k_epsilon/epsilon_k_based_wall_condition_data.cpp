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

// Application includes
#include "custom_utilities/fluid_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "epsilon_k_based_wall_condition_data.h"

namespace Kratos
{
namespace KEpsilonWallConditionData
{
const Variable<double>& EpsilonKBasedWallConditionData::GetScalarVariable()
{
    return TURBULENT_ENERGY_DISSIPATION_RATE;
}

void EpsilonKBasedWallConditionData::Check(
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA))
        << "TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";

    for (const auto& r_node : rCondition.GetGeometry()) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}

void EpsilonKBasedWallConditionData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mEpsilonSigma = rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    mCmu25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);

    KRATOS_CATCH("");
}

double EpsilonKBasedWallConditionData::CalculateWallFlux(
    const Vector& rShapeFunctions,
    const ScalarWallFluxConditionData::Parameters& rParameters)
{
    double tke;
    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rShapeFunctions,
        std::tie(tke, TURBULENT_KINETIC_ENERGY));

    const double u_tau = mCmu25 * std::sqrt(std::max(tke, 0.0));

    return (rParameters.mKinematicViscosity + rParameters.mWallTurbulentViscosity / mEpsilonSigma) * std::pow(u_tau, 5) /
           (rParameters.mKappa * std::pow(rParameters.mYPlus * rParameters.mKinematicViscosity, 2));
}

} // namespace KEpsilonWallConditionData

} // namespace Kratos