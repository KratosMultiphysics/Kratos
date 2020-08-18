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
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/variables.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
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

const Variable<double>& EpsilonUBasedWallConditionData::GetScalarRateVariable()
{
    return TURBULENT_ENERGY_DISSIPATION_RATE_2;
}

void EpsilonUBasedWallConditionData::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const int number_of_nodes = rGeometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const auto& r_node = rGeometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE_2, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}
GeometryData::IntegrationMethod EpsilonUBasedWallConditionData::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_1;
}

void EpsilonUBasedWallConditionData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mEpsilonSigma = rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    mKappa = rCurrentProcessInfo[WALL_VON_KARMAN];
    mBeta = rCurrentProcessInfo[WALL_SMOOTHNESS_BETA];
    mInvKappa = 1.0 / mKappa;

    KRATOS_ERROR_IF(!(this->GetGeometry().Has(RANS_Y_PLUS)))
        << "RANS_Y_PLUS value is not set at " << this->GetGeometry() << "\n";

    const double y_plus_limit = rCurrentProcessInfo[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];
    mYPlus = std::max(this->GetGeometry().GetValue(RANS_Y_PLUS), y_plus_limit);

    KRATOS_CATCH("");
}

bool EpsilonUBasedWallConditionData::IsWallFluxComputable() const
{
    return true;
}

double EpsilonUBasedWallConditionData::CalculateWallFlux(
    const Vector& rShapeFunctions) const
{
    using namespace RansCalculationUtilities;

    double nu, nu_t;
    array_1d<double, 3> velocity;

    EvaluateInPoint(this->GetGeometry(), rShapeFunctions,
                    std::tie(nu, KINEMATIC_VISCOSITY),
                    std::tie(nu_t, TURBULENT_VISCOSITY),
                    std::tie(velocity, VELOCITY));

    const double velocity_magnitude = norm_2(velocity);

    const double u_tau = velocity_magnitude / (mInvKappa * std::log(mYPlus) + mBeta);

    return (nu + nu_t / mEpsilonSigma) * std::pow(u_tau, 5) /
           (mKappa * std::pow(mYPlus * nu, 2));
}

} // namespace KEpsilonWallConditionData

} // namespace Kratos