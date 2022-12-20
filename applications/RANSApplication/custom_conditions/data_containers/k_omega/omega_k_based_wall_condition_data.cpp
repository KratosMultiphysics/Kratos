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
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "omega_k_based_wall_condition_data.h"

namespace Kratos
{
namespace KOmegaWallConditionData
{
const Variable<double>& OmegaKBasedWallConditionData::GetScalarVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

void OmegaKBasedWallConditionData::Check(
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = rCondition.GetGeometry();
    const auto& r_properties = rCondition.GetProperties();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(VON_KARMAN))
        << "VON_KARMAN is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(r_properties.Has(RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT))
        << "RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT is not found in condition properties [ Condition.Id() = "
        << rCondition.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    for (const auto& r_node : r_geometry) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}
GeometryData::IntegrationMethod OmegaKBasedWallConditionData::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
}

void OmegaKBasedWallConditionData::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mOmegaSigma = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA];
    mCmu25 = std::pow(rCurrentProcessInfo[TURBULENCE_RANS_C_MU], 0.25);
    mKappa = rCurrentProcessInfo[VON_KARMAN];

    KRATOS_ERROR_IF(!(this->GetGeometry().Has(RANS_Y_PLUS)))
        << "RANS_Y_PLUS value is not set at " << this->GetGeometry() << "\n";

    mDensity = this->GetElementProperties()[DENSITY];

    const double y_plus_limit = this->GetConditionProperties()[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];

    mYPlus = std::max(this->GetGeometry().GetValue(RANS_Y_PLUS), y_plus_limit);

    KRATOS_CATCH("");
}

bool OmegaKBasedWallConditionData::IsWallFluxComputable() const
{
    return true;
}

double OmegaKBasedWallConditionData::CalculateWallFlux(
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

    return (nu + mOmegaSigma * nu_t) * std::pow(u_tau, 3) /
           (mKappa * std::pow(mCmu25 * mYPlus * nu, 2));
}

} // namespace KOmegaWallConditionData

} // namespace Kratos