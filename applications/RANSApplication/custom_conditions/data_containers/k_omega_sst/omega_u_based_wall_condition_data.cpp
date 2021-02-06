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
#include "omega_u_based_wall_condition_data.h"

namespace Kratos
{
namespace KOmegaSSTWallConditionData
{
template<unsigned int TDim>
const Variable<double>& OmegaUBasedWallConditionData<TDim>::GetScalarVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

template<unsigned int TDim>
void OmegaUBasedWallConditionData<TDim>::Check(
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
    KRATOS_ERROR_IF_NOT(r_properties.Has(WALL_SMOOTHNESS_BETA))
        << "WALL_SMOOTHNESS_BETA is not found in condition properties [ Condition.Id() = "
        << rCondition.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    KRATOS_ERROR_IF_NOT(rCondition.Has(TURBULENT_VISCOSITY))
        << "TURBULENT_VISCOSITY value is not set at in condition with id " << rCondition.Id() << "\n";
    KRATOS_ERROR_IF_NOT(rCondition.Has(RANS_Y_PLUS))
        << "RANS_Y_PLUS value is not set at in condition with id " << rCondition.Id() << "\n";

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}
template<unsigned int TDim>
GeometryData::IntegrationMethod OmegaUBasedWallConditionData<TDim>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_1;
}

template<unsigned int TDim>
void OmegaUBasedWallConditionData<TDim>::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    mBetaStar = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    mSigmaOmega1 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1];
    mSigmaOmega2 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2];
    mKappa = rCurrentProcessInfo[VON_KARMAN];
    mCmu25 = std::pow(mBetaStar, 0.25);

    mDensity = this->GetElementProperties()[DENSITY];

    const auto& r_properties = this->GetConditionProperties();
    mBeta = r_properties[WALL_SMOOTHNESS_BETA];
    const double y_plus_limit = r_properties[RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT];

    mInvKappa = 1.0 / mKappa;
    const auto& r_geometry = this->GetGeometry();
    mYPlus = std::max(r_geometry.GetValue(RANS_Y_PLUS), y_plus_limit);
    mTurbulentViscosity = r_geometry.GetValue(TURBULENT_VISCOSITY);

    auto normal = r_geometry.GetValue(NORMAL);
    normal /= norm_2(normal);
    mWallHeight = inner_prod(r_geometry.Center() - r_geometry.GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry().Center(), normal);

    KRATOS_CATCH("");
}

template<unsigned int TDim>
bool OmegaUBasedWallConditionData<TDim>::IsWallFluxComputable() const
{
    return true;
}

template<unsigned int TDim>
double OmegaUBasedWallConditionData<TDim>::CalculateWallFlux(
    const Vector& rShapeFunctions)
{
    using namespace RansCalculationUtilities;

    auto& cl_parameters = this->GetConstitutiveLawParameters();
    cl_parameters.SetShapeFunctionsValues(rShapeFunctions);

    double nu, tke, omega;
    array_1d<double, 3> velocity;

    this->GetConstitutiveLaw().CalculateValue(cl_parameters, EFFECTIVE_VISCOSITY, nu);
    nu /= mDensity;

    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rShapeFunctions,
        std::tie(velocity, VELOCITY),
        std::tie(tke, TURBULENT_KINETIC_ENERGY),
        std::tie(omega, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));

    const double velocity_magnitude = norm_2(velocity);

    const double u_tau = velocity_magnitude / (mInvKappa * std::log(mYPlus) + mBeta);

    const double f1 = KOmegaSSTElementData::CalculateF1(
        tke, omega, nu, mWallHeight, mBetaStar, 1e-10, mSigmaOmega2);

    const double blended_sigma_omega =
        KOmegaSSTElementData::CalculateBlendedPhi(mSigmaOmega1, mSigmaOmega2, f1);

    return (nu + blended_sigma_omega * mTurbulentViscosity) *
           std::pow(u_tau, 3) / (mKappa * std::pow(mCmu25 * mYPlus * nu, 2));
}

// template instantiations
template class OmegaUBasedWallConditionData<2>;
template class OmegaUBasedWallConditionData<3>;

} // namespace KOmegaSSTWallConditionData

} // namespace Kratos