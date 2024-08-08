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
#include "includes/variables.h"

// Application includes
#include "custom_elements/data_containers/k_epsilon/element_data_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "element_data_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "omega_element_data.h"

namespace Kratos
{
namespace KOmegaSSTElementData
{
template <unsigned int TDim>
const Variable<double>& OmegaElementData<TDim>::GetScalarVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim>
void OmegaElementData<TDim>::Check(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = rElement.GetGeometry();
    const auto& r_properties = rElement.GetProperties();
    const int number_of_nodes = r_geometry.PointsNumber();

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_BETA_1))
        << "TURBULENCE_RANS_BETA_1 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_BETA_2))
        << "TURBULENCE_RANS_BETA_2 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2))
        << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2 is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_C_MU))
        << "TURBULENCE_RANS_C_MU is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(VON_KARMAN))
        << "VON_KARMAN is not found in process info.\n";
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(TURBULENCE_RANS_A1))
        << "TURBULENCE_RANS_A1 is not found in process info.\n";

    KRATOS_ERROR_IF_NOT(r_properties.Has(DYNAMIC_VISCOSITY))
        << "DYNAMIC_VISCOSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";
    KRATOS_ERROR_IF_NOT(r_properties.Has(DENSITY))
        << "DENSITY is not found in element properties [ Element.Id() = "
        << rElement.Id() << ", Properties.Id() = " << r_properties.Id() << " ].\n";

    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = r_geometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);

        KRATOS_ERROR_IF_NOT(r_node.Has(TURBULENT_KINETIC_ENERGY))
            << "TURBULENT_KINETIC_ENERGY is not found in non-historical data "
               "value container of node with id "
            << r_node.Id() << ".\n";
        KRATOS_ERROR_IF_NOT(r_node.Has(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE))
            << "TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE is not found in non-historical data value container of node with id "
            << r_node.Id() << ".\n";

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim>
void OmegaElementData<TDim>::CalculateConstants(
    const ProcessInfo& rCurrentProcessInfo)
{
    mBeta1 = rCurrentProcessInfo[TURBULENCE_RANS_BETA_1];
    mBeta2 = rCurrentProcessInfo[TURBULENCE_RANS_BETA_2];
    mSigmaOmega1 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1];
    mSigmaOmega2 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2];
    mBetaStar = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    mKappa = rCurrentProcessInfo[VON_KARMAN];
    mA1 = rCurrentProcessInfo[TURBULENCE_RANS_A1];

    const auto& r_properties = this->GetProperties();
    mDensity = r_properties[DENSITY];
    mKinematicViscosity = r_properties[DYNAMIC_VISCOSITY] / mDensity;
}

template <unsigned int TDim>
void OmegaElementData<TDim>::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    const auto& r_geometry = this->GetGeometry();

    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rShapeFunctions, Step,
        std::tie(mTurbulentKineticEnergy, TURBULENT_KINETIC_ENERGY),
        std::tie(mTurbulentSpecificEnergyDissipationRate, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(mWallDistance, DISTANCE),
        std::tie(mEffectiveVelocity, VELOCITY));

    KRATOS_ERROR_IF(mWallDistance < 0.0) << "Wall distance is negative at " << r_geometry;

    FluidCalculationUtilities::EvaluateGradientInPoint(
        this->GetGeometry(), rShapeFunctionDerivatives,
        std::tie(mTurbulentKineticEnergyGradient, TURBULENT_KINETIC_ENERGY),
        std::tie(mTurbulentSpecificEnergyDissipationRateGradient, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(mVelocityGradient, VELOCITY));

    mCrossDiffusion = KOmegaSSTElementData::CalculateCrossDiffusionTerm<TDim>(
        mSigmaOmega2, mTurbulentSpecificEnergyDissipationRate,
        mTurbulentKineticEnergyGradient, mTurbulentSpecificEnergyDissipationRateGradient);

    mF1 = KOmegaSSTElementData::CalculateF1(
        mTurbulentKineticEnergy, mTurbulentSpecificEnergyDissipationRate,
        mKinematicViscosity, mWallDistance, mBetaStar, mCrossDiffusion, mSigmaOmega2);

    mBlendedSigmaOmega =
        KOmegaSSTElementData::CalculateBlendedPhi(mSigmaOmega1, mSigmaOmega2, mF1);
    mBlendedBeta = KOmegaSSTElementData::CalculateBlendedPhi(mBeta1, mBeta2, mF1);

    const double gamma_1 = KOmegaSSTElementData::CalculateGamma(
        mBeta1, mBetaStar, mSigmaOmega1, mKappa);
    const double gamma_2 = KOmegaSSTElementData::CalculateGamma(
        mBeta2, mBetaStar, mSigmaOmega2, mKappa);

    mBlendedGamma = KOmegaSSTElementData::CalculateBlendedPhi(gamma_1, gamma_2, mF1);

    mVelocityDivergence = CalculateMatrixTrace<TDim>(mVelocityGradient);

    double tke_old, omega_old;
    FluidCalculationUtilities::EvaluateNonHistoricalInPoint(
        this->GetGeometry(), rShapeFunctions,
        std::tie(tke_old, TURBULENT_KINETIC_ENERGY),
        std::tie(omega_old, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE));

    const double f_2 = KOmegaSSTElementData::CalculateF2(tke_old, omega_old, mKinematicViscosity, mWallDistance, mBetaStar);

    const BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient =
        (mVelocityGradient + trans(mVelocityGradient)) * 0.5;

    const double t = norm_frobenius(symmetric_velocity_gradient) * 1.414;

    mTurbulentKinematicViscosity = KOmegaSSTElementData::CalculateTurbulentKinematicViscosity(tke_old, omega_old, t, f_2, mA1);
    mTurbulentKinematicViscosity = std::max(mTurbulentKinematicViscosity, 1e-12);

    mEffectiveVelocity -= mTurbulentKineticEnergyGradient * (2.0 * (1-mF1) * mSigmaOmega2 / omega_old);

    mEffectiveKinematicViscosity = mKinematicViscosity + mTurbulentKinematicViscosity * mBlendedSigmaOmega;

    const double omega = std::max(mTurbulentSpecificEnergyDissipationRate, 1e-12);
    mReactionTerm = mBlendedBeta * omega;
    mReactionTerm -= (1.0 - mF1) * mCrossDiffusion / omega;
    mReactionTerm += mBlendedGamma * 2.0 * mVelocityDivergence / 3.0;
    mReactionTerm = std::max(mReactionTerm, 0.0);

    mSourceTerm = KEpsilonElementData::CalculateProductionTerm<TDim>(mVelocityGradient, 1.0) * (mBlendedGamma);

    KRATOS_CATCH("");
}

// template instantiations

template class OmegaElementData<2>;
template class OmegaElementData<3>;

} // namespace KOmegaSSTElementData

} // namespace Kratos