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
    const GeometryType& rGeometry,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const int number_of_nodes = rGeometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node) {
        const auto& r_node = rGeometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE, r_node);

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

    mDensity = this->GetProperties().GetValue(DENSITY);
}

template <unsigned int TDim>
void OmegaElementData<TDim>::CalculateGaussPointData(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const int Step)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    auto& cl_parameters = this->GetConstitutiveLawParameters();
    cl_parameters.SetShapeFunctionsValues(rShapeFunctions);

    this->GetConstitutiveLaw().CalculateValue(cl_parameters, EFFECTIVE_VISCOSITY, mKinematicViscosity);
    mKinematicViscosity /= mDensity;

    const auto& r_geometry = this->GetGeometry();

    FluidCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), rShapeFunctions, Step,
        std::tie(mTurbulentKineticEnergy, TURBULENT_KINETIC_ENERGY),
        std::tie(mTurbulentSpecificEnergyDissipationRate, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE),
        std::tie(mTurbulentKinematicViscosity, TURBULENT_VISCOSITY),
        std::tie(mWallDistance, DISTANCE),
        std::tie(mEffectiveVelocity, VELOCITY));

    KRATOS_ERROR_IF(mWallDistance < 0.0) << "Wall distance is negative at " << r_geometry;

    CalculateGradient(mTurbulentKineticEnergyGradient, r_geometry,
                      TURBULENT_KINETIC_ENERGY, rShapeFunctionDerivatives, Step);

    CalculateGradient(mTurbulentSpecificEnergyDissipationRateGradient,
                      r_geometry, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE,
                      rShapeFunctionDerivatives, Step);

    mCrossDiffusion = KOmegaSSTElementData::CalculateCrossDiffusionTerm(
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

    mVelocityDivergence = GetDivergence(r_geometry, VELOCITY, rShapeFunctionDerivatives);

    CalculateGradient<TDim>(mVelocityGradient, r_geometry, VELOCITY,
                            rShapeFunctionDerivatives, Step);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
array_1d<double, 3> OmegaElementData<TDim>::CalculateEffectiveVelocity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    return mEffectiveVelocity;
}

template <unsigned int TDim>
double OmegaElementData<TDim>::CalculateEffectiveKinematicViscosity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    return mKinematicViscosity + mTurbulentKinematicViscosity * mBlendedSigmaOmega;
}

template <unsigned int TDim>
double OmegaElementData<TDim>::CalculateReactionTerm(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    const double omega = std::max(mTurbulentSpecificEnergyDissipationRate, 1e-12);
    double value = mBlendedBeta * omega;
    value -= (1.0 - mF1) * mCrossDiffusion / omega;
    value += mBlendedGamma * 2.0 * mVelocityDivergence / 3.0;
    return std::max(value, 0.0);
}

template <unsigned int TDim>
double OmegaElementData<TDim>::CalculateSourceTerm(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives) const
{
    double production = 0.0;

    production = KEpsilonElementData::CalculateSourceTerm<TDim>(
        mVelocityGradient, mTurbulentKinematicViscosity);

    production *= (mBlendedGamma / mTurbulentKinematicViscosity);

    return production;
}

// template instantiations

template class OmegaElementData<2>;
template class OmegaElementData<3>;

} // namespace KOmegaSSTElementData

} // namespace Kratos