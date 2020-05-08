//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

// System includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/variables.h"

// Application includes
#include "custom_elements/element_data/evm_k_epsilon/evm_k_epsilon_element_data_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "evm_k_omega_sst_element_data_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "evm_k_omega_sst_omega_element_data.h"

namespace Kratos
{
namespace EvmKOmegaSSTElementDataUtilities
{
template <unsigned int TDim>
const Variable<double>& OmegaElementData<TDim>::GetScalarVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim>
const Variable<double>& OmegaElementData<TDim>::GetScalarRateVariable()
{
    return TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2;
}

template <unsigned int TDim>
const Variable<double>& OmegaElementData<TDim>::GetScalarRelaxedRateVariable()
{
    return RANS_AUXILIARY_VARIABLE_2;
}

template <unsigned int TDim>
void OmegaElementData<TDim>::Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const int number_of_nodes = rGeometry.PointsNumber();

    for (int i_node = 0; i_node < number_of_nodes; ++i_node)
    {
        const NodeType& r_node = rGeometry[i_node];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(KINEMATIC_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_VISCOSITY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(
            TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(
            TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_2, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}
template <unsigned int TDim>
GeometryData::IntegrationMethod OmegaElementData<TDim>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_2;
}

template <unsigned int TDim>
void OmegaElementData<TDim>::CalculateConstants(const ProcessInfo& rCurrentProcessInfo)
{
    mBeta1 = rCurrentProcessInfo[TURBULENCE_RANS_BETA_1];
    mBeta2 = rCurrentProcessInfo[TURBULENCE_RANS_BETA_2];
    mSigmaOmega1 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1];
    mSigmaOmega2 = rCurrentProcessInfo[TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2];
    mBetaStar = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    mKappa = rCurrentProcessInfo[WALL_VON_KARMAN];
}

template <unsigned int TDim>
void OmegaElementData<TDim>::CalculateGaussPointData(const Vector& rShapeFunctions,
                                                     const Matrix& rShapeFunctionDerivatives,
                                                     const int Step)
{
    KRATOS_TRY

    using namespace RansCalculationUtilities;

    const GeometryType& r_geometry = this->GetGeometry();

    mTurbulentKineticEnergy = std::max(
        EvaluateInPoint(r_geometry, TURBULENT_KINETIC_ENERGY, rShapeFunctions), 1e-12);
    mTurbulentSpecificEnergyDissipationRate = std::max(
        EvaluateInPoint(r_geometry, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, rShapeFunctions),
        1e-12);
    mKinematicViscosity = EvaluateInPoint(r_geometry, KINEMATIC_VISCOSITY, rShapeFunctions);
    mWallDistance = EvaluateInPoint(r_geometry, DISTANCE, rShapeFunctions);
    mTurbulentKinematicViscosity =
        EvaluateInPoint(r_geometry, TURBULENT_VISCOSITY, rShapeFunctions);

    CalculateGradient(mTurbulentKineticEnergyGradient, r_geometry,
                      TURBULENT_KINETIC_ENERGY, rShapeFunctionDerivatives, Step);

    CalculateGradient(mTurbulentSpecificEnergyDissipationRateGradient,
                      r_geometry, TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE,
                      rShapeFunctionDerivatives, Step);

    mCrossDiffusion = EvmKOmegaSSTElementDataUtilities::CalculateCrossDiffusionTerm(
        mSigmaOmega2, mTurbulentSpecificEnergyDissipationRate,
        mTurbulentKineticEnergyGradient, mTurbulentSpecificEnergyDissipationRateGradient);

    mF1 = EvmKOmegaSSTElementDataUtilities::CalculateF1(
        mTurbulentKineticEnergy, mTurbulentSpecificEnergyDissipationRate,
        mKinematicViscosity, mWallDistance, mBetaStar, mCrossDiffusion, mSigmaOmega2);

    mBlendedSigmaOmega = EvmKOmegaSSTElementDataUtilities::CalculateBlendedPhi(
        mSigmaOmega1, mSigmaOmega2, mF1);
    mBlendedBeta =
        EvmKOmegaSSTElementDataUtilities::CalculateBlendedPhi(mBeta1, mBeta2, mF1);

    const double gamma_1 = EvmKOmegaSSTElementDataUtilities::CalculateGamma(
        mBeta1, mBetaStar, mSigmaOmega1, mKappa);
    const double gamma_2 = EvmKOmegaSSTElementDataUtilities::CalculateGamma(
        mBeta2, mBetaStar, mSigmaOmega2, mKappa);

    mBlendedGamma =
        EvmKOmegaSSTElementDataUtilities::CalculateBlendedPhi(gamma_1, gamma_2, mF1);

    mVelocityDivergence = GetDivergence(r_geometry, VELOCITY, rShapeFunctionDerivatives);

    CalculateGradient<TDim>(mVelocityGradient, r_geometry, VELOCITY,
                            rShapeFunctionDerivatives, Step);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
array_1d<double, 3> OmegaElementData<TDim>::CalculateEffectiveVelocity(
    const Vector& rShapeFunctions, const Matrix& rShapeFunctionDerivatives) const
{
    const array_1d<double, 3>& r_velocity = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), VELOCITY, rShapeFunctions);

    return r_velocity;
}

template <unsigned int TDim>
double OmegaElementData<TDim>::CalculateEffectiveKinematicViscosity(
    const Vector& rShapeFunctions, const Matrix& rShapeFunctionDerivatives) const
{
    return mKinematicViscosity + mTurbulentKinematicViscosity * mBlendedSigmaOmega;
}

template <unsigned int TDim>
double OmegaElementData<TDim>::CalculateReactionTerm(const Vector& rShapeFunctions,
                                                     const Matrix& rShapeFunctionDerivatives) const
{
    double value = mBlendedBeta * mTurbulentSpecificEnergyDissipationRate;
    value -= (1.0 - mF1) * mCrossDiffusion / mTurbulentSpecificEnergyDissipationRate;
    value += mBlendedGamma * 2.0 * mVelocityDivergence / 3.0;
    return std::max(value, 0.0);
}

template <unsigned int TDim>
double OmegaElementData<TDim>::CalculateSourceTerm(const Vector& rShapeFunctions,
                                                   const Matrix& rShapeFunctionDerivatives) const
{
    double production = 0.0;

    production = EvmKEpsilonElementDataUtilities::CalculateSourceTerm<TDim>(
        mVelocityGradient, mTurbulentKinematicViscosity);

    production *= (mBlendedGamma / mTurbulentKinematicViscosity);

    return production;
}

// template instantiations

template class OmegaElementData<2>;
template class OmegaElementData<3>;

} // namespace EvmKOmegaSSTElementDataUtilities

} // namespace Kratos