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
#include "custom_utilities/rans_calculation_utilities.h"
#include "evm_k_epsilon_element_data_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "evm_k_epsilon_epsilon_element_data.h"

namespace Kratos
{
namespace EvmKEpsilonElementDataUtilities
{
template <unsigned int TDim>
const Variable<double>& EpsilonElementData<TDim>::GetScalarVariable()
{
    return TURBULENT_ENERGY_DISSIPATION_RATE;
}

template <unsigned int TDim>
const Variable<double>& EpsilonElementData<TDim>::GetScalarRateVariable()
{
    return TURBULENT_ENERGY_DISSIPATION_RATE_2;
}

template <unsigned int TDim>
const Variable<double>& EpsilonElementData<TDim>::GetScalarRelaxedRateVariable()
{
    return RANS_AUXILIARY_VARIABLE_2;
}

template <unsigned int TDim>
void EpsilonElementData<TDim>::Check(const GeometryType& rGeometry,
                                     const ProcessInfo& rCurrentProcessInfo)
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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_ENERGY_DISSIPATION_RATE_2, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_2, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_ENERGY_DISSIPATION_RATE, r_node);
    }

    KRATOS_CATCH("");
}
template <unsigned int TDim>
GeometryData::IntegrationMethod EpsilonElementData<TDim>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_2;
}

template <unsigned int TDim>
void EpsilonElementData<TDim>::CalculateGaussPointData(const Vector& rShapeFunctions,
                                                       const Matrix& rShapeFunctionDerivatives,
                                                       const ProcessInfo& rCurrentProcessInfo,
                                                       const int Step)
{
    KRATOS_TRY

    mC1 = rCurrentProcessInfo[TURBULENCE_RANS_C1];
    mC2 = rCurrentProcessInfo[TURBULENCE_RANS_C2];

    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    mKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), KINEMATIC_VISCOSITY, rShapeFunctions, Step);
    mTurbulentKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), TURBULENT_VISCOSITY, rShapeFunctions, Step);
    mTurbulentKineticEnergy = RansCalculationUtilities::EvaluateInPoint(
        this->GetGeometry(), TURBULENT_KINETIC_ENERGY, rShapeFunctions, Step);
    mGamma = EvmKEpsilonElementDataUtilities::CalculateGamma(
        c_mu, mTurbulentKineticEnergy, mTurbulentKinematicViscosity);

    mVelocityDivergence = RansCalculationUtilities::GetDivergence(
        this->GetGeometry(), VELOCITY, rShapeFunctionDerivatives);

    RansCalculationUtilities::CalculateGradient<TDim>(
        this->mVelocityGradient, this->GetGeometry(), VELOCITY,
        rShapeFunctionDerivatives, Step);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
double EpsilonElementData<TDim>::CalculateEffectiveKinematicViscosity(
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const double epsilon_sigma =
        rCurrentProcessInfo[TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA];
    return mKinematicViscosity + mTurbulentKinematicViscosity / epsilon_sigma;
}

template <unsigned int TDim>
double EpsilonElementData<TDim>::CalculateReactionTerm(const Vector& rShapeFunctions,
                                                       const Matrix& rShapeFunctionDerivatives,
                                                       const ProcessInfo& rCurrentProcessInfo) const
{
    return std::max(mC2 * mGamma + mC1 * 2.0 * mVelocityDivergence / 3.0, 0.0);
}

template <unsigned int TDim>
double EpsilonElementData<TDim>::CalculateSourceTerm(const Vector& rShapeFunctions,
                                                     const Matrix& rShapeFunctionDerivatives,
                                                     const ProcessInfo& rCurrentProcessInfo) const
{
    double production = 0.0;

    production = EvmKEpsilonElementDataUtilities::CalculateSourceTerm<TDim>(
        mVelocityGradient, mTurbulentKinematicViscosity);

    production *= (mC1 * mGamma);
    return production;
}

// template instantiations

template class EpsilonElementData<2>;
template class EpsilonElementData<3>;

} // namespace EvmKEpsilonElementDataUtilities

} // namespace Kratos