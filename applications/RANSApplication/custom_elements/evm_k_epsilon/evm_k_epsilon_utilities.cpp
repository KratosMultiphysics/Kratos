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
#include <cmath>

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "evm_k_epsilon_utilities.h"

namespace Kratos
{
namespace EvmKepsilonModelUtilities
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
void EpsilonElementData<TDim>::CalculateElementData(const GeometryType& rGeometry,
                                                    const Vector& rShapeFunctions,
                                                    const Matrix& rShapeFunctionDerivatives,
                                                    const ProcessInfo& rCurrentProcessInfo,
                                                    const int Step)
{
    KRATOS_TRY

    mC1 = rCurrentProcessInfo[TURBULENCE_RANS_C1];
    mC2 = rCurrentProcessInfo[TURBULENCE_RANS_C2];

    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];
    mKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        rGeometry, KINEMATIC_VISCOSITY, rShapeFunctions, Step);
    mTurbulentKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        rGeometry, TURBULENT_VISCOSITY, rShapeFunctions, Step);
    mTurbulentKineticEnergy = RansCalculationUtilities::EvaluateInPoint(
        rGeometry, TURBULENT_KINETIC_ENERGY, rShapeFunctions, Step);
    mGamma = EvmKepsilonModelUtilities::CalculateGamma(
        c_mu, 1.0, mTurbulentKineticEnergy, mTurbulentKinematicViscosity);

    mVelocityDivergence = RansCalculationUtilities::GetDivergence(
        rGeometry, VELOCITY, rShapeFunctionDerivatives);

    RansCalculationUtilities::CalculateGradient<TDim>(
        this->mVelocityGradient, rGeometry, VELOCITY, rShapeFunctionDerivatives, Step);

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
    return RansCalculationUtilities::SoftPositive(
        mC2 * mGamma + mC1 * 2.0 * mVelocityDivergence / 3.0);
}

template <unsigned int TDim>
double EpsilonElementData<TDim>::CalculateSourceTerm(const Vector& rShapeFunctions,
                                                     const Matrix& rShapeFunctionDerivatives,
                                                     const ProcessInfo& rCurrentProcessInfo) const
{
    double production = 0.0;

    production = EvmKepsilonModelUtilities::CalculateSourceTerm<TDim>(
        mVelocityGradient, mTurbulentKinematicViscosity, mTurbulentKineticEnergy);

    production *= (mC1 * mGamma);
    return production;
}

template <unsigned int TDim>
const Variable<double>& KElementData<TDim>::GetScalarVariable()
{
    return TURBULENT_KINETIC_ENERGY;
}

template <unsigned int TDim>
const Variable<double>& KElementData<TDim>::GetScalarRateVariable()
{
    return TURBULENT_KINETIC_ENERGY_RATE;
}

template <unsigned int TDim>
const Variable<double>& KElementData<TDim>::GetScalarRelaxedRateVariable()
{
    return RANS_AUXILIARY_VARIABLE_1;
}

template <unsigned int TDim>
void KElementData<TDim>::Check(const GeometryType& rGeometry, const ProcessInfo& rCurrentProcessInfo)
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
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(TURBULENT_KINETIC_ENERGY_RATE, r_node);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(RANS_AUXILIARY_VARIABLE_1, r_node);

        KRATOS_CHECK_DOF_IN_NODE(TURBULENT_KINETIC_ENERGY, r_node);
    }

    KRATOS_CATCH("");
}
template <unsigned int TDim>
GeometryData::IntegrationMethod KElementData<TDim>::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_2;
}

template <unsigned int TDim>
void KElementData<TDim>::CalculateElementData(const GeometryType& rGeometry,
                                              const Vector& rShapeFunctions,
                                              const Matrix& rShapeFunctionDerivatives,
                                              const ProcessInfo& rCurrentProcessInfo,
                                              const int Step)
{
    KRATOS_TRY

    const double c_mu = rCurrentProcessInfo[TURBULENCE_RANS_C_MU];

    mTurbulentKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        rGeometry, TURBULENT_VISCOSITY, rShapeFunctions);
    mTurbulentKineticEnergy = RansCalculationUtilities::EvaluateInPoint(
        rGeometry, TURBULENT_KINETIC_ENERGY, rShapeFunctions);
    mKinematicViscosity = RansCalculationUtilities::EvaluateInPoint(
        rGeometry, KINEMATIC_VISCOSITY, rShapeFunctions);
    mGamma = EvmKepsilonModelUtilities::CalculateGamma(
        c_mu, 1.0, mTurbulentKineticEnergy, mTurbulentKinematicViscosity);

    mVelocityDivergence = RansCalculationUtilities::GetDivergence(
        rGeometry, VELOCITY, rShapeFunctionDerivatives);

    RansCalculationUtilities::CalculateGradient<TDim>(
        this->mVelocityGradient, rGeometry, VELOCITY, rShapeFunctionDerivatives, Step);

    KRATOS_CATCH("");
}

template <unsigned int TDim>
double KElementData<TDim>::CalculateEffectiveKinematicViscosity(const Vector& rShapeFunctions,
                                                                const Matrix& rShapeFunctionDerivatives,
                                                                const ProcessInfo& rCurrentProcessInfo) const
{
    const double tke_sigma = rCurrentProcessInfo[TURBULENT_KINETIC_ENERGY_SIGMA];
    return mKinematicViscosity + mTurbulentKinematicViscosity / tke_sigma;
}

template <unsigned int TDim>
double KElementData<TDim>::CalculateReactionTerm(const Vector& rShapeFunctions,
                                                 const Matrix& rShapeFunctionDerivatives,
                                                 const ProcessInfo& rCurrentProcessInfo) const
{
    return RansCalculationUtilities::SoftPositive(mGamma + (2.0 / 3.0) * mVelocityDivergence);
}

template <unsigned int TDim>
double KElementData<TDim>::CalculateSourceTerm(const Vector& rShapeFunctions,
                                               const Matrix& rShapeFunctionDerivatives,
                                               const ProcessInfo& rCurrentProcessInfo) const
{
    double production = 0.0;

    production = EvmKepsilonModelUtilities::CalculateSourceTerm<TDim>(
        mVelocityGradient, mTurbulentKinematicViscosity, mTurbulentKineticEnergy);

    return production;
}

double CalculateTurbulentViscosity(const double C_mu,
                                   const double turbulent_kinetic_energy,
                                   const double turbulent_energy_dissipation_rate,
                                   const double f_mu)
{
    return C_mu * f_mu * std::pow(turbulent_kinetic_energy, 2) / turbulent_energy_dissipation_rate;
}

double CalculateFmu(const double y_plus)
{
    return 1.0 - std::exp(-0.0115 * y_plus);
}

double CalculateF2(const double turbulent_kinetic_energy,
                   const double kinematic_viscosity,
                   const double turbulent_energy_dissipation_rate)
{
    if (turbulent_energy_dissipation_rate == 0.0)
        return 1.0;
    else
    {
        const double Re_t = std::pow(turbulent_kinetic_energy, 2) /
                            (kinematic_viscosity * turbulent_energy_dissipation_rate);
        const double f2 = 1.0 - 0.22 * std::exp(-1.0 * std::pow(Re_t * (1.0 / 6.0), 2));

        return f2;
    }
}

template <unsigned int TDim>
double CalculateSourceTerm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                           const double turbulent_kinematic_viscosity,
                           const double turbulent_kinetic_energy)
{
    const double velocity_divergence =
        RansCalculationUtilities::CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient;
    noalias(symmetric_velocity_gradient) = rVelocityGradient + trans(rVelocityGradient);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;

    noalias(reynolds_stress_tensor) =
        turbulent_kinematic_viscosity *
        (symmetric_velocity_gradient - (2.0 / 3.0) * velocity_divergence * identity);

    double source = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            source += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);

    return source;
}

double CalculateGamma(const double C_mu,
                      const double f_mu,
                      const double turbulent_kinetic_energy,
                      const double turbulent_kinematic_viscosity)
{
    return RansCalculationUtilities::SoftPositive(
        C_mu * f_mu * turbulent_kinetic_energy / turbulent_kinematic_viscosity);
}

template double CalculateSourceTerm<2>(const BoundedMatrix<double, 2, 2>&, const double, const double);
template double CalculateSourceTerm<3>(const BoundedMatrix<double, 3, 3>&, const double, const double);

template class EpsilonElementData<2>;
template class EpsilonElementData<3>;

template class KElementData<2>;
template class KElementData<3>;

} // namespace EvmKepsilonModelUtilities

///@}

} // namespace Kratos
