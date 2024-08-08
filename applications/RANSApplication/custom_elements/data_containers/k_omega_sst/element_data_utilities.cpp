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
#include <cmath>

// Project includes

// Application includes

// Include base h
#include "element_data_utilities.h"

namespace Kratos
{
namespace KOmegaSSTElementData
{
double CalculateTanh(const double value)
{
    // const double inv_exp2x = std::exp(-2.0 * value);
    // return (1.0 - inv_exp2x) / (1 + inv_exp2x);
    return std::tanh(value);
}

double CalculateBlendedPhi(
    const double Phi1,
    const double Phi2,
    const double F1)
{
    return F1 * Phi1 + (1.0 - F1) * Phi2;
}

template<unsigned int TDim>
double CalculateCrossDiffusionTerm(
    const double SigmaTurbulentSpecificEnergyDissipationRate2,
    const double TurbulentSpecificEnergyDissipationRate,
    const array_1d<double, TDim>& rTurbulentKineticEnergyGradient,
    const array_1d<double, TDim>& rTurbulentSpecificEnergyDissipationRateGradient)
{
    KRATOS_TRY

    double value = inner_prod(rTurbulentKineticEnergyGradient,
                              rTurbulentSpecificEnergyDissipationRateGradient);
    value *= (2.0 * SigmaTurbulentSpecificEnergyDissipationRate2 /
              std::max(TurbulentSpecificEnergyDissipationRate, 1e-10));
    return value;

    KRATOS_CATCH("");
}

double CalculateF1(
    const double TurbulentKineticEnergy,
    const double TurbulentSpecificEnergyDissipationRate,
    const double KinematicViscosity,
    const double WallDistance,
    const double BetaStar,
    const double CrossDiffusion,
    const double SigmaTurbulentSpecificEnergyDissipationRate2)
{
    KRATOS_TRY

    const double y = std::max(WallDistance, 1e-10);
    const double y_2 = std::pow(y, 2);

    const double tke = std::max(TurbulentKineticEnergy, 0.0);
    const double omega = std::max(TurbulentSpecificEnergyDissipationRate, 1e-10);

    const double arg1 = CalculateArg1(BetaStar, tke, omega, y);
    const double arg2 = CalculateArg2(KinematicViscosity, omega, y_2);
    const double arg3 = CalculateArg3(SigmaTurbulentSpecificEnergyDissipationRate2, tke, CrossDiffusion, y_2);

    const double arg4 = std::max(arg1, arg2);

    const double arg5 = std::min(arg4, arg3);

    const double arg6 = std::min(arg5, 10.0);

    return CalculateTanh(std::pow(arg6, 4));

    KRATOS_CATCH("");
}

double CalculateF2(
    const double TurbulentKineticEnergy,
    const double TurbulentSpecificEnergyDissipationRate,
    const double KinematicViscosity,
    const double WallDistance,
    const double BetaStar)
{
    KRATOS_TRY

    const double y = std::max(WallDistance, 1e-10);
    const double y_2 = std::pow(y, 2);

    const double tke = std::max(TurbulentKineticEnergy, 0.0);
    const double omega = std::max(TurbulentSpecificEnergyDissipationRate, 1e-10);

    const double arg1 = CalculateArg1(BetaStar, tke, omega, y);
    const double arg2 = CalculateArg2(KinematicViscosity, omega, y_2);

    const double arg4 = std::max(2.0 * arg1, arg2);

    const double arg5 = std::min(arg4, 100.0);

    return CalculateTanh(std::pow(arg5, 2));

    KRATOS_CATCH("");
}

double CalculateTurbulentKinematicViscosity(
    const double TurbulentKineticEnergy,
    const double TurbulentSpecificEnergyDissipationRate,
    const double VorticityNorm,
    const double F2,
    const double A1)
{
    KRATOS_TRY

    return A1 * TurbulentKineticEnergy /
           std::max(std::max(A1 * TurbulentSpecificEnergyDissipationRate, VorticityNorm * F2), 1e-10);

    KRATOS_CATCH("");
}

double CalculateGamma(
    const double Beta,
    const double BetaStar,
    const double Sigma,
    const double Kappa)
{
    return Beta / BetaStar - Sigma * std::pow(Kappa, 2) / std::sqrt(BetaStar);
}

template<unsigned int TDim>
double CalculateVorticityNorm(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient)
{
    const BoundedMatrix<double, TDim, TDim> symmetric_velocity_gradient =
        (rVelocityGradient + trans(rVelocityGradient)) * 0.5;
    return norm_frobenius(symmetric_velocity_gradient) * 1.414;
}

double CalculateArg1(
    const double BetaStar,
    const double TurbulentKineticEnergy,
    const double TurbulentSpecificEnergyDissipationRate,
    const double WallDistance)
{
    return std::sqrt(TurbulentKineticEnergy) /
           std::max(BetaStar * TurbulentSpecificEnergyDissipationRate * WallDistance, 1e-10);
}

double CalculateArg2(
    const double KinematicViscosity,
    const double TurbulentSpecificEnergyDissipationRate,
    const double WallDistanceSquare)
{
    return 500.0 * KinematicViscosity / std::max(WallDistanceSquare * TurbulentSpecificEnergyDissipationRate, 1e-10);
}

double CalculateArg3(
    const double SigmaTurbulentSpecificEnergyDissipationRate2,
    const double TurbulentKineticEnergy,
    const double CrossDiffusion,
    const double WallDistanceSquare)
{
    return 4.0 * SigmaTurbulentSpecificEnergyDissipationRate2 * TurbulentKineticEnergy /
           std::max(std::max(CrossDiffusion, 1e-10) * WallDistanceSquare, 1e-10);
}

// template instantiations

template double CalculateCrossDiffusionTerm<2>(const double, const double, const array_1d<double, 2>&, const array_1d<double, 2>&);
template double CalculateCrossDiffusionTerm<3>(const double, const double, const array_1d<double, 3>&, const array_1d<double, 3>&);

template double CalculateVorticityNorm<2>(const BoundedMatrix<double, 2, 2>&);
template double CalculateVorticityNorm<3>(const BoundedMatrix<double, 3, 3>&);

} // namespace KOmegaSSTElementData

} // namespace Kratos
