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
    const array_1d<double, TDim>& rTurbulentSpecificEnergyDissipationRate)
{
    KRATOS_TRY

    double value = inner_prod(rTurbulentKineticEnergyGradient,
                              rTurbulentSpecificEnergyDissipationRate);
    value *= (2.0 * SigmaTurbulentSpecificEnergyDissipationRate2 /
              TurbulentSpecificEnergyDissipationRate);
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

    const double y = std::max(WallDistance, 1e-12);
    const double y_2 = std::pow(y, 2);

    const double tke = std::max(TurbulentKineticEnergy, 0.0);
    const double omega = std::max(TurbulentSpecificEnergyDissipationRate, 1e-12);

    double arg1 = std::max(std::sqrt(tke) / (BetaStar * omega * y),
                           500.0 * KinematicViscosity / (y_2 * omega));

    arg1 = std::min(std::min(arg1, 4.0 * SigmaTurbulentSpecificEnergyDissipationRate2 *
                                       tke / (std::max(CrossDiffusion, 1e-12) * y_2)),
                    10.0);

    return std::tanh(std::pow(arg1, 4));

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

    const double y = std::max(WallDistance, 1e-12);
    const double y_2 = std::pow(y, 2);

    const double tke = std::max(TurbulentKineticEnergy, 0.0);
    const double omega = std::max(TurbulentSpecificEnergyDissipationRate, 1e-12);

    const double arg2 =
        std::min(std::max(2.0 * std::sqrt(tke) / (BetaStar * omega * y),
                          500.0 * KinematicViscosity / (y_2 * omega)),
                 100.0);

    return std::tanh(std::pow(arg2, 2));

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
           std::max(A1 * TurbulentSpecificEnergyDissipationRate, VorticityNorm * F2);

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

// template instantiations

template double CalculateCrossDiffusionTerm<2>(const double, const double, const array_1d<double, 2>&, const array_1d<double, 2>&);
template double CalculateCrossDiffusionTerm<3>(const double, const double, const array_1d<double, 3>&, const array_1d<double, 3>&);

} // namespace KOmegaSSTElementData

} // namespace Kratos
