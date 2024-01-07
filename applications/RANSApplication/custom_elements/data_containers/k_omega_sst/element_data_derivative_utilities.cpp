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
#include "element_data_utilities.h"

// Include base h
#include "element_data_derivative_utilities.h"

namespace Kratos
{
namespace KOmegaSSTElementData
{

template<unsigned int TDim>
double AdjointUtilities<TDim>::CalculateTanhDerivative(
    const double Value,
    const double ValueDerivative)
{
    const double exp_2value = std::exp(-2.0 * Value);
    return 4.0 * ValueDerivative * exp_2value / std::pow(exp_2value + 1, 2);
}

template<unsigned int TDim>
double AdjointUtilities<TDim>::CalculateBlendedPhiDerivative(
    const double Phi1,
    const double Phi1Derivative,
    const double Phi2,
    const double Phi2Derivative,
    const double F1,
    const double F1Derivative)
{
    double value = 0.0;

    value += F1Derivative * Phi1;
    value += F1 * Phi1Derivative;
    value -= F1Derivative * Phi2;
    value += (1.0 - F1) * Phi2Derivative;

    return value;
}

template<unsigned int TDim>
double AdjointUtilities<TDim>::CalculateCrossDiffusionTermDerivative(
    const double SigmaTurbulentSpecificEnergyDissipationRate2,
    const double TurbulentSpecificEnergyDissipationRate,
    const double TurbulentSpecificEnergyDissipationRateDerivative,
    const ArrayD& rTurbulentKineticEnergyGradient,
    const ArrayD& rTurbulentKineticEnergyGradientDerivative,
    const ArrayD& rTurbulentSpecificEnergyDissipationRateGradient,
    const ArrayD& rTurbulentSpecificEnergyDissipationRateGradientDerivative)
{
    double value = 0.0;

    double omega;
    double omega_detivative;
    if (TurbulentSpecificEnergyDissipationRate <= 1e-12) {
        omega = 1e-12;
        omega_detivative = 0.0;
    } else {
        omega = TurbulentSpecificEnergyDissipationRate;
        omega_detivative = TurbulentSpecificEnergyDissipationRateDerivative;
    }

    value += inner_prod(rTurbulentKineticEnergyGradientDerivative, rTurbulentSpecificEnergyDissipationRateGradient);
    value += inner_prod(rTurbulentKineticEnergyGradient, rTurbulentSpecificEnergyDissipationRateGradientDerivative);

    value = value / omega -
            inner_prod(rTurbulentKineticEnergyGradient,
                       rTurbulentSpecificEnergyDissipationRateGradient) *
                omega_detivative /
                std::pow(omega, 2);

    value *= 2.0 * SigmaTurbulentSpecificEnergyDissipationRate2;

    return value;
}

template<unsigned int TDim>
double AdjointUtilities<TDim>::CalculateF1Derivative(
    const double TurbulentKineticEnergy,
    const double TurbulentKineticEnergyDerivative,
    const double TurbulentSpecificEnergyDissipationRate,
    const double TurbulentSpecificEnergyDissipationRateDerivative,
    const double KinematicViscosity,
    const double WallDistance,
    const double WallDistanceDerivative,
    const double BetaStar,
    const double CrossDiffusion,
    const double CrossDiffusionDerivative,
    const double SigmaTurbulentSpecificEnergyDissipationRate2)
{
    KRATOS_TRY

    const double y = std::max(WallDistance, 1e-12);
    const double y_derivative = (WallDistance > 1e-12) ? WallDistanceDerivative : 0.0;

    const double y_2 = std::pow(y, 2);
    const double y_2_derivative = 2.0 * y * y_derivative;

    const double tke = std::max(TurbulentKineticEnergy, 1e-12);
    const double tke_derivative = (TurbulentKineticEnergy > 0.0) ? TurbulentKineticEnergyDerivative : 0.0;

    const double omega = std::max(TurbulentSpecificEnergyDissipationRate, 1e-12);
    const double omega_derivative = (TurbulentSpecificEnergyDissipationRate > 1e-12) ? TurbulentSpecificEnergyDissipationRateDerivative : 0.0;

    const double arg1 = CalculateArg1(BetaStar, tke, omega, y);
    const double arg1_derivative = CalculateArg1Derivative(BetaStar, tke, tke_derivative, omega, omega_derivative, y, y_derivative);

    const double arg2 = CalculateArg2(KinematicViscosity, omega, y_2);
    const double arg2_derivative = CalculateArg2Derivative(KinematicViscosity, omega, omega_derivative, y_2, y_2_derivative);

    const double arg3 = CalculateArg3(SigmaTurbulentSpecificEnergyDissipationRate2, tke, CrossDiffusion, y_2);
    const double arg3_derivative = CalculateArg3Derivative(SigmaTurbulentSpecificEnergyDissipationRate2, tke, tke_derivative, CrossDiffusion, CrossDiffusionDerivative, y_2, y_2_derivative);

    const double arg4 = std::max(arg1, arg2);
    const double arg4_derivative = (arg1 > arg2) ? arg1_derivative : arg2_derivative;

    const double arg5 = std::min(arg4, arg3);
    const double arg5_derivative = (arg4 < arg3) ? arg4_derivative : arg3_derivative;

    const double arg6 = std::min(arg5, 10.0);
    const double arg6_derivative = (arg5 < 10.0) ? arg5_derivative : 0.0;

    const double value = std::pow(arg6, 4);
    const double value_derivative = 4.0 * std::pow(arg6, 3) * arg6_derivative;

    return CalculateTanhDerivative(value, value_derivative);

    KRATOS_CATCH("");
}

template<unsigned int TDim>
double AdjointUtilities<TDim>::CalculateF2Derivative(
        const double TurbulentKineticEnergy,
        const double TurbulentKineticEnergyDerivative,
        const double TurbulentSpecificEnergyDissipationRate,
        const double TurbulentSpecificEnergyDissipationRateDerivative,
        const double KinematicViscosity,
        const double WallDistance,
        const double WallDistanceDerivative,
        const double BetaStar)
{
    KRATOS_TRY

    const double y = std::max(WallDistance, 1e-12);
    const double y_derivative = (WallDistance > 1e-12) ? WallDistanceDerivative : 0.0;

    const double y_2 = std::pow(y, 2);
    const double y_2_derivative = 2.0 * y * y_derivative;

    const double tke = std::max(TurbulentKineticEnergy, 1e-12);
    const double tke_derivative = (TurbulentKineticEnergy > 0.0) ? TurbulentKineticEnergyDerivative : 0.0;

    const double omega = std::max(TurbulentSpecificEnergyDissipationRate, 1e-12);
    const double omega_derivative = (TurbulentSpecificEnergyDissipationRate > 1e-12) ? TurbulentSpecificEnergyDissipationRateDerivative : 0.0;

    const double arg1 = 2.0 * CalculateArg1(BetaStar, tke, omega, y);
    const double arg1_derivative = 2.0 * CalculateArg1Derivative(BetaStar, tke, tke_derivative, omega, omega_derivative, y, y_derivative);

    const double arg2 = CalculateArg2(KinematicViscosity, omega, y_2);
    const double arg2_derivative = CalculateArg2Derivative(KinematicViscosity, omega, omega_derivative, y_2, y_2_derivative);

    const double arg4 = std::max(arg1, arg2);
    const double arg4_derivative = (arg1 > arg2) ? arg1_derivative : arg2_derivative;

    const double arg5 = std::min(arg4, 100.0);
    const double arg5_derivative = (arg4 < 100.0) ? arg4_derivative : 0.0;

    const double arg6 = std::pow(arg5, 2);
    const double arg6_derivative = 2.0 * arg5 * arg5_derivative;

    return CalculateTanhDerivative(arg6, arg6_derivative);

    KRATOS_CATCH("");
}

template<unsigned int TDim>
double AdjointUtilities<TDim>::CalculateTurbulentKinematicViscosityDerivative(
    const double TurbulentKineticEnergy,
    const double TurbulentKineticEnergyDerivative,
    const double TurbulentSpecificEnergyDissipationRate,
    const double TurbulentSpecificEnergyDissipationRateDerivative,
    const double VorticityNorm,
    const double VorticityNormDerivative,
    const double F2,
    const double F2Derivative,
    const double A1)
{
    KRATOS_TRY

    const double arg1 = A1 * TurbulentSpecificEnergyDissipationRate;
    const double arg1_derivative = A1 * TurbulentSpecificEnergyDissipationRateDerivative;

    const double arg2 = VorticityNorm * F2;
    const double arg2_derivative = VorticityNormDerivative * F2 + VorticityNorm * F2Derivative;

    const double denominator = std::max(arg1, arg2);
    const double denominator_derivative = (arg1 > arg2) ? arg1_derivative : arg2_derivative;

    double value = 0.0;

    value += A1 * TurbulentKineticEnergyDerivative / denominator;
    value -= A1 * TurbulentKineticEnergy * denominator_derivative / std::pow(denominator, 2);

    return value;

    KRATOS_CATCH("");
}

template<unsigned int TDim>
double AdjointUtilities<TDim>::CalculateArg1Derivative(
    const double BetaStar,
    const double TurbulentKineticEnergy,
    const double TurbulentKineticEnergyDerivative,
    const double TurbulentSpecificEnergyDissipationRate,
    const double TurbulentSpecificEnergyDissipationRateDerivative,
    const double WallDistance,
    const double WallDistanceDerivative)
{
    double value = 0.0;

    value += 0.5 * TurbulentKineticEnergyDerivative / ( std::sqrt(TurbulentKineticEnergy) * BetaStar * TurbulentSpecificEnergyDissipationRate * WallDistance);
    value -= std::sqrt(TurbulentKineticEnergy) * TurbulentSpecificEnergyDissipationRateDerivative / (BetaStar * std::pow(TurbulentSpecificEnergyDissipationRate, 2) * WallDistance);
    value -= std::sqrt(TurbulentKineticEnergy) * WallDistanceDerivative / (BetaStar * TurbulentSpecificEnergyDissipationRate * std::pow(WallDistance, 2));

    return value;
}

template<unsigned int TDim>
double AdjointUtilities<TDim>::CalculateArg2Derivative(
    const double KinematicViscosity,
    const double TurbulentSpecificEnergyDissipationRate,
    const double TurbulentSpecificEnergyDissipationRateDerivative,
    const double WallDistanceSquare,
    const double WallDistanceSquareDerivative)
{
    double value = 0.0;

    value -= WallDistanceSquareDerivative / (std::pow(WallDistanceSquare, 2) * TurbulentSpecificEnergyDissipationRate);
    value -= TurbulentSpecificEnergyDissipationRateDerivative / (WallDistanceSquare * std::pow(TurbulentSpecificEnergyDissipationRate, 2));

    return value * 500.0 * KinematicViscosity;
}

template<unsigned int TDim>
double AdjointUtilities<TDim>::CalculateArg3Derivative(
    const double SigmaTurbulentSpecificEnergyDissipationRate2,
    const double TurbulentKineticEnergy,
    const double TurbulentKineticEnergyDerivative,
    const double CrossDiffusion,
    const double CrossDiffusionDerivative,
    const double WallDistanceSquare,
    const double WallDistanceSquareDerivative)
{
    double value = 0.0;

    value += TurbulentKineticEnergyDerivative / (std::max(CrossDiffusion, 1e-12) * WallDistanceSquare);
    value -= (CrossDiffusion > 1e-12) ? TurbulentKineticEnergy * CrossDiffusionDerivative / (std::pow(CrossDiffusion, 2) * WallDistanceSquare) : 0.0;
    value -= TurbulentKineticEnergy * WallDistanceSquareDerivative / (std::max(CrossDiffusion, 1e-12) * std::pow(WallDistanceSquare, 2));

    return 4.0 * SigmaTurbulentSpecificEnergyDissipationRate2 * value;
}

// template instantiations

template class AdjointUtilities<2>;
template class AdjointUtilities<3>;

} // namespace KOmegaSSTElementData

} // namespace Kratos
