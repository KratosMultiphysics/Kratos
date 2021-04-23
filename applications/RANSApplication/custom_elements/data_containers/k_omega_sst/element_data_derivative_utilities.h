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

#if !defined(KRATOS_K_OMEGA_SST_ELEMENT_DATA_DERIVATIVE_UTILITIES_H_INCLUDED)
#define KRATOS_K_OMEGA_SST_ELEMENT_DATA_DERIVATIVE_UTILITIES_H_INCLUDED

// System includes

// Project includes
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name  Functions
///@{

namespace KOmegaSSTElementData
{

///@name Classes
///{

template<unsigned int TDim>
class AdjointUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using ArrayD = array_1d<double, TDim>;

    ///@}
    ///@name Static Operations
    ///@{

    static double CalculateBlendedPhiDerivative(
        const double Phi1,
        const double Phi1Derivative,
        const double Phi2,
        const double Phi2Derivative,
        const double F1,
        const double F1Derivative);

    static double CalculateCrossDiffusionTermDerivative(
        const double SigmaTurbulentSpecificEnergyDissipationRate2,
        const double TurbulentSpecificEnergyDissipationRate,
        const double TurbulentSpecificEnergyDissipationRateDerivative,
        const ArrayD& rTurbulentKineticEnergyGradient,
        const ArrayD& rTurbulentKineticEnergyGradientDerivative,
        const ArrayD& rTurbulentSpecificEnergyDissipationRateGradient,
        const ArrayD& rTurbulentSpecificEnergyDissipationRateGradientDerivative);

    static double CalculateF1Derivative(
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
        const double SigmaTurbulentSpecificEnergyDissipationRate2);

    static double CalculateF2Derivative(
        const double TurbulentKineticEnergy,
        const double TurbulentKineticEnergyDerivative,
        const double TurbulentSpecificEnergyDissipationRate,
        const double TurbulentSpecificEnergyDissipationRateDerivative,
        const double KinematicViscosity,
        const double WallDistance,
        const double WallDistanceDerivative,
        const double BetaStar);

    static double CalculateTurbulentKinematicViscosityDerivative(
        const double TurbulentKineticEnergy,
        const double TurbulentKineticEnergyDerivative,
        const double TurbulentSpecificEnergyDissipationRate,
        const double TurbulentSpecificEnergyDissipationRateDerivative,
        const double VorticityNorm,
        const double VorticityNormDerivative,
        const double F2,
        const double F2Derivative,
        const double A1);

    static double CalculateArg1Derivative(
        const double BetaStar,
        const double TurbulentKineticEnergy,
        const double TurbulentKineticEnergyDerivative,
        const double TurbulentSpecificEnergyDissipationRate,
        const double TurbulentSpecificEnergyDissipationRateDerivative,
        const double WallDistance,
        const double WallDistanceDerivative);

    static double CalculateArg2Derivative(
        const double KinematicViscosity,
        const double TurbulentSpecificEnergyDissipationRate,
        const double TurbulentSpecificEnergyDissipationRateDerivative,
        const double WallDistanceSquare,
        const double WallDistanceSquareDerivative);

    static double CalculateArg3Derivative(
        const double SigmaTurbulentSpecificEnergyDissipationRate2,
        const double TurbulentKineticEnergy,
        const double TurbulentKineticEnergyDerivative,
        const double CrossDiffusion,
        const double CrossDiffusionDerivative,
        const double WallDistanceSquare,
        const double WallDistanceSquareDerivative);

    static double CalculateTanhDerivative(
        const double Value,
        const double ValueDerivative);

    ///@}
};

///@}

} // namespace KOmegaSSTElementData

///@}

} // namespace Kratos

#endif  // KRATOS_K_OMEGA_SST_ELEMENT_DATA_DERIVATIVE_UTILITIES_H_INCLUDED