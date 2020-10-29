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

#if !defined(KRATOS_K_OMEGA_SST_ELEMENT_DATA_UTILITIES_H_INCLUDED)
#define KRATOS_K_OMEGA_SST_ELEMENT_DATA_UTILITIES_H_INCLUDED

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
double CalculateBlendedPhi(
    const double Phi1,
    const double Phi2,
    const double F1);

double CalculateCrossDiffusionTerm(
    const double SigmaTurbulentSpecificEnergyDissipationRate2,
    const double TurbulentSpecificEnergyDissipationRate,
    const array_1d<double, 3>& rTurbulentKineticEnergyGradient,
    const array_1d<double, 3>& rTurbulentSpecificEnergyDissipationRate);

double CalculateF1(
    const double TurbulentKineticEnergy,
    const double TurbulentSpecificEnergyDissipationRate,
    const double KinematicViscosity,
    const double WallDistance,
    const double BetaStar,
    const double CrossDiffusion,
    const double SigmaTurbulentSpecificEnergyDissipationRate2);

double CalculateF2(
    const double TurbulentKineticEnergy,
    const double TurbulentSpecificEnergyDissipationRate,
    const double KinematicViscosity,
    const double WallDistance,
    const double BetaStar);

double CalculateTurbulentKinematicViscosity(
    const double TurbulentKineticEnergy,
    const double TurbulentSpecificEnergyDissipationRate,
    const double VorticityNorm,
    const double F2,
    const double A1);

template <unsigned int TDim>
array_1d<double, 3> CalculateVorticity(
    const BoundedMatrix<double, TDim, TDim>& rVelocityGradient);

double CalculateGamma(
    const double Beta,
    const double BetaStar,
    const double Sigma,
    const double Kappa);

} // namespace KOmegaSSTElementData

///@}

} // namespace Kratos

#endif