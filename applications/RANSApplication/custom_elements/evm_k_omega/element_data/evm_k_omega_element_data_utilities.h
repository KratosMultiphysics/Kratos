//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Dharmin Shah (https://github.com/sdharmin)
//                   Bence Rochlitz (https://github.com/bencerochlitz)
//
//  Supervised by:   Jordi Cotela (https://github.com/jcotela)
//                   Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_EVM_K_OMEGA_ELEMENT_DATA_UTILITIES_H_INCLUDED)
#define KRATOS_EVM_K_OMEGA_ELEMENT_DATA_UTILITIES_H_INCLUDED // ALL Epsilons renamed for omega

// System includes

// Project includes
#include "includes/ublas_interface.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

namespace EvmKOmegaElementDataUtilities
{
// change calculation in this one:
double CalculateTurbulentViscosity(const double TurbulentKineticEnergy,
                                   const double TurbulentSpecificEnergyDissipationRate);

// we need to calculate our own parameters, which are: Beta and Beta_Star --> parameters of the CDR equations.
template <unsigned int TDim>
double CalculateBeta(const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                     const double VelocityDivergence,
                     const double TurbulentSpecificDissipationRate,
                     const double BetaZero,
                     const double BetaStar); // define the input parameters

double CalculateSigmaD(const array_1d<double, 3>& rTurbulentKineticEnergyGradient,
                       const array_1d<double, 3>& rTurbulentSpecificEnergyDissipationRateGradient);

} // namespace EvmKOmegaElementDataUtilities

///@}

} // namespace Kratos

#endif