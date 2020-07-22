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

#if !defined(KRATOS_K_EPSILON_HIGH_RE_ELEMENT_DATA_UTILITIES_H_INCLUDED)
#define KRATOS_K_EPSILON_HIGH_RE_ELEMENT_DATA_UTILITIES_H_INCLUDED

// System includes

// Project includes
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name  Functions
///@{

namespace KEpsilonElementData
{
double CalculateTurbulentViscosity(
    const double Cmu,
    const double TurbulentKineticEnergy,
    const double TurbulentEnergyDissipationRate);

template <unsigned int TDim>
double CalculateSourceTerm(
    const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
    const double TurbulentKinematicViscosity);

double CalculateGamma(
    const double Cmu,
    const double TurbulentKineticEnergy,
    const double TurbulentKinematicViscosity);

} // namespace KEpsilonElementData

///@}

} // namespace Kratos

#endif