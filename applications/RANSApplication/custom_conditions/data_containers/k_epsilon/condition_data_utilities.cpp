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
#include "condition_data_utilities.h"

namespace Kratos
{
///@name  Classes
///@{

double KEpsilonConditionDataUtilities::CalculateWallFlux(
    const double KinematicViscosity,
    const double EpsilonSigma,
    const double UTau,
    const double Kappa,
    const double YPlus)
{
    const double wall_turbulent_viscosity = Kappa * KinematicViscosity * YPlus;
    return (KinematicViscosity + wall_turbulent_viscosity / EpsilonSigma) *
           std::pow(UTau, 5) / (Kappa * std::pow(YPlus * KinematicViscosity, 2));
}

double KEpsilonConditionDataUtilities::CalculateWallFluxDerivative(
    const double KinematicViscosity,
    const double EpsilonSigma,
    const double UTau,
    const double UTauDerivative,
    const double Kappa,
    const double YPlus,
    const double YPlusDerivative)
{
    const double wall_turbulent_viscosity = Kappa * KinematicViscosity * YPlus;
    const double wall_turbulent_viscosity_derivative = Kappa * KinematicViscosity * YPlusDerivative;

    const double effective_viscosity = KinematicViscosity + wall_turbulent_viscosity / EpsilonSigma;
    const double effective_viscosity_derivative = wall_turbulent_viscosity_derivative / EpsilonSigma;

    const double coeff = 1.0 / (Kappa * std::pow(YPlus * KinematicViscosity, 2));
    const double coeff_derivative = (-1.0 * std::pow(coeff, 2)) * (Kappa * std::pow(KinematicViscosity, 2) * 2.0 * YPlus * YPlusDerivative);

    const double u_tau_5 = std::pow(UTau, 5);

    double value = 0.0;

    value += effective_viscosity_derivative * u_tau_5 * coeff;
    value += effective_viscosity * 5.0 * std::pow(UTau, 4) * UTauDerivative * coeff;
    value += effective_viscosity * u_tau_5 * coeff_derivative;

    return value;
}

} // namespace Kratos
