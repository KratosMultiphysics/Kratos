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

double KOmegaConditionDataUtilities::CalculateWallFlux(
    const double KinematicViscosity,
    const double OmegaSigma,
    const double UTau,
    const double Cmu25,
    const double Kappa,
    const double YPlus)
{
    const double wall_turbulent_viscosity = Kappa * KinematicViscosity * YPlus;
    return (KinematicViscosity + OmegaSigma * wall_turbulent_viscosity) *
           std::pow(UTau, 3) /
           (Kappa * std::pow(Cmu25 * YPlus * KinematicViscosity, 2));
}

double KOmegaConditionDataUtilities::CalculateWallFluxDerivative(
    const double KinematicViscosity,
    const double OmegaSigma,
    const double OmegaSigmaDerivative,
    const double UTau,
    const double UTauDerivative,
    const double Cmu25,
    const double Kappa,
    const double YPlus,
    const double YPlusDerivative)
{
    const double wall_turbulent_viscosity = Kappa * KinematicViscosity * YPlus;
    const double wall_turbulent_viscosity_derivative = Kappa * KinematicViscosity * YPlusDerivative;

    const double effective_viscosity = KinematicViscosity + wall_turbulent_viscosity * OmegaSigma;
    const double effective_viscosity_derivative = wall_turbulent_viscosity_derivative * OmegaSigma + wall_turbulent_viscosity * OmegaSigmaDerivative;

    const double coeff = 1.0 / (Kappa * std::pow(Cmu25 * YPlus * KinematicViscosity, 2));
    const double coeff_derivative = (-1.0 * std::pow(coeff, 2)) * (Kappa * std::pow(Cmu25 * KinematicViscosity, 2) * 2.0 * YPlus * YPlusDerivative);

    const double u_tau_3 = std::pow(UTau, 3);

    double value = 0.0;

    value += effective_viscosity_derivative * u_tau_3 * coeff;
    value += effective_viscosity * 3.0 * std::pow(UTau, 2) * UTauDerivative * coeff;
    value += effective_viscosity * u_tau_3 * coeff_derivative;

    return value;
}

} // namespace Kratos
