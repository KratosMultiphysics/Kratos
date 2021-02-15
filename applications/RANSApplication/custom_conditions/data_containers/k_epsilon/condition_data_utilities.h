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

#if !defined(KRATOS_K_EPSILON_CONDITION_DATA_UTILITIES_H_INCLUDED)
#define KRATOS_K_EPSILON_CONDITION_DATA_UTILITIES_H_INCLUDED

// System includes

// Project includes

// Application includes

namespace Kratos
{
///@name  Classes
///@{

class KEpsilonConditionDataUtilities
{
public:
    ///@name Static Operations
    ///@{

    static double CalculateWallFlux(
        const double KinematicViscosity,
        const double EpsilonSigma,
        const double UTau,
        const double Kappa,
        const double YPlus);

    static double CalculateWallFluxDerivative(
        const double KinematicViscosity,
        const double EpsilonSigma,
        const double UTau,
        const double UTauDerivative,
        const double Kappa,
        const double YPlus,
        const double YPlusDerivative);

    ///@}
};

///@}

} // namespace Kratos

#endif // KRATOS_K_EPSILON_CONDITION_DATA_UTILITIES_H_INCLUDED