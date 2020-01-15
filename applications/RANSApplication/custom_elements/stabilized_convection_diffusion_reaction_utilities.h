//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//

#if !defined(KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_UTILITIES_H_INCLUDED)
#define KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_UTILITIES_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "includes/ublas_interface.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

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

namespace StabilizedConvectionDiffusionReactionUtilities
{
inline double CalculatePsiOne(const double VelocityNorm, const double Tau, const double DynamicReaction)
{
    return VelocityNorm + Tau * VelocityNorm * DynamicReaction;
}

inline double CalculatePsiTwo(const double DynamicReaction, const double Tau, const double ElementLength)
{
    return (DynamicReaction + Tau * DynamicReaction * std::abs(DynamicReaction)) *
           std::pow(ElementLength, 2) * (1.0 / 6.0);
}

inline void CalculateStabilizationTau(double& rTau,
                                      double& rElementLength,
                                      const array_1d<double, 3>& rVelocity,
                                      const Matrix& rContravariantMetricTensor,
                                      const double Reaction,
                                      const double EffectiveKinematicViscosity,
                                      const double Alpha,
                                      const double Gamma,
                                      const double DeltaTime,
                                      const double DynamicTau)
{
    unsigned int dim = rContravariantMetricTensor.size2();
    const Vector velocity = RansCalculationUtilities::GetVector(rVelocity, dim);
    Vector temp(dim);
    noalias(temp) = prod(rContravariantMetricTensor, velocity);
    const double velocity_norm = norm_2(rVelocity);

    if (velocity_norm > 0.0)
    {
        rElementLength = 2.0 * velocity_norm / std::sqrt(inner_prod(velocity, temp));
    }
    else
    {
        rElementLength = 0.0;
        for (unsigned int i = 0; i < dim; ++i)
            for (unsigned int j = 0; j < dim; ++j)
                rElementLength += rContravariantMetricTensor(i, j);
        rElementLength = std::sqrt(1.0 / rElementLength) * 2.0;
    }

    const double stab_convection = std::pow(2.0 * norm_2(rVelocity) / rElementLength, 2);
    const double stab_diffusion = std::pow(
        12.0 * EffectiveKinematicViscosity / (rElementLength * rElementLength), 2);
    const double stab_dynamics =
        std::pow(DynamicTau * (1 - Alpha) / (Gamma * DeltaTime), 2);
    const double stab_reaction = std::pow(Reaction, 2);

    rTau = 1.0 / std::sqrt(stab_dynamics + stab_convection + stab_diffusion + stab_reaction);
}

inline void CalculateCrossWindDiffusionParameters(double& rChi,
                                                  double& rStreamLineDiffusionCoeff,
                                                  double& rCrossWindDiffusionCoeff,
                                                  const double VelocityMagnitude,
                                                  const double Tau,
                                                  const double EffectiveKinematicViscosity,
                                                  const double Reaction,
                                                  const double Alpha,
                                                  const double Gamma,
                                                  const double DeltaTime,
                                                  const double ElementLength,
                                                  const double DynamicTau)
{
    const double reaction_dynamics =
        Reaction + DynamicTau * (1 - Alpha) / (Gamma * DeltaTime);

    rChi = 2.0 / (std::abs(reaction_dynamics) * ElementLength + 2.0 * VelocityMagnitude);

    const double psi_one = CalculatePsiOne(VelocityMagnitude, Tau, reaction_dynamics);
    const double psi_two = CalculatePsiTwo(reaction_dynamics, Tau, ElementLength);

    double value =
        0.5 * std::abs(psi_one - Tau * VelocityMagnitude * reaction_dynamics) * ElementLength;
    value -= (EffectiveKinematicViscosity + Tau * std::pow(VelocityMagnitude, 2));
    value += psi_two;

    rStreamLineDiffusionCoeff = RansCalculationUtilities::SoftPositive(value);

    value = 0.5 * std::abs(psi_one) * ElementLength;
    value -= EffectiveKinematicViscosity;
    value += psi_two;

    rCrossWindDiffusionCoeff = RansCalculationUtilities::SoftPositive(value);
}
} // namespace StabilizedConvectionDiffusionReactionUtilities

///@}
///@name Kratos Classes
///@{
///@}

} // namespace Kratos

#endif