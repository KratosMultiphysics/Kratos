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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_STABILIZATION_UTILITIES_H_INCLUDED)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_STABILIZATION_UTILITIES_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "includes/ublas_interface.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"

namespace Kratos
{
///@name  Functions
///@{

namespace ConvectionDiffusionReactionStabilizationUtilities
{
inline double CalculatePsiOne(
    const double VelocityNorm,
    const double Tau,
    const double DynamicReaction)
{
    return VelocityNorm + Tau * VelocityNorm * DynamicReaction;
}

inline double CalculatePsiTwo(
    const double DynamicReaction,
    const double Tau,
    const double ElementLength)
{
    return (DynamicReaction + Tau * DynamicReaction * std::abs(DynamicReaction)) *
           std::pow(ElementLength, 2) * (1.0 / 6.0);
}

inline void CalculateStabilizationTau(
    double& rTau,
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
    const Vector& velocity = RansCalculationUtilities::GetVector(rVelocity, dim);
    Vector temp(dim);
    noalias(temp) = prod(rContravariantMetricTensor, velocity);
    const double velocity_norm = norm_2(rVelocity);

    if (velocity_norm > 0.0) {
        rElementLength = 2.0 * velocity_norm / std::sqrt(inner_prod(velocity, temp));
    } else {
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

inline double CalculateStabilizationTau(
    const double ElementLength,
    const double Velocity,
    const double Reaction,
    const double EffectiveKinematicViscosity,
    const double Alpha,
    const double Gamma,
    const double DeltaTime,
    const double DynamicTau)
{
    const double stab_convection = std::pow(2.0 * Velocity / ElementLength, 2);
    const double stab_diffusion = std::pow(
        12.0 * EffectiveKinematicViscosity / (ElementLength * ElementLength), 2);
    const double stab_dynamics =
        std::pow(DynamicTau * (1 - Alpha) / (Gamma * DeltaTime), 2);
    const double stab_reaction = std::pow(Reaction, 2);

    return 1.0 / std::sqrt(stab_dynamics + stab_convection + stab_diffusion + stab_reaction);
}

inline void CalculateCrossWindDiffusionParameters(
    double& rChi,
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

template <unsigned int TSize>
inline void CalculateDiscreteUpwindOperator(
    double& rScalarCoeff,
    BoundedMatrix<double, TSize, TSize>& rDiffusionMatrix,
    const BoundedMatrix<double, TSize, TSize>& rInputMatrix)
{
    rDiffusionMatrix.clear();

    for (unsigned int a = 0; a < TSize; ++a) {
        for (unsigned int b = a + 1; b < TSize; ++b) {
            rDiffusionMatrix(a, b) =
                -std::max(std::max(rInputMatrix(a, b), rInputMatrix(b, a)), 0.0);
            rDiffusionMatrix(b, a) = rDiffusionMatrix(a, b);
        }
    }

    for (unsigned int a = 0; a < TSize; ++a) {
        double row_sum = 0.0;
        for (unsigned int b = 0; b < TSize; ++b) {
            // all the diagonal terms are initialized with zero
            row_sum += rDiffusionMatrix(a, b);
        }
        rDiffusionMatrix(a, a) = -row_sum;
    }

    rScalarCoeff = norm_frobenius(rDiffusionMatrix);
}

inline double CalculatePositivityPreservingMatrix(
    const Matrix& rInputMatrix)
{
    double coefficient = 0.0;
    for (unsigned int a = 0; a < rInputMatrix.size1(); ++a) {
        double row_sum = 0.0;
        for (unsigned int b = 0; b < rInputMatrix.size2(); ++b) {
            row_sum += rInputMatrix(a, b);
        }
        coefficient = std::max(coefficient, -row_sum);
    }
    return coefficient;
}

inline void AddMassMatrixSUPGStabilizationGaussPointContributions(
    Matrix& rMassMatrix,
    const double AbsoluteReactionTerm,
    const double Tau,
    const Vector& rVelocityConvectiveTerms,
    const double GaussWeight,
    const Vector& rGaussShapeFunctions)
{
    const IndexType n = rMassMatrix.size1();

    for (IndexType i = 0; i < n; ++i) {
        for (IndexType j = 0; j < n; ++j) {
            rMassMatrix(i, j) += GaussWeight * Tau *
                                 (rVelocityConvectiveTerms[i] +
                                  AbsoluteReactionTerm * rGaussShapeFunctions[i]) *
                                 rGaussShapeFunctions[j];
        }
    }
}

inline void AddDampingMatrixSUPGStabilizationGaussPointContributions(
    Matrix& rDampingMatrix,
    const double ReactionTerm,
    const double Tau,
    const Vector& rVelocityConvectiveTerms,
    const double GaussWeight,
    const Vector& rGaussShapeFunctions)
{
    const IndexType n = rDampingMatrix.size1();
    const double s = std::abs(ReactionTerm);

    for (IndexType a = 0; a < n; ++a) {
        for (IndexType b = 0; b < n; ++b) {
            double value = 0.0;

            // Adding SUPG stabilization terms
            value += Tau * (rVelocityConvectiveTerms[a] + s * rGaussShapeFunctions[a]) *
                     rVelocityConvectiveTerms[b];
            value += Tau * (rVelocityConvectiveTerms[a] + s * rGaussShapeFunctions[a]) *
                     ReactionTerm * rGaussShapeFunctions[b];

            rDampingMatrix(a, b) += value * GaussWeight;
        }
    }
}

inline void AddSourceTermWithSUPGStabilizationGaussPointContributions(
    Vector& rRightHandSideVector,
    const double SourceTerm,
    const double AbsoluteReactionTerm,
    const double Tau,
    const Vector& rVelocityConvectiveTerms,
    const double GaussWeight,
    const Vector& rGaussShapeFunctions)
{
    for (IndexType a = 0; a < rRightHandSideVector.size(); ++a) {
        double value = 0.0;

        value += rGaussShapeFunctions[a] * SourceTerm;

        // Add supg stabilization terms
        value += (rVelocityConvectiveTerms[a] +
                  AbsoluteReactionTerm * rGaussShapeFunctions[a]) *
                 Tau * SourceTerm;

        rRightHandSideVector[a] += GaussWeight * value;
    }
}

} // namespace ConvectionDiffusionReactionStabilizationUtilities

///@}

} // namespace Kratos

#endif