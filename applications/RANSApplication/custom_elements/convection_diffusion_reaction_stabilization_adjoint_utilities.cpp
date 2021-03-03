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

// External includes

// Project includes

// Application includes

// Include base h
#include "convection_diffusion_reaction_stabilization_adjoint_utilities.h"

namespace Kratos
{
namespace ConvectionDiffusionReactionStabilizationUtilities
{
template <>
double AdjointUtilities<2, 3>::Derivatives::Shape::CalculateElementLengthDerivative(
    const double DetJDerivative,
    const double ElementLength)
{
    return 0.31830988618378353 * DetJDerivative / ElementLength;
}

template <>
double AdjointUtilities<3, 4>::Derivatives::Shape::CalculateElementLengthDerivative(
    const double DetJDerivative,
    const double ElementLength)
{
    return 0.4714045207910277 * DetJDerivative / std::pow(ElementLength, 2);
}

template <unsigned int TDim, unsigned int TNumNodes>
double AdjointUtilities<TDim, TNumNodes>::CalculateVectorNormDerivative(
    const double VectorNorm,
    const ArrayD& rVector,
    const ArrayD& rVectorDerivative)
{
    if (VectorNorm > 0.0) {
        return inner_prod(rVector, rVectorDerivative) / VectorNorm;
    } else {
        return 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double AdjointUtilities<TDim, TNumNodes>::CalculateStabilizationTauDerivative(
    const double StabilizationTau,
    const double EffectiveVelocityMagnitude,
    const double EffectiveKinematicViscosity,
    const double ReactionTerm,
    const double ElementLength,
    const double EffectiveVelocityMagnitudeDerivative,
    const double EffectiveKinematicViscosityDerivative,
    const double ReactionTermDerivative,
    const double ElementLengthDerivative)
{
    double value = 0.0;

    value += 4.0 * EffectiveVelocityMagnitude * EffectiveVelocityMagnitudeDerivative /  std::pow(ElementLength, 2);
    value -= 4.0 * std::pow(EffectiveVelocityMagnitude, 2) * ElementLengthDerivative / std::pow(ElementLength, 3);

    value += 144.0 * EffectiveKinematicViscosity * EffectiveKinematicViscosityDerivative / std::pow(ElementLength, 4);
    value -= 288.0 * std::pow(EffectiveKinematicViscosity, 2) * ElementLengthDerivative / std::pow(ElementLength, 5);

    value += ReactionTerm * ReactionTermDerivative;

    return value  * (-1.0 * std::pow(StabilizationTau, 3));
}

template <unsigned int TDim, unsigned int TNumNodes>
double AdjointUtilities<TDim, TNumNodes>::CalculateAbsoluteValueDerivative(
    const double Value,
    const double ValueDerivative)
{
    return ValueDerivative * ((Value >= 0.0) ? 1.0 : -1.0);
}

template <unsigned int TDim, unsigned int TNumNodes>
void AdjointUtilities<TDim, TNumNodes>::CalculateDiscreteUpwindOperatorDerivative(
    VectorN& rOutput,
    const VectorN& rNodalScalarValues,
    const MatrixNN& rInputMatrix,
    const MatrixNN& rInputMatrixDerivative)
{
    rOutput.clear();
    for (IndexType a = 0; a < TNumNodes; ++a) {
        for (IndexType b = 0; b < TNumNodes; ++b) {
            if (a != b && (rInputMatrix(a, b) > 0.0 || rInputMatrix(b, a) > 0.0)) {
                const double flux = rNodalScalarValues[b] - rNodalScalarValues[a];
                if (rInputMatrix(a, b) > rInputMatrix(b, a)) {
                    rOutput[a] -= rInputMatrixDerivative(a, b) * flux;
                } else {
                    rOutput[a] -= rInputMatrixDerivative(b, a) * flux;
                }
            }
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
double AdjointUtilities<TDim, TNumNodes>::CalculatePositivityPreservingCoefficientDerivative(
    const double PositivityPreservingMatrixCoefficient,
    const MatrixNN& rInputMatrix,
    const MatrixNN& rInputMatrixDerivative)
{
    double value = 0.0;

    if (PositivityPreservingMatrixCoefficient > 0.0) {
        VectorN derivative_row = ZeroVector(TNumNodes);
        double max_negative_row_sum = 0.0;
        for (IndexType a = 0; a < TNumNodes; ++a) {
            double row_sum = 0.0;
            for (IndexType b = 0; b < TNumNodes; ++b) {
                row_sum += rInputMatrix(a, b);
            }

            if (row_sum < max_negative_row_sum) {
                max_negative_row_sum = row_sum;
                noalias(derivative_row) = row(rInputMatrixDerivative, a);
            }
        }

        for (IndexType a = 0; a < TNumNodes; ++a) {
            value -= derivative_row[a];
        }
    }

    return value;
}

// template instantiations
template class AdjointUtilities<2, 3>;
template class AdjointUtilities<3, 4>;

} // namespace ConvectionDiffusionReactionStabilizationUtilities
} // namespace Kratos
