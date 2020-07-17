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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_STABILIZATION_ADJOINT_UTILITIES_H_INCLUDED)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_STABILIZATION_ADJOINT_UTILITIES_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name  Functions
///@{

namespace ConvectionDiffusionReactionStabilizationUtilities
{
template <unsigned int TDim, unsigned int TNumNodes>
class AdjointUtilities
{
public:
    using IndexType = unsigned int;
    using BoundedArray3D = array_1d<double, 3>;
    using BoundedVectorN = BoundedVector<double, TNumNodes>;
    using BoundedVectorND = BoundedVector<double, TNumNodes * TDim>;
    using BoundedMatrixND = BoundedMatrix<double, TNumNodes, TDim>;
    using BoundedMatrixNN = BoundedMatrix<double, TNumNodes, TNumNodes>;
    using BoundedMatrixNDD = BoundedMatrix<double, TNumNodes * TDim, TDim>;
    using BoundedMatrixNNN =
        BoundedVector<BoundedMatrix<double, TNumNodes, TNumNodes>, TNumNodes>;

    template <unsigned int TRowSize>
    static BoundedVector<double, TRowSize> MatrixVectorProduct(
        const BoundedMatrix<double, TRowSize, TDim>& rM, const BoundedArray3D& rV)
    {
        BoundedVector<double, TRowSize> output;
        for (IndexType i = 0; i < TRowSize; ++i)
        {
            output[i] = rM(i, 0) * rV[0];
            for (IndexType j = 1; j < TDim; ++j)
            {
                output[i] += rM(i, j) * rV[j];
            }
        }

        return output;
    }

    static BoundedMatrixND ConvertNDVectorToNDMatrix(const BoundedVectorND& rInput)
    {
        BoundedMatrixND output;
        for (IndexType c = 0; c < TNumNodes; ++c)
        {
            for (IndexType k = 0; k < TDim; ++k)
            {
                output(c, k) = rInput[c * TDim + k];
            }
        }

        return output;
    }

    static double CalculateElementLengthShapeDerivative(const double ElementLength,
                                                        const double DetJDerivative);

    static void CalculateEffectiveVelocityMagnitudeVelocityDerivative(
        BoundedMatrix<double, TNumNodes, TDim>& rOutput,
        const double EffectiveVelocityMagnitude,
        const array_1d<double, 3>& rEffectiveVelocity,
        const BoundedMatrixNDD& rEffectiveVelocityDerivative)
    {
        rOutput.clear();
        const BoundedVectorND& r_velocity_dot_velocity_derivative =
            MatrixVectorProduct<TNumNodes * TDim>(rEffectiveVelocityDerivative,
                                                  rEffectiveVelocity);

        if (EffectiveVelocityMagnitude > 0.0)
        {
            noalias(rOutput) = ConvertNDVectorToNDMatrix(r_velocity_dot_velocity_derivative);
            noalias(rOutput) = rOutput * (1 / EffectiveVelocityMagnitude);
        }
    }

    static double CalculateEffectiveVelocityMagnitudeShapeDerivative(
        const double EffectiveVelocityMagnitude,
        const array_1d<double, 3>& rEffectiveVelocity,
        const array_1d<double, 3>& rEffectiveVelocityDerivative)
    {
        if (EffectiveVelocityMagnitude > 0.0)
        {
            return inner_prod(rEffectiveVelocity, rEffectiveVelocityDerivative) /
                   EffectiveVelocityMagnitude;
        }
        else
        {
            return 0.0;
        }
    }

    static void CalculateAbsoluteResidualUnrelatedScalarDerivatives(
        BoundedVectorN& rOutput,
        const double Residual,
        const double ScalarValue,
        const BoundedMatrixND& rEffectiveVelocityDerivative,
        const BoundedArray3D& rScalarGradient,
        const BoundedVectorN& rReactionTermDerivative,
        const BoundedVectorN& rSourceTermDerivative)
    {
        const double coeff = (Residual > 0.0) ? 1.0 : -1.0;

        noalias(rOutput) = MatrixVectorProduct<TNumNodes>(
            rEffectiveVelocityDerivative, rScalarGradient);
        noalias(rOutput) += rReactionTermDerivative * ScalarValue;
        noalias(rOutput) -= rSourceTermDerivative;
        noalias(rOutput) = rOutput * coeff;
    }

    static void CalculateAbsoluteResidualRelatedScalarDerivativeAdditionalTerms(
        BoundedVectorN& rOutput,
        const double Residual,
        const double ReactionTerm,
        const BoundedVectorN& rRelaxedAccelerationDerivative,
        const BoundedArray3D& rEffectiveVelocity,
        const Vector& rShapeFunctionValues,
        const Matrix& rShapeFunctionDerivatives)
    {
        const double coeff = (Residual > 0.0) ? 1.0 : -1.0;

        noalias(rOutput) = rRelaxedAccelerationDerivative;
        noalias(rOutput) += MatrixVectorProduct<TNumNodes>(
            rShapeFunctionDerivatives, rEffectiveVelocity);
        noalias(rOutput) += rShapeFunctionValues * ReactionTerm;
        noalias(rOutput) = rOutput * coeff;
    }

    static void CalculateAbsoluteResidualVelocityDerivativeTerms(
        BoundedMatrixND& rOutput,
        const double Residual,
        const double ScalarValue,
        const BoundedArray3D& rScalarGradient,
        const BoundedMatrixND& rRelaxedAccelerationDerivative,
        const BoundedMatrixNDD& rEffectiveVelocityDerivative,
        const BoundedMatrixND& rReactionTermDerivative,
        const BoundedMatrixND& rSourceTermDerivative)
    {
        const double coeff = (Residual > 0.0) ? 1.0 : -1.0;
        noalias(rOutput) = rRelaxedAccelerationDerivative;
        noalias(rOutput) += rReactionTermDerivative * ScalarValue;
        noalias(rOutput) -= rSourceTermDerivative;

        const BoundedVectorND& r_effective_velocity_dot_scalar_gradient =
            MatrixVectorProduct<TNumNodes * TDim>(rEffectiveVelocityDerivative, rScalarGradient);
        noalias(rOutput) +=
            ConvertNDVectorToNDMatrix(r_effective_velocity_dot_scalar_gradient);

        noalias(rOutput) = rOutput * coeff;
    }

    static double CalculateAbsoluteResidualShapeDerivative(
        const double Residual,
        const double ScalarValue,
        const double RelaxedAccelerationDerivative,
        const double ReactionTermDerivative,
        const double SourceTermDerivative,
        const double EffectiveVelocityDerivativeDotScalarVariableGradient,
        const double EffectiveVelocityDotScalarVariableGradientDerivative)
    {
        const double coeff = (Residual > 0.0) ? 1.0 : -1.0;

        double value = 0.0;

        value += RelaxedAccelerationDerivative;
        value += EffectiveVelocityDerivativeDotScalarVariableGradient;
        value += EffectiveVelocityDotScalarVariableGradientDerivative;
        value += ReactionTermDerivative * ScalarValue;
        value -= SourceTermDerivative;

        return value * coeff;
    }

    static void CalculateTauScalarDerivatives(BoundedVectorN& rOutput,
                                              const double Tau,
                                              const double EffectiveKinematicViscosity,
                                              const double ElementLength,
                                              const double ReactionTerm,
                                              const BoundedVectorN& EffectiveKinematicViscosityDerivative,
                                              const BoundedVectorN& ReactionTermDerivative)
    {
        noalias(rOutput) =
            EffectiveKinematicViscosityDerivative *
            (144 * EffectiveKinematicViscosity / std::pow(ElementLength, 4));
        noalias(rOutput) += ReactionTermDerivative * ReactionTerm;
        noalias(rOutput) = rOutput * (-1.0 * std::pow(Tau, 3));
    }

    static void CalculateTauVelocityDerivatives(
        BoundedMatrixND& rOutput,
        const double Tau,
        const double EffectiveVelocityMagnitude,
        const double EffectiveKinematicViscosity,
        const double ElementLength,
        const double ReactionTerm,
        const BoundedMatrixND& rEffectiveVelocityMagnitudeDerivative,
        const BoundedMatrixND& rEffectiveKinematicViscosityDerivative,
        const BoundedMatrixND& rReactionTermDerivative)
    {
        noalias(rOutput) = rEffectiveVelocityMagnitudeDerivative *
                           (4 * EffectiveVelocityMagnitude / std::pow(ElementLength, 2));
        noalias(rOutput) +=
            rEffectiveKinematicViscosityDerivative *
            (144 * EffectiveKinematicViscosity / std::pow(ElementLength, 4));
        noalias(rOutput) += rReactionTermDerivative * ReactionTerm;
        noalias(rOutput) = rOutput * (-1.0 * std::pow(Tau, 3));
    }

    static double CalculateTauShapeDerivatives(const double Tau,
                                               const double EffectiveVelocityMagnitude,
                                               const double EffectiveKinematicViscosity,
                                               const double ElementLength,
                                               const double ReactionTerm,
                                               const double EffectiveVelocityMagnitudeDerivative,
                                               const double EffectiveKinematicViscosityDerivative,
                                               const double ReactionTermDerivative,
                                               const double DetJDerivative)

    {
        const double element_length_derivative =
            CalculateElementLengthShapeDerivative(ElementLength, DetJDerivative);

        double shape_sensitivity = 0.0;

        shape_sensitivity +=
            (4 * EffectiveVelocityMagnitude / std::pow(ElementLength, 2)) *
            (EffectiveVelocityMagnitudeDerivative -
             EffectiveVelocityMagnitude * element_length_derivative / ElementLength);
        shape_sensitivity +=
            (144 * EffectiveKinematicViscosity / std::pow(ElementLength, 4)) *
            (EffectiveKinematicViscosityDerivative -
             2 * EffectiveKinematicViscosity * element_length_derivative / ElementLength);
        shape_sensitivity += ReactionTerm * ReactionTermDerivative;
        shape_sensitivity *= -1.0 * std::pow(Tau, 3);

        return shape_sensitivity;
    }

    static void AddRFCBetaUnrelatedScalarDerivatives(BoundedVectorN& rOutput,
                                                     const double ScalarValue,
                                                     const double AbsoluteResidual,
                                                     const double Tau,
                                                     const BoundedVectorN& rAbsoluteResidualUnrelatedDerivatives,
                                                     const BoundedVectorN& rTauDerivatives)
    {
        if (ScalarValue > 0.0)
        {
            noalias(rOutput) += rAbsoluteResidualUnrelatedDerivatives * Tau;
            noalias(rOutput) += rTauDerivatives * AbsoluteResidual;
            noalias(rOutput) = rOutput * (1.0 / ScalarValue);
        }
    }

    static void AddRFCBetaRelatedScalarDerivativeAdditionalTerms(
        BoundedVectorN& rOutput,
        const double ScalarValue,
        const double AbsoluteResidual,
        const double Tau,
        const BoundedVectorN& rAbsoluteResidualRelatedDerivativeAdditionalTerms,
        const Vector& rShapeFunctionValues)
    {
        if (ScalarValue > 0.0)
        {
            noalias(rOutput) += rAbsoluteResidualRelatedDerivativeAdditionalTerms * Tau;
            noalias(rOutput) = rOutput * (1.0 / ScalarValue);
            noalias(rOutput) -= rShapeFunctionValues *
                                (AbsoluteResidual * Tau / std::pow(ScalarValue, 2));
        }
    }

    static void AddRFCBetaVelocityDerivative(BoundedMatrixND& rOutput,
                                             const double ScalarValue,
                                             const double AbsoluteResidual,
                                             const double Tau,
                                             const BoundedMatrixND& rAbsoluteResidualDerivative,
                                             const BoundedMatrixND& rTauDerivative)
    {
        if (ScalarValue > 0.0)
        {
            noalias(rOutput) += rAbsoluteResidualDerivative * Tau;
            noalias(rOutput) += rTauDerivative * AbsoluteResidual;
            noalias(rOutput) = rOutput * (1.0 / ScalarValue);
        }
    }

    static double CalculateRFCBetaShapeDerivative(const double ScalarValue,
                                                  const double AbsoluteResidual,
                                                  const double Tau,
                                                  const double AbsoluteResidualDerivative,
                                                  const double TauDerivative)
    {
        if (ScalarValue > 0.0)
        {
            double value = 0.0;

            value += AbsoluteResidualDerivative * Tau;
            value += TauDerivative * AbsoluteResidual;
            value /= ScalarValue;

            return value;
        }
        else
        {
            return 0.0;
        }
    }

    template <unsigned int TDerivativesSize>
    void static CalculateDiscreteUpwindOperatorResidualContributionDerivatives(
        BoundedMatrix<double, TDerivativesSize, TNumNodes>& rOutput,
        const BoundedVectorN& rNodalScalarValues,
        const BoundedMatrixNN& rInputMatrix,
        const BoundedVector<BoundedMatrix<double, TNumNodes, TNumNodes>, TDerivativesSize>& rInputMatrixDerivatives)
    {
        rOutput.clear();
        for (IndexType c = 0; c < TDerivativesSize; ++c)
        {
            const BoundedMatrixNN& derivatives_matrix = rInputMatrixDerivatives[c];
            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                for (IndexType b = 0; b < TNumNodes; ++b)
                {
                    if (a != b && (rInputMatrix(a, b) > 0.0 || rInputMatrix(b, a) > 0.0))
                    {
                        const double flux =
                            rNodalScalarValues[b] - rNodalScalarValues[a];
                        if (rInputMatrix(a, b) > rInputMatrix(b, a))
                        {
                            rOutput(c, a) -= derivatives_matrix(a, b) * flux;
                        }
                        else
                        {
                            rOutput(c, a) -= derivatives_matrix(b, a) * flux;
                        }
                    }
                }
            }
        }
    }

    template <unsigned int TDerivativesSize>
    void static CalculatePositivityPreservingMatrixResidualContributionDerivatives(
        BoundedMatrix<double, TDerivativesSize, TNumNodes>& rOutput,
        const double PositivityPreservingMatrixCoefficient,
        const BoundedVectorN& rNodalScalarValues,
        const BoundedMatrixNN& rInputMatrix,
        const BoundedVector<BoundedMatrix<double, TNumNodes, TNumNodes>, TDerivativesSize>& rInputMatrixDerivatives)
    {
        rOutput.clear();

        if (PositivityPreservingMatrixCoefficient > 0.0)
        {
            int row_index = -1;
            double max_negative_row_sum = 0.0;
            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                double row_sum = 0.0;
                for (IndexType b = 0; b < TNumNodes; ++b)
                {
                    row_sum -= rInputMatrix(a, b);
                }

                if (row_sum > max_negative_row_sum)
                {
                    max_negative_row_sum = row_sum;
                    row_index = a;
                }
            }

            for (IndexType c = 0; c < TDerivativesSize; ++c)
            {
                const BoundedMatrixNN& derivatives_matrix = rInputMatrixDerivatives[c];
                for (IndexType a = 0; a < TNumNodes; ++a)
                {
                    const double nodal_scalar_value = rNodalScalarValues[a];
                    for (IndexType b = 0; b < TNumNodes; ++b)
                    {
                        rOutput(c, a) -= derivatives_matrix(row_index, b) * nodal_scalar_value;
                    }
                }
            }
        }
    }
};
} // namespace ConvectionDiffusionReactionStabilizationUtilities
} // namespace Kratos
///@}

#endif