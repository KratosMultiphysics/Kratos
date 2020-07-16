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

    static BoundedVectorN MatrixVectorProduct(const BoundedMatrixND& rM,
                                              const BoundedArray3D& rV)
    {
        BoundedVectorN output;
        for (IndexType i = 0; i < TNumNodes; ++i)
        {
            output[i] = rM(i, 0) * rV[0];
            for (IndexType j = 1; j < TDim; ++j)
            {
                output[i] += rM(i, j) * rV[j];
            }
        }

        return output;
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

        noalias(rOutput) =
            MatrixVectorProduct(rEffectiveVelocityDerivative, rScalarGradient);
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
        noalias(rOutput) +=
            MatrixVectorProduct(rShapeFunctionDerivatives, rEffectiveVelocity);
        noalias(rOutput) += rShapeFunctionValues * ReactionTerm;
        noalias(rOutput) = rOutput * coeff;
    }

    static void CalculateAbsoluteResidualVelocityDerivativeTerms(
        BoundedVectorND& rOutput,
        const double Residual,
        const double ScalarValue,
        const BoundedArray3D& rScalarGradient,
        const BoundedMatrixNDD& rRelaxedAccelerationDerivative,
        const BoundedMatrixNDD& rEffectiveVelocityDerivative,
        const BoundedMatrixND& rReactionTermDerivative,
        const BoundedMatrixND& rSourceTermDerivative)
    {
        const double coeff = (Residual > 0.0) ? 1.0 : -1.0;
        noalias(rOutput) = rRelaxedAccelerationDerivative;
        noalias(rOutput) += rEffectiveVelocityDerivative;
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

    static void CalculateRFCBetaUnrelatedScalarDerivatives(BoundedVectorN& rOutput,
                                                           const double ScalarValue,
                                                           const double AbsoluteResidual,
                                                           const double Tau,
                                                           const BoundedVectorN& rAbsoluteResidualUnrelatedDerivatives,
                                                           const BoundedVectorN& rTauDerivatives)
    {
        if (ScalarValue > 0.0)
        {
            noalias(rOutput) = rAbsoluteResidualUnrelatedDerivatives * Tau;
            noalias(rOutput) += rTauDerivatives * AbsoluteResidual;
            noalias(rOutput) = rOutput * (1.0 / ScalarValue);
        }
        else
        {
            rOutput.clear();
        }
    }

    static void CalculateRFCBetaRelatedScalarDerivativeAdditionalTerms(
        BoundedVectorN& rOutput,
        const double ScalarValue,
        const double AbsoluteResidual,
        const double Tau,
        const BoundedVectorN& rAbsoluteResidualRelatedDerivativeAdditionalTerms,
        const Vector& rShapeFunctionValues)
    {
        if (ScalarValue > 0.0)
        {
            noalias(rOutput) = rAbsoluteResidualRelatedDerivativeAdditionalTerms * Tau;
            noalias(rOutput) = rOutput * (1.0 / ScalarValue);
            noalias(rOutput) -= rShapeFunctionValues *
                                (AbsoluteResidual * Tau / std::pow(ScalarValue, 2));
        }
        else
        {
            rOutput.clear();
        }
    }

    void static CalculateDiscreteUpwindOperatorResidualContributionScalarDerivatives(
        BoundedMatrixNN& rOutput,
        const BoundedVectorN& rNodalScalarValues,
        const BoundedMatrixNN& rInputMatrix,
        const BoundedMatrixNNN& rInputMatrixDerivatives)
    {
        rOutput.clear();
        for (IndexType c = 0; c < TNumNodes; ++c)
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

    void static CalculatePositivityPreservingMatrixResidualContributionScalarDerivatives(
        BoundedMatrixNN& rOutput,
        const double PositivityPreservingMatrixCoefficient,
        const BoundedVectorN& rNodalScalarValues,
        const BoundedMatrixNN& rInputMatrix,
        const BoundedMatrixNNN& rInputMatrixDerivatives)
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

            for (IndexType c = 0; c < TNumNodes; ++c)
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