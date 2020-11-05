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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_STABILIZATION_ADJOINT_UTILITIES_H_INCLUDED)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_STABILIZATION_ADJOINT_UTILITIES_H_INCLUDED

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

template<unsigned int TDim, unsigned int TNumNodes>
class AdjointUtilities
{
public:
    ///@name Public type definitions
    ///@{

    using VectorN = BoundedVector<double, TNumNodes>;
    using MatrixN = BoundedMatrix<double, TNumNodes, TNumNodes>;
    using Array3D = array_1d<double, 3>;

    template<std::size_t TDerivativesSize>
    using VectorD = BoundedVector<double, TDerivativesSize>;

    template<std::size_t TDerivativesSize>
    using MatrixD = BoundedMatrix<double, TDerivativesSize, 3>;

    ///@}
    ///@name Public static operations
    ///@{

    template<std::size_t TSize1, std::size_t TSize2>
    static void DidacticProduct(
        BoundedMatrix<double, TSize1, TSize2>& rOutput,
        const Vector& rA,
        const Vector& rB)
    {
        KRATOS_TRY

        KRATOS_DEBUG_ERROR_IF(rA.size() != TSize1)
            << "Dimensions mismatch in didactic product [ TSize1 = " << TSize1
            << ", rA.size() = " << rA.size() << " ].\n";

        KRATOS_DEBUG_ERROR_IF(rB.size() != TSize2)
            << "Dimensions mismatch in didactic product [ TSize2 = " << TSize2
            << ", rB.size() = " << rB.size() << " ].\n";

        for (IndexType i = 0; i < TSize1; ++i) {
            for (IndexType j = 0; j < TSize2; ++j) {
                rOutput(i, j) = rA[i] * rB[j];
            }
        }

        KRATOS_CATCH("");
    }

    template<std::size_t TSize1, std::size_t TSize2>
    static void MatrixMatrixProduct(
        BoundedMatrix<double, TSize1, TSize2>& rOutput,
        const Matrix& rA,
        const Matrix& rB)
    {
        KRATOS_TRY

        rOutput.clear();

        KRATOS_DEBUG_ERROR_IF(rB.size2() != TSize2)
            << "Dimensions mismatch in matrix matrix product. [ TSize2 = " << TSize2
            << ", rB.size2() = " << rB.size2() << " ].\n";

        const auto a_2 = std::min(rA.size2(), rB.size1());

        for (IndexType i = 0; i < TSize1; ++i) {
            for (IndexType j = 0; j < TSize2; ++j) {
                for (IndexType k = 0; k < a_2; ++k) {
                    rOutput(i, j) += rA(i, k) * rB(k, j);
                }
            }
        }

        KRATOS_CATCH("");
    }

    template<std::size_t TSize>
    static void MatrixVectorProduct(
        BoundedVector<double, TSize>& rOutput,
        const Matrix& rA,
        const Vector& rB)
    {
        KRATOS_TRY

        rOutput.clear();

        KRATOS_DEBUG_ERROR_IF(rA.size1() != TSize)
            << "Dimensions mismatch in matrix vector product. [ TSize = " << TSize
            << ", rA.size1() = " << rA.size1() << " ].\n";

        const auto a_2 = std::min(rA.size2(), rB.size());

        for (IndexType i = 0; i < TSize; ++i) {
            for (IndexType k = 0; k < a_2; ++k) {
                rOutput[i] += rA(i, k) * rB[k];
            }
        }

        KRATOS_CATCH("");
    }

    static double InnerProduct(const Vector& rA, const Vector& rB)
    {
        const auto rows = std::min(rA.size(), rB.size());
        double result = 0.0;
        for (IndexType i = 0; i < rows; ++i) {
            result += rA[i] * rB[i];
        }
        return result;
    }

    template<std::size_t TDerivativesSize>
    static void CalculateEffectiveVelocityMagnitudeDerivative(
        VectorD<TDerivativesSize>& rOutput,
        const double EffectiveVelocityMagnitude,
        const Array3D& rEffectiveVelocity,
        const MatrixD<TDerivativesSize>& rEffectiveVelocityDerivatives)
    {
        rOutput.clear();

        if (EffectiveVelocityMagnitude > 0.0) {
            noalias(rOutput) = prod(rEffectiveVelocityDerivatives, rEffectiveVelocity);
            noalias(rOutput) = rOutput * (1 / EffectiveVelocityMagnitude);
        }
    }

    template<std::size_t TDerivativesSize>
    static void CalculateStabilizationTauDerivatives(
        VectorD<TDerivativesSize>& rOutput,
        const double EffectiveVelocityMagnitude,
        const double EffectiveKinematicViscosity,
        const double ReactionTerm,
        const double ElementLength,
        const double StabilizationTau,
        const VectorD<TDerivativesSize>& rEffectiveKinematicViscosityDerivatives,
        const VectorD<TDerivativesSize>& rReactionTermDerivatives,
        const VectorD<TDerivativesSize>& rEffectiveVelocityMagnitudeDerivatives)
    {
        noalias(rOutput) = rEffectiveVelocityMagnitudeDerivatives * (4 * EffectiveVelocityMagnitude / std::pow(ElementLength, 2));
        noalias(rOutput) += rEffectiveKinematicViscosityDerivatives * (144 * EffectiveKinematicViscosity / std::pow(ElementLength, 4));
        noalias(rOutput) += rReactionTermDerivatives * ReactionTerm;
        noalias(rOutput) = rOutput * (-1.0 * std::pow(StabilizationTau, 3));
    }


    template<std::size_t TDerivativesSize>
    static typename std::enable_if<(TNumNodes == TDerivativesSize), void>::type
    CalculateResidualDerivatives(
        VectorD<TDerivativesSize>& rOutput,
        const double SelfWeight,
        const double ScalarValue,
        const double ReactionTerm,
        const Array3D& rEffectiveVelocity,
        const Array3D& rScalarGradient,
        const Vector& rGPShapeFunctions,
        const Matrix& rGPShapeFunctionDerivatives,
        const VectorD<TDerivativesSize>& rReactionTermDerivatives,
        const VectorD<TDerivativesSize>& rSourceTermDerivatives,
        const MatrixD<TDerivativesSize>& rEffectiveVelocityDerivative)
    {
        MatrixVectorProduct(rOutput, rGPShapeFunctionDerivatives, rEffectiveVelocity * SelfWeight);
        noalias(rOutput) += prod(rEffectiveVelocityDerivative, rScalarGradient);
        noalias(rOutput) += rReactionTermDerivatives * ScalarValue;
        noalias(rOutput) += rGPShapeFunctions * (ReactionTerm * SelfWeight);
        noalias(rOutput) -= rSourceTermDerivatives;
    }

    template<std::size_t TDerivativesSize>
    static typename std::enable_if<(TNumNodes != TDerivativesSize), void>::type
    CalculateResidualDerivatives(
        VectorD<TDerivativesSize>& rOutput,
        const double SelfWeight,
        const double ScalarValue,
        const double ReactionTerm,
        const Array3D& rEffectiveVelocity,
        const Array3D& rScalarGradient,
        const Vector& rGPShapeFunctions,
        const Matrix& rGPShapeFunctionDerivatives,
        const VectorD<TDerivativesSize>& rReactionTermDerivatives,
        const VectorD<TDerivativesSize>& rSourceTermDerivatives,
        const MatrixD<TDerivativesSize>& rEffectiveVelocityDerivative)
    {
        noalias(rOutput) = prod(rEffectiveVelocityDerivative, rScalarGradient);
        noalias(rOutput) += rReactionTermDerivatives * ScalarValue;
        noalias(rOutput) -= rSourceTermDerivatives;
    }

    template<std::size_t TDerivativesSize>
    static void CalculateAbsoluteScalarDerivatives(
        VectorD<TDerivativesSize>& rOutput,
        const double ScalarValue,
        const VectorD<TDerivativesSize>& rScalarDerivatives)
    {
        noalias(rOutput) = rScalarDerivatives * ((ScalarValue >= 0.0) ? 1.0 : -1.0);
    }

    template <std::size_t TDerivativesSize>
    void static CalculateDiscreteUpwindOperatorResidualContributionDerivatives(
        BoundedMatrix<double, TDerivativesSize, TNumNodes>& rOutput,
        const BoundedVector<double, TNumNodes>& rNodalScalarValues,
        const BoundedMatrix<double, TNumNodes, TNumNodes>& rInputMatrix,
        const BoundedVector<BoundedMatrix<double, TNumNodes, TNumNodes>, TDerivativesSize>& rInputMatrixDerivatives)
    {
        rOutput.clear();
        for (IndexType c = 0; c < TDerivativesSize; ++c) {
            const BoundedMatrix<double, TNumNodes, TNumNodes>& derivatives_matrix = rInputMatrixDerivatives[c];
            for (IndexType a = 0; a < TNumNodes; ++a) {
                for (IndexType b = 0; b < TNumNodes; ++b) {
                    if (a != b && (rInputMatrix(a, b) > 0.0 || rInputMatrix(b, a) > 0.0)) {
                        const double flux = rNodalScalarValues[b] - rNodalScalarValues[a];
                        if (rInputMatrix(a, b) > rInputMatrix(b, a)) {
                            rOutput(c, a) -= derivatives_matrix(a, b) * flux;
                        } else {
                            rOutput(c, a) -= derivatives_matrix(b, a) * flux;
                        }
                    }
                }
            }
        }
    }

    template <std::size_t TDerivativesSize>
    void static CalculatePositivityPreservingCoefficientDerivatives(
        BoundedVector<double, TDerivativesSize>& rOutput,
        const double PositivityPreservingMatrixCoefficient,
        const BoundedMatrix<double, TNumNodes, TNumNodes>& rInputMatrix,
        const BoundedVector<BoundedMatrix<double, TNumNodes, TNumNodes>, TDerivativesSize>& rInputMatrixDerivatives)
    {
        rOutput.clear();

        if (PositivityPreservingMatrixCoefficient > 0.0) {
            int row_index = -1;
            double max_negative_row_sum = 0.0;
            for (IndexType a = 0; a < TNumNodes; ++a) {
                double row_sum = 0.0;
                for (IndexType b = 0; b < TNumNodes; ++b) {
                    row_sum += rInputMatrix(a, b);
                }

                if (row_sum < max_negative_row_sum) {
                    max_negative_row_sum = row_sum;
                    row_index = a;
                }
            }

            for (IndexType c = 0; c < TDerivativesSize; ++c) {
                const BoundedMatrix<double, TNumNodes, TNumNodes>& derivatives_matrix = rInputMatrixDerivatives[c];
                for (IndexType a = 0; a < TNumNodes; ++a) {
                    rOutput[c] -= derivatives_matrix(row_index, a);
                }
            }
        }
    }

    ///@}
};

} // namespace ConvectionDiffusionReactionStabilizationUtilities
} // namespace Kratos

#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_STABILIZATION_ADJOINT_UTILITIES_H_INCLUDED