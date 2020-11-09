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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_ADJOINT_SENSITIVITY_DERIVATIVES_H_INCLUDED)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_ADJOINT_SENSITIVITY_DERIVATIVES_H_INCLUDED

// System includes
#include <vector>

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/time_discretization.h"
#include "utilities/geometrical_sensitivity_utility.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_elements/convection_diffusion_reaction_stabilization_adjoint_utilities.h"
#include "custom_elements/convection_diffusion_reaction_stabilization_utilities.h"

namespace Kratos
{
namespace ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointUtilities
{
template <unsigned int TDim, unsigned int TNumNodes>
class StabilizationSensitivityDerivatives
{
public:
    ///@name Public type definitions
    ///@{

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    template<std::size_t TSize>
    using BVector = BoundedVector<double, TSize>;

    template<std::size_t TSize1, std::size_t TSize2>
    using BMatrix = BoundedMatrix<double, TSize1, TSize2>;

    ///@}
    ///@name Public classes
    ///@{

    template <class TDerivativesType, unsigned int TEquationOffset>
    class ShapeDerivatives
    {
    public:
        ///@name Public type definitions
        ///@{

        static constexpr unsigned int TDerivativesSize = TDerivativesType::TDerivativesSize;

        using VectorDerivativesType = BVector<TDerivativesSize>;

        template<unsigned int TSize>
        using MatrixDerivativesType = BMatrix<TDerivativesSize, TSize>;

        using AdjointUtilities = ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<TDim, TNumNodes>;

        ///@}
        ///@name Life cycle
        ///@{

        ShapeDerivatives(
            const GeometryType& rGeometry)
            : mrData(rGeometry)
        {
        }

        ///@}
        ///@name Public operations
        ///@{

        void Initialize(
            Matrix& rOutput,
            const ProcessInfo& rProcessInfo)
        {
            mrData.CalculateConstants(rProcessInfo);

            for (IndexType c = 0; c < TDerivativesSize; ++c) {
                mPrimalDampingMatrixDerivatives[c].clear();
            }

            for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
                const auto& r_node = mrData.GetGeometry()[i_node];
                mPrimalVariableValues[i_node] = r_node.FastGetSolutionStepValue(TDerivativesType::GetScalarVariable());
                mRelaxedPrimalVariableValues[i_node] = r_node.FastGetSolutionStepValue(TDerivativesType::GetScalarVariable().GetTimeDerivative().GetTimeDerivative());
            }

            mScalarMultiplierDerivatives.clear();
            mDerivativeBlockSize = rOutput.size2() / TNumNodes;
            mDeltaTime = rProcessInfo[DELTA_TIME] * -1.0;
            mBossakAlpha = rProcessInfo[BOSSAK_ALPHA];
            mBossakGamma = TimeDiscretization::Bossak(mBossakAlpha, 0.25, 0.5).GetGamma();
            mDynamicTau = rProcessInfo[DYNAMIC_TAU];
            mElementLength = mrData.GetGeometry().Length();
            mStabilizationDiscreteDiffusionMatrixUserCoefficient = rProcessInfo[RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT];
            mStabilizationPositivityPreservingMatrixUserCoefficient = rProcessInfo[RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT];
            mNumberOfGaussPoints = 0.0;
        }

        void CalculateGaussPointData(
            const double GPWeight,
            const Vector& rGPShapeFunctions,
            const Matrix& rGPShapeFunctionDerivatives)
        {
            using namespace RansCalculationUtilities;
            mrData.CalculateGaussPointData(rGPShapeFunctions, rGPShapeFunctionDerivatives);

            mEffectiveVelocity = mrData.CalculateEffectiveVelocity(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            mEffectiveKinematicViscosity = mrData.CalculateEffectiveKinematicViscosity(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            mReactionTerm = mrData.CalculateReactionTerm(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            mSourceTerm = mrData.CalculateSourceTerm(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            mAbsoluteReactionTerm = std::abs(mReactionTerm);

            const auto& r_scalar_variable = TDerivativesType::GetScalarVariable();

            EvaluateInPoint(
                mrData.GetGeometry(), rGPShapeFunctions,
                std::tie(mScalarVariableValue, r_scalar_variable),
                std::tie(mRelaxedScalarRateVariableValue, r_scalar_variable.GetTimeDerivative().GetTimeDerivative()));

            CalculateGradient(mScalarVariableGradient, mrData.GetGeometry(),
                              r_scalar_variable, rGPShapeFunctionDerivatives);

            mVelocityDotScalarVariableGradient = inner_prod(mScalarVariableGradient, mEffectiveVelocity);
            mEffectiveVelocityMagnitude = norm_2(mEffectiveVelocity);

            mStabilizationTau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
                mElementLength, mEffectiveVelocityMagnitude, mReactionTerm,
                mEffectiveKinematicViscosity, mBossakAlpha, mBossakGamma, mDeltaTime, mDynamicTau);

            mResidual = mRelaxedScalarRateVariableValue;
            mResidual += mVelocityDotScalarVariableGradient;
            mResidual += mReactionTerm * mScalarVariableValue;
            mResidual -= mSourceTerm;
            mAbsoluteResidual = std::abs(mResidual);

            if (mScalarVariableValue > 0.0) {
                mScalarMultiplier += mAbsoluteResidual * mStabilizationTau / mScalarVariableValue;
            }

            GetConvectionOperator(mVelocityConvectiveTerms, mEffectiveVelocity, rGPShapeFunctionDerivatives);

            ConvectionDiffusionReactionStabilizationUtilities::AddPrimalDampingMatrixGaussPointContributions(
                mPrimalDampingMatrix, mEffectiveKinematicViscosity, mReactionTerm,
                mAbsoluteReactionTerm, mStabilizationTau, GPWeight, mVelocityConvectiveTerms,
                rGPShapeFunctions, rGPShapeFunctionDerivatives);

            noalias(mShapeDerivativesInnerProduct) = prod(rGPShapeFunctionDerivatives, trans(rGPShapeFunctionDerivatives));

            mNumberOfGaussPoints += 1.0;
        }

        void CalculateResidualDerivatives(
            Matrix& rMatrix,
            const ShapeParameter& rShapeParameters,
            const double GPWeight,
            const Vector& rGPShapeFunctions,
            const Matrix& rGPShapeFunctionDerivatives,
            const double GPWeightDerivative,
            const double DetJDerivative,
            const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rGPShapeFunctionDerivativesDerivative)
        {
            KRATOS_TRY

            using namespace RansCalculationUtilities;

            const auto& r_scalar_variable = TDerivativesType::GetScalarVariable();
            const unsigned int deriv_index = rShapeParameters.NodeIndex * TDim + rShapeParameters.Direction;

            const auto& effective_velocity_derivative = mrData.CalculateEffectiveVelocityDerivatives(rShapeParameters, rGPShapeFunctions, rGPShapeFunctionDerivatives, DetJDerivative, rGPShapeFunctionDerivativesDerivative);
            const double effective_kinematic_viscosity_derivative = mrData.CalculateEffectiveKinematicViscosityDerivatives(rShapeParameters, rGPShapeFunctions, rGPShapeFunctionDerivatives, DetJDerivative, rGPShapeFunctionDerivativesDerivative);
            const double reaction_term_derivative = mrData.CalculateReactionTermDerivatives(rShapeParameters, rGPShapeFunctions, rGPShapeFunctionDerivatives, DetJDerivative, rGPShapeFunctionDerivativesDerivative);
            const double absolute_reaction_term_derivative = (mReactionTerm >= 0.0 ? reaction_term_derivative : reaction_term_derivative * -1.0);
            const double source_term_derivative = mrData.CalculateSourceTermDerivatives(rShapeParameters, rGPShapeFunctions, rGPShapeFunctionDerivatives, DetJDerivative, rGPShapeFunctionDerivativesDerivative);

            const double effective_velocity_magnitude_derivative =
                AdjointUtilities::CalculateEffectiveVelocityMagnitudeDerivative(
                    mEffectiveVelocityMagnitude, mEffectiveVelocity, effective_velocity_derivative);

            const double tau_derivative = AdjointUtilities::CalculateTauShapeDerivatives(
                mStabilizationTau, mEffectiveVelocityMagnitude, mEffectiveKinematicViscosity,
                mElementLength, mReactionTerm, effective_velocity_magnitude_derivative,
                effective_kinematic_viscosity_derivative,
                reaction_term_derivative, DetJDerivative);

            BVector<TNumNodes> velocity_convective_terms_derivative;
            GetConvectionOperator(velocity_convective_terms_derivative, mEffectiveVelocity, rGPShapeFunctionDerivativesDerivative);

            BVector<TNumNodes> velocity_derivative_convective_terms;
            GetConvectionOperator(velocity_derivative_convective_terms, effective_velocity_derivative, rGPShapeFunctionDerivatives);

            array_1d<double, 3> scalar_variable_gradient_derivative;
            CalculateGradient(scalar_variable_gradient_derivative,
                              mrData.GetGeometry(), r_scalar_variable,
                              rGPShapeFunctionDerivativesDerivative);

            const double velocity_dot_scalar_variable_gradient_derivative = inner_prod(mEffectiveVelocity, scalar_variable_gradient_derivative);
            const double velocity_derivative_dot_scalar_variable_gradient = inner_prod(effective_velocity_derivative, mScalarVariableGradient);

            BVector<TNumNodes> dN_dX_dot_scalar_variable_gradient_derivative;
            AdjointUtilities::MatrixVectorProduct(dN_dX_dot_scalar_variable_gradient_derivative, rGPShapeFunctionDerivatives, scalar_variable_gradient_derivative);
            BVector<TNumNodes> dN_dX_derivative_dot_scalar_variable_gradient;
            AdjointUtilities::MatrixVectorProduct(dN_dX_derivative_dot_scalar_variable_gradient, rGPShapeFunctionDerivativesDerivative, mScalarVariableGradient);
            BVector<TNumNodes> dN_dX_dot_scalar_variable_gradeint;
            AdjointUtilities::MatrixVectorProduct(dN_dX_dot_scalar_variable_gradeint, rGPShapeFunctionDerivatives, mScalarVariableGradient);

            BMatrix<TNumNodes, TNumNodes> dNa_dNb_derivative = prod(rGPShapeFunctionDerivatives, trans(rGPShapeFunctionDerivativesDerivative));
            BMatrix<TNumNodes, TNumNodes>& r_primal_damping_matrix_derivatives = mPrimalDampingMatrixDerivatives[deriv_index];

            for (IndexType a = 0; a < TNumNodes; ++a)
            {
                const double tau_operator = mVelocityConvectiveTerms[a] + mAbsoluteReactionTerm * rGPShapeFunctions[a];
                const double tau_operator_derivative = velocity_convective_terms_derivative[a] + velocity_derivative_convective_terms[a] + absolute_reaction_term_derivative * rGPShapeFunctions[a];

                double value = 0.0;

                // adding RHS derivative contributions
                value += rGPShapeFunctions[a] * source_term_derivative * GPWeight;
                value += rGPShapeFunctions[a] * mSourceTerm * GPWeightDerivative;

                // adding RHS SUPG derivative contributions
                value += tau_derivative * tau_operator * mSourceTerm * GPWeight;
                value += mStabilizationTau * tau_operator_derivative * mSourceTerm * GPWeight;
                value += mStabilizationTau * tau_operator * source_term_derivative * GPWeight;
                value += mStabilizationTau * tau_operator * mSourceTerm * GPWeightDerivative;

                // adding LHS derivative contributions
                value -= rGPShapeFunctions[a] * velocity_dot_scalar_variable_gradient_derivative * GPWeight;
                value -= rGPShapeFunctions[a] * velocity_derivative_dot_scalar_variable_gradient * GPWeight;
                value -= rGPShapeFunctions[a] * mVelocityDotScalarVariableGradient * GPWeightDerivative;

                value -= rGPShapeFunctions[a] * reaction_term_derivative * mScalarVariableValue * GPWeight;
                value -= rGPShapeFunctions[a] * mReactionTerm * mScalarVariableValue * GPWeightDerivative;

                value -= effective_kinematic_viscosity_derivative * dN_dX_dot_scalar_variable_gradeint[a] * GPWeight;
                value -= mEffectiveKinematicViscosity * dN_dX_derivative_dot_scalar_variable_gradient[a] * GPWeight;
                value -= mEffectiveKinematicViscosity * dN_dX_dot_scalar_variable_gradient_derivative[a] * GPWeight;
                value -= mEffectiveKinematicViscosity * dN_dX_dot_scalar_variable_gradeint[a] * GPWeightDerivative;

                // adding LHS SUPG derivative contributions
                value -= tau_derivative * tau_operator * mVelocityDotScalarVariableGradient * GPWeight;
                value -= mStabilizationTau * tau_operator_derivative * mVelocityDotScalarVariableGradient * GPWeight;
                value -= mStabilizationTau * tau_operator * velocity_dot_scalar_variable_gradient_derivative * GPWeight;
                value -= mStabilizationTau * tau_operator * velocity_derivative_dot_scalar_variable_gradient * GPWeight;
                value -= mStabilizationTau * tau_operator * mVelocityDotScalarVariableGradient * GPWeightDerivative;

                value -= tau_derivative * tau_operator * mReactionTerm * mScalarVariableValue * GPWeight;
                value -= mStabilizationTau * tau_operator_derivative * mReactionTerm * mScalarVariableValue * GPWeight;
                value -= mStabilizationTau * tau_operator * reaction_term_derivative * mScalarVariableValue * GPWeight;
                value -= mStabilizationTau * tau_operator * mReactionTerm * mScalarVariableValue * GPWeightDerivative;

                // adding Mass derivative contrirbutions
                value -= GPWeightDerivative / TNumNodes * mRelaxedPrimalVariableValues[a];

                value -= tau_derivative * tau_operator * mRelaxedScalarRateVariableValue * GPWeight;
                value -= mStabilizationTau * tau_operator_derivative * mRelaxedScalarRateVariableValue * GPWeight;
                value -= mStabilizationTau * tau_operator * mRelaxedScalarRateVariableValue * GPWeightDerivative;

                rMatrix(deriv_index, a * mDerivativeBlockSize + TEquationOffset) += value;

                // compute primal damping matrix derivatives
                for (IndexType b = 0; b < TNumNodes; ++b) {
                    double value = 0.0;

                    value += rGPShapeFunctions[a] * mVelocityConvectiveTerms[b] * GPWeightDerivative;
                    value += rGPShapeFunctions[a] * velocity_convective_terms_derivative[b] * GPWeight;
                    value += rGPShapeFunctions[a] * velocity_derivative_convective_terms[b] * GPWeight;

                    value += rGPShapeFunctions[a] * mReactionTerm * rGPShapeFunctions[b] * GPWeightDerivative;
                    value += rGPShapeFunctions[a] * reaction_term_derivative * rGPShapeFunctions[b] * GPWeight;

                    value += mEffectiveKinematicViscosity * mShapeDerivativesInnerProduct(a, b) * GPWeightDerivative;
                    value += effective_kinematic_viscosity_derivative * mShapeDerivativesInnerProduct(a, b) * GPWeight;
                    value += mEffectiveKinematicViscosity * dNa_dNb_derivative(a, b) * GPWeight;
                    value += mEffectiveKinematicViscosity * dNa_dNb_derivative(b, a) * GPWeight;

                    value += mStabilizationTau * tau_operator * mVelocityConvectiveTerms[b] * GPWeightDerivative;
                    value += tau_derivative * tau_operator * mVelocityConvectiveTerms[b] * GPWeight;
                    value += mStabilizationTau * tau_operator_derivative * mVelocityConvectiveTerms[b] * GPWeight;
                    value += mStabilizationTau * tau_operator * velocity_convective_terms_derivative[b] * GPWeight;
                    value += mStabilizationTau * tau_operator * velocity_derivative_convective_terms[b] * GPWeight;

                    value += mStabilizationTau * tau_operator * mReactionTerm * rGPShapeFunctions[b] * GPWeightDerivative;
                    value += tau_derivative * tau_operator * mReactionTerm * rGPShapeFunctions[b] * GPWeight;
                    value += mStabilizationTau * tau_operator_derivative * mReactionTerm * rGPShapeFunctions[b] * GPWeight;
                    value += mStabilizationTau * tau_operator * reaction_term_derivative * rGPShapeFunctions[b] * GPWeight;

                    r_primal_damping_matrix_derivatives(a, b) += value;
                }
            }

            // adding scalar multiplier derivatives
            if (mScalarVariableValue > 0.0) {
                double residual_derivative = velocity_dot_scalar_variable_gradient_derivative;
                residual_derivative += velocity_derivative_dot_scalar_variable_gradient;
                residual_derivative += reaction_term_derivative * mScalarVariableValue;
                residual_derivative -= source_term_derivative;

                const double absolute_residual_derivative =
                    (mResidual > 0.0 ? residual_derivative : residual_derivative * -1.0);

                mScalarMultiplierDerivatives[deriv_index] +=
                    (absolute_residual_derivative * mStabilizationTau +
                     mAbsoluteResidual * tau_derivative) /
                    mScalarVariableValue;
            }

            KRATOS_CATCH("");
        }

        void Finalize(Matrix& rOutput, const ProcessInfo& rProcessInfo)
        {
            // compute primal quantities
            mScalarMultiplier /= mNumberOfGaussPoints;

            // compute discrete diffusion matrix
            BoundedMatrix<double, TNumNodes, TNumNodes> discrete_diffusion_matrix;
            double discrete_diffusion_matrix_coefficient;
            ConvectionDiffusionReactionStabilizationUtilities::CalculateDiscreteUpwindOperator<TNumNodes>(
                discrete_diffusion_matrix_coefficient,
                discrete_diffusion_matrix, mPrimalDampingMatrix);
            const BVector<TNumNodes>& discrete_diffusion_residual_contribution =
                prod(discrete_diffusion_matrix, mPrimalVariableValues);

            // compute positivity preserving coefficient
            mPositivityPreservingMatrixCoefficient =
                ConvectionDiffusionReactionStabilizationUtilities::CalculatePositivityPreservingMatrix(
                    mPrimalDampingMatrix);

            // compute adjoint quantities
            noalias(mScalarMultiplierDerivatives) = mScalarMultiplierDerivatives / mNumberOfGaussPoints;

            BMatrix<TDerivativesSize, TNumNodes> discrete_diffusion_matrix_residual_derivatives;
            AdjointUtilities::CalculateStabilizationDiscreteUpwindMatrixResidualFristDerivatives(
                discrete_diffusion_matrix_residual_derivatives, mScalarMultiplier,
                mStabilizationDiscreteDiffusionMatrixUserCoefficient, mScalarMultiplierDerivatives,
                mPrimalVariableValues, discrete_diffusion_residual_contribution,
                mPrimalDampingMatrix, mPrimalDampingMatrixDerivatives);

            BMatrix<TDerivativesSize, TNumNodes> positivity_preserving_matrix_residual_derivatives;
            AdjointUtilities::CalculateStabilizationPositivityPreservingMatrixResidualFristDerivatives(
                positivity_preserving_matrix_residual_derivatives,
                mScalarMultiplier, mPositivityPreservingMatrixCoefficient,
                mStabilizationPositivityPreservingMatrixUserCoefficient,
                mScalarMultiplierDerivatives, mPrimalVariableValues,
                mPrimalDampingMatrix, mPrimalDampingMatrixDerivatives);

            for (IndexType c = 0; c < TNumNodes; ++c) {
                for (IndexType k = 0; k < TDim; ++k) {
                    for (IndexType a = 0; a < TNumNodes; ++a) {
                        double& value = rOutput(
                            c * TDim + k, a * mDerivativeBlockSize + TEquationOffset);
                        value -= discrete_diffusion_matrix_residual_derivatives(
                            c * TDim + k, a);
                        value -= positivity_preserving_matrix_residual_derivatives(
                            c * TDim + k, a);
                    }
                }
            }
        }

        ///@}
    private:
        ///@name Private members
        ///@{

        TDerivativesType mrData;
        IndexType mDerivativeBlockSize;
        BoundedVector<BMatrix<TNumNodes, TNumNodes>, TDerivativesSize> mPrimalDampingMatrixDerivatives;
        VectorDerivativesType mScalarMultiplierDerivatives;

        // primal data storage
        array_1d<double, 3> mEffectiveVelocity;
        array_1d<double, 3> mScalarVariableGradient;
        double mEffectiveVelocityMagnitude;
        double mScalarVariableValue;
        double mEffectiveKinematicViscosity;
        double mReactionTerm;
        double mAbsoluteReactionTerm;
        double mSourceTerm;
        double mStabilizationTau;
        double mVelocityDotScalarVariableGradient;
        double mResidual;
        double mAbsoluteResidual;
        double mScalarMultiplier;
        double mElementLength;
        double mRelaxedScalarRateVariableValue;
        double mPositivityPreservingMatrixCoefficient;
        double mDeltaTime;
        double mBossakAlpha;
        double mBossakGamma;
        double mDynamicTau;
        double mStabilizationDiscreteDiffusionMatrixUserCoefficient;
        double mStabilizationPositivityPreservingMatrixUserCoefficient;
        double mNumberOfGaussPoints;

        BVector<TNumNodes> mPrimalVariableValues;
        BVector<TNumNodes> mRelaxedPrimalVariableValues;
        BVector<TNumNodes> mVelocityConvectiveTerms;
        BMatrix<TNumNodes, TNumNodes> mPrimalDampingMatrix;
        BMatrix<TNumNodes, TNumNodes> mShapeDerivativesInnerProduct;

        ///@}

    };

    ///@}
};

} // namespace ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointUtilities
} // namespace Kratos

#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_ADJOINT_SENSITIVITY_DERIVATIVES_H_INCLUDED