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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_ADJOINT_UTILITIES_H_INCLUDED)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_ADJOINT_UTILITIES_H_INCLUDED

// System includes
#include <vector>

// Project includes
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/time_discretization.h"

// Application includes
#include "custom_utilities/rans_calculation_utilities.h"
#include "custom_elements/convection_diffusion_reaction_stabilization_adjoint_utilities.h"
#include "custom_elements/convection_diffusion_reaction_stabilization_utilities.h"

namespace Kratos
{
namespace ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointUtilities
{
template <unsigned int TDim, unsigned int TNumNodes, class TPrimalEquationStateDerivativesType>
class StabilizationStateDerivatives
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
    ///@name Public forward delcarations
    ///@{

    class FirstDerivativesData;

    ///@}
    ///@name Public static operations
    ///@{

    static void Check(
        const GeometryType& rGeometry,
        const ProcessInfo& rProcessInfo)
    {
        KRATOS_TRY

        TPrimalEquationStateDerivativesType::Data::Check(rGeometry, rProcessInfo);

        KRATOS_CATCH("");
    }

    static const std::vector<const Variable<double>*> GetDofVariablesList()
    {
        KRATOS_TRY

        return {&TPrimalEquationStateDerivativesType::Data::GetAdjointScalarVariable()};

        KRATOS_CATCH("");
    }

    static GeometryData::IntegrationMethod GetIntegrationMethod()
    {
        return GeometryData::GI_GAUSS_2;
    }

    ///@}
    ///@name Public static operations
    ///@{

    ///@}
    ///@name Public classes
    ///@{

    template <class TDerivativesType, unsigned int TEquationOffset, unsigned int TDerivativeOffset>
    class FirstDerivatives
    {
    public:
        ///@name Public type definitions
        ///@{

        static constexpr unsigned int TDerivativesSize = TDerivativesType::TDerivativesSize;

        static constexpr unsigned int TDerivativesDim = TDerivativesSize / TNumNodes;

        static constexpr unsigned int TSelfWeight = (TEquationOffset == TDerivativeOffset);

        using VectorDerivativesType = BVector<TDerivativesSize>;

        template<unsigned int TSize>
        using MatrixDerivativesType = BMatrix<TDerivativesSize, TSize>;

        using AdjointUtilities = ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<TDim, TNumNodes>;

        ///@}
        ///@name Life cycle
        ///@{

        FirstDerivatives(
            const FirstDerivativesData& rData)
            : mrData(rData)
        {
        }

        ///@}
        ///@name Public operations
        ///@{

        void Initialize(
            Matrix& rOutput,
            const ProcessInfo& rProcessInfo)
        {
            for (IndexType c = 0; c < TDerivativesSize; ++c) {
                mPrimalDampingMatrixDerivatives[c].clear();
            }

            mScalarMultiplierDerivatives.clear();
            mDerivativeBlockSize = rOutput.size1() / TNumNodes;
        }

        void CalculateResidualDerivatives(
            Matrix& rOutput,
            const double GPWeight,
            const Vector& rGPShapeFunctions,
            const Matrix& rGPShapeFunctionDerivatives)
        {
            TDerivativesType derivatives(mrData.mEquationData);

            // calculate residual equation coefficient derivatives
            MatrixDerivativesType<3> velocity_derivatives;
            derivatives.CalculateEffectiveVelocityDerivatives(velocity_derivatives, rGPShapeFunctions, rGPShapeFunctionDerivatives);

            VectorDerivativesType viscosity_derivatives;
            derivatives.CalculateEffectiveKinematicViscosityDerivatives(viscosity_derivatives, rGPShapeFunctions, rGPShapeFunctionDerivatives);

            VectorDerivativesType reaction_term_derivatives;
            derivatives.CalculateReactionTermDerivatives(reaction_term_derivatives, rGPShapeFunctions, rGPShapeFunctionDerivatives);

            VectorDerivativesType source_term_derivatives;
            derivatives.CalculateSourceTermDerivatives(source_term_derivatives, rGPShapeFunctions, rGPShapeFunctionDerivatives);

            // calculate auxiliary derivatives
            VectorDerivativesType velocity_magnitude_derivatives;
            AdjointUtilities::CalculateEffectiveVelocityMagnitudeDerivative(
                velocity_magnitude_derivatives, mrData.mEffectiveVelocityMagnitude,
                mrData.mEffectiveVelocity, velocity_derivatives);

            VectorDerivativesType absolute_reaction_term_derivatives;
            AdjointUtilities::CalculateAbsoluteScalarDerivatives(
                absolute_reaction_term_derivatives, mrData.mReactionTerm,
                reaction_term_derivatives);

            VectorDerivativesType stabilization_tau_derivatives;
            AdjointUtilities::CalculateStabilizationTauDerivatives(
                stabilization_tau_derivatives, mrData.mEffectiveVelocityMagnitude,
                mrData.mEffectiveKinematicViscosity, mrData.mReactionTerm,
                mrData.mElementLength, mrData.mStabilizationTau, viscosity_derivatives,
                reaction_term_derivatives, velocity_magnitude_derivatives);

            VectorDerivativesType absolute_residual_derivatives;
            AdjointUtilities::CalculateResidualDerivatives(
                absolute_residual_derivatives, TSelfWeight, mrData.mScalarVariableValue,
                mrData.mReactionTerm, mrData.mEffectiveVelocity,
                mrData.mScalarVariableGradient, rGPShapeFunctions,
                rGPShapeFunctionDerivatives, reaction_term_derivatives,
                source_term_derivatives, velocity_derivatives);
            AdjointUtilities::CalculateAbsoluteScalarDerivatives(
                absolute_residual_derivatives, mrData.mResidual, absolute_residual_derivatives);

            AddScalarMultiplierFirstDerivatives(
                mScalarMultiplierDerivatives, TSelfWeight, mrData.mScalarVariableValue,
                mrData.mStabilizationTau, mrData.mAbsoluteResidual, absolute_residual_derivatives,
                stabilization_tau_derivatives, rGPShapeFunctions);

            MatrixDerivativesType<TNumNodes> velocity_convective_term_derivatives;
            AdjointUtilities::MatrixMatrixProduct(velocity_convective_term_derivatives, velocity_derivatives, trans(rGPShapeFunctionDerivatives));

            const Matrix& dNa_dNb = prod(rGPShapeFunctionDerivatives, trans(rGPShapeFunctionDerivatives));
            BVector<TNumNodes> dNa_i_dPhi_i;
            AdjointUtilities::MatrixVectorProduct(dNa_i_dPhi_i, rGPShapeFunctionDerivatives, mrData.mScalarVariableGradient);

            // calculating primal damping matrix derivatives
            for (IndexType c = 0; c < TNumNodes; ++c) {
                for (IndexType k = 0; k < TDerivativesDim; ++k) {

                    const auto deriv_index = c * TDerivativesDim + k;

                    BoundedMatrix<double, TNumNodes, TNumNodes>& r_derivatives_matrix =  mPrimalDampingMatrixDerivatives[deriv_index];
                    const double dU_ic_dPhi_i = inner_prod(mrData.mScalarVariableGradient, row(velocity_derivatives, deriv_index));

                    for (IndexType a = 0; a < TNumNodes; ++a) {
                        const double tau_operator = mrData.mVelocityConvectiveTerms[a] + mrData.mAbsoluteReactionTerm * rGPShapeFunctions[a];
                        const double tau_operator_derivative = velocity_convective_term_derivatives(deriv_index, a) + absolute_reaction_term_derivatives[deriv_index] * rGPShapeFunctions[a];

                        for (IndexType b = 0; b < TNumNodes; ++b) {
                            double value = 0.0;

                            value += rGPShapeFunctions[a] * velocity_convective_term_derivatives(deriv_index, b);
                            value += rGPShapeFunctions[a] * reaction_term_derivatives[deriv_index] * rGPShapeFunctions[b];
                            value += dNa_dNb(a, b) * viscosity_derivatives[deriv_index];

                            // adding derivatives of SUPG
                            value += tau_operator * stabilization_tau_derivatives[deriv_index] * mrData.mVelocityConvectiveTerms[b];
                            value += tau_operator * mrData.mStabilizationTau * velocity_convective_term_derivatives(deriv_index, b);
                            value += tau_operator_derivative * mrData.mStabilizationTau * mrData.mVelocityConvectiveTerms[b];

                            value += tau_operator * stabilization_tau_derivatives[deriv_index]  * mrData.mReactionTerm * rGPShapeFunctions[b];
                            value += tau_operator * mrData.mStabilizationTau  * reaction_term_derivatives[deriv_index] * rGPShapeFunctions[b];
                            value += tau_operator_derivative * mrData.mStabilizationTau  * mrData.mReactionTerm * rGPShapeFunctions[b];

                            r_derivatives_matrix(a, b) += GPWeight * value;
                        }

                        double value = 0.0;

                        // add RHS derivative contributions
                        value += rGPShapeFunctions[a] * source_term_derivatives[deriv_index];

                        // add RHS SUPG stabilization derivative terms
                        value += stabilization_tau_derivatives[deriv_index] * tau_operator * mrData.mSourceTerm;
                        value += mrData.mStabilizationTau * tau_operator_derivative * mrData.mSourceTerm;
                        value += mrData.mStabilizationTau * tau_operator * source_term_derivatives[deriv_index];

                        // add LHS derivative contributions
                        value -= rGPShapeFunctions[a] * dU_ic_dPhi_i;
                        value -= rGPShapeFunctions[a] * mrData.mVelocityConvectiveTerms[c] * TSelfWeight;
                        value -= rGPShapeFunctions[a] * reaction_term_derivatives[deriv_index] * mrData.mScalarVariableValue;
                        value -= rGPShapeFunctions[a] * mrData.mReactionTerm * rGPShapeFunctions[c] * TSelfWeight;
                        value -= viscosity_derivatives[deriv_index] * dNa_i_dPhi_i[a];
                        value -= dNa_dNb(c, a) * mrData.mEffectiveKinematicViscosity * TSelfWeight;

                        // add LHS SUPG stabilization derivative terms
                        value -= stabilization_tau_derivatives[deriv_index] * tau_operator * mrData.mVelocityDotScalarVariableGradient;
                        value -= mrData.mStabilizationTau * tau_operator_derivative * mrData.mVelocityDotScalarVariableGradient;
                        value -= mrData.mStabilizationTau * tau_operator * dU_ic_dPhi_i;
                        value -= mrData.mStabilizationTau * tau_operator * mrData.mVelocityConvectiveTerms[c] * TSelfWeight;

                        value -= stabilization_tau_derivatives[deriv_index] * tau_operator * mrData.mReactionTerm * mrData.mScalarVariableValue;
                        value -= mrData.mStabilizationTau * tau_operator_derivative * mrData.mReactionTerm * mrData.mScalarVariableValue;
                        value -= mrData.mStabilizationTau * tau_operator * reaction_term_derivatives[deriv_index] * mrData.mScalarVariableValue;
                        value -= mrData.mStabilizationTau * tau_operator * mrData.mReactionTerm * rGPShapeFunctions[c] * TSelfWeight;

                        // add mass matrix contributions
                        value -= stabilization_tau_derivatives[deriv_index] * tau_operator * mrData.mRelaxedScalarRateVariableValue;
                        value -= mrData.mStabilizationTau * tau_operator_derivative * mrData.mRelaxedScalarRateVariableValue;

                        rOutput(c * mDerivativeBlockSize + TDerivativeOffset + k, a * mDerivativeBlockSize + TEquationOffset) += value * GPWeight;
                    }
                }
            }
        }

        void Finalize(Matrix& rOutput, const ProcessInfo& rProcessInfo)
        {
            noalias(mScalarMultiplierDerivatives) = mScalarMultiplierDerivatives / mrData.mNumberOfGaussPoints;

            MatrixDerivativesType<TNumNodes> discrete_diffusion_matrix_residual_derivatives;
            CalculateStabilizationDiscreteUpwindMatrixResidualFristDerivatives(
                discrete_diffusion_matrix_residual_derivatives, TSelfWeight,
                mrData.mScalarMultiplier,
                mrData.mStabilizationDiscreteDiffusionMatrixUserCoefficient,
                mScalarMultiplierDerivatives, mrData.mPrimalVariableValues,
                mrData.mDiscreteDiffusionResidualContribution, mrData.mPrimalDampingMatrix,
                mrData.mDiscreteDiffusionMatrix, mPrimalDampingMatrixDerivatives);

            MatrixDerivativesType<TNumNodes> positivity_preserving_matrix_residual_derivatives;
            CalculateStabilizationPositivityPreservingMatrixResidualFristDerivatives(
                positivity_preserving_matrix_residual_derivatives, TSelfWeight,
                mrData.mScalarMultiplier, mrData.mPositivityPreservingMatrixCoefficient,
                mrData.mStabilizationPositivityPreservingMatrixUserCoefficient,
                mScalarMultiplierDerivatives, mrData.mPrimalVariableValues,
                mrData.mPrimalDampingMatrix, mPrimalDampingMatrixDerivatives);

            for (IndexType c = 0; c < TNumNodes; ++c) {
                for (IndexType k = 0; k < TDerivativesDim; ++k) {
                    for (IndexType a = 0; a < TNumNodes; ++a) {
                        double& value =
                            rOutput(c * mDerivativeBlockSize + TDerivativeOffset + k,
                                    a * mDerivativeBlockSize + TEquationOffset);
                        value -= discrete_diffusion_matrix_residual_derivatives(
                            c * TDerivativesDim + k, a);
                        value -= positivity_preserving_matrix_residual_derivatives(
                            c * TDerivativesDim + k, a);
                    }
                }
            }
        }

        ///@}
    private:
        ///@name Private members
        ///@{

        const FirstDerivativesData& mrData;
        IndexType mDerivativeBlockSize;
        BoundedVector<BMatrix<TNumNodes, TNumNodes>, TDerivativesSize> mPrimalDampingMatrixDerivatives;
        VectorDerivativesType mScalarMultiplierDerivatives;

        ///@}
        ///@name Private operations
        ///@{

        ///@}
    };

    template <unsigned int TEquationOffset>
    class SecondDerivatives
    {
    public:
        ///@name Life cycle
        ///@{

        SecondDerivatives(
            const GeometryType& rGeometry)
            : mEquationData(rGeometry)
        {
        }

        ///@}
        ///@name Public operations
        ///@{

        void Initialize(Matrix& rOutput, const ProcessInfo& rProcessInfo)
        {
            mEquationData.CalculateConstants(rProcessInfo);
            mPrimalDampingMatrix.clear();
            mScalarMultiplierDerivatives.clear();

            mDeltaTime = rProcessInfo[DELTA_TIME] * -1.0;
            mBossakAlpha = rProcessInfo[BOSSAK_ALPHA];
            mBossakGamma = TimeDiscretization::Bossak(mBossakAlpha, 0.25, 0.5).GetGamma();
            mDynamicTau = rProcessInfo[DYNAMIC_TAU];
            mElementLength = mEquationData.GetGeometry().Length();
            mStabilizationDiscreteDiffusionMatrixUserCoefficient = rProcessInfo[RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT];
            mStabilizationPositivityPreservingMatrixUserCoefficient = rProcessInfo[RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT];
            mDerivativeBlockSize = rOutput.size1() / TNumNodes;
            mNumberOfGaussPoints = 0.0;
        }

        void CalculateResidualDerivatives(
            Matrix& rMatrix,
            const double GPWeight,
            const Vector& rGPShapeFunctions,
            const Matrix& rGPShapeFunctionDerivatives)
        {
            using namespace RansCalculationUtilities;

            mEquationData.CalculateGaussPointData(rGPShapeFunctions, rGPShapeFunctionDerivatives);

            const auto& r_scalar_variable = TPrimalEquationStateDerivativesType::Data::GetScalarVariable();

            const array_1d<double, 3>& velocity = mEquationData.CalculateEffectiveVelocity(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            const double viscosity = mEquationData.CalculateEffectiveKinematicViscosity(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            const double reaction_term = mEquationData.CalculateReactionTerm(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            const double source_term = mEquationData.CalculateSourceTerm(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            const double absolute_reaction_term = std::abs(reaction_term);

            const double tau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
                mElementLength, norm_2(velocity), reaction_term,
                viscosity, mBossakAlpha, mBossakGamma, mDeltaTime, mDynamicTau);

            BVector<TNumNodes> velocity_convective_terms;
            GetConvectionOperator(velocity_convective_terms, velocity, rGPShapeFunctionDerivatives);

            AddPrimalDampingMatrixGaussPointContributions(
                mPrimalDampingMatrix, viscosity, reaction_term,
                absolute_reaction_term, tau, GPWeight, velocity_convective_terms,
                rGPShapeFunctions, rGPShapeFunctionDerivatives);

            double scalar_variable_value, relaxed_scalar_rate_variable_value;
            EvaluateInPoint(
                mEquationData.GetGeometry(), rGPShapeFunctions,
                std::tie(scalar_variable_value, r_scalar_variable),
                std::tie(relaxed_scalar_rate_variable_value, r_scalar_variable.GetTimeDerivative().GetTimeDerivative()));

            array_1d<double, 3> scalar_variable_gradient;
            CalculateGradient(scalar_variable_gradient, mEquationData.GetGeometry(),
                              r_scalar_variable, rGPShapeFunctionDerivatives);

            double residual = relaxed_scalar_rate_variable_value;
            residual += inner_prod(scalar_variable_gradient, velocity);
            residual += reaction_term * scalar_variable_value;
            residual -= source_term;

            if (scalar_variable_value > 0.0) {
                noalias(mScalarMultiplierDerivatives) +=
                    rGPShapeFunctions *
                    (tau * (residual >= 0.0 ? 1.0 : -1.0) / scalar_variable_value);
            }

            const double mass = GPWeight / TNumNodes;
            for (IndexType c = 0; c < TNumNodes; ++c) {
                rMatrix(c * mDerivativeBlockSize + TEquationOffset,
                        c * mDerivativeBlockSize + TEquationOffset) -= mass;
                for (IndexType a = 0; a < TNumNodes; ++a) {
                    rMatrix(c * mDerivativeBlockSize + TEquationOffset,
                            a * mDerivativeBlockSize + TEquationOffset) -=
                        GPWeight * tau *
                        (velocity_convective_terms[a] +
                         absolute_reaction_term * rGPShapeFunctions[a]) *
                        rGPShapeFunctions[c];
                }
            }

            mNumberOfGaussPoints += 1.0;
        }

        void Finalize(Matrix& rOutput, const ProcessInfo& rProcessInfo)
        {
            noalias(mScalarMultiplierDerivatives) = mScalarMultiplierDerivatives / mNumberOfGaussPoints;

            BVector<TNumNodes> primal_values;
            for (IndexType i_node = 0; i_node < mEquationData.GetGeometry().size(); ++i_node) {
                primal_values[i_node] = mEquationData.GetGeometry()[i_node].FastGetSolutionStepValue(
                    TPrimalEquationStateDerivativesType::Data::GetScalarVariable());
            }

            // calculating residual derivatives of discrete upwind operator
            BoundedMatrix<double, TNumNodes, TNumNodes> discrete_diffusion_matrix;
            double discrete_diffusion_matrix_coefficient;
            ConvectionDiffusionReactionStabilizationUtilities::CalculateDiscreteUpwindOperator<TNumNodes>(
                discrete_diffusion_matrix_coefficient,
                discrete_diffusion_matrix, mPrimalDampingMatrix);

            BVector<TNumNodes> residual;
            noalias(residual) =
                prod(discrete_diffusion_matrix,
                     primal_values * mStabilizationDiscreteDiffusionMatrixUserCoefficient);

            // calculating residual derivatives of positivity preserving matrix
            const double positivity_preserving_coefficient =
                ConvectionDiffusionReactionStabilizationUtilities::CalculatePositivityPreservingMatrix(
                    mPrimalDampingMatrix);
            noalias(residual) += primal_values * positivity_preserving_coefficient *
                                 mStabilizationPositivityPreservingMatrixUserCoefficient;

            for (IndexType c = 0; c < TNumNodes; ++c) {
                for (IndexType a = 0; a < TNumNodes; ++a) {
                    rOutput(c * mDerivativeBlockSize + TEquationOffset,
                            a * mDerivativeBlockSize + TEquationOffset) -=
                        mScalarMultiplierDerivatives[c] * residual[a];
                }
            }
        }

        ///@}
    private:
        ///@name Private members
        ///@{

        typename TPrimalEquationStateDerivativesType::Data mEquationData;

        BMatrix<TNumNodes, TNumNodes> mPrimalDampingMatrix;
        BVector<TNumNodes> mScalarMultiplierDerivatives;
        IndexType mDerivativeBlockSize;
        double mNumberOfGaussPoints;
        double mElementLength;
        double mDeltaTime;
        double mBossakAlpha;
        double mBossakGamma;
        double mDynamicTau;
        double mStabilizationDiscreteDiffusionMatrixUserCoefficient;
        double mStabilizationPositivityPreservingMatrixUserCoefficient;

        ///@}
    };

    class FirstDerivativesData
    {
    public:
        ///@name Life cycle
        ///@{

        FirstDerivativesData(
            const GeometryType& rGeometry)
            : mEquationData(rGeometry)
        {
        }

        ///@}
        ///@name Public operations
        ///@{

        void Initialize(
            const ProcessInfo& rProcessInfo)
        {
            mEquationData.CalculateConstants(rProcessInfo);
            mPrimalDampingMatrix.clear();

            for (IndexType i_node = 0; i_node < mEquationData.GetGeometry().size(); ++i_node) {
                mPrimalVariableValues[i_node] =
                    mEquationData.GetGeometry()[i_node].FastGetSolutionStepValue(
                        TPrimalEquationStateDerivativesType::Data::GetScalarVariable());
            }

            mDeltaTime = rProcessInfo[DELTA_TIME] * -1.0;
            mBossakAlpha = rProcessInfo[BOSSAK_ALPHA];
            mBossakGamma = TimeDiscretization::Bossak(mBossakAlpha, 0.25, 0.5).GetGamma();
            mDynamicTau = rProcessInfo[DYNAMIC_TAU];
            mElementLength = mEquationData.GetGeometry().Length();
            mScalarMultiplier = 0.0;
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
            mEquationData.CalculateGaussPointData(rGPShapeFunctions, rGPShapeFunctionDerivatives);

            mEffectiveVelocity = mEquationData.CalculateEffectiveVelocity(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            mEffectiveKinematicViscosity = mEquationData.CalculateEffectiveKinematicViscosity(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            mReactionTerm = mEquationData.CalculateReactionTerm(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            mSourceTerm = mEquationData.CalculateSourceTerm(rGPShapeFunctions, rGPShapeFunctionDerivatives);
            mAbsoluteReactionTerm = std::abs(mReactionTerm);

            const auto& r_scalar_variable = TPrimalEquationStateDerivativesType::Data::GetScalarVariable();

            EvaluateInPoint(
                mEquationData.GetGeometry(), rGPShapeFunctions,
                std::tie(mScalarVariableValue, r_scalar_variable),
                std::tie(mRelaxedScalarRateVariableValue, r_scalar_variable.GetTimeDerivative().GetTimeDerivative()));

            CalculateGradient(mScalarVariableGradient, mEquationData.GetGeometry(),
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

            AddPrimalDampingMatrixGaussPointContributions(
                mPrimalDampingMatrix, mEffectiveKinematicViscosity, mReactionTerm,
                mAbsoluteReactionTerm, mStabilizationTau, GPWeight, mVelocityConvectiveTerms,
                rGPShapeFunctions, rGPShapeFunctionDerivatives);

            mNumberOfGaussPoints += 1.0;
        }

        void Finalize(
            const ProcessInfo& rProcessInfo)
        {
            mScalarMultiplier /= mNumberOfGaussPoints;

            ConvectionDiffusionReactionStabilizationUtilities::CalculateDiscreteUpwindOperator<TNumNodes>(
                mDiscreteDiffusionMatrixCoefficient, mDiscreteDiffusionMatrix,
                mPrimalDampingMatrix);

            noalias(mDiscreteDiffusionResidualContribution) = prod(mDiscreteDiffusionMatrix, mPrimalVariableValues);

            mPositivityPreservingMatrixCoefficient =
                ConvectionDiffusionReactionStabilizationUtilities::CalculatePositivityPreservingMatrix(
                    mPrimalDampingMatrix);
        }

        ///@}
    private:
        ///@name Private members
        ///@{

        typename TPrimalEquationStateDerivativesType::Data mEquationData;

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
        double mDiscreteDiffusionMatrixCoefficient;
        double mPositivityPreservingMatrixCoefficient;
        double mDeltaTime;
        double mBossakAlpha;
        double mBossakGamma;
        double mDynamicTau;
        double mStabilizationDiscreteDiffusionMatrixUserCoefficient;
        double mStabilizationPositivityPreservingMatrixUserCoefficient;
        double mNumberOfGaussPoints;

        BVector<TNumNodes> mPrimalVariableValues;
        BVector<TNumNodes> mVelocityConvectiveTerms;
        BVector<TNumNodes> mDiscreteDiffusionResidualContribution;
        BMatrix<TNumNodes, TNumNodes> mPrimalDampingMatrix;
        BMatrix<TNumNodes, TNumNodes> mDiscreteDiffusionMatrix;

        ///@}
        ///@name Friend class declarations for private member data access
        ///@{

        template <class TDerivativesType, unsigned int TEquationOffset, unsigned int TDerivativeOffset>
        friend class FirstDerivatives;

        ///@}
    };

    ///@}
    ///@name Private static operations
    ///@{

    static void AddPrimalDampingMatrixGaussPointContributions(
        BMatrix<TNumNodes, TNumNodes>& rOutput,
        const double EffectiveKinematicViscosity,
        const double ReactionTerm,
        const double AbsoluteReactionTerm,
        const double StabilizationTau,
        const double GPWeight,
        const BVector<TNumNodes>& rVelocityConvectiveTerms,
        const Vector& rGPShapeFunctions,
        const Matrix& rGPShapeFunctionDerivatives)
    {
        for (IndexType a = 0; a < TNumNodes; ++a) {
            for (IndexType b = 0; b < TNumNodes; ++b) {
                const double dNa_dNb =
                    inner_prod(row(rGPShapeFunctionDerivatives, a),
                               row(rGPShapeFunctionDerivatives, b));
                double value = 0.0;

                value += rGPShapeFunctions[a] * rVelocityConvectiveTerms[b];
                value += rGPShapeFunctions[a] * ReactionTerm * rGPShapeFunctions[b];
                value += EffectiveKinematicViscosity * dNa_dNb;

                // Adding SUPG stabilization terms
                value += StabilizationTau *
                         (rVelocityConvectiveTerms[a] +
                          AbsoluteReactionTerm * rGPShapeFunctions[a]) *
                         rVelocityConvectiveTerms[b];
                value += StabilizationTau *
                         (rVelocityConvectiveTerms[a] +
                          AbsoluteReactionTerm * rGPShapeFunctions[a]) *
                         ReactionTerm * rGPShapeFunctions[b];

                rOutput(a, b) += GPWeight * value;
            }
        }
    }

    template<std::size_t TDerivativesSize>
    static typename std::enable_if<(TNumNodes != TDerivativesSize), void>::type
    AddScalarMultiplierFirstDerivatives(
        BVector<TDerivativesSize>& rOutput,
        const double SelfWeight,
        const double ScalarVariableValue,
        const double StabilizationTau,
        const double AbsoluteResidual,
        const BVector<TDerivativesSize>& rAbsoluteResidualDerivatives,
        const BVector<TDerivativesSize>& rStabilizationTauDerivatives,
        const Vector& rGPShapeFunctions)
    {
        if (ScalarVariableValue > 0.0) {
            BVector<TDerivativesSize> output;
            noalias(output) = rAbsoluteResidualDerivatives * StabilizationTau;
            noalias(output) += AbsoluteResidual * rStabilizationTauDerivatives;
            noalias(rOutput) += output * (1.0 / ScalarVariableValue);
        }
    }

    template<std::size_t TDerivativesSize>
    static typename std::enable_if<(TNumNodes == TDerivativesSize), void>::type
    AddScalarMultiplierFirstDerivatives(
        BVector<TDerivativesSize>& rOutput,
        const double SelfWeight,
        const double ScalarVariableValue,
        const double StabilizationTau,
        const double AbsoluteResidual,
        const BVector<TDerivativesSize>& rAbsoluteResidualDerivatives,
        const BVector<TDerivativesSize>& rStabilizationTauDerivatives,
        const Vector& rGPShapeFunctions)
    {
        if (ScalarVariableValue > 0.0) {
            BVector<TDerivativesSize> output;
            noalias(output) = rAbsoluteResidualDerivatives * StabilizationTau;
            noalias(output) += AbsoluteResidual * rStabilizationTauDerivatives;
            noalias(output) -= rGPShapeFunctions * (AbsoluteResidual * StabilizationTau * SelfWeight / ScalarVariableValue);
            noalias(rOutput) += output * (1.0 / ScalarVariableValue);
        }
    }

    template<std::size_t TDerivativesSize>
    static typename std::enable_if<(TNumNodes == TDerivativesSize), void>::type
    CalculateStabilizationDiscreteUpwindMatrixResidualFristDerivatives(
        BMatrix<TDerivativesSize, TNumNodes>& rOutput,
        const double SelfWeight,
        const double ScalarMultiplier,
        const double StabilizationDiscreteDiffusionUserCoefficient,
        const BVector<TDerivativesSize>& rScalarMultiplierDerivatives,
        const BVector<TNumNodes>& rNodalValues,
        const BVector<TNumNodes>& rDiscreteDiffusionValues,
        const BMatrix<TNumNodes, TNumNodes>& rInputMatrix,
        const BMatrix<TNumNodes, TNumNodes>& rDiscreteDiffusionMatrix,
        const BoundedVector<BMatrix<TNumNodes, TNumNodes>, TDerivativesSize>& rInputMatrixDerivatives)
    {
        using AdjointUtilities = ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<TDim, TNumNodes>;

        BMatrix<TDerivativesSize, TNumNodes> discrete_upwind_operator_residual_derivatives;
        AdjointUtilities::CalculateDiscreteUpwindOperatorResidualContributionDerivatives(
            discrete_upwind_operator_residual_derivatives, rNodalValues,
            rInputMatrix, rInputMatrixDerivatives);

        AdjointUtilities::DidacticProduct(rOutput, rScalarMultiplierDerivatives, rDiscreteDiffusionValues);
        noalias(rOutput) += discrete_upwind_operator_residual_derivatives * ScalarMultiplier;
        noalias(rOutput) += rDiscreteDiffusionMatrix * (ScalarMultiplier * SelfWeight);
        noalias(rOutput) = rOutput * StabilizationDiscreteDiffusionUserCoefficient;
    }

    template<std::size_t TDerivativesSize>
    static typename std::enable_if<(TNumNodes != TDerivativesSize), void>::type
    CalculateStabilizationDiscreteUpwindMatrixResidualFristDerivatives(
        BMatrix<TDerivativesSize, TNumNodes>& rOutput,
        const double SelfWeight,
        const double ScalarMultiplier,
        const double StabilizationDiscreteDiffusionUserCoefficient,
        const BVector<TDerivativesSize>& rScalarMultiplierDerivatives,
        const BVector<TNumNodes>& rNodalValues,
        const BVector<TNumNodes>& rDiscreteDiffusionValues,
        const BMatrix<TNumNodes, TNumNodes>& rInputMatrix,
        const BMatrix<TNumNodes, TNumNodes>& rDiscreteDiffusionMatrix,
        const BoundedVector<BMatrix<TNumNodes, TNumNodes>, TDerivativesSize>& rInputMatrixDerivatives)
    {
        using AdjointUtilities = ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<TDim, TNumNodes>;

        BMatrix<TDerivativesSize, TNumNodes> discrete_upwind_operator_residual_derivatives;
        AdjointUtilities::CalculateDiscreteUpwindOperatorResidualContributionDerivatives(
            discrete_upwind_operator_residual_derivatives, rNodalValues,
            rInputMatrix, rInputMatrixDerivatives);

        AdjointUtilities::DidacticProduct(rOutput, rScalarMultiplierDerivatives, rDiscreteDiffusionValues);
        noalias(rOutput) += discrete_upwind_operator_residual_derivatives * ScalarMultiplier;
        noalias(rOutput) = rOutput * StabilizationDiscreteDiffusionUserCoefficient;
    }

    template<std::size_t TDerivativesSize>
    static typename std::enable_if<(TNumNodes != TDerivativesSize), void>::type
    CalculateStabilizationPositivityPreservingMatrixResidualFristDerivatives(
        BMatrix<TDerivativesSize, TNumNodes>& rOutput,
        const double SelfWeight,
        const double ScalarMultiplier,
        const double PositivityPreservingMatrixCoefficient,
        const double StabilizationPositivityPreservingUserCoefficient,
        const BVector<TDerivativesSize>& rScalarMultiplierDerivatives,
        const BVector<TNumNodes>& rNodalScalarValues,
        const BMatrix<TNumNodes, TNumNodes>& rInputMatrix,
        const BoundedVector<BMatrix<TNumNodes, TNumNodes>, TDerivativesSize>& rInputMatrixDerivatives)
    {
        using AdjointUtilities = ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<TDim, TNumNodes>;
        BVector<TDerivativesSize> positivity_preserving_coefficient_derivatives;
        AdjointUtilities::CalculatePositivityPreservingCoefficientDerivatives(
            positivity_preserving_coefficient_derivatives, PositivityPreservingMatrixCoefficient, rInputMatrix,
            rInputMatrixDerivatives);

        AdjointUtilities::DidacticProduct(
            rOutput,
            positivity_preserving_coefficient_derivatives * ScalarMultiplier +
                rScalarMultiplierDerivatives * PositivityPreservingMatrixCoefficient,
            rNodalScalarValues);
        noalias(rOutput) = rOutput * StabilizationPositivityPreservingUserCoefficient;
    }

    template<std::size_t TDerivativesSize>
    static typename std::enable_if<(TNumNodes == TDerivativesSize), void>::type
    CalculateStabilizationPositivityPreservingMatrixResidualFristDerivatives(
        BMatrix<TDerivativesSize, TNumNodes>& rOutput,
        const double SelfWeight,
        const double ScalarMultiplier,
        const double PositivityPreservingMatrixCoefficient,
        const double StabilizationPositivityPreservingUserCoefficient,
        const BVector<TDerivativesSize>& rScalarMultiplierDerivatives,
        const BVector<TNumNodes>& rNodalScalarValues,
        const BMatrix<TNumNodes, TNumNodes>& rInputMatrix,
        const BoundedVector<BMatrix<TNumNodes, TNumNodes>, TDerivativesSize>& rInputMatrixDerivatives)
    {
        using AdjointUtilities =
            ConvectionDiffusionReactionStabilizationUtilities::AdjointUtilities<TDim, TNumNodes>;
        BVector<TDerivativesSize> positivity_preserving_coefficient_derivatives;
        AdjointUtilities::CalculatePositivityPreservingCoefficientDerivatives(
            positivity_preserving_coefficient_derivatives,
            PositivityPreservingMatrixCoefficient, rInputMatrix, rInputMatrixDerivatives);

        AdjointUtilities::DidacticProduct(
            rOutput,
            positivity_preserving_coefficient_derivatives * ScalarMultiplier +
                rScalarMultiplierDerivatives * PositivityPreservingMatrixCoefficient,
            rNodalScalarValues);
        noalias(rOutput) += IdentityMatrix(TNumNodes) * ScalarMultiplier *
                            PositivityPreservingMatrixCoefficient * SelfWeight;
        noalias(rOutput) = rOutput * StabilizationPositivityPreservingUserCoefficient;
    }

    ///@}
};

} // namespace ConvectionDiffusionReactionResidualBasedFluxCorrectedAdjointUtilities
} // namespace Kratos

#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_RESIDUAL_BASED_FLUX_CORRECTED_ADJOINT_UTILITIES_H_INCLUDED