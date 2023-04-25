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
#include "containers/variable.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/cfd_variables.h"
#include "includes/node.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "utilities/time_discretization.h"

// Application includes
#include "custom_constitutive/fluid_constitutive_law.h"
#include "custom_elements/convection_diffusion_reaction_stabilization_adjoint_utilities.h"
#include "custom_elements/convection_diffusion_reaction_stabilization_utilities.h"
#include "custom_utilities/fluid_calculation_utilities.h"
#include "custom_utilities/fluid_element_utilities.h"
#include "rans_application_variables.h"

// Derivative data type includes
// stabilization validaton
#include "custom_elements/data_containers/stabilization_validation/circular_convection_element_data_derivatives.h"
#include "custom_elements/data_containers/stabilization_validation/diffusion_element_data_derivatives.h"

// k-epsilon
#include "custom_elements/data_containers/k_epsilon/k_element_data_derivatives.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data_derivatives.h"

// k-omega
#include "custom_elements/data_containers/k_omega/k_element_data_derivatives.h"
#include "custom_elements/data_containers/k_omega/omega_element_data_derivatives.h"

// k-omega
#include "custom_elements/data_containers/k_omega_sst/k_element_data_derivatives.h"
#include "custom_elements/data_containers/k_omega_sst/omega_element_data_derivatives.h"

// Include base h
#include "convection_diffusion_reaction_residual_based_flux_corrected_derivatives.h"

namespace Kratos
{

/***************************************************************************************************/
/***************************************** Static Methods ******************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Check(
    const Element& rElement,
    const ProcessInfo& rProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(BOSSAK_ALPHA)) << "BOSSAK_ALPHA is not found in process info. \n";
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(DYNAMIC_TAU)) << "DYNAMIC_TAU is not found in process info. \n";
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT)) << "RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT is not found in process info. \n";
    KRATOS_ERROR_IF_NOT(rProcessInfo.Has(RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT)) << "RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT is not found in process info. \n";

    TElementDataType::Check(rElement, rProcessInfo);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
GeometryData::IntegrationMethod ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

/***************************************************************************************************/
/********************************************** Data ***********************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Data::Data(
    const Element& rElement,
    const ProcessInfo& rProcessInfo)
    : mrElement(rElement),
      mrElementData(rElement.GetGeometry(), rElement.GetProperties(), rProcessInfo)
{
    KRATOS_TRY

    mrElementData.CalculateConstants(rProcessInfo);

    mNumberOfGaussPoints = 0.0;
    mScalarMultiplier = 0.0;
    mPrimalDampingMatrix.clear();

    mDeltaTime = rProcessInfo[DELTA_TIME];
    KRATOS_ERROR_IF(mDeltaTime >= 0.0)
        << "Adjoints are computed in reverse time, therefore DELTA_TIME should "
           "be negative. [ DELTA_TIME = "
        << mDeltaTime << " ].\n";
    mDeltaTime *= -1.0;

    mBossakAlpha = rProcessInfo[BOSSAK_ALPHA];
    mBossakGamma = TimeDiscretization::Bossak(mBossakAlpha, 0.25, 0.5).GetGamma();
    mDynamicTau = rProcessInfo[DYNAMIC_TAU];
    mDiscreteUpwindOperatorCoefficient = rProcessInfo[RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT];
    mDiagonalPositivityPreservingCoefficient = rProcessInfo[RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT];

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Data::CalculateGaussPointData(
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const int Step)
{
    KRATOS_TRY

    mrElementData.CalculateGaussPointData(rN, rdNdX, Step);

    mVelocityMagnitude = norm_2(mrElementData.GetEffectiveVelocity());
    mAbsoluteReactionTerm = std::abs(mrElementData.GetReactionTerm());
    mElementLength = mrElementData.GetGeometry().Length();

    const auto& primal_scalar_variable = TElementDataType::GetScalarVariable();
    const auto& primal_relaxed_scalar_variable = primal_scalar_variable.GetTimeDerivative().GetTimeDerivative();

    FluidCalculationUtilities::EvaluateInPoint(
        mrElementData.GetGeometry(), rN, Step,
        std::tie(mPrimalVariableValue, primal_scalar_variable),
        std::tie(mPrimalRelaxedVariableRateValue, primal_relaxed_scalar_variable));

    FluidCalculationUtilities::EvaluateGradientInPoint(
        mrElementData.GetGeometry(), rdNdX, Step,
        std::tie(mPrimalVariableGradient, TElementDataType::GetScalarVariable()));

    mStabilizationTau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
        mElementLength, mVelocityMagnitude, mrElementData.GetReactionTerm(),
        mrElementData.GetEffectiveKinematicViscosity(), mBossakAlpha,
        mBossakGamma, mDeltaTime, mDynamicTau);

    mPrimalVariableGradientDotVelocity =
        inner_prod(mPrimalVariableGradient, mrElementData.GetEffectiveVelocity());

    noalias(mdNadNb) = prod(rdNdX, trans(rdNdX));
    noalias(mVelocityConvectiveTerms) = prod(rdNdX, mrElementData.GetEffectiveVelocity());
    noalias(mPrimalVariableGradientDotDnDx) = prod(rdNdX, mPrimalVariableGradient);

    mResidual = mPrimalRelaxedVariableRateValue;
    mResidual += mPrimalVariableGradientDotVelocity;
    mResidual += mrElementData.GetReactionTerm() * mPrimalVariableValue;
    mResidual -= mrElementData.GetSourceTerm();
    mAbsoluteResidual = std::abs(mResidual);

    if (mPrimalVariableValue > 0.0) {
        mScalarMultiplier += mAbsoluteResidual * mStabilizationTau / mPrimalVariableValue;
    }

    for (IndexType i = 0; i < TNumNodes; ++i) {
        const auto& r_node = mrElementData.GetGeometry()[i];
        mNodalPrimalVariableValues[i] = r_node.FastGetSolutionStepValue(primal_scalar_variable);
        mNodalPrimalRelaxedVariableRateValues[i] = r_node.FastGetSolutionStepValue(primal_relaxed_scalar_variable);
    }

    AddDampingMatrixContributions(mPrimalDampingMatrix, W, rN, rdNdX);

    mNumberOfGaussPoints += 1.0;

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Data::AddDampingMatrixContributions(
    MatrixNN& rDampingMatrix,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX) const
{
    KRATOS_TRY

    const auto& element_data = mrElementData;

    // compute primal equation coefficients
    const auto& viscosity = element_data.GetEffectiveKinematicViscosity();
    const auto& reaction_term = element_data.GetReactionTerm();

    for (IndexType a = 0; a < TNumNodes; ++a) {
        const double tau_operator = mVelocityConvectiveTerms[a] + mAbsoluteReactionTerm * rN[a];

        for (IndexType b = 0; b < TNumNodes; ++b) {

            double value = 0.0;

            value += rN[a] * mVelocityConvectiveTerms[b];
            value += rN[a] * reaction_term * rN[b];
            value += viscosity * mdNadNb(a, b);

            // adding LHS SUPG stabilization terms
            value += tau_operator * mStabilizationTau * mVelocityConvectiveTerms[b];
            value += tau_operator * mStabilizationTau * reaction_term * rN[b];

            rDampingMatrix(a, b) += value * W;
        }
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Data::CalculateDataAfterGaussPointPointLoop()
{
    using namespace ConvectionDiffusionReactionStabilizationUtilities;
    CalculateDiscreteUpwindOperator<TNumNodes>(mDiscreteDiffusionMatrix, mPrimalDampingMatrix);
    mDiagonalCoefficient = CalculatePositivityPreservingMatrix(mPrimalDampingMatrix);
    mScalarMultiplier /= mNumberOfGaussPoints;
}

/***************************************************************************************************/
/************************************* Residual Contributions **************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualsContributions::ResidualsContributions(
    Data& rData)
    : mrData(rData)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualsContributions::AddGaussPointResidualsContributions(
    VectorN& rResidual,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
    KRATOS_TRY

    const auto& element_data = mrData.mrElementData;

    // compute primal equation coefficients
    const auto& velocity = element_data.GetEffectiveVelocity();
    const auto& viscosity = element_data.GetEffectiveKinematicViscosity();
    const auto& reaction_term = element_data.GetReactionTerm();
    const auto& source_term = element_data.GetSourceTerm();

    for (IndexType a = 0; a < TNumNodes; ++a) {
        const double tau_operator = mrData.mVelocityConvectiveTerms[a] + mrData.mAbsoluteReactionTerm * rN[a];

        // compute the residual derivative without discrete upwind operator and positivity preserving
        // coefficient contributions to derivatives.
        double value = 0.0;

        // adding RHS contributions
        value += rN[a] * source_term;

        // adding RHS SUPG stabilization terms
        value += tau_operator * mrData.mStabilizationTau * source_term;

        // adding LHS terms
        value -= rN[a] * mrData.mPrimalVariableGradientDotVelocity;
        value -= rN[a] * reaction_term * mrData.mPrimalVariableValue;
        value -= viscosity * mrData.mPrimalVariableGradientDotDnDx[a];

        // adding LHS SUPG stabilization terms
        value -= tau_operator * mrData.mStabilizationTau * mrData.mPrimalVariableGradientDotVelocity;
        value -= tau_operator * mrData.mStabilizationTau * reaction_term * mrData.mPrimalVariableValue;

        // adding Mass terms
        value -= (1.0/TNumNodes) * mrData.mNodalPrimalRelaxedVariableRateValues[a];
        value -= tau_operator * mrData.mStabilizationTau * mrData.mPrimalRelaxedVariableRateValue;

        rResidual[a] += value * W;
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualsContributions::AddResidualsContributionsAfterGaussPointLoop(
    VectorN& rResidual)
{
    noalias(rResidual) -= prod(mrData.mDiscreteDiffusionMatrix, mrData.mNodalPrimalVariableValues) * (mrData.mDiscreteUpwindOperatorCoefficient * mrData.mScalarMultiplier);
    noalias(rResidual) -= mrData.mNodalPrimalVariableValues * (mrData.mDiagonalPositivityPreservingCoefficient * mrData.mDiagonalCoefficient * mrData.mScalarMultiplier);
}

/***************************************************************************************************/
/************************************** Variable Derivatives ***************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
template <class TDerivativesType>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::VariableDerivatives<
    TDerivativesType>::VariableDerivatives(Data& rData)
    : mrData(rData),
      mResidualWeightDerivativeContributions(rData)
{
    // clear the matrices holding derivatives
    for (IndexType i = 0; i < TDerivativesSize; ++i) {
        mPrimalDampingMatrixDerivatives[i].clear();
        mScalarMultiplierDerivatives[i] = 0.0;
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
template <class TDerivativesType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::VariableDerivatives<TDerivativesType>::CalculateGaussPointResidualsDerivativeContributions(
    VectorN& rResidualDerivative,
    const int NodeIndex,
    const int DirectionIndex,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative,
    const double MassTermsDerivativesWeight)
{
    KRATOS_TRY

    using namespace ConvectionDiffusionReactionStabilizationUtilities;

    using adjoint_utilities = AdjointUtilities<TDim, TNumNodes>;

    const auto& element_data = mrData.mrElementData;

    // create derivatives object
    const TDerivativesType derivative(element_data, NodeIndex, DirectionIndex, W, rN,
                                      rdNdX, WDerivative, DetJDerivative, rdNdXDerivative);

    // compute primal equation coefficients
    const auto& velocity = element_data.GetEffectiveVelocity();
    const auto& viscosity = element_data.GetEffectiveKinematicViscosity();
    const auto& reaction_term = element_data.GetReactionTerm();
    const auto& source_term = element_data.GetSourceTerm();

    // compute primal equation coefficient derivatives
    const auto& velocity_derivative = derivative.CalculateEffectiveVelocityDerivative();
    const auto& viscosity_derivative = derivative.CalculateEffectiveKinematicViscosityDerivative();
    const auto& reaction_term_derivative = derivative.CalculateReactionTermDerivative();
    const auto& source_term_derivative = derivative.CalculateSourceTermDerivative();

    // utility derivatives
    const double element_length_derivative = TDerivativesType::UtilitiesDerivative::CalculateElementLengthDerivative(DetJDerivative, mrData.mElementLength);

    // compute auxiliary values
    const VectorN velocity_derivative_dot_dn_dx = prod(rdNdX, velocity_derivative);
    const VectorN velocity_dot_dn_dx_derivative = prod(rdNdXDerivative, velocity);

    const double absolute_reaction_term_derivative = adjoint_utilities::CalculateAbsoluteValueDerivative(reaction_term, reaction_term_derivative);
    const double velocity_magnitude_derivative = adjoint_utilities::CalculateVectorNormDerivative(mrData.mVelocityMagnitude, velocity, velocity_derivative);
    const double stabilization_tau_derivative = adjoint_utilities::CalculateStabilizationTauDerivative(
        mrData.mStabilizationTau, mrData.mVelocityMagnitude, viscosity,
        reaction_term, mrData.mElementLength, velocity_magnitude_derivative,
        viscosity_derivative, reaction_term_derivative, element_length_derivative);

    ArrayD primal_variable_gradient_derivative;
    FluidCalculationUtilities::EvaluateGradientInPoint(
        element_data.GetGeometry(), rdNdXDerivative,
        std::tie(primal_variable_gradient_derivative, TDerivativesType::DataType::GetScalarVariable()));

    const VectorN primal_variable_gradient_derivative_dot_dn_dx = prod(rdNdX, primal_variable_gradient_derivative);
    const VectorN primal_variable_gradient_dot_dn_dx_derivative = prod(rdNdXDerivative, mrData.mPrimalVariableGradient);

    const double velocity_derivative_dot_primal_variable_gradient = inner_prod(velocity_derivative, mrData.mPrimalVariableGradient);
    const double velocity_dot_primal_variable_gradient_derivative = inner_prod(velocity, primal_variable_gradient_derivative);

    const MatrixNN dn_dx_dot_dn_dx_derivative = prod(rdNdX, trans(rdNdXDerivative));

    MatrixNN& r_primal_damping_matrix_derivative = mPrimalDampingMatrixDerivatives[TDerivativeDimension * NodeIndex + DirectionIndex];

    // compute scalar multiplier derivative
    double& r_scalar_multiplier_derivative = mScalarMultiplierDerivatives[TDerivativeDimension * NodeIndex + DirectionIndex];

    // compute residual derivative
    double residual_derivative = 0.0;

    residual_derivative += velocity_derivative_dot_primal_variable_gradient;
    residual_derivative += velocity_dot_primal_variable_gradient_derivative;
    residual_derivative += mrData.mVelocityConvectiveTerms[NodeIndex] * TSelfWeight;

    residual_derivative += reaction_term_derivative * mrData.mPrimalVariableValue;
    residual_derivative += reaction_term * rN[NodeIndex] * TSelfWeight;

    residual_derivative -= source_term_derivative;

    const double absolute_residual_derivative = adjoint_utilities::CalculateAbsoluteValueDerivative(mrData.mResidual, residual_derivative);

    if (mrData.mPrimalVariableValue > 0.0) {
        r_scalar_multiplier_derivative += absolute_residual_derivative * mrData.mStabilizationTau / mrData.mPrimalVariableValue;
        r_scalar_multiplier_derivative += mrData.mAbsoluteResidual * stabilization_tau_derivative / mrData.mPrimalVariableValue;
        r_scalar_multiplier_derivative -= mrData.mAbsoluteResidual * mrData.mStabilizationTau  * rN[NodeIndex] * TSelfWeight / std::pow(mrData.mPrimalVariableValue, 2);
    }

    for (IndexType a = 0; a < TNumNodes; ++a) {
        const double tau_operator = mrData.mVelocityConvectiveTerms[a] + mrData.mAbsoluteReactionTerm * rN[a];
        const double tau_operator_derivative = velocity_derivative_dot_dn_dx[a] + velocity_dot_dn_dx_derivative[a] + absolute_reaction_term_derivative * rN[a];

        // compute LHS matrix derivative for residual based flux corrected stabilization
        for (IndexType b = 0; b < TNumNodes; ++b) {
            double value = 0;

            // adding LHS derivative contributions
            value += rN[a] * velocity_derivative_dot_dn_dx[b];
            value += rN[a] * velocity_dot_dn_dx_derivative[b];

            value += rN[a] * reaction_term_derivative * rN[b];

            value += viscosity_derivative * mrData.mdNadNb(a, b);
            value += viscosity * dn_dx_dot_dn_dx_derivative(b, a);
            value += viscosity * dn_dx_dot_dn_dx_derivative(a, b);

            // adding LHS SUPG stabilization terms derivatives
            value += tau_operator_derivative * mrData.mStabilizationTau * mrData.mVelocityConvectiveTerms[b];
            value += tau_operator * stabilization_tau_derivative * mrData.mVelocityConvectiveTerms[b];
            value += tau_operator * mrData.mStabilizationTau * velocity_derivative_dot_dn_dx[b];
            value += tau_operator * mrData.mStabilizationTau * velocity_dot_dn_dx_derivative[b];

            value += tau_operator_derivative * mrData.mStabilizationTau * reaction_term * rN[b];
            value += tau_operator * stabilization_tau_derivative * reaction_term * rN[b];
            value += tau_operator * mrData.mStabilizationTau * reaction_term_derivative * rN[b];

            r_primal_damping_matrix_derivative(a, b) += value * W;
        }

        // compute the residual derivative without discrete upwind operator and positivity preserving
        // coefficient contributions to derivatives.
        double value = 0.0;

        // adding RHS derivative contributions
        value += rN[a] * source_term_derivative;

        // adding RHS SUPG stabilization terms derivatives
        value += tau_operator_derivative * mrData.mStabilizationTau * source_term;
        value += tau_operator * stabilization_tau_derivative * source_term;
        value += tau_operator * mrData.mStabilizationTau * source_term_derivative;

        // adding LHS derivative contributions
        value -= rN[a] * velocity_derivative_dot_primal_variable_gradient;
        value -= rN[a] * velocity_dot_primal_variable_gradient_derivative;

        value -= rN[a] * reaction_term_derivative * mrData.mPrimalVariableValue;

        value -= viscosity_derivative * mrData.mPrimalVariableGradientDotDnDx[a];
        value -= viscosity * primal_variable_gradient_derivative_dot_dn_dx[a];
        value -= viscosity * primal_variable_gradient_dot_dn_dx_derivative[a];

        // adding LHS SUPG stabilization term derivative contributions
        value -= tau_operator_derivative * mrData.mStabilizationTau * mrData.mPrimalVariableGradientDotVelocity;
        value -= tau_operator * stabilization_tau_derivative * mrData.mPrimalVariableGradientDotVelocity;
        value -= tau_operator * mrData.mStabilizationTau * velocity_derivative_dot_primal_variable_gradient;
        value -= tau_operator * mrData.mStabilizationTau * velocity_dot_primal_variable_gradient_derivative;

        value -= tau_operator_derivative * mrData.mStabilizationTau * reaction_term * mrData.mPrimalVariableValue;
        value -= tau_operator * stabilization_tau_derivative * reaction_term * mrData.mPrimalVariableValue;
        value -= tau_operator * mrData.mStabilizationTau * reaction_term_derivative * mrData.mPrimalVariableValue;

        // add mass term derivatives
        value -= tau_operator_derivative * mrData.mStabilizationTau * mrData.mPrimalRelaxedVariableRateValue * MassTermsDerivativesWeight;
        value -= tau_operator * stabilization_tau_derivative * mrData.mPrimalRelaxedVariableRateValue * MassTermsDerivativesWeight;

        rResidualDerivative[a] = value * W;
    }

    mrData.AddDampingMatrixContributions(r_primal_damping_matrix_derivative, WDerivative, rN, rdNdX);
    mResidualWeightDerivativeContributions.AddGaussPointResidualsContributions(rResidualDerivative, WDerivative, rN, rdNdX);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
template <class TDerivativesType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::VariableDerivatives<TDerivativesType>::CalculateResidualsDerivativeContributionsAfterGaussPointPointLoop(
    VectorN& rResidualDerivative,
    const int NodeIndex,
    const int DirectionIndex)
{
    using namespace ConvectionDiffusionReactionStabilizationUtilities;

    using adjoint_utilities = AdjointUtilities<TDim, TNumNodes>;

    rResidualDerivative.clear();

    // for some reason linking fails if multiplication by 1.0 is removed. Thus, it is kept.
    noalias(rResidualDerivative) -= column(mrData.mPrimalDampingMatrix, NodeIndex) * (TSelfWeight * 1.0);

    const auto& primal_damping_matrix_derivative = mPrimalDampingMatrixDerivatives[NodeIndex * TDerivativeDimension + DirectionIndex];
    const auto& scalar_multiplier_derivative = mScalarMultiplierDerivatives[NodeIndex * TDerivativeDimension + DirectionIndex] / mrData.mNumberOfGaussPoints;

    // add derivative contributions from discrete diffusion matrix
    VectorN discrete_diffusion_matrix_derivative;
    adjoint_utilities::CalculateDiscreteUpwindOperatorDerivative(
        discrete_diffusion_matrix_derivative, mrData.mNodalPrimalVariableValues,
        mrData.mPrimalDampingMatrix, primal_damping_matrix_derivative);

    noalias(rResidualDerivative) -= discrete_diffusion_matrix_derivative * (mrData.mDiscreteUpwindOperatorCoefficient * mrData.mScalarMultiplier);
    noalias(rResidualDerivative) -= prod(mrData.mDiscreteDiffusionMatrix, mrData.mNodalPrimalVariableValues) * (mrData.mDiscreteUpwindOperatorCoefficient * scalar_multiplier_derivative);
    noalias(rResidualDerivative) -= column(mrData.mDiscreteDiffusionMatrix, NodeIndex) * (mrData.mDiscreteUpwindOperatorCoefficient * mrData.mScalarMultiplier * TSelfWeight);

    // add derivative contributions from positivity preserving matrix
    const double diagonal_coefficient_derivative =
        adjoint_utilities::CalculatePositivityPreservingCoefficientDerivative(
            mrData.mDiagonalCoefficient, mrData.mPrimalDampingMatrix,
            primal_damping_matrix_derivative);

    noalias(rResidualDerivative) -= mrData.mNodalPrimalVariableValues * (mrData.mDiagonalPositivityPreservingCoefficient * mrData.mDiagonalCoefficient * scalar_multiplier_derivative);
    noalias(rResidualDerivative) -= mrData.mNodalPrimalVariableValues * (mrData.mDiagonalPositivityPreservingCoefficient * diagonal_coefficient_derivative * mrData.mScalarMultiplier);
    rResidualDerivative[NodeIndex] -= (mrData.mDiagonalPositivityPreservingCoefficient * mrData.mDiagonalCoefficient * mrData.mScalarMultiplier * TSelfWeight);
}

/***************************************************************************************************/
/*************************************** Second Derivatives ****************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivatives::SecondDerivatives(
    Data& rData)
    : mrData(rData)
{
    mScalarMultiplierDerivatives.clear();
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivatives::CalculateGaussPointResidualsDerivativeContributions(
    VectorN& rResidualDerivative,
    const int NodeIndex,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
    using namespace ConvectionDiffusionReactionStabilizationUtilities;

    using adjoint_utilities = AdjointUtilities<TDim, TNumNodes>;

    rResidualDerivative.clear();

    // adding mass term derivatives
    rResidualDerivative[NodeIndex] -= (W / TNumNodes);

    for (IndexType a = 0; a < TNumNodes; ++a) {
        const double tau_operator = mrData.mVelocityConvectiveTerms[a] + mrData.mAbsoluteReactionTerm * rN[a];
        rResidualDerivative[a] -= tau_operator * mrData.mStabilizationTau * rN[NodeIndex] * W;
    }

    const double absolute_residual_derivative = adjoint_utilities::CalculateAbsoluteValueDerivative(mrData.mResidual, rN[NodeIndex]);

    if (mrData.mPrimalVariableValue > 0.0) {
        mScalarMultiplierDerivatives[NodeIndex] += absolute_residual_derivative * mrData.mStabilizationTau / mrData.mPrimalVariableValue;
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivatives::CalculateResidualsDerivativeContributionsAfterGaussPointPointLoop(
    VectorN& rResidualDerivative,
    const int NodeIndex)
{
    rResidualDerivative.clear();

    const auto& scalar_multiplier_derivative = mScalarMultiplierDerivatives[NodeIndex] / mrData.mNumberOfGaussPoints;

    // add derivative contributions from discrete diffusion matrix
    noalias(rResidualDerivative) -= prod(mrData.mDiscreteDiffusionMatrix, mrData.mNodalPrimalVariableValues) * (mrData.mDiscreteUpwindOperatorCoefficient * scalar_multiplier_derivative);

    // add derivative contributions from positivity preserving matrix
    noalias(rResidualDerivative) -= mrData.mNodalPrimalVariableValues * (mrData.mDiagonalPositivityPreservingCoefficient * mrData.mDiagonalCoefficient * scalar_multiplier_derivative);
}

// template instantiations

// stabilization validation element data derivatives
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, StabilizationValidationElementData::CircularConvectionElementDataDerivatives::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, StabilizationValidationElementData::CircularConvectionElementDataDerivatives::Data>::VariableDerivatives<StabilizationValidationElementData::CircularConvectionElementDataDerivatives::PhiDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, StabilizationValidationElementData::CircularConvectionElementDataDerivatives::Data>::VariableDerivatives<StabilizationValidationElementData::CircularConvectionElementDataDerivatives::ShapeDerivative>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, StabilizationValidationElementData::DiffusionElementDataDerivatives::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, StabilizationValidationElementData::DiffusionElementDataDerivatives::Data>::VariableDerivatives<StabilizationValidationElementData::DiffusionElementDataDerivatives::PhiDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, StabilizationValidationElementData::DiffusionElementDataDerivatives::Data>::VariableDerivatives<StabilizationValidationElementData::DiffusionElementDataDerivatives::ShapeDerivative>;

// k-epsilon k element derivatives
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::KElementDataDerivatives<2, 3>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KEpsilonElementData::KElementDataDerivatives<2, 3>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KEpsilonElementData::KElementDataDerivatives<2, 3>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KEpsilonElementData::KElementDataDerivatives<2, 3>::EpsilonDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KEpsilonElementData::KElementDataDerivatives<2, 3>::ShapeDerivative>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::KElementDataDerivatives<3, 4>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KEpsilonElementData::KElementDataDerivatives<3, 4>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KEpsilonElementData::KElementDataDerivatives<3, 4>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KEpsilonElementData::KElementDataDerivatives<3, 4>::EpsilonDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KEpsilonElementData::KElementDataDerivatives<3, 4>::ShapeDerivative>;

// k-epsilon epsilon element derivatives
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::EpsilonDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KEpsilonElementData::EpsilonElementDataDerivatives<2, 3>::ShapeDerivative>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::EpsilonDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KEpsilonElementData::EpsilonElementDataDerivatives<3, 4>::ShapeDerivative>;

// k-omega k element derivatives
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaElementData::KElementDataDerivatives<2, 3>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaElementData::KElementDataDerivatives<2, 3>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaElementData::KElementDataDerivatives<2, 3>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaElementData::KElementDataDerivatives<2, 3>::OmegaDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaElementData::KElementDataDerivatives<2, 3>::ShapeDerivative>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaElementData::KElementDataDerivatives<3, 4>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaElementData::KElementDataDerivatives<3, 4>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaElementData::KElementDataDerivatives<3, 4>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaElementData::KElementDataDerivatives<3, 4>::OmegaDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaElementData::KElementDataDerivatives<3, 4>::ShapeDerivative>;

// k-omega omega element derivatives
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaElementData::OmegaElementDataDerivatives<2, 3>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaElementData::OmegaElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaElementData::OmegaElementDataDerivatives<2, 3>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaElementData::OmegaElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaElementData::OmegaElementDataDerivatives<2, 3>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaElementData::OmegaElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaElementData::OmegaElementDataDerivatives<2, 3>::OmegaDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaElementData::OmegaElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaElementData::OmegaElementDataDerivatives<2, 3>::ShapeDerivative>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaElementData::OmegaElementDataDerivatives<3, 4>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaElementData::OmegaElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaElementData::OmegaElementDataDerivatives<3, 4>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaElementData::OmegaElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaElementData::OmegaElementDataDerivatives<3, 4>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaElementData::OmegaElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaElementData::OmegaElementDataDerivatives<3, 4>::OmegaDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaElementData::OmegaElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaElementData::OmegaElementDataDerivatives<3, 4>::ShapeDerivative>;

// k-omega-sst k element derivatives
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaSSTElementData::KElementDataDerivatives<2, 3>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaSSTElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaSSTElementData::KElementDataDerivatives<2, 3>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaSSTElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaSSTElementData::KElementDataDerivatives<2, 3>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaSSTElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaSSTElementData::KElementDataDerivatives<2, 3>::OmegaDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaSSTElementData::KElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaSSTElementData::KElementDataDerivatives<2, 3>::ShapeDerivative>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaSSTElementData::KElementDataDerivatives<3, 4>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaSSTElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaSSTElementData::KElementDataDerivatives<3, 4>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaSSTElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaSSTElementData::KElementDataDerivatives<3, 4>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaSSTElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaSSTElementData::KElementDataDerivatives<3, 4>::OmegaDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaSSTElementData::KElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaSSTElementData::KElementDataDerivatives<3, 4>::ShapeDerivative>;

// k-omega-sst omega element derivatives
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::OmegaDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<2, 3, KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::Data>::VariableDerivatives<KOmegaSSTElementData::OmegaElementDataDerivatives<2, 3>::ShapeDerivative>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaSSTElementData::OmegaElementDataDerivatives<3, 4>::Data>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaSSTElementData::OmegaElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaSSTElementData::OmegaElementDataDerivatives<3, 4>::UDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaSSTElementData::OmegaElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaSSTElementData::OmegaElementDataDerivatives<3, 4>::KDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaSSTElementData::OmegaElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaSSTElementData::OmegaElementDataDerivatives<3, 4>::OmegaDerivative>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<3, 4, KOmegaSSTElementData::OmegaElementDataDerivatives<3, 4>::Data>::VariableDerivatives<KOmegaSSTElementData::OmegaElementDataDerivatives<3, 4>::ShapeDerivative>;

} // namespace Kratos
