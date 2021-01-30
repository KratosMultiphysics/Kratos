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

// Derivative data type includes
#include "custom_elements/data_containers/k_epsilon/k_element_data_derivatives.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data_derivatives.h"

// Include base h
#include "convection_diffusion_reaction_residual_based_flux_corrected_derivatives.h"

namespace Kratos
{

/***************************************************************************************************/
/***************************************** Static Methods ******************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
int ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Check(
    const GeometryType& rGeometry,
    const ProcessInfo& rProcessInfo)
{
    return 0;
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
GeometryData::IntegrationMethod ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::GetIntegrationMethod()
{
    return GeometryData::IntegrationMethod::GI_GAUSS_1;
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

    mDeltaTime = rProcessInfo[DELTA_TIME];
    KRATOS_ERROR_IF(mDeltaTime > 0.0)
        << "Adjoints are computed in reverse time, therefore DELTA_TIME should "
           "be negative. [ DELTA_TIME = "
        << mDeltaTime << " ].\n";
    mDeltaTime *= -1.0;

    mBossakAlpha = rProcessInfo[BOSSAK_ALPHA];
    mBossakGamma = TimeDiscretization::Bossak(mBossakAlpha, 0.25, 0.5).GetGamma();
    mDynamicTau = rProcessInfo[DYNAMIC_TAU];

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::Data::CalculateGaussPointData(
    const Vector& rN,
    const Matrix& rdNdX,
    const int Step)
{
    KRATOS_TRY

    mrElementData.CalculateGaussPointData(rN, rdNdX, Step);

    mVelocityMagnitude = norm_2(mrElementData.GetEffectiveVelocity());
    mAbsoluteReactionTerm = std::abs(mrElementData.GetReactionTerm());
    mElementLength = mrElementData.GetGeometry().Length();

    FluidCalculationUtilities::EvaluateInPoint(
        mrElementData.GetGeometry(), rN, Step,
        std::tie(mPrimalVariableValue, TElementDataType::GetScalarVariable()));

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

    KRATOS_CATCH("");
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
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualsContributions::Initialize(
    const Vector& rResidual,
    const ProcessInfo& rProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualsContributions::AddResidualsContributions(
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

        rResidual[a] += value * W;
    }

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualsContributions::AddDampingMatrixContributions(
    MatrixNN& rDampingMatrix,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
    // KRATOS_TRY

    // const auto& element_data = mrData.mrElementData;

    // // compute primal equation coefficients
    // const auto& velocity = element_data.GetEffectiveVelocity();
    // const auto& viscosity = element_data.GetEffectiveKinematicViscosity();
    // const auto& reaction_term = element_data.GetReactionTerm();
    // const auto& source_term = element_data.GetSourceTerm();

    // for (IndexType a = 0; a < TNumNodes; ++a) {
    //     const double tau_operator = mrData.mVelocityConvectiveTerms[a] + mrData.mAbsoluteReactionTerm * rN[a];

    //     for (IndexType b = 0; b < TNumNodes; ++b) {

    //         double value = 0.0;

    //         value += rN[a] * mrData.mPrimalVariableGradientDotVelocity;
    //         value -= rN[a] * reaction_term * mrData.mPrimalVariableValue;

    //         value -= viscosity * mrData.mPrimalVariableGradientDotDnDx[a];

    //         rResidual[a] += value * W;
    //     }
    // }

    // KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::ResidualsContributions::Finalize(
    VectorN& rResidual,
    const ProcessInfo& rProcessInfo)
{
    rResidual.clear();
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
    // clear the matrices holding primal damping matrix derivatives
    for (IndexType i = 0; i < TDerivativesSize; ++i) {
        mPrimalMatrixDerivatives[i].clear();
    }
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
template <class TDerivativesType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::VariableDerivatives<TDerivativesType>::Initialize(
    const Vector& rResidualDerivative,
    const ProcessInfo& rProcessInfo)
{
    mResidualWeightDerivativeContributions.Initialize(rResidualDerivative, rProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
template <class TDerivativesType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::VariableDerivatives<TDerivativesType>::CalculateResidualsDerivativeContributions(
    VectorN& rResidualDerivative,
    const int NodeIndex,
    const int DirectionIndex,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX,
    const double WDerivative,
    const double DetJDerivative,
    const Matrix& rdNdXDerivative)
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

    const double velocity_derivative_dot_primal_variable_gradient = inner_prod(velocity_derivative, mrData.mPrimalVariableGradient);
    const double velocity_dot_primal_variable_gradient_derivative = inner_prod(velocity, primal_variable_gradient_derivative);

    const VectorN primal_variable_gradient_derivative_dot_dn_dx = prod(rdNdX, primal_variable_gradient_derivative);
    const VectorN primal_variable_gradient_dot_dn_dx_derivative = prod(rdNdXDerivative, mrData.mPrimalVariableGradient);

    const MatrixNN dn_dx_dot_dn_dx_derivative = prod(rdNdX, trans(rdNdXDerivative));

    MatrixNN& r_primal_matrix_derivative = mPrimalMatrixDerivatives[TDerivativeDimension * NodeIndex + DirectionIndex];

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

            r_primal_matrix_derivative(a, b) += value * W;
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
        value -= rN[a] * mrData.mVelocityConvectiveTerms[NodeIndex] * TSelfWeight;

        value -= rN[a] * reaction_term_derivative * mrData.mPrimalVariableValue;
        value -= rN[a] * reaction_term * rN[NodeIndex] * TSelfWeight;

        value -= viscosity_derivative * mrData.mPrimalVariableGradientDotDnDx[a];
        value -= viscosity * primal_variable_gradient_derivative_dot_dn_dx[a];
        value -= viscosity * primal_variable_gradient_dot_dn_dx_derivative[a];
        value -= viscosity * mrData.mdNadNb(a, NodeIndex) * TSelfWeight;

        // adding LHS SUPG stabilization terms derivatives
        value -= tau_operator_derivative * mrData.mStabilizationTau * mrData.mPrimalVariableGradientDotVelocity;
        value -= tau_operator * stabilization_tau_derivative * mrData.mPrimalVariableGradientDotVelocity;
        value -= tau_operator * mrData.mStabilizationTau * velocity_derivative_dot_primal_variable_gradient;
        value -= tau_operator * mrData.mStabilizationTau * velocity_dot_primal_variable_gradient_derivative;
        value -= tau_operator * mrData.mStabilizationTau * mrData.mVelocityConvectiveTerms[NodeIndex] * TSelfWeight;

        value -= tau_operator_derivative * mrData.mStabilizationTau * reaction_term * mrData.mPrimalVariableValue;
        value -= tau_operator * stabilization_tau_derivative * reaction_term * mrData.mPrimalVariableValue;
        value -= tau_operator * mrData.mStabilizationTau * reaction_term_derivative * mrData.mPrimalVariableValue;
        value -= tau_operator * mrData.mStabilizationTau * reaction_term * rN[NodeIndex] * TSelfWeight;

        rResidualDerivative[a] = value * W;
    }

    mResidualWeightDerivativeContributions.AddDampingMatrixContributions(r_primal_matrix_derivative, WDerivative, rN, rdNdX);
    mResidualWeightDerivativeContributions.AddResidualsContributions(rResidualDerivative, WDerivative, rN, rdNdX);

    KRATOS_CATCH("");
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
template <class TDerivativesType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::VariableDerivatives<TDerivativesType>::Finalize(
    VectorN& rResidualDerivative,
    const ProcessInfo& rProcessInfo)
{
    rResidualDerivative.clear();
}

/***************************************************************************************************/
/*************************************** Second Derivatives ****************************************/
/***************************************************************************************************/

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivatives::SecondDerivatives(
    Data& rData)
    : mrData(rData)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivatives::Initialize(
    const Vector& rResidualDerivative,
    const ProcessInfo& rProcessInfo)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivatives::CalculateResidualsDerivativeContributions(
    VectorN& rResidualDerivative,
    const double W,
    const Vector& rN,
    const Matrix& rdNdX)
{
}

template <unsigned int TDim, unsigned int TNumNodes, class TElementDataType>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedDerivatives<TDim, TNumNodes, TElementDataType>::SecondDerivatives::Finalize(
    VectorN& rResidualDerivative,
    const ProcessInfo& rProcessInfo)
{
}

// template instantiations

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

} // namespace Kratos
