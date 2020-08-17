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
#include <cmath>

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/geometry_utilities.h"
#include "utilities/time_discretization.h"

// Application includes
#include "convection_diffusion_reaction_stabilization_utilities.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_omega/k_element_data.h"
#include "custom_elements/data_containers/k_omega/omega_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/k_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/omega_element_data.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "rans_application_variables.h"

// Include base h
#include "convection_diffusion_reaction_residual_based_flux_corrected_element.h"

namespace Kratos
{
template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Check sizes and initialize matrix
    if (rLeftHandSideMatrix.size1() != TNumNodes || rLeftHandSideMatrix.size2() != TNumNodes) {
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes, TNumNodes);

    // Calculate RHS
    this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rRightHandSideVector.size() != TNumNodes) {
        rRightHandSideVector.resize(TNumNodes, false);
    }

    noalias(rRightHandSideVector) = ZeroVector(TNumNodes);

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const IndexType num_gauss_points = gauss_weights.size();

    const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
    const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
    const double bossak_gamma =
        TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
    const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
    const double element_length = this->GetGeometry().Length();

    const auto& r_geometry = this->GetGeometry();
    TConvectionDiffusionReactionData r_current_data(r_geometry);

    r_current_data.CalculateConstants(rCurrentProcessInfo);

    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector gauss_shape_functions = row(shape_functions, g);

        r_current_data.CalculateGaussPointData(gauss_shape_functions, r_shape_derivatives);
        const array_1d<double, 3>& velocity = r_current_data.CalculateEffectiveVelocity(
            gauss_shape_functions, r_shape_derivatives);
        const double effective_kinematic_viscosity =
            r_current_data.CalculateEffectiveKinematicViscosity(
                gauss_shape_functions, r_shape_derivatives);

        const double reaction = r_current_data.CalculateReactionTerm(
            gauss_shape_functions, r_shape_derivatives);
        const double source = r_current_data.CalculateSourceTerm(
            gauss_shape_functions, r_shape_derivatives);

        const double tau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
            element_length, norm_2(velocity), reaction, effective_kinematic_viscosity,
            bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

        BoundedVector<double, TNumNodes> velocity_convective_terms;
        this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);

        ConvectionDiffusionReactionStabilizationUtilities::AddSourceTermWithSUPGStabilizationGaussPointContributions(
            rRightHandSideVector, source, std::abs(reaction), tau,
            velocity_convective_terms, gauss_weights[g], gauss_shape_functions);
    }

    KRATOS_CATCH("");
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateLocalVelocityContribution(
    MatrixType& rDampingMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);

    // Now calculate an additional contribution to the residual: r -= rDampingMatrix * (u,p)
    BoundedVector<double, TNumNodes> U;
    this->GetValuesArray(U);
    noalias(rRightHandSideVector) -= prod(rDampingMatrix, U);
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rMassMatrix.size1() != TNumNodes || rMassMatrix.size2() != TNumNodes) {
        rMassMatrix.resize(TNumNodes, TNumNodes, false);
    }

    rMassMatrix.clear();

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const IndexType num_gauss_points = gauss_weights.size();

    const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
    const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
    const double bossak_gamma =
        TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
    const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
    const double element_length = this->GetGeometry().Length();

    const auto& r_geometry = this->GetGeometry();
    TConvectionDiffusionReactionData r_current_data(r_geometry);

    r_current_data.CalculateConstants(rCurrentProcessInfo);

    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector gauss_shape_functions = row(shape_functions, g);

        const double mass = gauss_weights[g] * (1.0 / TNumNodes);
        this->AddLumpedMassMatrix(rMassMatrix, mass);

        r_current_data.CalculateGaussPointData(gauss_shape_functions, r_shape_derivatives);
        const array_1d<double, 3>& velocity = r_current_data.CalculateEffectiveVelocity(
            gauss_shape_functions, r_shape_derivatives);
        BoundedVector<double, TNumNodes> velocity_convective_terms;
        this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);
        const double effective_kinematic_viscosity =
            r_current_data.CalculateEffectiveKinematicViscosity(
                gauss_shape_functions, r_shape_derivatives);

        const double reaction = r_current_data.CalculateReactionTerm(
            gauss_shape_functions, r_shape_derivatives);

        const double tau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
            element_length, norm_2(velocity), reaction, effective_kinematic_viscosity,
            bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

        ConvectionDiffusionReactionStabilizationUtilities::AddMassMatrixSUPGStabilizationGaussPointContributions(
            rMassMatrix, std::abs(reaction), tau, velocity_convective_terms,
            gauss_weights[g], gauss_shape_functions);
    }

    KRATOS_CATCH("");
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    const double scalar_multiplier =
        this->CalculatePrimalDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
    this->SetValue(ERROR_OVERALL, scalar_multiplier);

    double local_matrix_norm = norm_frobenius(rDampingMatrix);
    local_matrix_norm = (local_matrix_norm > 0.0 ? local_matrix_norm : 1.0);

    const double discrete_upwind_operator_coefficient =
        rCurrentProcessInfo[RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT];
    const double diagonal_positivity_preserving_coefficient =
        rCurrentProcessInfo[RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT];

    BoundedMatrix<double, TNumNodes, TNumNodes> discrete_diffusion_matrix;
    double matrix_norm;
    ConvectionDiffusionReactionStabilizationUtilities::CalculateDiscreteUpwindOperator<TNumNodes>(
        matrix_norm, discrete_diffusion_matrix, rDampingMatrix);

    double diagonal_coefficient =
        ConvectionDiffusionReactionStabilizationUtilities::CalculatePositivityPreservingMatrix(
            rDampingMatrix);

    diagonal_coefficient *= diagonal_positivity_preserving_coefficient * scalar_multiplier;

    noalias(rDampingMatrix) += discrete_diffusion_matrix *
                             (discrete_upwind_operator_coefficient * scalar_multiplier);
    noalias(rDampingMatrix) += IdentityMatrix(TNumNodes) * (diagonal_coefficient);

    this->SetValue(RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT,
                   discrete_upwind_operator_coefficient * matrix_norm *
                       scalar_multiplier / local_matrix_norm);
    this->SetValue(RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT,
                   diagonal_coefficient / local_matrix_norm);

}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
double ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculatePrimalDampingMatrix(
    Matrix& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rDampingMatrix.size1() != TNumNodes || rDampingMatrix.size2() != TNumNodes) {
        rDampingMatrix.resize(TNumNodes, TNumNodes, false);
    }

    rDampingMatrix.clear();

    // Get Shape function data
    Vector gauss_weights;
    Matrix shape_functions;
    ShapeFunctionDerivativesArrayType shape_derivatives;
    this->CalculateGeometryData(gauss_weights, shape_functions, shape_derivatives);
    const IndexType num_gauss_points = gauss_weights.size();

    const double delta_time = this->GetDeltaTime(rCurrentProcessInfo);
    const double bossak_alpha = rCurrentProcessInfo[BOSSAK_ALPHA];
    const double bossak_gamma =
        TimeDiscretization::Bossak(bossak_alpha, 0.25, 0.5).GetGamma();
    const double dynamic_tau = rCurrentProcessInfo[DYNAMIC_TAU];
    const double element_length = this->GetGeometry().Length();

    array_1d<double, 3> variable_gradient;
    const Variable<double>& primal_variable =
        TConvectionDiffusionReactionData::GetScalarVariable();
    const Variable<double>& relaxed_primal_rate_variable =
        TConvectionDiffusionReactionData::GetScalarRelaxedRateVariable();

    const auto& r_geometry = this->GetGeometry();
    TConvectionDiffusionReactionData r_current_data(r_geometry);
    double variable_value, relaxed_variable_acceleration;

    r_current_data.CalculateConstants(rCurrentProcessInfo);

    double scalar_multiplier = 0.0;
    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);

        r_current_data.CalculateGaussPointData(gauss_shape_functions, r_shape_derivatives);
        const array_1d<double, 3>& velocity = r_current_data.CalculateEffectiveVelocity(
            gauss_shape_functions, r_shape_derivatives);
        BoundedVector<double, TNumNodes> velocity_convective_terms;
        this->GetConvectionOperator(velocity_convective_terms, velocity, r_shape_derivatives);
        const double velocity_magnitude = norm_2(velocity);

        const double effective_kinematic_viscosity =
            r_current_data.CalculateEffectiveKinematicViscosity(
                gauss_shape_functions, r_shape_derivatives);

        const double reaction = r_current_data.CalculateReactionTerm(
            gauss_shape_functions, r_shape_derivatives);

        const double tau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
            element_length, velocity_magnitude, reaction, effective_kinematic_viscosity,
            bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

        const double source = r_current_data.CalculateSourceTerm(
            gauss_shape_functions, r_shape_derivatives);
        this->CalculateGradient(variable_gradient, primal_variable, r_shape_derivatives);
        const double velocity_dot_variable_gradient =
            inner_prod(velocity, variable_gradient);

        RansCalculationUtilities::EvaluateInPoint(
            r_geometry, gauss_shape_functions,
            std::tie(variable_value, primal_variable),
            std::tie(relaxed_variable_acceleration, relaxed_primal_rate_variable));

        double residual = relaxed_variable_acceleration;
        residual += velocity_dot_variable_gradient;
        residual += reaction * variable_value;
        residual -= source;

        if (variable_value > 0.0) {
            scalar_multiplier += std::abs(residual) * tau / variable_value;
        }

        const Matrix& dNa_dNb = prod(r_shape_derivatives, trans(r_shape_derivatives));

        this->AddDampingMatrixGaussPointContributions(
            rDampingMatrix, reaction, effective_kinematic_viscosity,
            velocity_convective_terms, gauss_weights[g], gauss_shape_functions, dNa_dNb);

        ConvectionDiffusionReactionStabilizationUtilities::AddDampingMatrixSUPGStabilizationGaussPointContributions(
            rDampingMatrix, reaction, tau, velocity_convective_terms,
            gauss_weights[g], gauss_shape_functions);
    }

    r_current_data.UpdateElementDataValueContainer(*this);

    return scalar_multiplier / static_cast<double>(num_gauss_points);

    KRATOS_CATCH("");
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
GeometryData::IntegrationMethod ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetIntegrationMethod() const
{
    return GeometryData::GI_GAUSS_1;
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
double ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetScalarVariableGradientNorm(
    const Matrix& rShapeFunctionDerivatives,
    const int Step) const
{
    KRATOS_TRY;

    array_1d<double, 3> scalar_variable_gradient;
    this->CalculateGradient(scalar_variable_gradient,
                            TConvectionDiffusionReactionData::GetScalarVariable(),
                            rShapeFunctionDerivatives, Step);
    return norm_2(scalar_variable_gradient);

    KRATOS_CATCH("");
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateGradient(
    BoundedMatrix<double, TDim, TDim>& rOutput,
    const Variable<array_1d<double, 3>>& rVariable,
    const Matrix& rShapeDerivatives,
    const int Step) const
{
    const auto& r_geometry = this->GetGeometry();

    RansCalculationUtilities::CalculateGradient<TDim>(
        rOutput, r_geometry, rVariable, rShapeDerivatives, Step);
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
void ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::CalculateGradient(
    array_1d<double, 3>& rOutput,
    const Variable<double>& rVariable,
    const Matrix& rShapeDerivatives,
    const int Step) const
{
    const auto& r_geometry = this->GetGeometry();
    RansCalculationUtilities::CalculateGradient(rOutput, r_geometry, rVariable,
                                                rShapeDerivatives, Step);
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
double ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetDeltaTime(
    const ProcessInfo& rProcessInfo) const
{
    return rProcessInfo[DELTA_TIME];
}

// template instantiations
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KEpsilonElementData::KElementData<2>>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KEpsilonElementData::EpsilonElementData<2>>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KOmegaElementData::KElementData<2>>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KOmegaElementData::OmegaElementData<2>>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KOmegaSSTElementData::KElementData<2>>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, KOmegaSSTElementData::OmegaElementData<2>>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KEpsilonElementData::KElementData<3>>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KEpsilonElementData::EpsilonElementData<3>>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KOmegaElementData::KElementData<3>>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KOmegaElementData::OmegaElementData<3>>;

template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KOmegaSSTElementData::KElementData<3>>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<3, 4, KOmegaSSTElementData::OmegaElementData<3>>;

} // namespace Kratos.
