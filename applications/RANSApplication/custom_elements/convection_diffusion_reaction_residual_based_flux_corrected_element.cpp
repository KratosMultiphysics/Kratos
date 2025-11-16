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
#include "custom_elements/data_containers/stabilization_validation/circular_convection_element_data.h"
#include "custom_elements/data_containers/stabilization_validation/body_force_governed_cdr_element_data.h"
#include "custom_elements/data_containers/stabilization_validation/diffusion_element_data.h"
#include "custom_elements/data_containers/k_epsilon/epsilon_element_data.h"
#include "custom_elements/data_containers/k_epsilon/k_element_data.h"
#include "custom_elements/data_containers/k_omega/k_element_data.h"
#include "custom_elements/data_containers/k_omega/omega_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/k_element_data.h"
#include "custom_elements/data_containers/k_omega_sst/omega_element_data.h"
#include "custom_utilities/fluid_calculation_utilities.h"
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
    TConvectionDiffusionReactionData r_current_data(r_geometry, this->GetProperties(), rCurrentProcessInfo);

    r_current_data.CalculateConstants(rCurrentProcessInfo);

    BoundedVector<double, TNumNodes> velocity_convective_terms;
    ArrayD  mesh_velocity; // added for Chimera_RANS simulation

    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector gauss_shape_functions = row(shape_functions, g);

        // Get MESH_VELOCITY for Chimera_RANS simulation
        FluidCalculationUtilities::EvaluateInPoint(
            r_geometry, gauss_shape_functions,
            std::tie(mesh_velocity, MESH_VELOCITY));

        r_current_data.CalculateGaussPointData(gauss_shape_functions, r_shape_derivatives);
        const auto& fluid_velocity = r_current_data.GetEffectiveVelocity();
        const ArrayD& velocity = fluid_velocity - mesh_velocity;  // convective or effective velocity
        const double effective_kinematic_viscosity = r_current_data.GetEffectiveKinematicViscosity();
        const double reaction = r_current_data.GetReactionTerm();
        const double source = r_current_data.GetSourceTerm();

        const double tau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
            element_length, norm_2(velocity), reaction, effective_kinematic_viscosity,
            bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

        noalias(velocity_convective_terms) = prod(r_shape_derivatives, velocity);

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
    TConvectionDiffusionReactionData r_current_data(r_geometry, this->GetProperties(), rCurrentProcessInfo);

    r_current_data.CalculateConstants(rCurrentProcessInfo);

    BoundedVector<double, TNumNodes> velocity_convective_terms;
    ArrayD  mesh_velocity; // added for Chimera_RANS simulation

    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector gauss_shape_functions = row(shape_functions, g);

        const double mass = gauss_weights[g] * (1.0 / TNumNodes);
        this->AddLumpedMassMatrix(rMassMatrix, mass);

        // Get MESH_VELOCITY for Chimera_RANS simulation
        FluidCalculationUtilities::EvaluateInPoint(
            r_geometry, gauss_shape_functions,
            std::tie(mesh_velocity, MESH_VELOCITY));

        r_current_data.CalculateGaussPointData(gauss_shape_functions, r_shape_derivatives);
        const auto& fluid_velocity = r_current_data.GetEffectiveVelocity();
        const ArrayD& velocity = fluid_velocity - mesh_velocity;  // convective or effective velocity
        const double effective_kinematic_viscosity = r_current_data.GetEffectiveKinematicViscosity();
        const double reaction = r_current_data.GetReactionTerm();

        noalias(velocity_convective_terms) = prod(r_shape_derivatives, velocity);

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
    const double residual_factor = this->CalculatePrimalDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
    double diagonal_coefficient = ConvectionDiffusionReactionStabilizationUtilities::CalculatePositivityPreservingMatrix(rDampingMatrix);
    const double acceleration_factor = rCurrentProcessInfo[RANS_RESIDUAL_BASED_FLUX_CORRECTED_ACCELERATON_FACTOR];
    const double scaling_factor = acceleration_factor * std::log(residual_factor + 1);
    const double domain_size = this->GetGeometry().DomainSize();

    const double discrete_upwind_operator_coefficient =
        rCurrentProcessInfo[RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT];
    const double diagonal_positivity_preserving_coefficient =
        rCurrentProcessInfo[RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT];

    BoundedMatrix<double, TNumNodes, TNumNodes> stabilization_diffusion_matrix;
    ConvectionDiffusionReactionStabilizationUtilities::CalculateDiscreteUpwindOperator<TNumNodes>(stabilization_diffusion_matrix, rDampingMatrix);

    this->SetValue(ERROR_OVERALL, scaling_factor);
    this->SetValue(LAMBDA, residual_factor);
    this->SetValue(RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, norm_frobenius(stabilization_diffusion_matrix) / domain_size);
    this->SetValue(RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, diagonal_coefficient / domain_size);

    diagonal_coefficient *= diagonal_positivity_preserving_coefficient * residual_factor;

    noalias(rDampingMatrix) += stabilization_diffusion_matrix *
                             (discrete_upwind_operator_coefficient * residual_factor);
    noalias(rDampingMatrix) += IdentityMatrix(TNumNodes) * (diagonal_coefficient);

    // noalias(rDampingMatrix) += (stabilization_diffusion_matrix + IdentityMatrix(TNumNodes) * diagonal_coefficient) * residual_factor;
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

    array_1d<double, TDim> variable_gradient;

    const Variable<double>& primal_variable =
        TConvectionDiffusionReactionData::GetScalarVariable();
    const Variable<double>& relaxed_primal_rate_variable =
        primal_variable.GetTimeDerivative().GetTimeDerivative();

    const auto& r_geometry = this->GetGeometry();
    TConvectionDiffusionReactionData r_current_data(r_geometry, this->GetProperties(), rCurrentProcessInfo);
    double variable_value, relaxed_variable_acceleration;

    r_current_data.CalculateConstants(rCurrentProcessInfo);

    BoundedVector<double, TNumNodes> velocity_convective_terms;
    ArrayD  mesh_velocity; // added for Chimera_RANS simulation

    double scalar_multiplier = 0.0;
    for (IndexType g = 0; g < num_gauss_points; ++g) {
        const Matrix& r_shape_derivatives = shape_derivatives[g];
        const Vector& gauss_shape_functions = row(shape_functions, g);

        // Get MESH_VELOCITY for Chimera_RANS simulation
        FluidCalculationUtilities::EvaluateInPoint(
            r_geometry, gauss_shape_functions,
            std::tie(mesh_velocity, MESH_VELOCITY));

        r_current_data.CalculateGaussPointData(gauss_shape_functions, r_shape_derivatives);
        const auto& fluid_velocity = r_current_data.GetEffectiveVelocity();
        const ArrayD& velocity = fluid_velocity - mesh_velocity;  // convective or effective velocity
        const double effective_kinematic_viscosity = r_current_data.GetEffectiveKinematicViscosity();
        const double reaction = r_current_data.GetReactionTerm();
        const double source = r_current_data.GetSourceTerm();

        noalias(velocity_convective_terms) = prod(r_shape_derivatives, velocity);
        const double velocity_magnitude = norm_2(velocity);

        const double tau = ConvectionDiffusionReactionStabilizationUtilities::CalculateStabilizationTau(
            element_length, velocity_magnitude, reaction, effective_kinematic_viscosity,
            bossak_alpha, bossak_gamma, delta_time, dynamic_tau);

        FluidCalculationUtilities::EvaluateGradientInPoint(
            r_geometry, r_shape_derivatives,
            std::tie(variable_gradient, primal_variable));

        const double velocity_dot_variable_gradient =
            inner_prod(velocity, variable_gradient);

        FluidCalculationUtilities::EvaluateInPoint(
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

    return scalar_multiplier / static_cast<double>(num_gauss_points);

    KRATOS_CATCH("");
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
GeometryData::IntegrationMethod ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

template <IndexType TDim, IndexType TNumNodes, class TConvectionDiffusionReactionData>
double ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<TDim, TNumNodes, TConvectionDiffusionReactionData>::GetDeltaTime(
    const ProcessInfo& rProcessInfo) const
{
    return rProcessInfo[DELTA_TIME];
}

// template instantiations
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, StabilizationValidationElementData::CircularConvectionElementData>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, StabilizationValidationElementData::BodyForceGovernedCDRElementData>;
template class ConvectionDiffusionReactionResidualBasedFluxCorrectedElement<2, 3, StabilizationValidationElementData::DiffusionElementData>;

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
