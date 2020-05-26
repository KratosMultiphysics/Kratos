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

#if !defined(KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_UTILITIES)
#define KRATOS_STABILIZED_CONVECTION_DIFFUSION_REACTION_ADJOINT_UTILITIES

#include "custom_utilities/rans_calculation_utilities.h"
#include "includes/define.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

namespace StabilizedConvectionDiffusionReactionAdjointUtilities
{
/**
 * @brief Calculates stabilization tau scalar derivatives
 *
 * \[
 *  -\tau_\phi^3\left[144\frac{\nu_\phi}{h^4_2}\left(\nu_{\phi,w}\right)^c + s_\phi\left(s_{\phi,w}\right)^c\right]
 * \]
 *
 * Where $w$ is the derivative variable
 *
 * @param rOutput                                        Scalar derivatives for each node w.r.t. $w$
 * @param Tau                                            Stabilization tau
 * @param EffectiveKinematicViscosity                    Effective kinematic viscosity $\nu_\phi$
 * @param Reaction                                       Reaction coefficient $s_\phi$
 * @param ElementLength                                  Element length $h_2$
 * @param rEffectiveKinematicViscosityScalarDerivatives  Scalar derivatives of effective kinematic viscosity $\left(\nu_{\phi,w}\right)$
 * @param rReactionScalarDerivatives                     Reaction scalar derivatives $\left(s_{\phi,w}\right)$
 */
template <std::size_t TNumNodes>
inline void CalculateStabilizationTauScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double Tau,
    const double EffectiveKinematicViscosity,
    const double Reaction,
    const double ElementLength,
    const BoundedVector<double, TNumNodes>& rEffectiveKinematicViscosityScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rReactionScalarDerivatives)
{
    noalias(rOutput) =
        (rEffectiveKinematicViscosityScalarDerivatives *
             (144 * EffectiveKinematicViscosity / std::pow(ElementLength, 4)) +
         rReactionScalarDerivatives * (Reaction)) *
        (-1.0 * std::pow(Tau, 3));
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateStabilizationTauVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double Tau,
    const double EffectiveKinematicViscosity,
    const double Reaction,
    const double ElementLength,
    const array_1d<double, 3>& rVelocity,
    const BoundedMatrix<double, TDim, TDim>& rContravariantMetricTensor,
    const BoundedMatrix<double, TNumNodes, TDim>& rEffectiveKinematicViscosityVelocityDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rReactionVelocityDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rElementLengthDerivatives,
    const Vector& rGaussShapeFunctions)
{
    Vector contravariant_metric_velocity(TDim);
    const Vector& velocity = RansCalculationUtilities::GetVector<TDim>(rVelocity);

    noalias(contravariant_metric_velocity) =
        prod(rContravariantMetricTensor, velocity) +
        prod(trans(rContravariantMetricTensor), velocity);

    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
        for (std::size_t i_dim = 0; i_dim < TDim; ++i_dim)
            rOutput(i_node, i_dim) = 0.5 * rGaussShapeFunctions[i_node] *
                                     contravariant_metric_velocity[i_dim];

    noalias(rOutput) +=
        rEffectiveKinematicViscosityVelocityDerivatives *
        (144.0 * EffectiveKinematicViscosity / std::pow(ElementLength, 4));
    noalias(rOutput) -=
        rElementLengthDerivatives *
        (288.0 * std::pow(EffectiveKinematicViscosity, 2) / std::pow(ElementLength, 5));
    noalias(rOutput) += rReactionVelocityDerivatives * (Reaction);
    noalias(rOutput) = rOutput * (-1.0 * std::pow(Tau, 3));
}

inline double CalculateStabilizationTauShapeSensitivity(const double Tau,
                                                        const double VelocityMagnitude,
                                                        const double ElementLength,
                                                        const double ElementLengthDeriv,
                                                        const double EffectiveKinematicViscosity,
                                                        const double EffectiveKinematicViscosityDeriv,
                                                        const double Reaction,
                                                        const double ReactionDeriv)
{
    double shape_sensitivity = 0.0;

    shape_sensitivity += 4.0 * std::pow(VelocityMagnitude, 2) *
                         ElementLengthDeriv / std::pow(ElementLength, 3);
    shape_sensitivity -= 144.0 * EffectiveKinematicViscosity *
                         EffectiveKinematicViscosityDeriv / std ::pow(ElementLength, 4);
    shape_sensitivity += 288.0 * std::pow(EffectiveKinematicViscosity, 2) *
                         ElementLengthDeriv / std::pow(ElementLength, 5);
    shape_sensitivity -= Reaction * ReactionDeriv;

    shape_sensitivity *= std::pow(Tau, 3);

    return shape_sensitivity;
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateGaussSensitivities(BoundedMatrix<double, TNumNodes, TDim>& rGaussSensitivities,
                                        const BoundedMatrix<double, TNumNodes, TDim>& rNodalSensitivities,
                                        const Vector& rGaussShapeFunctions)
{
    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
    {
        for (std::size_t i_dim = 0; i_dim < TDim; ++i_dim)
        {
            rGaussSensitivities(i_node, i_dim) =
                rGaussShapeFunctions[i_node] * rNodalSensitivities(i_node, i_dim);
        }
    }
}

template <std::size_t TNumNodes>
inline void CalculateGaussSensitivities(BoundedVector<double, TNumNodes>& rGaussSensitivities,
                                        const BoundedVector<double, TNumNodes>& rNodalSensitivities,
                                        const Vector& rGaussShapeFunctions)
{
    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
    {
        rGaussSensitivities[i_node] =
            rGaussShapeFunctions[i_node] * rNodalSensitivities[i_node];
    }
}

inline double CalculateScalarProduct(const Vector& rVector1, const array_1d<double, 3>& rVector2)
{
    double result = 0.0;
    for (std::size_t i_dim = 0; i_dim < rVector1.size(); ++i_dim)
        result += rVector1[i_dim] * rVector2[i_dim];
    return result;
}

template <std::size_t TNumNodes>
inline void CalculateAbsoluteScalarValueScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double scalar_value,
    const BoundedVector<double, TNumNodes>& rScalarValueDerivatives)
{
    if (scalar_value >= 0.0)
    {
        noalias(rOutput) = rScalarValueDerivatives;
    }
    else
    {
        noalias(rOutput) = rScalarValueDerivatives * (-1.0);
    }
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateAbsoluteScalarValueVectorDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double scalar_value,
    const BoundedMatrix<double, TNumNodes, TDim>& rScalarValueDerivatives)
{
    if (scalar_value >= 0)
    {
        noalias(rOutput) = rScalarValueDerivatives;
    }
    else
    {
        noalias(rOutput) = rScalarValueDerivatives * -1.0;
    }
}

template <std::size_t TNumNodes>
inline void CalculateAbsoluteScalarGradientScalarDerivative(
    BoundedVector<double, TNumNodes>& rOutput,
    const array_1d<double, 3> rScalarGradient,
    const Matrix& rShapeFunctionDerivatives)
{
    const double scalar_gradient_norm = norm_2(rScalarGradient);

    if (scalar_gradient_norm > 0.0)
    {
        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
        {
            const Vector& shape_function_gradient = row(rShapeFunctionDerivatives, i_node);
            rOutput[i_node] =
                CalculateScalarProduct(shape_function_gradient, rScalarGradient) / scalar_gradient_norm;
        }
    }
    else
    {
        rOutput.clear();
    }
}

inline double CalculateAbsoluteScalarGradientShapeSensitivity(
    const array_1d<double, 3>& rScalarGradient,
    const Matrix& rShapeFunctionDerivShapeSensitivity,
    const Vector& rNodalScalarValues)
{
    const double scalar_gradient_norm = norm_2(rScalarGradient);

    if (scalar_gradient_norm > 0.0)
    {
        Vector scalar_gradient_shape_sensitivity(
            rShapeFunctionDerivShapeSensitivity.size2());
        noalias(scalar_gradient_shape_sensitivity) =
            prod(trans(rShapeFunctionDerivShapeSensitivity), rNodalScalarValues);

        Vector scalar_gradient(rShapeFunctionDerivShapeSensitivity.size2());
        for (std::size_t i = 0; i < rShapeFunctionDerivShapeSensitivity.size2(); ++i)
            scalar_gradient[i] = rScalarGradient[i];
        return inner_prod(scalar_gradient_shape_sensitivity, scalar_gradient) / scalar_gradient_norm;
    }
    else
    {
        return 0.0;
    }
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateVelocityMagnitudeVelocityDerivative(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double VelocityMagnitude,
    const array_1d<double, 3>& rVelocity,
    const Vector& rGaussShapeFunctions)
{
    if (VelocityMagnitude > 0.0)
    {
        for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
            for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
                rOutput(i_node, i_dim) =
                    rVelocity[i_dim] * rGaussShapeFunctions[i_node] / VelocityMagnitude;
    }
    else
    {
        rOutput.clear();
    }
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateElementLengthH2VelocityDerivative(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double VelocityMagnitude,
    const array_1d<double, 3>& rVelocity,
    const BoundedMatrix<double, TNumNodes, TDim>& rVelocityMagnitudeVelocityDerivatives,
    const BoundedMatrix<double, TDim, TDim>& rContravariantMetricTensor,
    const Vector& rGaussShapeFunctions)
{
    if (VelocityMagnitude > 0.0)
    {
        const Vector& velocity = RansCalculationUtilities::GetVector<TDim>(rVelocity);

        const double sqrt_u_e_u =
            std::sqrt(inner_prod(velocity, prod(rContravariantMetricTensor, velocity)));

        Vector contravariant_metric_velocity(TDim);
        noalias(contravariant_metric_velocity) =
            prod(rContravariantMetricTensor, velocity) +
            prod(trans(rContravariantMetricTensor), velocity);

        for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
            for (std::size_t i_dim = 0; i_dim < TDim; ++i_dim)
                rOutput(i_node, i_dim) = rGaussShapeFunctions[i_node] *
                                         contravariant_metric_velocity[i_dim];

        noalias(rOutput) =
            rOutput * (-1.0 * VelocityMagnitude / std::pow(sqrt_u_e_u, 3));
        noalias(rOutput) += (rVelocityMagnitudeVelocityDerivatives) * (2.0 / sqrt_u_e_u);
    }
    else
    {
        rOutput.clear();
    }
}

template <std::size_t TDim>
inline double CalculateElementLengthH2ShapeSensitivity(
    const double VelocityMagnitude,
    const array_1d<double, 3>& rVelocity,
    const BoundedMatrix<double, TDim, TDim>& rContravariantMetricTensor,
    const BoundedMatrix<double, TDim, TDim>& rContravariantMetricTensorShapeSensitivity)
{
    if (VelocityMagnitude > 0.0)
    {
        const Vector& velocity = RansCalculationUtilities::GetVector<TDim>(rVelocity);

        const double u_e_u = std::pow(
            inner_prod(velocity, prod(rContravariantMetricTensor, velocity)), 1.5);

        return -VelocityMagnitude *
               (inner_prod(velocity, prod(rContravariantMetricTensorShapeSensitivity, velocity))) /
               u_e_u;
    }
    else
    {
        double sensitivity = 0.0;
        double element_length = 0.0;
        for (unsigned int i = 0; i < TDim; ++i)
            for (unsigned int j = 0; j < TDim; ++j)
            {
                sensitivity += rContravariantMetricTensorShapeSensitivity(i, j);
                element_length += rContravariantMetricTensor(i, j);
            }
        element_length = std::sqrt(1.0 / element_length) * 2.0;

        return sensitivity * std::pow(element_length, 3) * (-1.0 / 8.0);
    }
}

template <std::size_t TNumNodes>
inline void CalculatePsiOneScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double velocity_norm,
    const double reaction_tilde,
    const double tau,
    const BoundedVector<double, TNumNodes>& rTauScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rAbsoluteReactionTildeScalarDerivatives)
{
    const double absolute_reaction_tilde = std::abs(reaction_tilde);

    noalias(rOutput) = rTauScalarDerivatives * (velocity_norm * absolute_reaction_tilde);
    noalias(rOutput) += rAbsoluteReactionTildeScalarDerivatives * (tau * velocity_norm);
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculatePsiOneVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double velocity_norm,
    const double reaction_tilde,
    const double tau,
    const BoundedMatrix<double, TNumNodes, TDim>& rTauDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rAbsoluteReactionTildeDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rVelocityMagnitudeDerivatives)
{
    noalias(rOutput) = rVelocityMagnitudeDerivatives +
                       rTauDerivatives * (velocity_norm * reaction_tilde) +
                       rVelocityMagnitudeDerivatives * (tau * reaction_tilde) +
                       rAbsoluteReactionTildeDerivatives * (tau * velocity_norm);
}

inline double CalculatePsiOneShapeSensitivity(const double tau,
                                              const double tau_deriv,
                                              const double velocity_magnitude,
                                              const double reaction,
                                              const double reaction_deriv,
                                              const double bossak_alpha,
                                              const double bossak_gamma,
                                              const double delta_time,
                                              const double DynamicTau)
{
    const double reaction_dynamics =
        reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
    const double abs_reaction_dynamics = std::abs(reaction_dynamics);

    if (abs_reaction_dynamics > 0.0)
    {
        return tau_deriv * velocity_magnitude * abs_reaction_dynamics +
               tau * velocity_magnitude * reaction_dynamics * reaction_deriv / abs_reaction_dynamics;
    }
    else
    {
        return 0.0;
    }
}

template <std::size_t TNumNodes>
inline void CalculatePsiTwoScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double element_length,
    const double tau,
    const double reaction_tilde,
    const BoundedVector<double, TNumNodes>& rTauScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rReactionTildeDerivatives,
    const BoundedVector<double, TNumNodes>& rAbsoluteReactionTildeScalarDerivatives)
{
    const double absolute_reaction_tilde = std::abs(reaction_tilde);

    noalias(rOutput) = rReactionTildeDerivatives;
    noalias(rOutput) += rTauScalarDerivatives * (reaction_tilde * absolute_reaction_tilde);
    noalias(rOutput) += rReactionTildeDerivatives * (tau * absolute_reaction_tilde);
    noalias(rOutput) += rAbsoluteReactionTildeScalarDerivatives * (tau * reaction_tilde);
    noalias(rOutput) = rOutput * (std::pow(element_length, 2) / 6.0);
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculatePsiTwoVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double reaction_tilde,
    const double tau,
    const double element_length,
    const BoundedMatrix<double, TNumNodes, TDim>& rTauDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rReactionTildeDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rAbsoluteReactionTildeDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rElementLengthDerivatives)
{
    const double abs_reaction_tilde = std::abs(reaction_tilde);

    noalias(rOutput) =
        (rReactionTildeDerivatives + rTauDerivatives * (reaction_tilde * abs_reaction_tilde) +
         rReactionTildeDerivatives * (tau * abs_reaction_tilde) +
         rAbsoluteReactionTildeDerivatives * (tau * reaction_tilde)) *
            std::pow(element_length, 2) / 6.0 +
        rElementLengthDerivatives *
            (element_length *
             (reaction_tilde + tau * reaction_tilde * abs_reaction_tilde) / 3.0);
}

inline double CalculatePsiTwoShapeSensitivity(const double psi_two,
                                              const double element_length,
                                              const double element_length_deriv,
                                              const double reaction,
                                              const double reaction_deriv,
                                              const double tau,
                                              const double tau_deriv,
                                              const double bossak_alpha,
                                              const double bossak_gamma,
                                              const double delta_time,
                                              const double DynamicTau)
{
    double shape_sensitivity = 0.0;

    const double reaction_dynamics =
        reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
    const double abs_reaction_dynamics = std::abs(reaction_dynamics);

    shape_sensitivity += reaction_deriv;
    if (abs_reaction_dynamics > 0.0)
    {
        shape_sensitivity += tau_deriv * reaction_dynamics * abs_reaction_dynamics;
        shape_sensitivity += tau * reaction_deriv * abs_reaction_dynamics;
        shape_sensitivity += tau * reaction_dynamics * reaction_dynamics *
                             reaction_deriv / abs_reaction_dynamics;
    }

    shape_sensitivity *= std::pow(element_length, 2) / 6.0;

    shape_sensitivity += 2.0 * psi_two * element_length_deriv / element_length;

    return shape_sensitivity;
}

template <std::size_t TNumNodes>
inline void CalculateChiScalarDerivatives(BoundedVector<double, TNumNodes>& rOutput,
                                          const double Chi,
                                          const double ElementLength,
                                          const double BossakAlpha,
                                          const double BossakGamma,
                                          const double DeltaTime,
                                          const double Reaction,
                                          const double DynamicTau,
                                          const BoundedVector<double, TNumNodes>& rReactionScalarDerivatives)
{
    const double reaction_tilde =
        Reaction + DynamicTau * (1 - BossakAlpha) / (BossakGamma * DeltaTime);

    CalculateAbsoluteScalarValueScalarDerivatives(rOutput, reaction_tilde,
                                                  rReactionScalarDerivatives);
    noalias(rOutput) = rOutput * (-0.5 * std::pow(Chi, 2) * ElementLength);
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateChiVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double Chi,
    const double ElementLength,
    const double BossakAlpha,
    const double BossakGamma,
    const double DeltaTime,
    const double Reaction,
    const double DynamicTau,
    const BoundedMatrix<double, TNumNodes, TDim>& rReactionDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rVelocityMagnitudeDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rElementLengthDerivatives)
{
    const double reaction_tilde =
        Reaction + DynamicTau * (1 - BossakAlpha) / (BossakGamma * DeltaTime);
    const double abs_reaction_tilde = std::abs(reaction_tilde);

    CalculateAbsoluteScalarValueVectorDerivatives(rOutput, reaction_tilde, rReactionDerivatives);

    noalias(rOutput) = (rOutput * ElementLength + rElementLengthDerivatives * abs_reaction_tilde +
                        rVelocityMagnitudeDerivatives * 2.0) *
                       (-0.5 * std::pow(Chi, 2));
}

inline double CalculateChiShapeSensitivity(const double chi,
                                           const double reaction,
                                           const double reaction_deriv,
                                           const double element_length,
                                           const double element_length_deriv,
                                           const double bossak_alpha,
                                           const double bossak_gamma,
                                           const double delta_time,
                                           const double DynamicTau)
{
    const double reaction_tilde =
        reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
    const double abs_reaction_tilde = std::abs(reaction_tilde);

    if (abs_reaction_tilde > 0.0)
    {
        return -0.5 * std::pow(chi, 2) *
               (abs_reaction_tilde * element_length_deriv +
                reaction_tilde * element_length * reaction_deriv / abs_reaction_tilde);
    }
    else
    {
        return 0.0;
    }
}

template <std::size_t TNumNodes>
inline void CalculateResidualScalarDerivative(
    BoundedVector<double, TNumNodes>& rOutput,
    const double scalar_value,
    const double reaction,
    const array_1d<double, 3>& rVelocity,
    const BoundedVector<double, TNumNodes>& rReactionScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rSourceScalarDerivatives,
    const Vector& rShapeFunctions,
    const Matrix& rShapeFunctionDerivatives,
    const Variable<double>& rPrimalVariable,
    const Variable<double>& rDerivativeVariable)
{
    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node)
    {
        const Vector& shape_function_gradient = row(rShapeFunctionDerivatives, i_node);
        double value = 0.0;

        value += scalar_value * rReactionScalarDerivatives[i_node];
        value -= rSourceScalarDerivatives[i_node];

        if (rPrimalVariable == rDerivativeVariable)
        {
            value += reaction * rShapeFunctions[i_node];
            value += CalculateScalarProduct(shape_function_gradient, rVelocity);
        }

        rOutput[i_node] = value;
    }
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateResidualVelocityDerivative(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double primal_variable_value,
    const array_1d<double, 3>& rPrimalVariableGradient,
    const BoundedMatrix<double, TNumNodes, TDim>& rReactionDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rSourceDerivatives,
    const Vector& rGaussShapeFunctions)
{
    for (unsigned int i_node = 0; i_node < TNumNodes; ++i_node)
        for (unsigned int i_dim = 0; i_dim < TDim; ++i_dim)
            rOutput(i_node, i_dim) =
                rGaussShapeFunctions[i_node] * rPrimalVariableGradient[i_dim];

    noalias(rOutput) = rOutput + rReactionDerivatives * primal_variable_value - rSourceDerivatives;
}

inline double CalculateResidualShapeSensitivity(const double residual,
                                                const array_1d<double, 3>& rVelocity,
                                                const Matrix& rShapeFunctionDerivShapeSensitivity,
                                                const double scalar_value,
                                                const Vector& rNodalScalarValues,
                                                const double reaction_deriv,
                                                const double source_deriv)
{
    const double abs_residual = std::abs(residual);

    if (abs_residual > 0.0)
    {
        Vector r_velocity(rShapeFunctionDerivShapeSensitivity.size2());
        for (std::size_t i = 0; i < rShapeFunctionDerivShapeSensitivity.size2(); ++i)
            r_velocity[i] = rVelocity[i];
        // const Vector& r_velocity = RansCalculationUtilities::GetVector<TDim>(rVelocity);
        Vector primal_variable_gradient_shape_sensitivity(
            rShapeFunctionDerivShapeSensitivity.size2());
        noalias(primal_variable_gradient_shape_sensitivity) =
            prod(trans(rShapeFunctionDerivShapeSensitivity), rNodalScalarValues);

        return residual *
               (inner_prod(r_velocity, primal_variable_gradient_shape_sensitivity) +
                reaction_deriv * scalar_value - source_deriv) /
               abs_residual;
    }
    else
    {
        return 0.0;
    }
}

template <std::size_t TNumNodes>
inline void CalculatePositivityPreservationCoefficientScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double chi,
    const double residual,
    const double scalar_gradient_norm,
    const double velocity_norm_square,
    const BoundedVector<double, TNumNodes>& rChiScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rAbsoluteResidualScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rAbsoluteScalarGradientScalarDerivative,
    const Variable<double>& rPrimalVariable,
    const Variable<double>& rDerivativeVariable)
{
    const double abs_residual = std::abs(residual);

    noalias(rOutput) = rAbsoluteResidualScalarDerivatives *
                       (chi / (velocity_norm_square * scalar_gradient_norm));
    noalias(rOutput) += rChiScalarDerivatives *
                        (abs_residual / (velocity_norm_square * scalar_gradient_norm));

    if (rPrimalVariable == rDerivativeVariable)
        noalias(rOutput) -=
            rAbsoluteScalarGradientScalarDerivative *
            (chi * abs_residual / (std::pow(scalar_gradient_norm, 2) * velocity_norm_square));
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculatePositivityPreservationCoefficientVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double absolute_residual,
    const double primal_variable_gradient_norm,
    const double velocity_magnitude,
    const double chi,
    const BoundedMatrix<double, TNumNodes, TDim>& rChiDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rAbsoluteResidualDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rVelocityMagnitudeDerivatives)
{
    const double velocity_magnitude_square = std::pow(velocity_magnitude, 2);

    noalias(rOutput) =
        (rVelocityMagnitudeDerivatives * (-2.0 * chi / velocity_magnitude) + rChiDerivatives) *
            (absolute_residual / (velocity_magnitude_square * primal_variable_gradient_norm)) +
        rAbsoluteResidualDerivatives *
            (chi / (primal_variable_gradient_norm * velocity_magnitude_square));
}

inline double CalculatePositivityPreservationCoefficientShapeSensitivity(
    const double chi,
    const double chi_deriv,
    const double abs_residual,
    const double abs_residual_deriv,
    const double velocity_magnitude_square,
    const double scalar_gradient_norm,
    const double scalar_gradient_norm_deriv)
{
    return chi_deriv * abs_residual / (scalar_gradient_norm * velocity_magnitude_square) +
           chi * abs_residual_deriv / (scalar_gradient_norm * velocity_magnitude_square) -
           chi * abs_residual * scalar_gradient_norm_deriv /
               (std::pow(scalar_gradient_norm, 2) * velocity_magnitude_square);
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateCrossWindDiffusionCoeffVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double cross_wind_diffusion_coefficient,
    const double psi_one,
    const double element_length,
    const BoundedMatrix<double, TNumNodes, TDim>& rPsiOneDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rPsiTwoDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rEffectiveKinematicViscosityDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rElementLengthDerivatives)
{
    if (cross_wind_diffusion_coefficient > 0.0)
    {
        const double abs_psi_one = std::abs(psi_one);

        noalias(rOutput) = rPsiTwoDerivatives - rEffectiveKinematicViscosityDerivatives;

        if (abs_psi_one > 0.0)
        {
            noalias(rOutput) +=
                rPsiOneDerivatives * (0.5 * psi_one * element_length / abs_psi_one) +
                rElementLengthDerivatives * (0.5 * abs_psi_one);
        }
    }
    else
    {
        rOutput.clear();
    }
}

inline double CalculateCrossWindDiffusionCoeffShapeSensitivity(
    const double cross_wind_diffusion_coefficient,
    const double psi_one,
    const double psi_one_deriv,
    const double element_length,
    const double element_length_deriv,
    const double effective_kinematic_viscosity_deriv,
    const double psi_two_deriv)
{
    if (cross_wind_diffusion_coefficient > 0.0)
    {
        const double abs_psi_one = std::abs(psi_one);

        double value = psi_two_deriv - effective_kinematic_viscosity_deriv;

        if (abs_psi_one > 0.0)
        {
            value += 0.5 * psi_one * element_length * psi_one_deriv / abs_psi_one +
                     0.5 * abs_psi_one * element_length_deriv;
        }
        return value;
    }
    else
    {
        return 0.0;
    }
}

template <std::size_t TNumNodes>
inline void CalculateCrossWindDiffusionCoeffScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double cross_wind_diffusion_coefficient,
    const double psi_one,
    const double psi_two,
    const double element_length,
    const double effective_kinematic_viscosity,
    const BoundedVector<double, TNumNodes>& rPsiOneScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rPsiTwoScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rEffectiveKinematicViscosityScalarDerivatives)
{
    if (cross_wind_diffusion_coefficient > 0.0)
    {
        noalias(rOutput) = rPsiTwoScalarDerivatives;
        noalias(rOutput) -= rEffectiveKinematicViscosityScalarDerivatives;

        const double abs_psi_one = std::abs(psi_one);
        if (abs_psi_one > 0.0)
        {
            noalias(rOutput) += rPsiOneScalarDerivatives *
                                (0.5 * psi_one * element_length / abs_psi_one);
        }
    }
    else
    {
        rOutput.clear();
    }
}

template <std::size_t TNumNodes>
inline void CalculateStreamLineDiffusionCoeffScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const double streamline_diffusion_coefficient,
    const double element_length,
    const double tau,
    const double velocity_norm,
    const double reaction_tilde,
    const double psi_one,
    const double psi_two,
    const BoundedVector<double, TNumNodes>& rPsiOneScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rPsiTwoScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rTauScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rReactionTildeScalarDerivatives,
    const BoundedVector<double, TNumNodes>& rEffectiveViscosityScalarDerivatives)
{
    if (streamline_diffusion_coefficient > 0.0)
    {
        noalias(rOutput) = rPsiOneScalarDerivatives;
        noalias(rOutput) -= rTauScalarDerivatives * (velocity_norm * reaction_tilde);
        noalias(rOutput) -= rReactionTildeScalarDerivatives * (tau * velocity_norm);

        const double coeff = psi_one - tau * velocity_norm * reaction_tilde;
        const double abs_coeff = std::abs(coeff);

        if (abs_coeff > 0.0)
        {
            noalias(rOutput) = rOutput * (0.5 * element_length * coeff / abs_coeff);
        }
        else
        {
            rOutput.clear();
        }

        noalias(rOutput) += rPsiTwoScalarDerivatives;
        noalias(rOutput) -= rEffectiveViscosityScalarDerivatives;
        noalias(rOutput) -= rTauScalarDerivatives * std::pow(velocity_norm, 2);
    }
    else
    {
        rOutput.clear();
    }
}

template <std::size_t TDim, std::size_t TNumNodes>
inline void CalculateStreamLineDiffusionCoeffVelocityDerivatives(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double streamline_diffusion_coefficient,
    const double element_length,
    const double tau,
    const double velocity_norm,
    const double reaction_tilde,
    const double psi_one,
    const double psi_two,
    const BoundedMatrix<double, TNumNodes, TDim>& rVelocityMagnitudeDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rPsiOneDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rPsiTwoDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rTauDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rReactionTildeDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rEffectiveViscosityDerivatives,
    const BoundedMatrix<double, TNumNodes, TDim>& rElementLengthDerivatives)
{
    if (streamline_diffusion_coefficient > 0.0)
    {
        noalias(rOutput) =
            (rPsiOneDerivatives - rTauDerivatives * (velocity_norm * reaction_tilde) -
             rVelocityMagnitudeDerivatives * (tau * reaction_tilde) -
             rReactionTildeDerivatives * (tau * velocity_norm));

        const double coeff = psi_one - tau * velocity_norm * reaction_tilde;
        const double abs_coeff = std::abs(coeff);
        if (abs_coeff > 0.0)
        {
            noalias(rOutput) = rOutput * (0.5 * element_length * coeff / abs_coeff);
        }
        else
        {
            rOutput.clear();
        }

        noalias(rOutput) += rElementLengthDerivatives * (0.5 * abs_coeff);

        noalias(rOutput) += rPsiTwoDerivatives;
        noalias(rOutput) -= rEffectiveViscosityDerivatives;
        noalias(rOutput) -= rTauDerivatives * std::pow(velocity_norm, 2);
        noalias(rOutput) -= rVelocityMagnitudeDerivatives * (2.0 * tau * velocity_norm);
    }
    else
    {
        rOutput.clear();
    }
}

inline double CalculateStreamLineDiffusionCoeffShapeSensitivity(
    const double streamline_diffusion_coefficient,
    const double psi_one,
    const double psi_one_deriv,
    const double tau,
    const double tau_deriv,
    const double velocity_magnitude,
    const double reaction,
    const double reaction_deriv,
    const double element_length,
    const double element_length_deriv,
    const double effective_kinematic_viscosity_deriv,
    const double psi_two_deriv,
    const double bossak_alpha,
    const double bossak_gamma,
    const double delta_time,
    const double DynamicTau)
{
    if (streamline_diffusion_coefficient > 0.0)
    {
        const double reaction_dynamics =
            reaction + DynamicTau * (1 - bossak_alpha) / (bossak_gamma * delta_time);
        const double coeff = psi_one - tau * velocity_magnitude * reaction_dynamics;
        const double abs_coeff = std::abs(coeff);
        double shape_sensitivity = 0.0;

        shape_sensitivity += psi_one_deriv - tau_deriv * velocity_magnitude * reaction_dynamics -
                             tau * velocity_magnitude * reaction_deriv;
        if (abs_coeff > 0.0)
        {
            shape_sensitivity *= 0.5 * coeff * element_length / abs_coeff;
        }
        else
        {
            shape_sensitivity = 0.0;
        }

        shape_sensitivity += 0.5 * abs_coeff * element_length_deriv;
        shape_sensitivity -= effective_kinematic_viscosity_deriv +
                             tau_deriv * std::pow(velocity_magnitude, 2);
        shape_sensitivity += psi_two_deriv;

        return shape_sensitivity;
    }
    else
    {
        return 0.0;
    }
}

} // namespace StabilizedConvectionDiffusionReactionAdjointUtilities

} // namespace Kratos
#endif
