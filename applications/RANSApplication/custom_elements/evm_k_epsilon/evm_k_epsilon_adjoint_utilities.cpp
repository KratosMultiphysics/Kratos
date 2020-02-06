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

// System includes

// Project includes
#include "evm_k_epsilon_adjoint_utilities.h"
#include "custom_elements/stabilized_convection_diffusion_reaction_adjoint_utilities.h"
#include "custom_utilities/rans_calculation_utilities.h"
#include "evm_k_epsilon_utilities.h"

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

namespace EvmKepsilonModelAdjointUtilities
{
template <unsigned int TDim, unsigned int TNumNodes>
void CalculateProductionVelocitySensitivities(
    BoundedMatrix<double, TNumNodes, TDim>& rOutput,
    const double NuT,
    const BoundedMatrix<double, TNumNodes, TDim>& rNuTVelocityDerivatives,
    const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
    const Matrix& rShapeDerivatives)
{
    rOutput.clear();

    double velocity_divergence = 0.0;
    velocity_divergence =
        RansCalculationUtilities::CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;
    noalias(reynolds_stress_tensor) = rVelocityGradient + trans(rVelocityGradient) -
                                      (2.0 / 3.0) * velocity_divergence * identity;

    double value = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            value += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);

    noalias(rOutput) = rNuTVelocityDerivatives * value;

    for (std::size_t c = 0; c < TNumNodes; ++c)
    {
        for (std::size_t k = 0; k < TDim; ++k)
        {
            value = 0.0;

            for (std::size_t j = 0; j < TDim; ++j)
            {
                value += rShapeDerivatives(c, j) * rVelocityGradient(k, j);
                value += rShapeDerivatives(c, j) * rVelocityGradient(j, k);
                value -= rShapeDerivatives(c, k) * (2.0 / 3.0) * rVelocityGradient(j, j);
                value += rShapeDerivatives(c, j) * reynolds_stress_tensor(k, j);
            }

            rOutput(c, k) += NuT * value;
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
void CalculateProductionShapeSensitivities(
    double& rOutput,
    const double turbulent_kinematic_viscosity,
    const double turbulent_kinematic_viscosity_derivative,
    const Matrix& rNodalVelocity,
    const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
    const Matrix& rShapeDerivatives,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_DxDerivatives)
{
    rOutput = 0.0;

    double velocity_divergence = 0.0;
    velocity_divergence =
        RansCalculationUtilities::CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;
    noalias(reynolds_stress_tensor) = rVelocityGradient + trans(rVelocityGradient) -
                                      (2.0 / 3.0) * velocity_divergence * identity;

    double value = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            value += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);

    const double negative_two_thirds = -2.0 / 3.0;

    for (std::size_t a = 0; a < TNumNodes; ++a)
    {
        const Vector& r_dna_dx_derivative = row(rDN_DxDerivatives, a);
        const Vector& r_velocity_a = row(rNodalVelocity, a);
        const Vector& r_dna_dx = row(rShapeDerivatives, a);
        for (std::size_t b = 0; b < TNumNodes; ++b)
        {
            const Vector& r_dnb_dx_derivative = row(rDN_DxDerivatives, b);
            const Vector& r_dnb_dx = row(rShapeDerivatives, b);
            const Vector& r_velocity_b = row(rNodalVelocity, b);

            const double uai_ubi = inner_prod(r_velocity_a, r_velocity_b);

            rOutput += inner_prod(r_dna_dx_derivative, r_dnb_dx) * uai_ubi;
            rOutput += inner_prod(r_dna_dx, r_dnb_dx_derivative) * uai_ubi;
            rOutput += inner_prod(r_velocity_a, r_dnb_dx_derivative) *
                       inner_prod(r_velocity_b, r_dna_dx);
            rOutput += inner_prod(r_velocity_a, r_dnb_dx) *
                       inner_prod(r_velocity_b, r_dna_dx_derivative);
            rOutput += negative_two_thirds * inner_prod(r_velocity_a, r_dna_dx_derivative) *
                       inner_prod(r_velocity_b, r_dnb_dx);
            rOutput += negative_two_thirds * inner_prod(r_velocity_a, r_dna_dx) *
                       inner_prod(r_velocity_b, r_dnb_dx_derivative);
        }
    }

    rOutput *= turbulent_kinematic_viscosity;
    rOutput += turbulent_kinematic_viscosity_derivative * value;
}

template <unsigned int TDim, unsigned int TNumNodes>
void CalculateProductionScalarSensitivities(BoundedVector<double, TNumNodes>& rOutput,
                                            const BoundedVector<double, TNumNodes>& rNuTScalarDerivatives,
                                            const BoundedMatrix<double, TDim, TDim>& rVelocityGradient)
{
    double velocity_divergence = 0.0;
    velocity_divergence =
        RansCalculationUtilities::CalculateMatrixTrace<TDim>(rVelocityGradient);
    identity_matrix<double> identity(TDim);

    BoundedMatrix<double, TDim, TDim> reynolds_stress_tensor;
    noalias(reynolds_stress_tensor) = rVelocityGradient + trans(rVelocityGradient) -
                                      (2.0 / 3.0) * velocity_divergence * identity;

    double value = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
        for (unsigned int j = 0; j < TDim; ++j)
            value += reynolds_stress_tensor(i, j) * rVelocityGradient(i, j);

    noalias(rOutput) = rNuTScalarDerivatives * value;
}

template <unsigned int TDim, unsigned int TNumNodes>
void CalculateThetaVelocitySensitivity(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                       const double c_mu,
                                       const double f_mu,
                                       const double tke,
                                       const double nu_t,
                                       const BoundedMatrix<double, TNumNodes, TDim>& rFmuSensitivities,
                                       const BoundedMatrix<double, TNumNodes, TDim>& rNuTSensitivities)
{
    rOutput.clear();
    const double gamma = EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);
    noalias(rOutput) += rFmuSensitivities * gamma / f_mu;
    noalias(rOutput) -= rNuTSensitivities * (gamma / nu_t);
}

template <unsigned int TNumNodes>
void CalculateThetaTKESensitivity(BoundedVector<double, TNumNodes>& rOutput,
                                  const double c_mu,
                                  const double f_mu,
                                  const double tke,
                                  const double nu_t,
                                  const BoundedVector<double, TNumNodes>& rNuTGaussSensitivities,
                                  const Vector& rGaussShapeFunctions)
{
    rOutput.clear();

    const double gamma = EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);
    noalias(rOutput) += rGaussShapeFunctions * gamma / tke;
    noalias(rOutput) -= rNuTGaussSensitivities * (gamma / nu_t);
}

template <unsigned int TNumNodes>
void CalculateThetaEpsilonSensitivity(BoundedVector<double, TNumNodes>& rOutput,
                                      const double c_mu,
                                      const double f_mu,
                                      const double tke,
                                      const double nu_t,
                                      const BoundedVector<double, TNumNodes>& rNuTSensitivities)
{
    const double gamma = EvmKepsilonModelUtilities::CalculateGamma(c_mu, f_mu, tke, nu_t);
    noalias(rOutput) = rNuTSensitivities * (-1.0 * gamma / nu_t);
}

template <unsigned int TNumNodes>
void CalculateEffectiveKinematicViscosityScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const BoundedVector<double, TNumNodes>& rNutSensitivities,
    const double Sigma,
    const Vector& rGaussShapeFunctions)
{
    StabilizedConvectionDiffusionReactionAdjointUtilities::CalculateGaussSensitivities(
        rOutput, rNutSensitivities, rGaussShapeFunctions);

    noalias(rOutput) = rOutput / Sigma;
}

// template instantiations
template void CalculateEffectiveKinematicViscosityScalarDerivatives<3>(
    BoundedVector<double, 3>&, const BoundedVector<double, 3>&, const double, const Vector&);
template void CalculateEffectiveKinematicViscosityScalarDerivatives<4>(
    BoundedVector<double, 4>&, const BoundedVector<double, 4>&, const double, const Vector&);
template void CalculateEffectiveKinematicViscosityScalarDerivatives<8>(
    BoundedVector<double, 8>&, const BoundedVector<double, 8>&, const double, const Vector&);

template void CalculateThetaVelocitySensitivity<2, 3>(BoundedMatrix<double, 3, 2>&,
                                                      const double,
                                                      const double,
                                                      const double,
                                                      const double,
                                                      const BoundedMatrix<double, 3, 2>&,
                                                      const BoundedMatrix<double, 3, 2>&);
template void CalculateThetaVelocitySensitivity<2, 4>(BoundedMatrix<double, 4, 2>&,
                                                      const double,
                                                      const double,
                                                      const double,
                                                      const double,
                                                      const BoundedMatrix<double, 4, 2>&,
                                                      const BoundedMatrix<double, 4, 2>&);
template void CalculateThetaVelocitySensitivity<3, 4>(BoundedMatrix<double, 4, 3>&,
                                                      const double,
                                                      const double,
                                                      const double,
                                                      const double,
                                                      const BoundedMatrix<double, 4, 3>&,
                                                      const BoundedMatrix<double, 4, 3>&);
template void CalculateThetaVelocitySensitivity<3, 8>(BoundedMatrix<double, 8, 3>&,
                                                      const double,
                                                      const double,
                                                      const double,
                                                      const double,
                                                      const BoundedMatrix<double, 8, 3>&,
                                                      const BoundedMatrix<double, 8, 3>&);

template void CalculateThetaTKESensitivity<3>(BoundedVector<double, 3>&,
                                              const double,
                                              const double,
                                              const double,
                                              const double,
                                              const BoundedVector<double, 3>&,
                                              const Vector&);
template void CalculateThetaTKESensitivity<4>(BoundedVector<double, 4>&,
                                              const double,
                                              const double,
                                              const double,
                                              const double,
                                              const BoundedVector<double, 4>&,
                                              const Vector&);
template void CalculateThetaTKESensitivity<8>(BoundedVector<double, 8>&,
                                              const double,
                                              const double,
                                              const double,
                                              const double,
                                              const BoundedVector<double, 8>&,
                                              const Vector&);
template void CalculateThetaEpsilonSensitivity<3>(BoundedVector<double, 3>&,
                                                  const double,
                                                  const double,
                                                  const double,
                                                  const double,
                                                  const BoundedVector<double, 3>&);
template void CalculateThetaEpsilonSensitivity<4>(BoundedVector<double, 4>&,
                                                  const double,
                                                  const double,
                                                  const double,
                                                  const double,
                                                  const BoundedVector<double, 4>&);
template void CalculateThetaEpsilonSensitivity<8>(BoundedVector<double, 8>&,
                                                  const double,
                                                  const double,
                                                  const double,
                                                  const double,
                                                  const BoundedVector<double, 8>&);

template void CalculateProductionVelocitySensitivities<2, 3>(
    BoundedMatrix<double, 3, 2>&,
    const double,
    const BoundedMatrix<double, 3, 2>&,
    const BoundedMatrix<double, 2, 2>&,
    const Matrix&);
template void CalculateProductionVelocitySensitivities<2, 4>(
    BoundedMatrix<double, 4, 2>&,
    const double,
    const BoundedMatrix<double, 4, 2>&,
    const BoundedMatrix<double, 2, 2>&,
    const Matrix&);
template void CalculateProductionVelocitySensitivities<3, 4>(
    BoundedMatrix<double, 4, 3>&,
    const double,
    const BoundedMatrix<double, 4, 3>&,
    const BoundedMatrix<double, 3, 3>&,
    const Matrix&);
template void CalculateProductionVelocitySensitivities<3, 8>(
    BoundedMatrix<double, 8, 3>&,
    const double,
    const BoundedMatrix<double, 8, 3>&,
    const BoundedMatrix<double, 3, 3>&,
    const Matrix&);
template void CalculateProductionShapeSensitivities<2, 3>(
    double&,
    const double,
    const double,
    const Matrix&,
    const BoundedMatrix<double, 2, 2>&,
    const Matrix&,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType&);
template void CalculateProductionShapeSensitivities<2, 4>(
    double&,
    const double,
    const double,
    const Matrix&,
    const BoundedMatrix<double, 2, 2>&,
    const Matrix&,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType&);
template void CalculateProductionShapeSensitivities<3, 4>(
    double&,
    const double,
    const double,
    const Matrix&,
    const BoundedMatrix<double, 3, 3>&,
    const Matrix&,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType&);
template void CalculateProductionShapeSensitivities<3, 8>(
    double&,
    const double,
    const double,
    const Matrix&,
    const BoundedMatrix<double, 3, 3>&,
    const Matrix&,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType&);
template void CalculateProductionScalarSensitivities<2, 3>(
    BoundedVector<double, 3>&,
    const BoundedVector<double, 3>&,
    const BoundedMatrix<double, 2, 2>&);
template void CalculateProductionScalarSensitivities<2, 4>(
    BoundedVector<double, 4>&,
    const BoundedVector<double, 4>&,
    const BoundedMatrix<double, 2, 2>&);

template void CalculateProductionScalarSensitivities<3, 4>(
    BoundedVector<double, 4>&,
    const BoundedVector<double, 4>&,
    const BoundedMatrix<double, 3, 3>&);

template void CalculateProductionScalarSensitivities<3, 8>(
    BoundedVector<double, 8>&,
    const BoundedVector<double, 8>&,
    const BoundedMatrix<double, 3, 3>&);

} // namespace EvmKepsilonModelAdjointUtilities

///@}

} // namespace Kratos
