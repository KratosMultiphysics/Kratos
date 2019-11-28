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

#if !defined(KRATOS_RANS_EVM_K_EPSILON_ADJOINT_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_EVM_K_EPSILON_ADJOINT_UTILITIES_H_INCLUDED

// System includes

// Project includes
#include "includes/ublas_interface.h"
#include "utilities/geometrical_sensitivity_utility.h"

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
    const Matrix& rShapeDerivatives);

template <unsigned int TDim, unsigned int TNumNodes>
void CalculateProductionShapeSensitivities(
    double& rOutput,
    const double turbulent_kinematic_viscosity,
    const double turbulent_kinematic_viscosity_derivative,
    const Matrix& rNodalVelocity,
    const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
    const Matrix& rShapeDerivatives,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_DxDerivatives);

template <unsigned int TDim, unsigned int TNumNodes>
void CalculateProductionScalarSensitivities(BoundedVector<double, TNumNodes>& rOutput,
                                            const BoundedVector<double, TNumNodes>& rNuTScalarDerivatives,
                                            const BoundedMatrix<double, TDim, TDim>& rVelocityGradient);

template <unsigned int TDim, unsigned int TNumNodes>
void CalculateThetaVelocitySensitivity(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                       const double c_mu,
                                       const double f_mu,
                                       const double tke,
                                       const double nu_t,
                                       const BoundedMatrix<double, TNumNodes, TDim>& rFmuSensitivities,
                                       const BoundedMatrix<double, TNumNodes, TDim>& rNuTSensitivities);

template <unsigned int TNumNodes>
void CalculateThetaTKESensitivity(BoundedVector<double, TNumNodes>& rOutput,
                                  const double c_mu,
                                  const double f_mu,
                                  const double tke,
                                  const double nu_t,
                                  const BoundedVector<double, TNumNodes>& rNuTSensitivities,
                                  const Vector& rGaussShapeFunctions);

template <unsigned int TNumNodes>
void CalculateThetaEpsilonSensitivity(BoundedVector<double, TNumNodes>& rOutput,
                                      const double c_mu,
                                      const double f_mu,
                                      const double tke,
                                      const double nu_t,
                                      const BoundedVector<double, TNumNodes>& rNuTSensitivities);

template <unsigned int TNumNodes>
void CalculateEffectiveKinematicViscosityScalarDerivatives(
    BoundedVector<double, TNumNodes>& rOutput,
    const BoundedVector<double, TNumNodes>& rNutSensitivities,
    const double Sigma,
    const Vector& rGaussShapeFunctions);

} // namespace EvmKepsilonModelAdjointUtilities

///@}

} // namespace Kratos

#endif