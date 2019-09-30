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
void CalculateGaussSensitivities(Matrix& rGaussSensitivities,
                                 const Matrix& rNodalSensitivities,
                                 const Vector& rGaussShapeFunctions);

void CalculateGaussSensitivities(Vector& rGaussSensitivities,
                                 const Vector& rNodalSensitivities,
                                 const Vector& rGaussShapeFunctions);

void CalculateNodalFmuVectorSensitivities(Matrix& rFmuNodalSensitivities,
                                          const Vector& nodal_y_plus,
                                          const Matrix& rYPlusNodalSensitivities);

void CalculateGaussFmuVectorSensitivities(Matrix& rFmuGaussSensitivities,
                                          const double y_plus,
                                          const Matrix& rYPlusNodalSensitivities,
                                          const Vector& rGaussShapeFunctions);

void CalculateNodalTurbulentViscosityVectorSensitivities(
    Matrix& rTurbulentViscosityNodalSensitivities,
    const double c_mu,
    const Vector& nodal_turbulent_kinetic_energy,
    const Vector& nodal_turbulent_energy_dissipation_rate,
    const Matrix& rFmuNodalSensitivities);

void CalculateNodalTurbulentViscosityTKESensitivities(Vector& rTurbulentViscosityNodalSensitivities,
                                                      const double c_mu,
                                                      const Vector& nodal_turbulent_kinetic_energy,
                                                      const Vector& nodal_turbulent_energy_dissipation_rate,
                                                      const Vector& nodal_f_mu);

void CalculateNodalTurbulentViscosityEpsilonSensitivities(
    Vector& rTurbulentViscosityNodalSensitivities,
    const double c_mu,
    const Vector& nodal_turbulent_kinetic_energy,
    const Vector& nodal_turbulent_energy_dissipation_rate,
    const Vector& nodal_f_mu);

template <unsigned int TDim>
void CalculateProductionVelocitySensitivities(Matrix& rOutput,
                                              const double NuT,
                                              const Matrix& rNuTVelocityDerivatives,
                                              const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                                              const Matrix& rShapeDerivatives);

template <unsigned int TDim>
void CalculateProductionShapeSensitivities(
    double& rOutput,
    const double turbulent_kinematic_viscosity,
    const double turbulent_kinematic_viscosity_derivative,
    const Matrix& rNodalVelocity,
    const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
    const Matrix& rShapeDerivatives,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_DxDerivatives);

template <unsigned int TDim>
void CalculateProductionScalarSensitivities(Vector& rOutput,
                                            const Vector& rNuTScalarDerivatives,
                                            const BoundedMatrix<double, TDim, TDim>& rVelocityGradient);

void CalculateThetaVelocitySensitivity(Matrix& rOutput,
                                       const double c_mu,
                                       const double f_mu,
                                       const double tke,
                                       const double nu_t,
                                       const Matrix& rFmuSensitivities,
                                       const Matrix& rNuTSensitivities);
void CalculateThetaTKESensitivity(Vector& rOutput,
                                  const double c_mu,
                                  const double f_mu,
                                  const double tke,
                                  const double nu_t,
                                  const Vector& rNuTSensitivities,
                                  const Vector& rGaussShapeFunctions);

void CalculateThetaEpsilonSensitivity(Vector& rOutput,
                                      const double c_mu,
                                      const double f_mu,
                                      const double tke,
                                      const double nu_t,
                                      const Vector& rNuTSensitivities);
void CalculateTurbulentReynoldsNumberVelocitySensitivity(Matrix& rOutput,
                                                         const double tke,
                                                         const double epsilon,
                                                         const double nu,
                                                         const Matrix& rNuTSensitivities);

void CalculateTurbulentReynoldsNumberTKESensitivity(Vector& rOutput,
                                                    const double tke,
                                                    const double epsilon,
                                                    const double nu,
                                                    const Vector& rGaussShapeFunctions);

void CalculateTurbulentReynoldsNumberEpsilonSensitivity(Vector& rOutput,
                                                        const double tke,
                                                        const double epsilon,
                                                        const double nu,
                                                        const Vector& rGaussShapeFunctions);

void CalculateF2VelocitySensitivity(Matrix& rOutput,
                                    const double tke,
                                    const double epsilon,
                                    const double nu_t,
                                    const Matrix& rNuTSensitivities);

void CalculateF2ScalarSensitivity(Vector& rOutput,
                                  const double epsilon,
                                  const double Re_t,
                                  const Vector& rReTSensitivities);

} // namespace EvmKepsilonModelAdjointUtilities

///@}

} // namespace Kratos

#endif