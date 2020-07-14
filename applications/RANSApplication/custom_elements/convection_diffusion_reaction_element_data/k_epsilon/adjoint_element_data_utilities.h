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

#if !defined(KRATOS_RANS_K_EPSILON_ADJOINT_ELEMENT_DATA_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_K_EPSILON_ADJOINT_ELEMENT_DATA_UTILITIES_H_INCLUDED

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

namespace KEpsilonElementData
{
template <unsigned int TNumNodes>
void CalculateNodalNutTKESensitivity(BoundedVector<double, TNumNodes>& rOutput,
                                     const double Cmu,
                                     const BoundedVector<double, TNumNodes>& rNodalTurbulentKineticEnergy,
                                     const BoundedVector<double, TNumNodes>& rNodalTurbulentEnergyDissipationRate);

template <unsigned int TNumNodes>
void CalculateNodalNutEpsilonSensitivity(
    BoundedVector<double, TNumNodes>& rOutput,
    const double Cmu,
    const BoundedVector<double, TNumNodes>& rNodalTurbulentKineticEnergy,
    const BoundedVector<double, TNumNodes>& rNodalTurbulentEnergyDissipationRate);

template <unsigned int TNumNodes>
void CalculateGaussGammaTkeSensitivity(BoundedVector<double, TNumNodes>& rOutput,
                                       const double Cmu,
                                       const double Gamma,
                                       const double TurbulentKinematicViscosity,
                                       const BoundedVector<double, TNumNodes>& rGaussNutTKESensitivities,
                                       const Vector& rShapeFunctions);

template <unsigned int TNumNodes>
void CalculateGaussGammaEpsilonSensitivity(BoundedVector<double, TNumNodes>& rOutput,
                                           const double Gamma,
                                           const double TurbulentKinematicViscosity,
                                           const BoundedVector<double, TNumNodes>& rGaussNutEpsilonSensitivities);

template <unsigned int TDim, unsigned int TNumNodes>
void CalculateProductionVelocitySensitivities(BoundedMatrix<double, TNumNodes, TDim>& rOutput,
                                              const double TurbulentKinematicViscosity,
                                              const double ProductionTerm,
                                              const BoundedMatrix<double, TDim, TDim>& rVelocityGradient,
                                              const Matrix& rShapeDerivatives);

template <unsigned int TDim, unsigned int TNumNodes>
double CalculateProductionShapeSensitivities(
    const double TurbulentKinematicViscosity,
    const double TurbulentKinematicViscosityDerivative,
    const double ProductionTerm,
    const BoundedMatrix<double, TNumNodes, TDim>& rNodalVelocity,
    const Matrix& rShapeDerivatives,
    const GeometricalSensitivityUtility::ShapeFunctionsGradientType& rDN_DxDerivatives);

template <unsigned int TNumNodes>
void CalculateProductionScalarSensitivities(BoundedVector<double, TNumNodes>& rOutput,
                                            const double TurbulentKinematicViscosity,
                                            const double ProductionTerm,
                                            const BoundedVector<double, TNumNodes>& rGaussNutScalarDerivatives);

} // namespace KEpsilonElementData

///@}

} // namespace Kratos

#endif