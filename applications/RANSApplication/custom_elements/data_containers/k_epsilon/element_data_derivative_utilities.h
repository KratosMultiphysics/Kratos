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

#if !defined(KRATOS_RANS_K_EPSILON_ADJOINT_ELEMENT_DATA_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_K_EPSILON_ADJOINT_ELEMENT_DATA_UTILITIES_H_INCLUDED

// System includes

// Project includes
#include "containers/array_1d.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
///@name  Functions
///@{

namespace KEpsilonElementData
{

template<unsigned int TDim, unsigned int TNumNodes>
class AdjointUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using IndexType = std::size_t;

    using ArrayD = array_1d<double, 3>;

    using MatrixDD = BoundedMatrix<double, TDim, TDim>;

    using MatrixND = BoundedMatrix<double, TNumNodes, TDim>;

    ///@}
    ///@name Static Operations
    ///@{

    static double CalculateProductionVelocityDerivative(
        const IndexType NodeIndex,
        const IndexType DirectionIndex,
        const double ProductionTerm,
        const double TurbulentKinematicViscosity,
        const double TurbulentKinematicViscosityDerivative,
        const MatrixDD& rVelocityGradient,
        const Matrix& rdNdX);

    static double CalculateProductionScalarDerivative(
        const double TurbulentKinematicViscosity,
        const double ProductionTerm,
        const double GaussNutScalarDerivative);

    static double CalculateProductionShapeDerivative(
        const double TurbulentKinematicViscosity,
        const double TurbulentKinematicViscosityDerivative,
        const double ProductionTerm,
        const MatrixND& rNodalVelocity,
        const Matrix& rdNdX,
        const Matrix& rdNdXDerivative);

    static double CalculateGammaKDerivative(
        const IndexType NodeIndex,
        const double Cmu,
        const double Gamma,
        const double TurbulentKinematicViscosity,
        const double GaussNutTKEDerivative,
        const Vector& rShapeFunctions);

    static double CalculateGammaEpsilonDerivative(
        const double Gamma,
        const double TurbulentKinematicViscosity,
        const double GaussNutEpsilonDerivative);

    ///@}
};

} // namespace KEpsilonElementData

///@}

} // namespace Kratos
#endif // KRATOS_RANS_K_EPSILON_ADJOINT_ELEMENT_DATA_UTILITIES_H_INCLUDED