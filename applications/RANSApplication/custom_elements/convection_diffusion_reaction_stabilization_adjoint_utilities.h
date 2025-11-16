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

#if !defined(KRATOS_CONVECTION_DIFFUSION_REACTION_STABILIZATION_ADJOINT_UTILITIES_H_INCLUDED)
#define KRATOS_CONVECTION_DIFFUSION_REACTION_STABILIZATION_ADJOINT_UTILITIES_H_INCLUDED

// System includes
#include <cmath>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

// Application includes


namespace Kratos
{
///@name  Functions
///@{

namespace ConvectionDiffusionReactionStabilizationUtilities
{

template<unsigned int TDim, unsigned int TNumNodes>
class AdjointUtilities
{
public:
    ///@name Public Type Definitions
    ///@{

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using ArrayD = array_1d<double, TDim>;

    using VectorN = BoundedVector<double, TNumNodes>;

    using MatrixNN = BoundedMatrix<double, TNumNodes, TNumNodes>;

    using IndexType = std::size_t;

    ///@}
    ///@name Classes
    ///@{

    class Derivatives
    {
    public:
        ///@name Classes
        ///@{

        class NonRelatedVariable
        {
        public:
            ///@name Static Operations
            ///@{

            static double CalculateElementLengthDerivative(
                const double DetJDerivative,
                const double ElementLength)
            {
                return 0.0;
            }

            ///@}
        };

        class Shape
        {
        public:
            ///@name Static Operations
            ///@{

            static double CalculateElementLengthDerivative(
                const double DetJDerivative,
                const double ElementLength);

            ///@}
        };

        ///@}
    };

    ///@}
    ///@name Static Operations
    ///@{

    static double CalculateVectorNormDerivative(
        const double VectorNorm,
        const ArrayD& rVector,
        const ArrayD& rVectorDerivative);

    static double CalculateStabilizationTauDerivative(
        const double StabilizationTau,
        const double EffectiveVelocityMagnitude,
        const double EffectiveKinematicViscosity,
        const double ReactionTerm,
        const double ElementLength,
        const double EffectiveVelocityMagnitudeDerivative,
        const double EffectiveKinematicViscosityDerivative,
        const double ReactionTermDerivative,
        const double ElementLengthDerivative);

    static double CalculateAbsoluteValueDerivative(
        const double Value,
        const double ValueDerivative);

    static void CalculateDiscreteUpwindOperatorDerivative(
        VectorN& rOutput,
        const VectorN& rNodalScalarValues,
        const MatrixNN& rInputMatrix,
        const MatrixNN& rInputMatrixDerivative);

    static double CalculatePositivityPreservingCoefficientDerivative(
        const double PositivityPreservingMatrixCoefficient,
        const MatrixNN& rInputMatrix,
        const MatrixNN& rInputMatrixDerivative);

    ///@}
};

} // namespace ConvectionDiffusionReactionStabilizationUtilities
} // namespace Kratos

#endif // KRATOS_CONVECTION_DIFFUSION_REACTION_STABILIZATION_ADJOINT_UTILITIES_H_INCLUDED