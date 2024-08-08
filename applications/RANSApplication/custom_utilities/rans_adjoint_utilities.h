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

#if !defined(KRATOS_RANS_APPLICATION_ADJOINT_UTILITIES_H_INCLUDED)
#define KRATOS_RANS_APPLICATION_ADJOINT_UTILITIES_H_INCLUDED

// System includes
#include <cmath>
#include <tuple>

// Project includes
#include "includes/node.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/model_part.h"
#include "utilities/geometrical_sensitivity_utility.h"

namespace Kratos
{
///@name Kratos Globals
///@{

class RansAdjointUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using NodeType = Node<3>;

    using GeometryType = Geometry<NodeType>;

    using ElementType = ModelPart::ElementType;

    using ConditionType = ModelPart::ConditionType;

    using IndexType = std::size_t;

    ///@}
    ///@name Static Operations
    ///@{

    static double CalculateVectorNormDerivative(
        const double VectorNorm,
        const array_1d<double, 3>& rVector,
        const array_1d<double, 3>& rVectorDerivative);

    static array_1d<double, 3> CalculateUnitVectorDerivative(
        const double VectorMagnitude,
        const array_1d<double, 3>& rUnitVector,
        const array_1d<double, 3>& rVectorDerivative);

    static double CalculateWallHeightConditionDerivative(
        const GeometryType& rConditionGeometry,
        const GeometryType& rParentElementGeometry,
        const IndexType DirectionIndex,
        const array_1d<double, 3>& rUnitNormal,
        const array_1d<double, 3>& rUnitNormalDerivative);

    static double CalculateWallHeightParentElementDerivative(
        const GeometryType& rConditionGeometry,
        const GeometryType& rParentElementGeometry,
        const IndexType DirectionIndex,
        const array_1d<double, 3>& rUnitNormal,
        const array_1d<double, 3>& rUnitNormalDerivative);

    static void CalculateYPlusAndUtauDerivative(
        double& rYPlusDerivative,
        double& rUTauDerivative,
        const double YPlus,
        const double WallVelocity,
        const double WallVelocityDerivative,
        const double WallHeight,
        const double WallHeightDerivative,
        const double KinematicViscosity,
        const double Kappa,
        const double Beta,
        const double YPlusLimit);

    static void CopyAdjointSolutionToNonHistorical(ModelPart& rModelPart);

    static void RescaleAdjointSolution(ModelPart& rModelPart);

    static void RescaleShapeSensitivity(ModelPart& rModelPart);

    static void CalculateTransientReponseFunctionInterpolationError(
        ModelPart& rModelPart,
        const double Gamma,
        const double DeltaTime);

    ///@}
    ///@name Classes
    ///@{

    template<unsigned int TDim, unsigned int TNumNodes>
    class GeometricalDerivatives
    {
    public:
        ///@name Static Operations
        ///@{

        static double DomainSizeDerivative(
            const GeometryType& rGeometry,
            const IndexType NodeIndex,
            const IndexType DirectionIndex);

        ///@}
    };

    ///@}
};

} // namespace Kratos

#endif // KRATOS_RANS_APPLICATION_ADJOINT_UTILITIES_H_INCLUDED