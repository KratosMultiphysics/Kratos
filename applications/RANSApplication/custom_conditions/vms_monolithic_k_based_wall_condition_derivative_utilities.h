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

#if !defined(KRATOS_VMS_MONOLITHIC_K_BASED_WALL_CONDITION_DERIVATIVE_UTILITIES_H)
#define KRATOS_VMS_MONOLITHIC_K_BASED_WALL_CONDITION_DERIVATIVE_UTILITIES_H

// System includes

// External includes

// Project includes
#include "containers/array_1d.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"
#include "includes/node.h"
#include "includes/ublas_interface.h"

// Application includes

namespace Kratos
{
///@name Kratos Classes
///@{


template <unsigned int TDim>
class VMSMonolithicKBasedWallConditionDerivativeUtilities
{
public:
    ///@name Public Type Definitions
    ///@{

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using IndexType = std::size_t;

    using ArrayD = array_1d<double, 3>;

    ///@}
    ///@name Classes
    ///@{

    class NonRelatedDerivative
    {
    public:
        ///@name TypeDefinitions
        ///@{

        static constexpr IndexType TDerivativeDimension = 1;

        static constexpr IndexType VelocityDerivativeFactor = 0.0;

        static constexpr IndexType TurbulentKineticEnergyDerivativeFactor = 0.0;

        ///@}
        ///@name Static Operations
        ///@{

        static ArrayD CalculateWallVelocityDerivative(
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const Vector& rN) { return ZeroVector(3); }

        static double CalculateWallHeightConditionDerivative(
            const GeometryType& rConditionGeometry,
            const GeometryType& rParentElementGeometry,
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const double NormalMagnitude,
            const ArrayD& rUnitNormal) { return 0.0; }

        static double CalculateWallHeightParentElementDerivative(
            const GeometryType& rConditionGeometry,
            const GeometryType& rParentElementGeometry,
            const IndexType DirectionIndex,
            const ArrayD& rUnitNormal) { return 0.0; }

        ///@}
    };

    class VelocityDerivative : public NonRelatedDerivative
    {
    public:
        ///@name TypeDefinitions
        ///@{

        static constexpr IndexType TDerivativeDimension = TDim;

        static constexpr IndexType VelocityDerivativeFactor = 1.0;

        static constexpr IndexType TurbulentKineticEnergyDerivativeFactor = 0.0;

        ///@}
        ///@name Static Operations
        ///@{

        static ArrayD CalculateWallVelocityDerivative(
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const Vector& rN);

        using NonRelatedDerivative::CalculateWallHeightConditionDerivative;

        using NonRelatedDerivative::CalculateWallHeightParentElementDerivative;

        ///@}
    };

    class KDerivative : public NonRelatedDerivative
    {
    public:
        ///@name TypeDefinitions
        ///@{

        static constexpr IndexType TDerivativeDimension = 1;

        static constexpr IndexType VelocityDerivativeFactor = 0.0;

        static constexpr IndexType TurbulentKineticEnergyDerivativeFactor = 1.0;

        ///@}
        ///@name Static Operations
        ///@{

        using NonRelatedDerivative::CalculateWallVelocityDerivative;

        using NonRelatedDerivative::CalculateWallHeightConditionDerivative;

        using NonRelatedDerivative::CalculateWallHeightParentElementDerivative;

        ///@}
    };



    class ShapeDerivative : public NonRelatedDerivative
    {
    public:
        ///@name TypeDefinitions
        ///@{

        static constexpr IndexType TDerivativeDimension = TDim;

        static constexpr IndexType VelocityDerivativeFactor = 0.0;

        static constexpr IndexType TurbulentKineticEnergyDerivativeFactor = 0.0;

        ///@}
        ///@name Static Operations
        ///@{

        using NonRelatedDerivative::CalculateWallVelocityDerivative;

        static double CalculateWallHeightConditionDerivative(
            const GeometryType& rConditionGeometry,
            const GeometryType& rParentElementGeometry,
            const IndexType NodeIndex,
            const IndexType DirectionIndex,
            const double NormalMagnitude,
            const ArrayD& rUnitNormal);

        static double CalculateWallHeightParentElementDerivative(
            const GeometryType& rConditionGeometry,
            const GeometryType& rParentElementGeometry,
            const IndexType DirectionIndex,
            const ArrayD& rUnitNormal);

        ///@}
    };

    ///@}
};

} // namespace Kratos

#endif // KRATOS_VMS_MONOLITHIC_K_BASED_WALL_CONDITION_DERIVATIVE_UTILITIES_H