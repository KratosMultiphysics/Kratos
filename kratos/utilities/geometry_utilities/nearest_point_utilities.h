//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "utilities/geometrical_projection_utilities.h"

namespace Kratos {

///@addtogroup KratosCore
///@{

///@name Kratos Globals
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class NearestPointUtilities
 * @ingroup KratosCore
 * @brief Tools to calculate the nearest point in different geometries
 * @details These tools are generic enough to be used in different contexts while used in the geometries
 * @author Pooyan Dadvand
 */
class KRATOS_API(KRATOS_CORE) NearestPointUtilities
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NearestPointUtilities
    KRATOS_CLASS_POINTER_DEFINITION(NearestPointUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NearestPointUtilities() = delete;

    /// Copy constructor.
    NearestPointUtilities(NearestPointUtilities const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    NearestPointUtilities& operator=(NearestPointUtilities const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Finds the nearest point to the given point on a line segment.
     * @details It first projects the point into the line. If the projected point is inside the segment boundary
     * it returns the projected point. If not it returns the nearest end point of the line.
     * @tparam Type of the Point
     * @tparam TGeometryType The type of the line. Assumes to have [] access and IsInside method
     * @param rPoint The query point which we want to get nearest point to it on the line
     * @param rLine The line in which we want to find the nearest point to rPoint
     * @return The nearest point to rPoint
     */
    template<class TPointType, class TGeometryType>
    static Point LineNearestPoint(
        const TPointType& rPoint,
        const TGeometryType& rLine
        )
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rLine.size() == 2) << "This function only accepts Line2D2 as input" << std::endl;
        return LineNearestPoint(rPoint, rLine[0], rLine[1]);
    }

    /**
     * @brief Finds the nearest point to the given point on a line segment (given 2 points).
     * @details It first projects the point into the line. If the projected point is inside the segment boundary
     * it returns the projected point. If not it returns the nearest end point of the line.
     * @tparam Type of the Point
     * @tparam TGeometryType The type of the line. Assumes to have [] access and IsInside method
     * @param rPoint The query point which we want to get nearest point to it on the line
     * @param rLinePointA The first point of the line
     * @param rLinePointB The second point of the line
     * @return The nearest point to rPoint
     */
    static Point LineNearestPoint(
        const array_1d<double, 3>& rPoint,
        const array_1d<double, 3>& rLinePointA,
        const array_1d<double, 3>& rLinePointB
        );

    /**
     * @brief Computes the nearest point on a triangle to a given point in 3D space.
     * @details This method calculates the nearest point on a triangle defined by three points
     * (`rTrianglePoint0`, `rTrianglePoint1`, and `rTrianglePoint2`) to a given point
     * `rPoint` in 3D space. If the nearest point lies within the triangle, it returns
     * the projection. If the nearest point lies on one of the triangle's edges, it
     * returns the nearest point on that edge.
     * Based on https://www.shadertoy.com/view/ttfGWl
     * @param rPoint The point from which the nearest point on the triangle is to be found.
     * @param rTrianglePoint0 The first vertex of the triangle.
     * @param rTrianglePoint1 The second vertex of the triangle.
     * @param rTrianglePoint2 The third vertex of the triangle.
     * @return Point The nearest point on the triangle to the given point.
     */
    static Point TriangleNearestPoint(
        const array_1d<double, 3>& rPoint,
        const array_1d<double, 3>& rTrianglePoint0,
        const array_1d<double, 3>& rTrianglePoint1,
        const array_1d<double, 3>& rTrianglePoint2
        );

    /**
     * @brief Finds the nearest point to the given point on a triangle.
     * @details It first projects the point into the triangle surface. If the projected point is inside the triangle
     * it returns the projected point. If not it returns the nearest point on the edges of the triangle.
     * Dividing the plane of the triangle in 7 zones and find the nearest reflecting those zones
     *
     *             \       /
     *              \  6  /
     *               \   /
     *                \ /
     *                 /\
     *                /  \
     *          2    /    \     3
     *              /   1  \
     *             /        \
     *    _______ /__________\_____________
     *           /            \
     *     5    /       4      \    7
     *         /                \
     *        /                  \
     * @tparam Type of the Point
     * @tparam TGeometryType The type of the triangle. Assumes to have [] access and IsInside method
     * @param rPoint The query point which we want to get nearest point to it on the line
     * @param rTriangle The triangle in which we want to find the nearest point to rPoint
     * @return The nearest point to rPoint
     */
    template<class TPointType, class TGeometryType>
    static Point TriangleNearestPoint(
        const TPointType& rPoint,
        const TGeometryType& rTriangle
        )
    {
        constexpr double Tolerance = 1e-12;
        Point result;
        array_1d<double,3> local_coordinates = ZeroVector(3);
        const Point center = rTriangle.Center();
        const array_1d<double, 3> normal = rTriangle.UnitNormal(local_coordinates);
        double distance = 0.0;
        const Point point_projected = GeometricalProjectionUtilities::FastProject( center, rPoint, normal, distance);
        rTriangle.PointLocalCoordinates(local_coordinates, point_projected);

        if(local_coordinates[0] < -Tolerance) { // case 2,5,6
            if(local_coordinates[1] < -Tolerance) { // case 5
                result = rTriangle[0];
            } else if ((local_coordinates[0] + local_coordinates[1]) > (1.0+Tolerance)) { // case 6
                result = rTriangle[2];
            } else {
                result = LineNearestPoint(rPoint, rTriangle.GetPoint(0), rTriangle.GetPoint(2));
            }
        } else if(local_coordinates[1] < -Tolerance) { // case 4,7 (case 5 is already covered in previous if)
            if ((local_coordinates[0] + local_coordinates[1]) > (1.0+Tolerance)) { // case 7
                result = rTriangle[1];
            } else { // case 4
                result = LineNearestPoint(rPoint, rTriangle.GetPoint(0), rTriangle.GetPoint(1));
            }
        } else if ((local_coordinates[0] + local_coordinates[1]) > (1.0+Tolerance)) { // case 3
            result = LineNearestPoint(rPoint, rTriangle.GetPoint(1), rTriangle.GetPoint(2));
        } else {  // inside
            result = point_projected;
        }

        return result;
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
}; // Class NearestPointUtilities

///@}

///@name Type Definitions
///@{

///@}

///@} addtogroup block

} // namespace Kratos