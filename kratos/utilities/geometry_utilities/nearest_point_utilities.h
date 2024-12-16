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
///@name  Enum's
///@{

/**
 * @brief Enum class that defines the possible locations of the nearest point to a triangle
 * @details This enum class defines the possible locations of the nearest point to a triangle
 * INSIDE_TRIANGLE: The nearest point is inside the triangle
 * ON_TRIANGLE_EDGE_01: The nearest point is on the edge 01 of the triangle
 * ON_TRIANGLE_EDGE_12: The nearest point is on the edge 12 of the triangle
 * ON_TRIANGLE_EDGE_20: The nearest point is on the edge 20 of the triangle
 * ON_TRIANGLE_VERTEX_0: The nearest point is on the vertex 0 of the triangle
 * ON_TRIANGLE_VERTEX_1: The nearest point is on the vertex 1 of the triangle
 * ON_TRIANGLE_VERTEX_2: The nearest point is on the vertex 2 of the triangle
 */
enum class TriangleNearestPointLocation
{
    INSIDE_TRIANGLE,
    ON_TRIANGLE_EDGE_01,
    ON_TRIANGLE_EDGE_12,
    ON_TRIANGLE_EDGE_20,
    ON_TRIANGLE_VERTEX_0,
    ON_TRIANGLE_VERTEX_1,
    ON_TRIANGLE_VERTEX_2
};

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
     * it returns the projected point. If not it returns the nearest end point of the line. Implementation based on shadertoy.com/view/ttfGWl.
     * @tparam Type of the Point
     * @tparam TGeometryType The type of the line. Assumes to have [] access and IsInside method
     * @param rPoint The query point which we want to get nearest point to it on the line
     * @param rLine The line in which we want to find the nearest point to rPoint
     * @return The nearest point to rPoint
     */
    template<class TGeometryType>
    static Point LineNearestPoint(
        const array_1d<double, 3>& rPoint,
        const TGeometryType& rLine,
        const double Tolerance = 1.0e-14
        )
    {
        KRATOS_DEBUG_ERROR_IF_NOT(rLine.size() == 2) << "This function only accepts Line2D2 as input" << std::endl;
        return LineNearestPoint(rPoint, rLine[0], rLine[1], Tolerance);
    }

    /**
     * @brief Finds the nearest point to the given point on a line segment (given 2 points).
     * @details It first projects the point into the line. If the projected point is inside the segment boundary
     * it returns the projected point. If not it returns the nearest end point of the line. Implementation based on shadertoy.com/view/ttfGWl.
     * @tparam Type of the Point
     * @tparam TGeometryType The type of the line. Assumes to have [] access and IsInside method
     * @param rPoint The query point which we want to get nearest point to it on the line
     * @param rLinePointA The first point of the line
     * @param rLinePointB The second point of the line
     * @param Tolerance The tolerance to consider the point inside the line
     * @return The nearest point to rPoint
     */
    static Point LineNearestPoint(
        const array_1d<double, 3>& rPoint,
        const array_1d<double, 3>& rLinePointA,
        const array_1d<double, 3>& rLinePointB,
        const double Tolerance = 1.0e-14
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
     * @param rResult Point The nearest point on the triangle to the given point.
     * @return TriangleNearestPointLocation The location of the nearest point on the triangle.
     */
    static TriangleNearestPointLocation TriangleNearestPoint(
        const array_1d<double, 3>& rPoint,
        const array_1d<double, 3>& rTrianglePoint0,
        const array_1d<double, 3>& rTrianglePoint1,
        const array_1d<double, 3>& rTrianglePoint2,
        Point& rResult
        );

    /**
     * @brief Computes the nearest point on a triangle to a given point in 3D space.
     * @details This method calculates the nearest point on a triangle defined by three points
     * (`rTrianglePoint0`, `rTrianglePoint1`, and `rTrianglePoint2`) to a given point
     * `rPoint` in 3D space. If the nearest point lies within the triangle, it returns
     * the projection. If the nearest point lies on one of the triangle's edges, it
     * returns the nearest point on that edge.
     * Based on https://www.shadertoy.com/view/ttfGWl
     * @tparam Type of the Point
     * @tparam TGeometryType The type of the triangle.
     * @param rPoint The query point which we want to get nearest point to it on the line
     * @param rTriangle The triangle in which we want to find the nearest point to rPoint
     * @return The nearest point to rPoint
     */
    template<class TGeometryType>
    static Point TriangleNearestPoint(
        const array_1d<double, 3>& rPoint,
        const TGeometryType& rTriangle
        )
    {
        return TriangleNearestPoint(rPoint, rTriangle[0], rTriangle[1], rTriangle[2]);
    }

    /**
     * @brief Computes the nearest point on a triangle to a given point in 3D space.
     * @details This method calculates the nearest point on a triangle defined by three points
     * (`rTrianglePoint0`, `rTrianglePoint1`, and `rTrianglePoint2`) to a given point
     * `rPoint` in 3D space. If the nearest point lies within the triangle, it returns
     * the projection. If the nearest point lies on one of the triangle's edges, it
     * returns the nearest point on that edge.
     * Based on https://www.shadertoy.com/view/ttfGWl
     * @tparam Type of the Point
     * @tparam TGeometryType The type of the triangle.
     * @param rPoint The query point which we want to get nearest point to it on the line
     * @param rTriangle The triangle in which we want to find the nearest point to rPoint
     * @param rResult Point The nearest point on the triangle to the given point.
     * @return TriangleNearestPointLocation The location of the nearest point on the triangle.
     */
    template<class TGeometryType>
    static TriangleNearestPointLocation TriangleNearestPoint(
        const array_1d<double, 3>& rPoint,
        const TGeometryType& rTriangle,
        Point& rResult
        )
    {
        return TriangleNearestPoint(rPoint, rTriangle[0], rTriangle[1], rTriangle[2], rResult);
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