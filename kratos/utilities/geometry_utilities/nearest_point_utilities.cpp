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

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities/nearest_point_utilities.h"

namespace Kratos
{

Point NearestPointUtilities::LineNearestPoint(
    const array_1d<double, 3>& rPoint,
    const array_1d<double, 3>& rLinePointA,
    const array_1d<double, 3>& rLinePointB,
    const double Tolerance
    )
{
    // Project point globally into the line
    const array_1d<double, 3> ab = rLinePointB - rLinePointA;
    const double factor = (inner_prod(rLinePointB, rPoint) - inner_prod(rLinePointA, rPoint) - inner_prod(rLinePointB, rLinePointA) + inner_prod(rLinePointA, rLinePointA)) / inner_prod(ab, ab);
    Point result(rLinePointA + factor * ab);

    // Compute square_length of the line
    double lx = rLinePointA[0] - rLinePointB[0];
    double ly = rLinePointA[1] - rLinePointB[1];
    double lz = rLinePointA[2] - rLinePointB[2];
    const double square_length = lx * lx + ly * ly + lz * lz;

    // Compute square_length of the projection (A)
    lx = result[0] - rLinePointA[0];
    ly = result[1] - rLinePointA[1];
    lz = result[2] - rLinePointA[2];
    const double square_length_1 = lx * lx + ly * ly + lz * lz;

    // Compute square_length of the projection (B)
    lx = result[0] - rLinePointB[0];
    ly = result[1] - rLinePointB[1];
    lz = result[2] - rLinePointB[2];
    const double square_length_2 = lx * lx + ly * ly + lz * lz;

    // Project point locally into the line
    array_1d<double,3> projected_local;
    if (square_length_1 <= (square_length + Tolerance) && square_length_2 <= (square_length + Tolerance)) {
        return result;
    } else {
        // If the projected point is outside the line, return the nearest end point
        if (square_length_1 < square_length_2) {
            noalias(result.Coordinates()) = rLinePointA;
        } else {
            noalias(result.Coordinates()) = rLinePointB;
        }
        return result;
    }
}

/***********************************************************************************/
/***********************************************************************************/

TriangleNearestPointLocation NearestPointUtilities::TriangleNearestPoint(
    const array_1d<double, 3>& rPoint,
    const array_1d<double, 3>& rTrianglePoint0,
    const array_1d<double, 3>& rTrianglePoint1,
    const array_1d<double, 3>& rTrianglePoint2,
    Point& rResult
    )
{
    // Winding-order-invariant algorithm based on Voronoi region classification.
    // Reference: "Real-Time Collision Detection" by Christer Ericson, Chapter 5.1.5

    const array_1d<double, 3> ab = rTrianglePoint1 - rTrianglePoint0;
    const array_1d<double, 3> ac = rTrianglePoint2 - rTrianglePoint0;
    const array_1d<double, 3> ap = rPoint - rTrianglePoint0;

    // Check if P in vertex region outside A (vertex 0)
    const double d1 = MathUtils<double>::Dot3(ab, ap);
    const double d2 = MathUtils<double>::Dot3(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) {
        noalias(rResult.Coordinates()) = rTrianglePoint0;
        return TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_0;
    }

    // Check if P in vertex region outside B (vertex 1)
    const array_1d<double, 3> bp = rPoint - rTrianglePoint1;
    const double d3 = MathUtils<double>::Dot3(ab, bp);
    const double d4 = MathUtils<double>::Dot3(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) {
        noalias(rResult.Coordinates()) = rTrianglePoint1;
        return TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_1;
    }

    // Check if P in edge region of AB (edge 01), if so return projection of P onto AB
    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const double v = d1 / (d1 - d3);
        noalias(rResult.Coordinates()) = rTrianglePoint0 + v * ab;
        return TriangleNearestPointLocation::ON_TRIANGLE_EDGE_01;
    }

    // Check if P in vertex region outside C (vertex 2)
    const array_1d<double, 3> cp = rPoint - rTrianglePoint2;
    const double d5 = MathUtils<double>::Dot3(ab, cp);
    const double d6 = MathUtils<double>::Dot3(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) {
        noalias(rResult.Coordinates()) = rTrianglePoint2;
        return TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_2;
    }

    // Check if P in edge region of AC (edge 02), if so return projection of P onto AC
    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const double w = d2 / (d2 - d6);
        noalias(rResult.Coordinates()) = rTrianglePoint0 + w * ac;
        return TriangleNearestPointLocation::ON_TRIANGLE_EDGE_20;
    }

    // Check if P in edge region of BC (edge 12), if so return projection of P onto BC
    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        noalias(rResult.Coordinates()) = rTrianglePoint1 + w * (rTrianglePoint2 - rTrianglePoint1);
        return TriangleNearestPointLocation::ON_TRIANGLE_EDGE_12;
    }

    // P inside face region. Compute Q through its barycentric coordinates (u, v, w)
    const double denom = 1.0 / (va + vb + vc);
    const double v = vb * denom;
    const double w = vc * denom;
    noalias(rResult.Coordinates()) = rTrianglePoint0 + ab * v + ac * w;
    return TriangleNearestPointLocation::INSIDE_TRIANGLE;
}

/***********************************************************************************/
/***********************************************************************************/

Point NearestPointUtilities::TriangleNearestPoint(
    const array_1d<double, 3>& rPoint,
    const array_1d<double, 3>& rTrianglePoint0,
    const array_1d<double, 3>& rTrianglePoint1,
    const array_1d<double, 3>& rTrianglePoint2
    )
{
    // Delegate to the overload that also returns location info
    Point result;
    TriangleNearestPoint(rPoint, rTrianglePoint0, rTrianglePoint1, rTrianglePoint2, result);
    return result;
}

} // namespace Kratos