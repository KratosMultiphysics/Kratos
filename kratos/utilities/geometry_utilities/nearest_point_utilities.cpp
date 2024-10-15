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
    // Compute side vectors and normal
    const array_1d<double, 3> v10 = rTrianglePoint1 - rTrianglePoint0;
    const array_1d<double, 3> p0 = rPoint - rTrianglePoint0;
    const array_1d<double, 3> v21 = rTrianglePoint2 - rTrianglePoint1;
    const array_1d<double, 3> p1 = rPoint - rTrianglePoint1;
    const array_1d<double, 3> v02 = rTrianglePoint0 - rTrianglePoint2;
    const array_1d<double, 3> p2 = rPoint - rTrianglePoint2;
    array_1d<double, 3> normal;
    MathUtils<double>::CrossProduct(normal, v10, v02);

    // Compute the cross product
    array_1d<double, 3> auxiliary_cross_product;

    // Check first size
    ClampScenario clamp_scenario = ClampScenario::WITHIN_BOUNDS;
    MathUtils<double>::CrossProduct(auxiliary_cross_product, v10,normal);
    if( MathUtils<double>::Dot3(auxiliary_cross_product, p0) < 0.0 ) {
        const double clamp_value = MathUtils<double>::Clamp( MathUtils<double>::Dot3(p0,v10)/MathUtils<double>::Dot3(v10, v10), 0.0, 1.0, clamp_scenario);
        noalias(rResult.Coordinates()) =  rTrianglePoint0 + v10 * clamp_value;
        switch (clamp_scenario) {
            case ClampScenario::WITHIN_BOUNDS:
                return TriangleNearestPointLocation::ON_TRIANGLE_EDGE_01;
            case ClampScenario::BELOW_MINIMUM:
                return TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_0;
            case ClampScenario::ABOVE_MAXIMUM:
                return TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_1;
            default:
                KRATOS_ERROR << "Invalid clamp scenario" << std::endl;
        }
    }

    // Check second side
    MathUtils<double>::CrossProduct(auxiliary_cross_product, v21, normal);
    if( MathUtils<double>::Dot3(auxiliary_cross_product, p1) < 0.0 ) {
        const double clamp_value = MathUtils<double>::Clamp( MathUtils<double>::Dot3(p1,v21)/MathUtils<double>::Dot3(v21, v21), 0.0, 1.0, clamp_scenario);
        noalias(rResult.Coordinates()) =  rTrianglePoint1 + v21 * clamp_value;
        switch (clamp_scenario) {
            case ClampScenario::WITHIN_BOUNDS:
                return TriangleNearestPointLocation::ON_TRIANGLE_EDGE_12;
            case ClampScenario::BELOW_MINIMUM:
                return TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_1;
            case ClampScenario::ABOVE_MAXIMUM:
                return TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_2;
            default:
                KRATOS_ERROR << "Invalid clamp scenario" << std::endl;
        }
    }

    // Check third side
    MathUtils<double>::CrossProduct(auxiliary_cross_product, v02, normal);
    if( MathUtils<double>::Dot3(auxiliary_cross_product, p2) < 0.0 ) {
        const double clamp_value = MathUtils<double>::Clamp( MathUtils<double>::Dot3(p2,v02)/MathUtils<double>::Dot3(v02, v02), 0.0, 1.0, clamp_scenario);
        noalias(rResult.Coordinates()) = rTrianglePoint2 + v02 * clamp_value;
        switch (clamp_scenario) {
            case ClampScenario::WITHIN_BOUNDS:
                return TriangleNearestPointLocation::ON_TRIANGLE_EDGE_20;
            case ClampScenario::BELOW_MINIMUM:
                return TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_2;
            case ClampScenario::ABOVE_MAXIMUM:
                return TriangleNearestPointLocation::ON_TRIANGLE_VERTEX_0;
            default:
                KRATOS_ERROR << "Invalid clamp scenario" << std::endl;
        }
    }

    // Compute the projection
    noalias(rResult.Coordinates()) = rPoint - normal * MathUtils<double>::Dot3(normal,p0)/MathUtils<double>::Dot3(normal, normal);
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
    // Compute side vectors and normal
    const array_1d<double, 3> v10 = rTrianglePoint1 - rTrianglePoint0;
    const array_1d<double, 3> p0 = rPoint - rTrianglePoint0;
    const array_1d<double, 3> v21 = rTrianglePoint2 - rTrianglePoint1;
    const array_1d<double, 3> p1 = rPoint - rTrianglePoint1;
    const array_1d<double, 3> v02 = rTrianglePoint0 - rTrianglePoint2;
    const array_1d<double, 3> p2 = rPoint - rTrianglePoint2;
    array_1d<double, 3> normal;
    MathUtils<double>::CrossProduct(normal, v10, v02);

    // Define the result
    Point result;

    // Compute the cross product
    array_1d<double, 3> auxiliary_cross_product;

    // Check first size
    MathUtils<double>::CrossProduct(auxiliary_cross_product, v10,normal);
    if( MathUtils<double>::Dot3(auxiliary_cross_product, p0) < 0.0 ) {
        noalias(result.Coordinates()) =  rTrianglePoint0 + v10 * MathUtils<double>::Clamp( MathUtils<double>::Dot3(p0,v10)/MathUtils<double>::Dot3(v10, v10), 0.0, 1.0 );
        return result;
    }

    // Check second side
    MathUtils<double>::CrossProduct(auxiliary_cross_product, v21, normal);
    if( MathUtils<double>::Dot3(auxiliary_cross_product, p1) < 0.0 ) {
        noalias(result.Coordinates()) =  rTrianglePoint1 + v21 * MathUtils<double>::Clamp( MathUtils<double>::Dot3(p1,v21)/MathUtils<double>::Dot3(v21, v21), 0.0, 1.0 );
        return result;
    }

    // Check third side
    MathUtils<double>::CrossProduct(auxiliary_cross_product, v02, normal);
    if( MathUtils<double>::Dot3(auxiliary_cross_product, p2) < 0.0 ) {
        noalias(result.Coordinates()) = rTrianglePoint2 + v02 * MathUtils<double>::Clamp( MathUtils<double>::Dot3(p2,v02)/MathUtils<double>::Dot3(v02, v02), 0.0, 1.0 );
        return result;
    }

    // Compute the projection
    noalias(result.Coordinates()) = rPoint - normal * MathUtils<double>::Dot3(normal,p0)/MathUtils<double>::Dot3(normal, normal);
    return result;
}

} // namespace Kratos