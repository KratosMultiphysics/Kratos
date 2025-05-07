//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/intersection_utilities.h"

namespace Kratos
{

int IntersectionUtilities::Compute2PlaneIntersection(
    const array_1d<double, 3>& rPlanePoint1,
    const array_1d<double, 3>& rNormal1,
    const array_1d<double, 3>& rPlanePoint2,
    const array_1d<double, 3>& rNormal2,
    array_1d<double, 3>& rPoint,
    array_1d<double, 3>& rTangent,
    const double Tolerance
    )
{
    // Compute the direction of the intersection line as the cross product of the normals.
    rTangent = MathUtils<double>::CrossProduct(rNormal1, rNormal2);

    // Compute the squared norm of the tangent vector.
    const double norm_tangent_sq = norm_2(rTangent);

    // Define a tolerance for determining when the tangent is effectively zero.
    if (norm_tangent_sq < Tolerance) {
        // The planes are parallel or coincident; a unique intersection line does not exist.
        return 0;
    }

    // Compute the plane constants for the plane equations:
    // For plane 1: n1 · x = d1, where d1 = rNormal1 · rPlanePoint1.
    // For plane 2: n2 · x = d2, where d2 = rNormal2 · rPlanePoint2.
    double d1 = MathUtils<double>::Dot(rNormal1, rPlanePoint1);
    double d2 = MathUtils<double>::Dot(rNormal2, rPlanePoint2);

    // Compute the intermediate vector v = d1 * rNormal2 - d2 * rNormal1.
    array_1d<double, 3> v;
    for (int i = 0; i < 3; i++) {
        v[i] = d1 * rNormal2[i] - d2 * rNormal1[i];
    }

    // Compute a point on the intersection line using:
    // rPoint = ( v x rTangent ) / ||rTangent||^2.
    const array_1d<double, 3> tmp = MathUtils<double>::CrossProduct(v, rTangent);
    for (int i = 0; i < 3; i++) {
        rPoint[i] = tmp[i] / norm_tangent_sq;
    }

    return 1;
}

/***********************************************************************************/
/***********************************************************************************/

int IntersectionUtilities::Compute3PlaneIntersection(
    const array_1d<double, 3>& rPlanePoint1,
    const array_1d<double, 3>& rNormal1,
    const array_1d<double, 3>& rPlanePoint2,
    const array_1d<double, 3>& rNormal2,
    const array_1d<double, 3>& rPlanePoint3,
    const array_1d<double, 3>& rNormal3,
    array_1d<double, 3>& rPoint1,
    array_1d<double, 3>& rVector1,
    array_1d<double, 3>& rPoint2,
    array_1d<double, 3>& rVector2,
    const double Tolerance
    )
{
    // First, compute the intersection line between plane 1 and plane 2.
    const int ret = Compute2PlaneIntersection(rPlanePoint1, rNormal1,
                                              rPlanePoint2, rNormal2,
                                              rPoint2, rVector2,
                                              Tolerance);

    // Plane 1 and 2 are parallel (or nearly so).
    if (ret == 0) {
        // Check whether rNormal1 and rNormal3 are also parallel.
        const double dot_n1_n3 = std::abs(MathUtils<double>::Dot(rNormal1, rNormal3));
        const double norm1 = std::sqrt(norm_2(rNormal1));
        const double norm3 = std::sqrt(norm_2(rNormal3));
        if (dot_n1_n3 / (norm1 * norm3) > 1.0 - Tolerance) {
            // All three planes are (nearly) parallel.
            return 0;
        } else {
            // Ambiguous configuration: two planes are parallel so no unique intersection point exists.
            Compute2PlaneIntersection(rPlanePoint1, rNormal1,
                rPlanePoint3, rNormal3,
                rPoint1, rVector1,
                Tolerance);
            Compute2PlaneIntersection(rPlanePoint2, rNormal2,
                rPlanePoint3, rNormal3,
                rPoint2, rVector2,
                Tolerance);
            return 2;
        }
    }

    // Now, intersect the line (rPoint2 + t * rVector2) with plane 3.
    const double d3 = MathUtils<double>::Dot(rNormal3, rPlanePoint3);
    const double denom = MathUtils<double>::Dot(rNormal3, rVector2);

    // If the line is parallel to plane 3.
    if (std::abs(denom) < Tolerance) {
        const double lineValue = MathUtils<double>::Dot(rNormal3, rPoint2);
        if (std::abs(lineValue - d3) < Tolerance) {
            // The entire line lies in plane 3.
            Compute2PlaneIntersection(rPlanePoint3, rNormal3,
                rPlanePoint1, rNormal1,
                rPoint1, rVector1,
                Tolerance);
            return 1;
        } else {
            // The line does not intersect plane 3 in a unique point.
            Compute2PlaneIntersection(rPlanePoint3, rNormal3,
                rPlanePoint1, rNormal1,
                rPoint1, rVector1,
                Tolerance);
            Compute2PlaneIntersection(rPlanePoint2, rNormal2,
                rPlanePoint1, rNormal1,
                rPoint2, rVector2,
                Tolerance);
            return 2;
        }
    }

    // Compute the parameter t such that: rPoint2 + t * rVector2 lies in plane 3.
    const double t = (d3 - MathUtils<double>::Dot(rNormal3, rPoint2)) / denom;
    for (int i = 0; i < 3; i++) {
        rPoint1[i] = rPoint2[i] + t * rVector2[i];
    }

    return 3;
}

}