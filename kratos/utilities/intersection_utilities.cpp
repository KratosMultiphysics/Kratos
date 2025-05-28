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
#include "utilities/math_utils.h"
#include "utilities/intersection_utilities.h"

namespace Kratos
{

int IntersectionUtilities::ComputeTwoPlaneIntersection(
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

int IntersectionUtilities::ComputeThreePlaneIntersection(
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
    const int ret = ComputeTwoPlaneIntersection(rPlanePoint1, rNormal1,
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
            ComputeTwoPlaneIntersection(rPlanePoint1, rNormal1,
                rPlanePoint3, rNormal3,
                rPoint1, rVector1,
                Tolerance);
            ComputeTwoPlaneIntersection(rPlanePoint2, rNormal2,
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
            ComputeTwoPlaneIntersection(rPlanePoint3, rNormal3,
                rPlanePoint1, rNormal1,
                rPoint1, rVector1,
                Tolerance);
            return 1;
        } else {
            // The line does not intersect plane 3 in a unique point.
            ComputeTwoPlaneIntersection(rPlanePoint3, rNormal3,
                rPlanePoint1, rNormal1,
                rPoint1, rVector1,
                Tolerance);
            ComputeTwoPlaneIntersection(rPlanePoint2, rNormal2,
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

/***********************************************************************************/
/***********************************************************************************/

int IntersectionUtilities::ComputeTriangleLineIntersection(
    const array_1d<double,3>& rTrianglePoint1,
    const array_1d<double,3>& rTrianglePoint2,
    const array_1d<double,3>& rTrianglePoint3,
    const array_1d<double,3>& rLinePoint1,
    const array_1d<double,3>& rLinePoint2,
    array_1d<double,3>& rIntersectionPoint,
    const double Epsilon) 
{
    // This is the adaption of the implementation provided in:
    // http://www.softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()
    // Based on Tomas Möller & Ben Trumbore (1997) Fast, Minimum Storage Ray-Triangle Intersection, Journal of Graphics Tools, 2:1, 21-28, DOI: 10.1080/10867651.1997.10487468 

    // Get triangle edge vectors and plane normal
    const array_1d<double,3> u = rTrianglePoint2 - rTrianglePoint1;
    const array_1d<double,3> v = rTrianglePoint3 - rTrianglePoint1;
    array_1d<double,3> n;
    MathUtils<double>::CrossProduct<array_1d<double,3>,array_1d<double,3>,array_1d<double,3>>(n,u,v);

    // Check if the triangle is degenerate (do not deal with this case)
    if (MathUtils<double>::Norm3(n) < Epsilon) {
        return -1;
    }

    const array_1d<double,3> dir = rLinePoint2 - rLinePoint1; // Edge direction vector
    const array_1d<double,3> w_0 = rLinePoint1 - rTrianglePoint1;
    const double a = -inner_prod(n,w_0);
    const double b = inner_prod(n,dir);

    // Check if the ray is parallel to the triangle plane
    if (std::abs(b) < Epsilon) {
        if (a == 0.0) {
            return 2;    // Edge lies in the triangle plane
        } else {
            return 0;    // Edge does not lie in the triangle plane
        }
    }

    // If the edge is not parallel, compute the intersection point
    const double r = a / b;
    if (r < 0.0) {
        return 0;    // Edge goes away from triangle
    } else if (r > 1.0) {
        return 0;    // Edge goes away from triangle
    }

    rIntersectionPoint = rLinePoint1 + r*dir;

    // Check if the intersection point is inside the triangle
    if (PointInTriangle(rTrianglePoint1, rTrianglePoint2, rTrianglePoint3, rIntersectionPoint)) {
        return 1;
    }
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

bool IntersectionUtilities::TriangleLineIntersection2D(
    const array_1d<double,3>& rVert0,
    const array_1d<double,3>& rVert1,
    const array_1d<double,3>& rVert2,
    const array_1d<double,3>& rPoint0,
    const array_1d<double,3>& rPoint1)
{
    // Check the intersection of each edge against the intersecting object
    array_1d<double,3> int_point;
    if (ComputeLineLineIntersection(rVert0, rVert1, rPoint0, rPoint1, int_point)) return true;
    if (ComputeLineLineIntersection(rVert1, rVert2, rPoint0, rPoint1, int_point)) return true;
    if (ComputeLineLineIntersection(rVert2, rVert0, rPoint0, rPoint1, int_point)) return true;

    // Let check second geometry is inside the first one.
    if (PointInTriangle(rVert0, rVert1, rVert2, rPoint0)) return true;

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool IntersectionUtilities::PointInTriangle(
    const array_1d<double,3>& rVert0,
    const array_1d<double,3>& rVert1,
    const array_1d<double,3>& rVert2,
    const array_1d<double,3>& rPoint,
    const double Tolerance)
{
    const array_1d<double,3> u = rVert1 - rVert0;
    const array_1d<double,3> v = rVert2 - rVert0;
    const array_1d<double,3> w = rPoint - rVert0;
    
    const double uu = inner_prod(u, u);
    const double uv = inner_prod(u, v);
    const double vv = inner_prod(v, v);
    const double wu = inner_prod(w, u);
    const double wv = inner_prod(w, v);
    const double denom = uv * uv - uu * vv;

    const double xi  = (uv * wv - vv * wu) / denom;
    const double eta = (uv * wu - uu * wv) / denom;

    if (xi < -Tolerance) return false;
    if (eta < -Tolerance) return false;
    if (xi + eta > 1.0 + Tolerance) return false;
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

int IntersectionUtilities::ComputeLineLineIntersection(
    const array_1d<double,3>& rLine1Point0,
    const array_1d<double,3>& rLine1Point1,
    const array_1d<double,3>& rLine2Point0,
    const array_1d<double,3>& rLine2Point1,
    array_1d<double,3>& rIntersectionPoint,
    const double Epsilon)
{
    const array_1d<double,3> r = rLine1Point1 - rLine1Point0;
    const array_1d<double,3> s = rLine2Point1 - rLine2Point0;
    const array_1d<double,3> q_p = rLine2Point0 - rLine1Point0;        // q - p

    const double aux_1 = CrossProd2D(r,s);
    const double aux_2 = CrossProd2D(q_p,r);
    const double aux_3 = CrossProd2D(q_p,s);

    // Check intersection cases
    // Check that the lines are not collinear
    if (std::abs(aux_1) < Epsilon && std::abs(aux_2) < Epsilon){
        const double aux_4 = inner_prod(r,r);
        const double aux_5 = inner_prod(s,r);
        const double t_0 = inner_prod(q_p,r)/aux_4;
        const double t_1 = t_0 + aux_5/aux_4;
        if (aux_5 < 0.0){
            if (t_1 >= 0.0 && t_0 <= 1.0){
                return 2;    // The lines are collinear and overlapping
            }
        } else {
            if (t_0 >= 0.0 && t_1 <= 1.0){
                return 2;    // The lines are collinear and overlapping
            }
        }
    // Check that the lines are parallel and non-intersecting
    } else if (std::abs(aux_1) < Epsilon && std::abs(aux_2) > Epsilon){
        return 0;
    // Check that the lines intersect
    } else if (std::abs(aux_1) > Epsilon){
        const double u = aux_2/aux_1;
        const double t = aux_3/aux_1;
        if (((u >= 0.0) && (u <= 1.0)) && ((t >= 0.0) && (t <= 1.0))){
            rIntersectionPoint = rLine2Point0 + u*s;
            // Check if the intersection occurs in one of the end points
            if (u < Epsilon || (1.0 - u) < Epsilon) {
                return 3;
            } else {
                return 1;
            }
        }
    }
    // Otherwise, the lines are non-parallel but do not intersect
    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

int IntersectionUtilities::ComputePlaneLineIntersection(
    const array_1d<double,3>& rPlaneBasePoint,
    const array_1d<double,3>& rPlaneNormal,
    const array_1d<double,3>& rLinePoint1,
    const array_1d<double,3>& rLinePoint2,
    array_1d<double,3>& rIntersectionPoint,
    const double Epsilon)
{
    // This is the adaption of the implementation provided in:
    // http://www.softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm#intersect_RayTriangle()
    // (Intersection of a Segment with a Plane)

    // Get direction vector of edge
    const array_1d<double,3> line_dir = rLinePoint2 - rLinePoint1;

    // Check if the segment is parallel to the plane or even coincides with it
    const double a = inner_prod(rPlaneNormal,( rPlaneBasePoint - rLinePoint1 ));
    const double b = inner_prod(rPlaneNormal,line_dir);
    if (std::abs(b) < Epsilon){
        if (std::abs(a) < Epsilon){
            return 2;    // Segment lies in the plane
        } else {
            return 0;    // Segment does not lie in the plane, but is parallel to it
        }
    }

    // Compute the intersection point and check if it is inside the bounds of the segment
    const double r = a / b;
    if (r < 0.0 - Epsilon){
        return 0;    // Intersection point lies outside the bounds of the segment
    } else if (r > 1.0 + Epsilon) {
        return 0;    // Intersection point lies outside the bounds of the segment
    }
    rIntersectionPoint = rLinePoint1 + r * line_dir;

    return 1;
}

/***********************************************************************************/
/***********************************************************************************/

int IntersectionUtilities::ComputeLineBoxIntersection(
    const array_1d<double,3>& rBoxPoint0,
    const array_1d<double,3>& rBoxPoint1,
    const array_1d<double,3>& rLinePoint0,
    const array_1d<double,3>& rLinePoint1)
{
    array_1d<double,3> intersection_point = ZeroVector(3);

    if (rLinePoint1[0] < rBoxPoint0[0] && rLinePoint0[0] < rBoxPoint0[0]) return false;
    if (rLinePoint1[0] > rBoxPoint1[0] && rLinePoint0[0] > rBoxPoint1[0]) return false;
    if (rLinePoint1[1] < rBoxPoint0[1] && rLinePoint0[1] < rBoxPoint0[1]) return false;
    if (rLinePoint1[1] > rBoxPoint1[1] && rLinePoint0[1] > rBoxPoint1[1]) return false;
    if (rLinePoint1[2] < rBoxPoint0[2] && rLinePoint0[2] < rBoxPoint0[2]) return false;
    if (rLinePoint1[2] > rBoxPoint1[2] && rLinePoint0[2] > rBoxPoint1[2]) return false;
    if (rLinePoint0[0] > rBoxPoint0[0] && rLinePoint0[0] < rBoxPoint1[0] &&
        rLinePoint0[1] > rBoxPoint0[1] && rLinePoint0[1] < rBoxPoint1[1] &&
        rLinePoint0[2] > rBoxPoint0[2] && rLinePoint0[2] < rBoxPoint1[2]) {
        return true;
    }
    if ((GetLineBoxIntersection(rLinePoint0[0]-rBoxPoint0[0], rLinePoint1[0]-rBoxPoint0[0], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 1 )) ||
        (GetLineBoxIntersection(rLinePoint0[1]-rBoxPoint0[1], rLinePoint1[1]-rBoxPoint0[1], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 2 )) ||
        (GetLineBoxIntersection(rLinePoint0[2]-rBoxPoint0[2], rLinePoint1[2]-rBoxPoint0[2], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 3 )) ||
        (GetLineBoxIntersection(rLinePoint0[0]-rBoxPoint1[0], rLinePoint1[0]-rBoxPoint1[0], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 1 )) ||
        (GetLineBoxIntersection(rLinePoint0[1]-rBoxPoint1[1], rLinePoint1[1]-rBoxPoint1[1], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 2 )) ||
        (GetLineBoxIntersection(rLinePoint0[2]-rBoxPoint1[2], rLinePoint1[2]-rBoxPoint1[2], rLinePoint0, rLinePoint1, intersection_point) && InBox(intersection_point, rBoxPoint0, rBoxPoint1, 3 ))){
        return true;
    }

    return false;
}

}  /* namespace Kratos.*/
