//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

// System includes
#include <algorithm>

// External includes

// Project includes
#include "geometries/triangle_object.h"

namespace Kratos {

double TriangleObject::LengthSquared(
    const array_1d<double, 3>& rVector
) 
{
    return rVector[0] * rVector[0] + rVector[1] * rVector[1] + rVector[2] * rVector[2];
}

array_1d<double, 3> TriangleObject::GetNormal(
    const TriangleObject::GeometryType& rTriangle
) 
{
    array_1d<double, 3> normal;
    
    MathUtils<double>::CrossProduct(normal, rTriangle[1] - rTriangle[0], rTriangle[2] - rTriangle[0]);
    
    return normal;
}

void TriangleObject::ProjectTriangle(
    const TriangleObject::GeometryType& rTriangle, 
    const array_1d<double, 3>& axis, 
    double& min_val, 
    double& max_val
) 
{
    double p0 = MathUtils<double>::Dot3(rTriangle[0], axis);
    double p1 = MathUtils<double>::Dot3(rTriangle[1], axis);
    double p2 = MathUtils<double>::Dot3(rTriangle[2], axis);
    
    min_val = std::min({p0, p1, p2});
    max_val = std::max({p0, p1, p2});
}

// Tests if a given axis separates the two triangles
bool TriangleObject::IsSeparatingAxis(
    const TriangleObject::GeometryType& t1, 
    const TriangleObject::GeometryType& t2, 
    const array_1d<double, 3>& axis, 
    double tolerance
) 
{  
    // If the cross product yielded a zero vector, it's an invalid axis
    if (TriangleObject::LengthSquared(axis) < TriangleObject::EPSILON * TriangleObject::EPSILON) return false;
    
    double min1, max1, min2, max2;

    TriangleObject::ProjectTriangle(t1, axis, min1, max1);
    TriangleObject::ProjectTriangle(t2, axis, min2, max2);
    
    // If the intervals do not overlap, this is a separating axis
    if (max1 < min2 - tolerance || max2 < min1 - tolerance) {
        return true; 
    }
    
    return false;
}

bool TriangleObject::Intersects(
    const GeometryType& rTriangleA, 
    const GeometryType& rTriangleB, 
    double tolerance
) 
{
    auto n1 = TriangleObject::GetNormal(rTriangleA);
    auto n2 = TriangleObject::GetNormal(rTriangleB);

    // Test Face Normals (Checks if one triangle is entirely above/below the other)
    if (TriangleObject::IsSeparatingAxis(rTriangleA, rTriangleB, n1, tolerance)) return false;
    if (TriangleObject::IsSeparatingAxis(rTriangleA, rTriangleB, n2, tolerance)) return false;

    // Edges of both triangles
    std::vector<array_1d<double, 3>> e1 = { rTriangleA[1] - rTriangleA[0], rTriangleA[2] - rTriangleA[1], rTriangleA[0] - rTriangleA[2] };
    std::vector<array_1d<double, 3>> e2 = { rTriangleB[1] - rTriangleB[0], rTriangleB[2] - rTriangleB[1], rTriangleB[0] - rTriangleB[2] };

    // Test Edge-to-Edge Cross Products (9 axes)
    // This catches edge-to-edge and point-to-face penetrations
    array_1d<double, 3> axis;

    for (std::size_t i = 0; i < 3; i++) {
        for (std::size_t j = 0; j < 3; j++) {
            MathUtils<double>::CrossProduct(axis, e1[i], e2[j]);
            if (TriangleObject::IsSeparatingAxis(rTriangleA, rTriangleB, axis, tolerance)) return false;
        }
    }

    // Test Coplanar Triangles (Overlapping Faces)
    // If normals are parallel, the cross products of edges above will yield zero-vectors.
    // We must test 2D separating axes (the perpendiculars of the edges within the plane).
    array_1d<double, 3> n1_cross_n2;
    MathUtils<double>::CrossProduct(n1_cross_n2, n1, n2);
    
    // Triangles are parallel/coplanar. Test 6 edge normals lying in the shared plane.
    if (TriangleObject::LengthSquared(n1_cross_n2) < TriangleObject::EPSILON * TriangleObject::EPSILON) {
        for (std::size_t i = 0; i < 3; i++) {
            MathUtils<double>::CrossProduct(axis, e1[i], n1);
            if (TriangleObject::IsSeparatingAxis(rTriangleA, rTriangleB, axis, tolerance)) return false;
            MathUtils<double>::CrossProduct(axis, e2[i], n1);
            if (TriangleObject::IsSeparatingAxis(rTriangleA, rTriangleB, axis, tolerance)) return false;
        }
    }

    // If no separating axis is found across all checks, the triangles must intersect.
    return true;
}

} // namespace Kratos