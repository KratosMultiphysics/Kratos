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

// External includes

// Project includes

#pragma once

namespace Kratos {

class TriangleObject {
public:
    using GeometryType = Geometry<Node>;

    // Epsilon for floating point inaccuracies
    static constexpr double EPSILON = std::numeric_limits<double>::epsilon(); // Adjust multiplier as needed

    static double LengthSquared(const array_1d<double, 3>& rVector);

    static array_1d<double, 3> GetNormal(const GeometryType& rTriangle);

    // Helper function to project a triangle onto an axis and return the min/max bounds
    static void ProjectTriangle(const GeometryType& rTriangle, const array_1d<double, 3>& axis, double& min_val, double& max_val);

    // Tests if a given axis separates the two triangles
    static bool IsSeparatingAxis(const GeometryType& t1, const GeometryType& t2, const array_1d<double, 3>& axis, double tolerance);

    static bool Intersects(const GeometryType& rTriangleA, const GeometryType& rTriangleB, double tolerance = 1e-4);

};

} // namespace Kratos