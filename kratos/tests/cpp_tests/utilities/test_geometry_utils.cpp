//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:     BSD License
//           Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// Project includes
#include "testing/testing.h"
#include "includes/checks.h"
#include "includes/gid_io.h"
#include "utilities/geometry_utilities.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsPointDistanceToTriangleOutOfPlane, KratosCoreFastSuite)
    {
        // Triangle plane points
        Point triangle_point_1(-1.0,-1.0, 0.0);
        Point triangle_point_2( 1.0,-1.0, 0.0);
        Point triangle_point_3(-1.0, 1.0, 0.0);

        // Point out of plane to compute the distance to
        Point distance_point(0.357143, -0.214286, 0.0714286);

        
        const double dist = GeometryUtils::PointDistanceToTriangle3D(triangle_point_1, triangle_point_2, triangle_point_3, distance_point);

        KRATOS_CHECK_NEAR(dist, 0.123718, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(GeometryUtilsPointDistanceToTriangleInPlane, KratosCoreFastSuite)
    {
        // Triangle plane points
        Point triangle_point_1( 1.0,-1.0, 0.0);
        Point triangle_point_2( 1.0, 1.0, 0.0);
        Point triangle_point_3(-1.0, 1.0, 0.0);

        // Point over the plane to compute the distance to
        Point distance_point(0.357143, -0.214286, 0.0714286);
        
        const double dist = GeometryUtils::PointDistanceToTriangle3D(triangle_point_1, triangle_point_2, triangle_point_3, distance_point);

        KRATOS_CHECK_NEAR(dist, distance_point.Z(), 1e-6);
    }

}  // namespace Testing.
}  // namespace Kratos.
