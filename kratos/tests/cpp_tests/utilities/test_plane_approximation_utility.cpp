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
#include "geometries/point.h"
#include "includes/checks.h"
#include "testing/testing.h"
#include "utilities/plane_approximation_utility.h"

namespace Kratos {
namespace Testing {

    KRATOS_TEST_CASE_IN_SUITE(PlaneApproximationUtilityTriangle, KratosCoreFastSuite)
    {
        // Set of points to be approximated by the plane
        Point point_1(0.0, 0.0, 0.0);
        Point point_2(1.0, 0.0, 0.0);
        Point point_3(0.5, 1.0, 0.0);

        // Fill the array of coordinates
        std::vector< array_1d<double,3> > points_coordinates;
        points_coordinates.push_back(point_1.Coordinates());
        points_coordinates.push_back(point_2.Coordinates());
        points_coordinates.push_back(point_3.Coordinates());

        // Compute the plane approximation
        array_1d<double,3> base_point_coords, normal;
        PlaneApproximationUtility<2>::ComputePlaneApproximation(points_coordinates, base_point_coords, normal);

        // Compute and check the obtained plane distance values
        const double d_1 = inner_prod(normal , point_1.Coordinates() - base_point_coords);
        const double d_2 = inner_prod(normal , point_2.Coordinates() - base_point_coords);
        const double d_3 = inner_prod(normal , point_3.Coordinates() - base_point_coords);

        KRATOS_CHECK_NEAR(d_1, -0.5, 1e-6);
        KRATOS_CHECK_NEAR(d_2,  0.5, 1e-6);
        KRATOS_CHECK_NEAR(d_3,  0.0, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(PlaneApproximationUtilityDiamond, KratosCoreFastSuite)
    {
        // Set of points to be approximated by the plane
        Point point_1(0.0, 0.0, 0.0);
        Point point_2(1.0, 0.0, 0.0);
        Point point_3(0.5, 1.0, 0.0);
        Point point_4(1.5, 1.0, 0.0);

        // Fill the array of coordinates
        std::vector< array_1d<double,3> > points_coordinates;
        points_coordinates.push_back(point_1.Coordinates());
        points_coordinates.push_back(point_2.Coordinates());
        points_coordinates.push_back(point_3.Coordinates());
        points_coordinates.push_back(point_4.Coordinates());

        // Compute the plane approximation
        array_1d<double,3> base_point_coords, normal;
        PlaneApproximationUtility<2>::ComputePlaneApproximation(points_coordinates, base_point_coords, normal);

        // Compute and check the obtained plane distance values
        const double d_1 = inner_prod(normal , point_1.Coordinates() - base_point_coords);
        const double d_2 = inner_prod(normal , point_2.Coordinates() - base_point_coords);
        const double d_3 = inner_prod(normal , point_3.Coordinates() - base_point_coords);
        const double d_4 = inner_prod(normal , point_4.Coordinates() - base_point_coords);

        KRATOS_CHECK_NEAR(d_1, 0.0674564, 1e-6);
        KRATOS_CHECK_NEAR(d_2, -0.547956, 1e-6);
        KRATOS_CHECK_NEAR(d_3, 0.547956, 1e-6);
        KRATOS_CHECK_NEAR(d_4, -0.0674564, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(PlaneApproximationUtilityCircleTriangleIntersection, KratosCoreFastSuite)
    {
        // Set of points to be approximated by the plane
        Point point_1(0.9, 0.0, 0.0);
        Point point_2(0.0, 0.9, 0.0);
        Point point_3(0.8937, 0.1063, 0.0);
        Point point_4(0.1063, 0.8937, 0.0);

        // Fill the array of coordinates
        std::vector< array_1d<double,3> > points_coordinates;
        points_coordinates.push_back(point_1.Coordinates());
        points_coordinates.push_back(point_2.Coordinates());
        points_coordinates.push_back(point_3.Coordinates());
        points_coordinates.push_back(point_4.Coordinates());

        // Compute the plane approximation
        array_1d<double,3> base_point_coords, normal;
        PlaneApproximationUtility<2>::ComputePlaneApproximation(points_coordinates, base_point_coords, normal);

        // Compute and check the obtained plane distance values
        const double d_1 = inner_prod(normal , point_1.Coordinates() - base_point_coords);
        const double d_2 = inner_prod(normal , point_2.Coordinates() - base_point_coords);
        const double d_3 = inner_prod(normal , point_3.Coordinates() - base_point_coords);
        const double d_4 = inner_prod(normal , point_4.Coordinates() - base_point_coords);

        KRATOS_CHECK_NEAR(d_1, -0.0353553, 1e-6);
        KRATOS_CHECK_NEAR(d_2, -0.0353553, 1e-6);
        KRATOS_CHECK_NEAR(d_3, 0.0353553, 1e-6);
        KRATOS_CHECK_NEAR(d_4, 0.0353553, 1e-6);
    }

    KRATOS_TEST_CASE_IN_SUITE(PlaneApproximationUtilityAlignedIntersection, KratosCoreFastSuite)
    {
        // Set of points to be approximated by the plane
        // Note that they are selected such that, a 2D case is solved in 3D
        Point point_1(0.0, 0.0, 0.0);
        Point point_2(1.0, 0.0, 0.0);
        Point point_3(0.0, 1.0, 0.0);

        // Fill the array of coordinates
        std::vector< array_1d<double,3> > points_coordinates;
        points_coordinates.push_back(point_1.Coordinates());
        points_coordinates.push_back(point_2.Coordinates());
        points_coordinates.push_back(point_3.Coordinates());

        // Compute the plane approximation
        array_1d<double,3> base_point_coords, normal;
        PlaneApproximationUtility<3>::ComputePlaneApproximation(points_coordinates, base_point_coords, normal);

        // Compute and check the obtained plane values
        KRATOS_CHECK_NEAR(base_point_coords[0], 1.0/3.0, 1e-6);
        KRATOS_CHECK_NEAR(base_point_coords[1], 1.0/3.0, 1e-6);
        KRATOS_CHECK_NEAR(base_point_coords[2], 0.0, 1e-6);
        KRATOS_CHECK_NEAR(normal[0], 0.0, 1e-6);
        KRATOS_CHECK_NEAR(normal[1], 0.0, 1e-6);
        KRATOS_CHECK_NEAR(normal[2], 1.0, 1e-6);
    }

}  // namespace Testing.
}  // namespace Kratos.
