//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/hexahedra_3d_27.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos::Testing
{
using Hexa27GeometryType = Hexahedra3D27<Node>;
using Hexa27GeometryPtrType = Hexa27GeometryType::Pointer;

/** Generates a sample Hexahedra3D27.
 * Generates a hexahedra with center on the origin with positive volume and side 1.
 * @return  Pointer to a Hexahedra3D27
 */
Hexa27GeometryPtrType GenerateCanonicalHexahedra3D27()
{
    auto p0  = GeneratePoint<Node>(0.0, 0.0, 0.0);
    auto p1  = GeneratePoint<Node>(1.0, 0.0, 0.0);
    auto p2  = GeneratePoint<Node>(1.0, 1.0, 0.0);
    auto p3  = GeneratePoint<Node>(0.0, 1.0, 0.0);
    auto p4  = GeneratePoint<Node>(0.0, 0.0, 1.0);
    auto p5  = GeneratePoint<Node>(1.0, 0.0, 1.0);
    auto p6  = GeneratePoint<Node>(1.0, 1.0, 1.0);
    auto p7  = GeneratePoint<Node>(0.0, 1.0, 1.0);
    auto p8  = GeneratePoint<Node>(0.5, 0.0, 0.0);
    auto p9  = GeneratePoint<Node>(1.0, 0.5, 0.0);
    auto p10 = GeneratePoint<Node>(0.5, 1.0, 0.0);
    auto p11 = GeneratePoint<Node>(0.0, 0.5, 0.0);
    auto p12 = GeneratePoint<Node>(0.0, 0.0, 0.5);
    auto p13 = GeneratePoint<Node>(1.0, 0.0, 0.5);
    auto p14 = GeneratePoint<Node>(1.0, 1.0, 0.5);
    auto p15 = GeneratePoint<Node>(0.0, 1.0, 0.5);
    auto p16 = GeneratePoint<Node>(0.5, 0.0, 1.0);
    auto p17 = GeneratePoint<Node>(1.0, 0.5, 1.0);
    auto p18 = GeneratePoint<Node>(0.5, 1.0, 1.0);
    auto p19 = GeneratePoint<Node>(0.0, 0.5, 1.0);
    auto p20 = GeneratePoint<Node>(0.5, 0.5, 0.0);
    auto p21 = GeneratePoint<Node>(0.5, 0.0, 0.5);
    auto p22 = GeneratePoint<Node>(1.0, 0.5, 0.5);
    auto p23 = GeneratePoint<Node>(0.5, 1.0, 0.5);
    auto p24 = GeneratePoint<Node>(0.0, 0.5, 0.5);
    auto p25 = GeneratePoint<Node>(0.5, 0.5, 1.0);
    auto p26 = GeneratePoint<Node>(0.5, 0.5, 0.5);
    return Hexa27GeometryPtrType(new Hexa27GeometryType(p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21, p22, p23, p24, p25, p26));
}

/** Checks if the number of edges is correct.
 * Checks if the number of edges is correct.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27EdgesNumber, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D27();

    KRATOS_EXPECT_EQ(geom->EdgesNumber(), 12);
}

/** Checks if the number of faces is correct.
 * Checks if the number of faces is correct.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27FacesNumber, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D27();

    KRATOS_EXPECT_EQ(geom->FacesNumber(), 6);
}

KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27Jacobian, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D27();
    Matrix J = ZeroMatrix(3,3), expected_J = 0.5*IdentityMatrix(3, 3);
    geom->Jacobian(J, array_1d<double, 3>(3, 0.0));

    KRATOS_EXPECT_MATRIX_NEAR(J, expected_J, TOLERANCE);

    KRATOS_EXPECT_NEAR(geom->DeterminantOfJacobian(array_1d<double,3>(3, 0.0)), 0.125, TOLERANCE);

    Hexa27GeometryType::JacobiansType integration_point_jacobians;
    geom->Jacobian(integration_point_jacobians);

    KRATOS_EXPECT_EQ(integration_point_jacobians.size(), 27);
    for (const auto& jacobian: integration_point_jacobians) {
        KRATOS_EXPECT_NEAR(MathUtils<double>::Det3(jacobian), 0.125, 1.0);
    }
}

/** Checks if the characteristic length of the hexahedra is calculated correctly.
 * Checks if the characteristic length of the hexahedra is calculated correctly.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27Length, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D27();

    KRATOS_EXPECT_NEAR(geom->Length(), 0.353553, TOLERANCE);
}

/** Checks if the area of the hexahedra is calculated correctly.
 * Checks if the area of the hexahedra is calculated correctly.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27Area, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D27();

    KRATOS_EXPECT_NEAR(geom->Area(), 1.0, TOLERANCE);
}

/** Checks if the volume of the hexahedra is calculated correctly.
 * Checks if the volume of the hexahedra is calculated correctly.
 * For hexahedra 3D8 'volume()' call defaults to 'area()'
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27Volume, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D27();

    KRATOS_EXPECT_NEAR(geom->Volume(), 1.0, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D27();
    TestAllShapeFunctionsLocalGradients(*geom);
}

/**
 * This test checks the HasIntersection method for bounding boxes
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27BoxIntersection, KratosCoreGeometriesFastSuite)
{
    auto hexahedron = GenerateCanonicalHexahedra3D27();

    //hexahedron inside the box
    KRATOS_EXPECT_TRUE(hexahedron->HasIntersection(Point(-0.1,-0.1,-0.1), Point(1.1,1.1,1.1)));

    //hexahedron contains the box
    KRATOS_EXPECT_TRUE(hexahedron->HasIntersection(Point(.25,.25,.25), Point(.75,.75,.75)));

    //hexahedron intersects the box
    KRATOS_EXPECT_TRUE(hexahedron->HasIntersection(Point(.25,.25,.25), Point(1.5,1.5,1.5)));

    //hexahedron not intersects the box
    KRATOS_EXPECT_FALSE(hexahedron->HasIntersection(Point(1.1,1.1,1.1), Point(2.1,2.1,2.1)));
}

KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27AverageEdgeLength, KratosCoreGeometriesFastSuite)
{
    auto hexahedron = GenerateCanonicalHexahedra3D27();

    KRATOS_EXPECT_NEAR(hexahedron->AverageEdgeLength(), 1.0, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D27PointsLocalCoordinates, KratosCoreGeometriesFastSuite)
{
    // Test local PointsLocalCoordinates by checking that the provided local coordinates
    // can be used to interpolate the correct global coordinate for each point
    auto hexahedron = GenerateCanonicalHexahedra3D27();
    Matrix points_local_coordinates = ZeroMatrix(27, 3);
    hexahedron->PointsLocalCoordinates(points_local_coordinates);
    array_1d<double, 3> local_coordinates(3, 0.0);
    array_1d<double, 3> global_coordinates(3, 0.0);

    const auto& r_points = hexahedron->Points();

    for (std::size_t i_point = 0; i_point < 27; i_point++) {
        local_coordinates[0] = points_local_coordinates(i_point, 0);
        local_coordinates[1] = points_local_coordinates(i_point, 1);
        local_coordinates[2] = points_local_coordinates(i_point, 2);

        global_coordinates = hexahedron->GlobalCoordinates(global_coordinates, local_coordinates);

        KRATOS_EXPECT_NEAR(global_coordinates[0], r_points[i_point].X(), TOLERANCE);
        KRATOS_EXPECT_NEAR(global_coordinates[1], r_points[i_point].Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR(global_coordinates[2], r_points[i_point].Z(), TOLERANCE);
    }
}

}  // namespace Kratos::Testing.
