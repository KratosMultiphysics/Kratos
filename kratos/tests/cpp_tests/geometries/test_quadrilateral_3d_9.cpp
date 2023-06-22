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
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/quadrilateral_3d_9.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

// Utility includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos::Testing {
namespace {
    /// Test utility functions

    /** Generates a sample Quadrilateral3D9.
    * Generates a right quadrilateral with origin in the origin and leg size 1.
    * @return  Pointer to a Quadrilateral3D9
    */
    template<class TPointType>
    typename Quadrilateral3D9<TPointType>::Pointer GenerateFlatQuadrilateral3D9()
    {
        return typename Quadrilateral3D9<TPointType>::Pointer(new Quadrilateral3D9<TPointType>(
        GeneratePoint<TPointType>( 0.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 1.0, 0.0),
        GeneratePoint<TPointType>( 0.0, 1.0, 0.0),
        GeneratePoint<TPointType>( 0.5, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 0.5, 0.0),
        GeneratePoint<TPointType>( 0.5, 1.0, 0.0),
        GeneratePoint<TPointType>( 0.0, 0.5, 0.0),
        GeneratePoint<TPointType>( 0.5, 0.5, 0.0)
        ));
    }

    template<class TPointType>
    void CheckSamePoint(const TPointType& rThisPoint, const TPointType& rThatPoint) {
        KRATOS_CHECK_NEAR(rThisPoint.X(), rThatPoint.X(), TOLERANCE);
        KRATOS_CHECK_NEAR(rThisPoint.Y(), rThatPoint.Y(), TOLERANCE);
        KRATOS_CHECK_NEAR(rThisPoint.Z(), rThatPoint.Z(), TOLERANCE);
    }
}

/// Tests

/** Checks if the number of edges is correct.
* Checks if the number of edges is correct.
*/
KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D9EdgesNumber, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateFlatQuadrilateral3D9<Node>();

    KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 4);
}

/** Checks if the number of edges is correct.
* Checks if the number of edges is correct.
*/
KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D9Edges, KratosCoreGeometriesFastSuite)
{
    auto p_geom = GenerateFlatQuadrilateral3D9<Node>();

    const auto& r_edges = p_geom->GenerateEdges();
    CheckSamePoint(r_edges[0][0], p_geom->GetPoint(0));
    CheckSamePoint(r_edges[0][1], p_geom->GetPoint(1));
    CheckSamePoint(r_edges[0][2], p_geom->GetPoint(4));

    CheckSamePoint(r_edges[1][0], p_geom->GetPoint(1));
    CheckSamePoint(r_edges[1][1], p_geom->GetPoint(2));
    CheckSamePoint(r_edges[1][2], p_geom->GetPoint(5));

    CheckSamePoint(r_edges[2][0], p_geom->GetPoint(2));
    CheckSamePoint(r_edges[2][1], p_geom->GetPoint(3));
    CheckSamePoint(r_edges[2][2], p_geom->GetPoint(6));

    CheckSamePoint(r_edges[3][0], p_geom->GetPoint(3));
    CheckSamePoint(r_edges[3][1], p_geom->GetPoint(0));
    CheckSamePoint(r_edges[3][2], p_geom->GetPoint(7));
}

/** Checks if the area of the quadrilateral is calculated correctly.
* Checks if the area of the quadrilateral is calculated correctly.
*/
KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D9Area, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateFlatQuadrilateral3D9<Node>();
    KRATOS_CHECK_NEAR(geom->Area(), 1.0, TOLERANCE);
}

/** Tests the PointLocalCoordinates for Quadrilateral3D9.
 * Tests the PointLocalCoordinates for Quadrilateral3D9.
 */
KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D9PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateFlatQuadrilateral3D9<Node>();

    Point TestPointA(1.0, 1.0, 0.0);
    Point TestPointB(0.5, 0.5, 0.0);
    Point TestResultA(0.0, 0.0, 0.0);
    Point TestResultB(0.0, 0.0, 0.0);

    geom->PointLocalCoordinates(TestResultA, TestPointA);
    geom->PointLocalCoordinates(TestResultB, TestPointB);

    // Test transformation in the edge
    KRATOS_CHECK_NEAR(TestResultA[0], 1.0, TOLERANCE);
    KRATOS_CHECK_NEAR(TestResultA[1], 1.0, TOLERANCE);
    KRATOS_CHECK_NEAR(TestResultA[2], 0.0, TOLERANCE);

    // Test transformation in the center
    KRATOS_CHECK_NEAR(TestResultB[0], 0.0, TOLERANCE);
    KRATOS_CHECK_NEAR(TestResultB[1], 0.0, TOLERANCE);
    KRATOS_CHECK_NEAR(TestResultB[2], 0.0, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D9ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateFlatQuadrilateral3D9<Node>();
    CrossCheckShapeFunctionsValues(*geom);
}

KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D9ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateFlatQuadrilateral3D9<Node>();
    TestAllShapeFunctionsLocalGradients(*geom);
}

KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D9Normal, KratosCoreGeometriesFastSuite) {
    Geometry<Point>::Pointer p_geom = Kratos::make_shared<Quadrilateral3D9<Point>>(
        Kratos::make_shared<Point>(0.9, 0.25, -0.15),
        Kratos::make_shared<Point>(1.0, 0.25, -0.15),
        Kratos::make_shared<Point>(1.0, 0.25, 0.05),
        Kratos::make_shared<Point>(0.9, 0.25, 0.05),
        Kratos::make_shared<Point>(0.95, 0.25, -0.15),
        Kratos::make_shared<Point>(1.0, 0.25, -0.05),
        Kratos::make_shared<Point>(0.95, 0.25, 0.05),
        Kratos::make_shared<Point>(0.9, 0.25, -0.05),
        Kratos::make_shared<Point>(0.95, 0.25, -0.05)
    );

    array_1d<double, 3> local_coordinates{1.0, 1.0, 0.0};
    auto normal = p_geom->Normal(local_coordinates);
    auto unit_normal = p_geom->UnitNormal(local_coordinates);

    array_1d<double, 3> expected_unit_normal{0.0, -1.0, 0.0};
    KRATOS_CHECK_VECTOR_NEAR(unit_normal, expected_unit_normal, TOLERANCE);

    KRATOS_CHECK_NEAR(MathUtils<double>::Norm3(normal), MathUtils<double>::Dot3(unit_normal, normal), TOLERANCE);
}

}
