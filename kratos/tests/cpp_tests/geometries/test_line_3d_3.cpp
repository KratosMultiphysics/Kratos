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
#include "geometries/line_3d_3.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos::Testing {

namespace {
    /// Factory functions

    /** Generates a point type sample Line3D3N
    * @return  Pointer to a Line3D3N
    */
    Line3D3<Point>::Pointer GeneratePointsUnitXDirectionLine3D3() {
        return Kratos::make_shared<Line3D3<Point>>(
        Kratos::make_shared<Point>(0.0, 0.0, 0.0),
        Kratos::make_shared<Point>(1.0, 0.0, 0.0),
        Kratos::make_shared<Point>(0.5, 0.0, 0.0)
        );
    }

    // /** Generates a point type sample Line3D3N
    // * @return  Pointer to a Line3D3N
    // */
    // Line3D3<Point>::Pointer GeneratePointsUnitYDirectionLine3D3() {
    //     return Kratos::make_shared<Line3D3<Point>>(
    //     Kratos::make_shared<Point>(0.0, 0.0, 0.0),
    //     Kratos::make_shared<Point>(0.0, 1.0, 0.0),
    //     Kratos::make_shared<Point>(0.0, 0.5, 0.0)
    //     );
    // }

    /** Generates a point type sample Line3D3N
    * @return  Pointer to a Line3D3N
    */
    Line3D3<Point>::Pointer GeneratePointsDiagonalLine3D3() {
        return Kratos::make_shared<Line3D3<Point>>(
        Kratos::make_shared<Point>(0.0, 0.0, 0.0),
        Kratos::make_shared<Point>(1.0, 1.0, 0.0),
        Kratos::make_shared<Point>(0.5, 0.5, 0.0)
        );
    }

    /** Generates a point type sample Line3D3N
    * @return  Pointer to a Line3D3N
    */
    Line3D3<Point>::Pointer GeneratePointsParabolaLine3D3() {
        return Kratos::make_shared<Line3D3<Point>>(
        Kratos::make_shared<Point>(0.0, 0.0, 0.0),
        Kratos::make_shared<Point>(1.0, 0.0, 0.0),
        Kratos::make_shared<Point>(0.5, 0.5, 0.0)
        );
    }

    // /** Generates a point type sample Line3D3N.
    // * @return  Pointer to a Line3D3N
    // */
    // Line3D3<Point>::Pointer GenerateLine3D3WithPoints(Point::Pointer pPointOne, Point::Pointer pPointTwo, Point::Pointer pPointThree ) {
    //     return Kratos::make_shared<Line3D3<Point>>(pPointOne, pPointTwo, pPointThree);
    // }
}

/** Checks if the number of edges is correct.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3EdgesNumber, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    KRATOS_EXPECT_EQ(p_geometry->EdgesNumber(), 1);
}

/** Checks if the edges are correct.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3Edges, KratosCoreGeometriesFastSuite) {
    auto p_geom = GeneratePointsUnitXDirectionLine3D3();

    const auto& r_edges = p_geom->GenerateEdges();
    ASSERT_EQ(r_edges.size(), 1);
    for (std::size_t i = 0; i < r_edges.front().PointsNumber(); ++i) {
        KRATOS_EXPECT_NEAR(r_edges.front()[i].X(), (p_geom->pGetPoint(i))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR(r_edges.front()[i].Y(), (p_geom->pGetPoint(i))->Y(), TOLERANCE);
    }
}

/** Checks if the number of faces is correct.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3FacesNumber, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    KRATOS_EXPECT_EQ(p_geometry->FacesNumber(), 0);
}

/** Checks if the length of the line is calculated correctly.
* Checks if the length of the line is calculated correctly.
*/
KRATOS_TEST_CASE_IN_SUITE(LengthLine3D3, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsDiagonalLine3D3();

    KRATOS_EXPECT_NEAR(p_geometry->Length(), std::sqrt(2.0), TOLERANCE);

    p_geometry = GeneratePointsParabolaLine3D3();

    KRATOS_EXPECT_NEAR(p_geometry->Length(), 1.46884, 1.0e-5); // NOTE: Analytic 1.47894
}

/** Checks if the bounding box of the line is calculated correctly.
* Checks if the bounding box of the line is calculated correctly.
*/
KRATOS_TEST_CASE_IN_SUITE(BoundingBoxLine3D3, KratosCoreGeometriesFastSuite) {
    auto p_geom = GeneratePointsDiagonalLine3D3();

    Point low_point, high_point;
    p_geom->BoundingBox(low_point, high_point);

    KRATOS_EXPECT_NEAR(low_point.X(), (p_geom->pGetPoint(0))->X(), TOLERANCE);
    KRATOS_EXPECT_NEAR(low_point.Y(), (p_geom->pGetPoint(0))->Y(), TOLERANCE);
    KRATOS_EXPECT_NEAR(high_point.X(), (p_geom->pGetPoint(1))->X(), TOLERANCE);
    KRATOS_EXPECT_NEAR(high_point.Y(), (p_geom->pGetPoint(1))->Y(), TOLERANCE);

    p_geom = GeneratePointsParabolaLine3D3();

    p_geom->BoundingBox(low_point, high_point);

    KRATOS_EXPECT_NEAR(low_point.X(), (p_geom->pGetPoint(0))->X(), TOLERANCE);
    KRATOS_EXPECT_NEAR(low_point.Y(), (p_geom->pGetPoint(0))->Y(), TOLERANCE);
    KRATOS_EXPECT_NEAR(high_point.X(), (p_geom->pGetPoint(1))->X(), TOLERANCE);
    KRATOS_EXPECT_NEAR(high_point.Y(), (p_geom->pGetPoint(2))->Y(), TOLERANCE);
}

/** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
* Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    const double expected_jacobian = 0.5;

    Vector jacobians_determinants;
    p_geometry->DeterminantOfJacobian( jacobians_determinants, GeometryData::IntegrationMethod::GI_GAUSS_1 );

    for (unsigned int i=0; i<jacobians_determinants.size(); ++i) {
        KRATOS_EXPECT_NEAR(jacobians_determinants[i], expected_jacobian, TOLERANCE);
    }
}

/** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
* Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3DeterminantOfJacobianArray2, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    const double expected_jacobian = 0.5;

    Vector jacobians_determinants;
    p_geometry->DeterminantOfJacobian( jacobians_determinants, GeometryData::IntegrationMethod::GI_GAUSS_2 );

    for (unsigned int i=0; i<jacobians_determinants.size(); ++i) {
        KRATOS_EXPECT_NEAR(jacobians_determinants[i], expected_jacobian, TOLERANCE);
    }
}

/** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
* Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3DeterminantOfJacobianArray3, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    const double expected_jacobian = 0.5;

    Vector jacobians_determinants;
    p_geometry->DeterminantOfJacobian( jacobians_determinants, GeometryData::IntegrationMethod::GI_GAUSS_3 );

    for (unsigned int i=0; i<jacobians_determinants.size(); ++i) {
        KRATOS_EXPECT_NEAR(jacobians_determinants[i], expected_jacobian, TOLERANCE);
    }
}

/** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
* Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3DeterminantOfJacobianArray4, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    const double expected_jacobian = 0.5;

    Vector jacobians_determinants;
    p_geometry->DeterminantOfJacobian( jacobians_determinants, GeometryData::IntegrationMethod::GI_GAUSS_4 );

    for (unsigned int i=0; i<jacobians_determinants.size(); ++i) {
        KRATOS_EXPECT_NEAR(jacobians_determinants[i], expected_jacobian, TOLERANCE);
    }
}

/** Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
* Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3DeterminantOfJacobianArray5, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    const double expected_jacobian = 0.5;

    Vector jacobians_determinants;
    p_geometry->DeterminantOfJacobian( jacobians_determinants, GeometryData::IntegrationMethod::GI_GAUSS_5 );

    for (unsigned int i=0; i<jacobians_determinants.size(); ++i) {
        KRATOS_EXPECT_NEAR(jacobians_determinants[i], expected_jacobian, TOLERANCE);
    }
}

/** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
* Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    const double expected_jacobian = 0.5;

    double jacobian_determinant = p_geometry->DeterminantOfJacobian( 0, GeometryData::IntegrationMethod::GI_GAUSS_1 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);
}

/** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
* Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3DeterminantOfJacobianIndex2, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    double jacobian_determinant = 0.0;
    const double expected_jacobian = 0.5;

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 0, GeometryData::IntegrationMethod::GI_GAUSS_2 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_2 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);
}

/** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
* Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3DeterminantOfJacobianIndex3, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    double jacobian_determinant = 0.0;
    const double expected_jacobian = 0.5;

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 0, GeometryData::IntegrationMethod::GI_GAUSS_3 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_3 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_3 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);
}

/** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
* Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3DeterminantOfJacobianIndex4, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    double jacobian_determinant = 0.0;
    const double expected_jacobian = 0.5;

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 0, GeometryData::IntegrationMethod::GI_GAUSS_4 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_4 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_4 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 3, GeometryData::IntegrationMethod::GI_GAUSS_4 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);
}

/** Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
* Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Line3D3DeterminantOfJacobianIndex5, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsUnitXDirectionLine3D3();
    double jacobian_determinant = 0.0;
    const double expected_jacobian = 0.5;

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 0, GeometryData::IntegrationMethod::GI_GAUSS_5 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_5 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_5 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 3, GeometryData::IntegrationMethod::GI_GAUSS_5 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);

    jacobian_determinant = p_geometry->DeterminantOfJacobian( 4, GeometryData::IntegrationMethod::GI_GAUSS_5 );
    KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Line3D3ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsParabolaLine3D3();
    auto& r_geom = *p_geometry;
    auto p_p_geom_nodes = Kratos::make_shared<Line3D3<Node>>(
    Kratos::make_intrusive<Node>(1, r_geom[0].X(), r_geom[0].Y(), r_geom[0].Z()),
    Kratos::make_intrusive<Node>(2, r_geom[1].X(), r_geom[1].Y(), r_geom[1].Z()),
    Kratos::make_intrusive<Node>(3, r_geom[2].X(), r_geom[2].Y(), r_geom[2].Z())
    );
    CrossCheckShapeFunctionsValues(*p_p_geom_nodes);
}

KRATOS_TEST_CASE_IN_SUITE(Line3D3ShapeFunctionsValuesMatrix, KratosCoreGeometriesFastSuite) {

    Geometry<Point>::Pointer p_geom = GeneratePointsDiagonalLine3D3();

    const Matrix N_values_geom = p_geom->ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

    KRATOS_EXPECT_NEAR(N_values_geom(0, 0), 0.455342, TOLERANCE);
    KRATOS_EXPECT_NEAR(N_values_geom(0, 1), -0.122008, TOLERANCE);
    KRATOS_EXPECT_NEAR(N_values_geom(0, 2), 0.666667, TOLERANCE);
    KRATOS_EXPECT_NEAR(N_values_geom(1, 0), -0.122008, TOLERANCE);
    KRATOS_EXPECT_NEAR(N_values_geom(1, 1), 0.455342, TOLERANCE);
    KRATOS_EXPECT_NEAR(N_values_geom(1, 2), 0.666667, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Line3D3ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
    auto p_geometry = GeneratePointsParabolaLine3D3();
    auto& r_geom = *p_geometry;
    auto p_p_geom_nodes = Kratos::make_shared<Line3D3<Node>>(
    Kratos::make_intrusive<Node>(1, r_geom[0].X(), r_geom[0].Y(), r_geom[0].Z()),
    Kratos::make_intrusive<Node>(2, r_geom[1].X(), r_geom[1].Y(), r_geom[1].Z()),
    Kratos::make_intrusive<Node>(3, r_geom[2].X(), r_geom[2].Y(), r_geom[2].Z())
    );
    TestAllShapeFunctionsLocalGradients(*p_p_geom_nodes);
}

KRATOS_TEST_CASE_IN_SUITE(Line3D3PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
    auto p0 = Kratos::make_shared<Point>(0.0, 0.0, 0.0);
    auto p1 = Kratos::make_shared<Point>(1.0, 1.0, 1.0);
    auto p2 = Kratos::make_shared<Point>(-0.2, 0.5, 0.7);
    auto geometry = Line3D3<Point>(p0, p1, p2);

    array_1d<double, 3> result(3, 0.0);
    array_1d<double, 3> expected{0.7, 0.0, 0.0};
    Vector shapes;
    geometry.ShapeFunctionsValues(shapes, expected);
    Point p(shapes[0]*(*p0) + shapes[1]*(*p1) + shapes[2]*(*p2));

    geometry.PointLocalCoordinates(result, p);
    KRATOS_EXPECT_VECTOR_NEAR(result, expected, TOLERANCE);

    array_1d<double, 3> local0{-1.0, 0.0, 0.0};
    geometry.PointLocalCoordinates(result, *p0);
    KRATOS_EXPECT_VECTOR_NEAR(result, local0, TOLERANCE);

    array_1d<double, 3> local1{1.0, 0.0, 0.0};
    geometry.PointLocalCoordinates(result, *p1);
    KRATOS_EXPECT_VECTOR_NEAR(result, local1, TOLERANCE);

    array_1d<double, 3> local2{0.0, 0.0, 0.0};
    geometry.PointLocalCoordinates(result, *p2);
    KRATOS_EXPECT_VECTOR_NEAR(result, local2, TOLERANCE);

    // This one should delegate to line3D2
    auto pl = Kratos::make_shared<Point>(0.5, 0.5, 0.5);
    auto px = Kratos::make_shared<Point>(0.75, 0.75, 0.75);
    auto flat_geometry = Line3D3<Point>(p0, p1, pl);
    flat_geometry.PointLocalCoordinates(result, *px);
    array_1d<double,3> expected_flat{0.5, 0.0, 0.0};
    KRATOS_EXPECT_VECTOR_NEAR(result, expected_flat, TOLERANCE);
}

} // namespace Kratos::Testing.
