//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohamed Nabi
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/line_2d_4.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos {
namespace Testing {

    /// Factory functions

    /** Generates a point type sample Line2D3N
    * @return  Pointer to a Line2D4N
    */
    Line2D4<Point>::Pointer GeneratePointsUnitXDirectionLine2D4() {
        return Kratos::make_shared<Line2D4<Point>>(
        Kratos::make_shared<Point>(0.0    , 0.0, 0.0),
        Kratos::make_shared<Point>(1.0    , 0.0, 0.0),
        Kratos::make_shared<Point>(1.0/3.0, 0.0, 0.0),
        Kratos::make_shared<Point>(2.0/3.0, 0.0, 0.0)
        );
    }

    /** Generates a point type sample Line2D4N
    * @return  Pointer to a Line2D4N
    */
    Line2D4<Point>::Pointer GeneratePointsUnitYDirectionLine2D4() {
        return Kratos::make_shared<Line2D4<Point>>(
        Kratos::make_shared<Point>(0.0, 0.0    , 0.0),
        Kratos::make_shared<Point>(0.0, 1.0    , 0.0),
        Kratos::make_shared<Point>(0.0, 1.0/3.0, 0.0),
        Kratos::make_shared<Point>(0.0, 2.0/3.0, 0.0)
        );
    }

    /** Generates a point type sample Line2D4N
    * @return  Pointer to a Line2D4N
    */
    Line2D4<Point>::Pointer GeneratePointsDiagonalLine2D4() {
        return Kratos::make_shared<Line2D4<Point>>(
        Kratos::make_shared<Point>(0.0    , 0.0    , 0.0),
        Kratos::make_shared<Point>(1.0    , 1.0    , 0.0),
        Kratos::make_shared<Point>(1.0/3.0, 1.0/3.0, 0.0),
        Kratos::make_shared<Point>(2.0/3.0, 2.0/3.0, 0.0)
        );
    }

    /** Generates a point type sample Line2D4N
    * @return  Pointer to a Line2D4N
    */
    Line2D4<Point>::Pointer GeneratePointsParabolaLine2D4() {
        return Kratos::make_shared<Line2D4<Point>>(
        Kratos::make_shared<Point>(0.0    , 0.0    , 0.0),
        Kratos::make_shared<Point>(1.0    , 0.0    , 0.0),
        Kratos::make_shared<Point>(1.0/3.0, 4.0/9.0, 0.0),
        Kratos::make_shared<Point>(2.0/3.0, 4.0/9.0, 0.0)
        );
    }

    /** Generates a point type sample Line2D4N.
    * @return  Pointer to a Line2D4N
    */
    Line2D4<Point>::Pointer GenerateLine2D4WithPoints(Point::Pointer pPoint01, Point::Pointer pPoint02,
                                                      Point::Pointer pPoint03, Point::Pointer pPoint04) {
        return Kratos::make_shared<Line2D4<Point>>(pPoint01, pPoint02, pPoint03, pPoint04);
    }

    /** Checks if the number of edges is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D4EdgesNumber, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
        KRATOS_EXPECT_EQ(p_geometry->EdgesNumber(), 1);
    }

    /** Checks if the edges are correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D4Edges, KratosCoreGeometriesFastSuite) {
        auto p_geom = GeneratePointsUnitXDirectionLine2D4();

        const auto& r_edges = p_geom->GenerateEdges();
        ASSERT_EQ(r_edges.size(), 1);
        for (std::size_t i = 0; i < r_edges.front().PointsNumber(); ++i) {
            KRATOS_EXPECT_NEAR(r_edges.front()[i].X(), (p_geom->pGetPoint(i))->X(), TOLERANCE);
            KRATOS_EXPECT_NEAR(r_edges.front()[i].Y(), (p_geom->pGetPoint(i))->Y(), TOLERANCE);
        }
    }

    /** Checks if the number of faces is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D4FacesNumber, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
        KRATOS_EXPECT_EQ(p_geometry->FacesNumber(), 0);
    }

    /** Checks if the length of the line is calculated correctly.
    * Checks if the length of the line is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(LengthLine2D4, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsDiagonalLine2D4();

        KRATOS_EXPECT_NEAR(p_geometry->Length(), std::sqrt(2.0), TOLERANCE);

        p_geometry = GeneratePointsParabolaLine2D4();

        KRATOS_EXPECT_NEAR(p_geometry->Length(), 1.48139, 1.0e-5); // NOTE: Analytic 1.478942858
    }

    /** Checks if the bounding box of the line is calculated correctly.
    * Checks if the bounding box of the line is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(BoundingBoxLine2D4, KratosCoreGeometriesFastSuite) {
        auto p_geom = GeneratePointsDiagonalLine2D4();

        Point low_point, high_point;
        p_geom->BoundingBox(low_point, high_point);

        KRATOS_EXPECT_NEAR(low_point.X(), (p_geom->pGetPoint(0))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR(low_point.Y(), (p_geom->pGetPoint(0))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR(high_point.X(), (p_geom->pGetPoint(1))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR(high_point.Y(), (p_geom->pGetPoint(1))->Y(), TOLERANCE);

        p_geom = GeneratePointsParabolaLine2D4();

        p_geom->BoundingBox(low_point, high_point);

        KRATOS_EXPECT_NEAR(low_point.X(), (p_geom->pGetPoint(0))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR(low_point.Y(), (p_geom->pGetPoint(0))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR(high_point.X(), (p_geom->pGetPoint(1))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR(high_point.Y(), (p_geom->pGetPoint(2))->Y(), TOLERANCE);
    }

    /** Checks the inside test for a given point respect to the line
    * Checks the inside test for a given point respect to the line
    * It performs 4 tests:
    * A Point inside the line: Expected result TRUE
    * A Point outside the line: Expected result FALSE
    * A Point over a vertex 1 of the line: Expected result TRUE
    * A Point over an vertex 2 of the line: Expected result TRUE
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D4IsInside, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsDiagonalLine2D4();

        Point PointInside(1.0/3.0, 1.0/3.0, 0.0);
        Point PointOutside(1.66, 0.66, 0.0);
        Point PointInVertex(0.0, 0.0, 0.0);
        Point PointInEdge(1.0, 1.0, 0.0);

        Point LocalCoords;

        KRATOS_EXPECT_TRUE(p_geometry->IsInside(PointInside, LocalCoords, EPSILON));
        KRATOS_EXPECT_FALSE(p_geometry->IsInside(PointOutside, LocalCoords, EPSILON));
        KRATOS_EXPECT_TRUE(p_geometry->IsInside(PointInVertex, LocalCoords, EPSILON));
        KRATOS_EXPECT_TRUE(p_geometry->IsInside(PointInEdge, LocalCoords, EPSILON));
    }

    /** Checks the point local coordinates for a given point respect to the
    * line. The baricentre of the triangle is selected due to its known
    * solution.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D4PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsDiagonalLine2D4();

        // Compute the global coordinates of the baricentre
        const Geometry<Point>::PointsArrayType p_geom_pts = p_geometry->Points();
        Point centre = Point{p_geom_pts[0] + p_geom_pts[1] + p_geom_pts[2] + p_geom_pts[3]};
        centre *= 1.0/4.0;

        // Compute the centre local coordinates
        array_1d<double, 3> centre_local_coords;
        p_geometry->PointLocalCoordinates(centre_local_coords, centre);

        KRATOS_EXPECT_NEAR(centre_local_coords(0), 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(centre_local_coords(1), 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(centre_local_coords(2), 0.0, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D4DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
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
    KRATOS_TEST_CASE_IN_SUITE(Line2D4DeterminantOfJacobianArray2, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
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
    KRATOS_TEST_CASE_IN_SUITE(Line2D4DeterminantOfJacobianArray3, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
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
    KRATOS_TEST_CASE_IN_SUITE(Line2D4DeterminantOfJacobianArray4, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
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
    KRATOS_TEST_CASE_IN_SUITE(Line2D4DeterminantOfJacobianArray5, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
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
    KRATOS_TEST_CASE_IN_SUITE(Line2D4DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
        const double expected_jacobian = 0.5;

        double jacobian_determinant = p_geometry->DeterminantOfJacobian( 0, GeometryData::IntegrationMethod::GI_GAUSS_1 );
        KRATOS_EXPECT_NEAR(jacobian_determinant, expected_jacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D4DeterminantOfJacobianIndex2, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
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
    KRATOS_TEST_CASE_IN_SUITE(Line2D4DeterminantOfJacobianIndex3, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
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
    KRATOS_TEST_CASE_IN_SUITE(Line2D4DeterminantOfJacobianIndex4, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
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

    /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D4DeterminantOfJacobianIndex5, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsUnitXDirectionLine2D4();
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

    KRATOS_TEST_CASE_IN_SUITE(Line2D4ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsParabolaLine2D4();
        auto& r_geom = *p_geometry;
        auto p_p_geom_nodes = Kratos::make_shared<Line2D4<Node>>(
        Kratos::make_intrusive<Node>(1, r_geom[0].X(), r_geom[0].Y(), r_geom[0].Z()),
        Kratos::make_intrusive<Node>(2, r_geom[1].X(), r_geom[1].Y(), r_geom[1].Z()),
        Kratos::make_intrusive<Node>(3, r_geom[2].X(), r_geom[2].Y(), r_geom[2].Z()),
        Kratos::make_intrusive<Node>(4, r_geom[3].X(), r_geom[3].Y(), r_geom[3].Z())
        );
        CrossCheckShapeFunctionsValues(*p_p_geom_nodes);
    }

    KRATOS_TEST_CASE_IN_SUITE(Line2D4ShapeFunctionsValuesMatrix, KratosCoreGeometriesFastSuite) {

        Geometry<Point>::Pointer p_geom = GeneratePointsDiagonalLine2D4();

        const Matrix N_values_geom = p_geom->ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

        KRATOS_EXPECT_NEAR(N_values_geom(0, 0),  0.19716878380, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_geom(0, 1),  0.05283121638, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_geom(0, 2),  1.02451905300, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_geom(0, 3), -0.27451905300, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_geom(1, 0),  0.05283121638, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_geom(1, 1),  0.19716878380, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_geom(1, 2), -0.27451905300, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_geom(1, 3),  1.02451905300, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Line2D4ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto p_geometry = GeneratePointsParabolaLine2D4();
        auto& r_geom = *p_geometry;
        auto p_p_geom_nodes = Kratos::make_shared<Line2D4<Node>>(
        Kratos::make_intrusive<Node>(1, r_geom[0].X(), r_geom[0].Y(), r_geom[0].Z()),
        Kratos::make_intrusive<Node>(2, r_geom[1].X(), r_geom[1].Y(), r_geom[1].Z()),
        Kratos::make_intrusive<Node>(3, r_geom[2].X(), r_geom[2].Y(), r_geom[2].Z()),
        Kratos::make_intrusive<Node>(4, r_geom[3].X(), r_geom[3].Y(), r_geom[3].Z())
        );
        TestAllShapeFunctionsLocalGradients(*p_p_geom_nodes);
    }

} // namespace Testing.
} // namespace Kratos.
