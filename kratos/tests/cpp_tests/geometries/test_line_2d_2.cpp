//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/line_2d_2.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos {
namespace Testing {

    /// Factory functions

    /** Generates a point type sample Line2D2N
    * @return  Pointer to a Line2D2N
    */
    Line2D2<Point>::Pointer GeneratePointsUnitXDirectionLine2D2() {
        return Kratos::make_shared<Line2D2<Point>>(
        Kratos::make_shared<Point>(0.0, 0.0, 0.0),
        Kratos::make_shared<Point>(1.0, 0.0, 0.0)
        );
    }

    /** Generates a point type sample Line2D2N
    * @return  Pointer to a Line2D2N
    */
    Line2D2<Point>::Pointer GeneratePointsUnitYDirectionLine2D2() {
        return Kratos::make_shared<Line2D2<Point>>(
        Kratos::make_shared<Point>(0.0, 0.0, 0.0),
        Kratos::make_shared<Point>(0.0, 1.0, 0.0)
        );
    }

    /** Generates a point type sample Line2D2N
    * @return  Pointer to a Line2D2N
    */
    Line2D2<Point>::Pointer GeneratePointsDiagonalLine2D2() {
        return Kratos::make_shared<Line2D2<Point>>(
        Kratos::make_shared<Point>(0.0, 0.0, 0.0),
        Kratos::make_shared<Point>(1.0, 1.0, 0.0)
        );
    }

    /** Generates a point type sample Line2D2N.
    * @return  Pointer to a Line2D2N
    */
    Line2D2<Point>::Pointer GenerateLine2D2WithPoints(Point::Pointer rPointOne, Point::Pointer rPointTwo ) {
        return Kratos::make_shared<Line2D2<Point>>(rPointOne, rPointTwo);
    }

    /** Checks if the number of edges is correct.
    * Checks if the number of edges is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2EdgesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 1);
    }

    /** Checks if the edges are correct.
    * Checks if the edges are correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2Edges, KratosCoreGeometriesFastSuite) {
        auto p_geom = GeneratePointsUnitXDirectionLine2D2();

        const auto& r_edges = p_geom->GenerateEdges();
        KRATOS_CHECK_NEAR((r_edges[0])[0].X(), (p_geom->pGetPoint(0))->X(), TOLERANCE);
        KRATOS_CHECK_NEAR((r_edges[0])[0].Y(), (p_geom->pGetPoint(0))->Y(), TOLERANCE);
        KRATOS_CHECK_NEAR((r_edges[0])[1].X(), (p_geom->pGetPoint(1))->X(), TOLERANCE);
        KRATOS_CHECK_NEAR((r_edges[0])[1].Y(), (p_geom->pGetPoint(1))->Y(), TOLERANCE);
    }

    /** Checks if the number of faces is correct.
    * Checks if the number of faces is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2FacesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        KRATOS_CHECK_EQUAL(geom->FacesNumber(), 0);
    }

    /** Checks if the length of the line is calculated correctly.
    * Checks if the length of the line is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(LengthLine2D2, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsDiagonalLine2D2();

        KRATOS_CHECK_NEAR(geom->Length(), std::sqrt(2.0), TOLERANCE);
    }

    /** Checks if the bounding box of the line is calculated correctly.
    * Checks if the bounding box of the line is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(BoundingBoxLine2D2, KratosCoreGeometriesFastSuite) {
        auto p_geom = GeneratePointsDiagonalLine2D2();

        Point low_point, high_point;
        p_geom->BoundingBox(low_point, high_point);

        KRATOS_CHECK_NEAR(low_point.X(), (p_geom->pGetPoint(0))->X(), TOLERANCE);
        KRATOS_CHECK_NEAR(low_point.Y(), (p_geom->pGetPoint(0))->Y(), TOLERANCE);
        KRATOS_CHECK_NEAR(high_point.X(), (p_geom->pGetPoint(1))->X(), TOLERANCE);
        KRATOS_CHECK_NEAR(high_point.Y(), (p_geom->pGetPoint(1))->Y(), TOLERANCE);
    }

    /** Checks the ProjectionPoint test for a given point respect to the line
    * Checks the ProjectionPoint test for a given point respect to the line
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2ProjectionPoint, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsDiagonalLine2D2();

        Point point(0.5, 0.55, 0.0);

        Geometry<Point>::CoordinatesArrayType global_coords;
        Geometry<Point>::CoordinatesArrayType local_coords;

        geom->ProjectionPoint(point.Coordinates(), global_coords, local_coords);

        const Point point_to_project(point);
        Point point_projected;
        GeometricalProjectionUtilities::FastProjectOnLine2D(*geom, point_to_project, point_projected);

        KRATOS_CHECK_RELATIVE_NEAR(global_coords[0], point_projected[0], TOLERANCE);
        KRATOS_CHECK_RELATIVE_NEAR(global_coords[1], point_projected[1], TOLERANCE);
        KRATOS_CHECK_RELATIVE_NEAR(global_coords[2], point_projected[2], TOLERANCE);

        KRATOS_CHECK_RELATIVE_NEAR(local_coords[0], 0.05, TOLERANCE);
        KRATOS_CHECK_NEAR(local_coords[1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(local_coords[2], 0.0, TOLERANCE);
    }

    /** Checks the inside test for a given point respect to the line
    * Checks the inside test for a given point respect to the line
    * It performs 4 tests:
    * A Point inside the line: Expected result TRUE
    * A Point outside the line: Expected result FALSE
    * A Point over a vertex 1 of the line: Expected result TRUE
    * A Point over an vertex 2 of the line: Expected result TRUE
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2IsInside, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsDiagonalLine2D2();

        Point PointInside(0.5, 0.5, 0.0);
        Point PointOutside(1.66, 0.66, 0.0);
        Point PointInVertex(0.0, 0.0, 0.0);
        Point PointInEdge(1.0, 1.0, 0.0);

        Point LocalCoords;

        KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
        KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
        KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
        KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
    }

    /** Checks the point local coordinates for a given point respect to the
    * line. The baricentre of the triangle is selected due to its known
    * solution.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsDiagonalLine2D2();

        // Compute the global coordinates of the baricentre
        const Geometry<Point>::PointsArrayType geom_pts = geom->Points();
        Point centre = Point{geom_pts[0] + geom_pts[1]};
        centre *= 1.0/2.0;

        // Compute the centre local coordinates
        array_1d<double, 3> centre_local_coords;
        geom->PointLocalCoordinates(centre_local_coords, centre);

        KRATOS_CHECK_NEAR(centre_local_coords(0), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(centre_local_coords(1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(centre_local_coords(2), 0.0, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        const double ExpectedJacobian = 0.5;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_1 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i) {
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2DeterminantOfJacobianArray2, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        const double ExpectedJacobian = 0.5;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_2 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i) {
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2DeterminantOfJacobianArray3, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        const double ExpectedJacobian = 0.5;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_3 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i) {
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2DeterminantOfJacobianArray4, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        const double ExpectedJacobian = 0.5;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_4 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i) {
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2DeterminantOfJacobianArray5, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        const double ExpectedJacobian = 0.5;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_5 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i) {
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        const double ExpectedJacobian = 0.5;

        double JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_1 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2DeterminantOfJacobianIndex2, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 0.5;

        JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_2 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::GI_GAUSS_2 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2DeterminantOfJacobianIndex3, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 0.5;

        JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_3 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::GI_GAUSS_3 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 3, GeometryData::GI_GAUSS_3 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2DeterminantOfJacobianIndex4, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 0.5;

        JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_4 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::GI_GAUSS_4 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 3, GeometryData::GI_GAUSS_4 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 4, GeometryData::GI_GAUSS_4 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2DeterminantOfJacobianIndex5, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 0.5;

        JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_5 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::GI_GAUSS_5 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 3, GeometryData::GI_GAUSS_5 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 4, GeometryData::GI_GAUSS_5 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 5, GeometryData::GI_GAUSS_5 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Line2D2ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsDiagonalLine2D2();
        array_1d<double, 3> coord(3);
        coord[0] = 2.0/3.0;
        coord[1] = 2.0/3.0;
        coord[2] = 0.0;
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), 1.0/6.0, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 5.0/6.0, TOLERANCE);
        auto& r_geom = *geom;
        auto p_geom_nodes = Kratos::make_shared<Line2D2<Node<3>>>(
        Kratos::make_intrusive<Node<3>>(1, r_geom[0].X(), r_geom[0].Y(), r_geom[0].Z()),
        Kratos::make_intrusive<Node<3>>(2, r_geom[1].X(), r_geom[1].Y(), r_geom[1].Z())
        );
        CrossCheckShapeFunctionsValues(*p_geom_nodes);
    }

    KRATOS_TEST_CASE_IN_SUITE(Line2D2ShapeFunctionsValuesMatrix, KratosCoreGeometriesFastSuite) {

        Geometry<Point>::Pointer p_geom = GeneratePointsDiagonalLine2D2();
        Line2D2<Point>::Pointer p_line = GeneratePointsDiagonalLine2D2();

        const Matrix N_values_geom = p_geom->ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const Matrix N_values_line = p_line->ShapeFunctionsValues(GeometryData::GI_GAUSS_2);

        KRATOS_CHECK_NEAR(N_values_geom(0, 0), 0.788675, TOLERANCE);
        KRATOS_CHECK_NEAR(N_values_geom(0, 1), 0.211325, TOLERANCE);
        KRATOS_CHECK_NEAR(N_values_geom(1, 0), 0.211325, TOLERANCE);
        KRATOS_CHECK_NEAR(N_values_geom(1, 1), 0.788675, TOLERANCE);

        KRATOS_CHECK_NEAR(N_values_line(0, 0), 0.788675, TOLERANCE);
        KRATOS_CHECK_NEAR(N_values_line(0, 1), 0.211325, TOLERANCE);
        KRATOS_CHECK_NEAR(N_values_line(1, 0), 0.211325, TOLERANCE);
        KRATOS_CHECK_NEAR(N_values_line(1, 1), 0.788675, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Line2D2ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsDiagonalLine2D2();
        auto& r_geom = *geom;
        auto p_geom_nodes = Kratos::make_shared<Line2D2<Node<3>>>(
        Kratos::make_intrusive<Node<3>>(1, r_geom[0].X(), r_geom[0].Y(), r_geom[0].Z()),
        Kratos::make_intrusive<Node<3>>(2, r_geom[1].X(), r_geom[1].Y(), r_geom[1].Z())
        );
        TestAllShapeFunctionsLocalGradients(*p_geom_nodes);
    }

    /**
     * Test an overlaping box and line (line has only one node in the box) HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2XIntersectionBoxSingleNodeInside, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        Point point_1(0.5, -0.1, 0.0);
        Point point_2(1.5, 0.1, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }
    KRATOS_TEST_CASE_IN_SUITE(Line2D2YIntersectionBoxSingleNodeInside, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitYDirectionLine2D2();
        Point point_1(-0.1, 0.5, 0.0);
        Point point_2(0.1, 1.5, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

    /**
     * Test an overlaping box and line (line has both nodes in the box) HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2XIntersectionBoxTwoNodesInside, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        Point point_1(-0.5, -0.1, 0.0);
        Point point_2(1.5, 0.1, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }
    KRATOS_TEST_CASE_IN_SUITE(Line2D2YIntersectionBoxTwoNodesInside, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitYDirectionLine2D2();
        Point point_1(-0.1, -0.5, 0.0);
        Point point_2(0.1, 1.5, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

    /**
     * Test an intersection with another line
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2IntersectionWithAnotherLine, KratosCoreGeometriesFastSuite) {
        auto geom_1 = GeneratePointsUnitXDirectionLine2D2();
        Point::Pointer point_1 = Kratos::make_shared<Point>(0.5, 0.5, 0.0);
        Point::Pointer point_2 = Kratos::make_shared<Point>(0.5, -0.5, 0.0);
        auto geom_2 = GenerateLine2D2WithPoints(point_1, point_2);
        KRATOS_CHECK(geom_1->HasIntersection(*geom_2));
    }

    /**
     * Test an intersection with another parallel line
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2IntersectionWithAnotherParallelLine, KratosCoreGeometriesFastSuite) {
        auto geom_1 = GeneratePointsUnitXDirectionLine2D2();
        Point::Pointer point_1 = Kratos::make_shared<Point>(0.0, 0.5, 0.0);
        Point::Pointer point_2 = Kratos::make_shared<Point>(0.5, 0.5, 0.0);
        auto geom_2 = GenerateLine2D2WithPoints(point_1, point_2);
        KRATOS_CHECK_IS_FALSE(geom_1->HasIntersection(*geom_2));
    }

    /**
     * Test a box inside a line HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2IntersectionBoxInsideX, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        Point point_1(0.25, -0.1, 0.0);
        Point point_2(0.75, 0.1, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

    /**
     * Test a box inside a line HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2IntersectionBoxInsideY, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitYDirectionLine2D2();
        Point point_1(-0.1,0.25, 0.0);
        Point point_2(0.1, 0.75, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }
    /**
     * Test a non overlaping box HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Line2D2NoIntersectionBox, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine2D2();
        Point point_1(1, 1, 0.0);
        Point point_2(2, 2, 0.0);
        KRATOS_CHECK_IS_FALSE(geom->HasIntersection(point_1, point_2));
    }

} // namespace Testing.
} // namespace Kratos.
