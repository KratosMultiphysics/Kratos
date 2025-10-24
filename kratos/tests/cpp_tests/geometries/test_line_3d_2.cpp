//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/line_3d_2.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos {
namespace Testing {

    /// Factory functions

    /** Generates a point type sample Line3D2N
    * @return  Pointer to a Line3D2N
    */
    Line3D2<Point>::Pointer GeneratePointsUnitXDirectionLine3D2() {
        return Kratos::make_shared<Line3D2<Point>>(
        Kratos::make_shared<Point>(0.0, 0.0, 0.0),
        Kratos::make_shared<Point>(1.0, 0.0, 0.0)
        );
    }

    /** Generates a point type sample Line3D2N
    * @return  Pointer to a Line3D2N
    */
    Line3D2<Point>::Pointer GeneratePointsUnitYDirectionLine3D2() {
        return Kratos::make_shared<Line3D2<Point>>(
        Kratos::make_shared<Point>(0.0, 0.0, 0.0),
        Kratos::make_shared<Point>(0.0, 1.0, 0.0)
        );
    }

    /** Generates a point type sample Line3D2N
    * @return  Pointer to a Line3D2N
    */
    Line3D2<Point>::Pointer GeneratePointsDiagonalLine3D2() {
        return Kratos::make_shared<Line3D2<Point>>(
        Kratos::make_shared<Point>(0.0, 0.0, 0.0),
        Kratos::make_shared<Point>(1.0, 1.0, 1.0)
        );
    }

    /** Generates a point type sample Line3D2N.
    * @return  Pointer to a Line3D2N
    */
    Line3D2<Point>::Pointer GenerateLine3D2WithPoints(Point::Pointer rPointOne, Point::Pointer rPointTwo ) {
        return Kratos::make_shared<Line3D2<Point>>(rPointOne, rPointTwo);
    }

    /** Checks if the number of edges is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2EdgesNumber, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        KRATOS_EXPECT_EQ(p_geom->EdgesNumber(), 1);
    }

    /** Checks if the edges are correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2Edges, KratosCoreGeometriesFastSuite) {
        auto p_geom = GeneratePointsUnitXDirectionLine3D2();

        const auto& r_edges = p_geom->GenerateEdges();
        KRATOS_EXPECT_NEAR((r_edges[0])[0].X(), (p_geom->pGetPoint(0))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[0])[0].Y(), (p_geom->pGetPoint(0))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[0])[1].X(), (p_geom->pGetPoint(1))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[0])[1].Y(), (p_geom->pGetPoint(1))->Y(), TOLERANCE);
    }

    /** Checks if the number of faces is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2FacesNumber, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        KRATOS_EXPECT_EQ(p_geom->FacesNumber(), 0);
    }

    /** Checks if the length of the line is calculated correctly.
    * Checks if the length of the line is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(LengthLine3D2, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsDiagonalLine3D2();

        KRATOS_EXPECT_NEAR(p_geom->Length(), std::sqrt(3.0), TOLERANCE);
    }

    /** Checks the inside test for a given point respect to the line
    * Checks the inside test for a given point respect to the line
    * It performs 4 tests:
    * A Point inside the line: Expected result TRUE
    * A Point outside the line: Expected result FALSE
    * A Point over a vertex 1 of the line: Expected result TRUE
    * A Point over an vertex 2 of the line: Expected result TRUE
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2IsInside, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsDiagonalLine3D2();

        Point PointInside(0.5, 0.5, 0.5);
        Point PointOutside(1.66, 0.66, 0.0);
        Point PointInVertex(0.0, 0.0, 0.0);
        Point PointInEdge(1.0, 1.0, 1.0);

        Point LocalCoords;

        KRATOS_EXPECT_TRUE(p_geom->IsInside(PointInside, LocalCoords, EPSILON));
        KRATOS_EXPECT_FALSE(p_geom->IsInside(PointOutside, LocalCoords, EPSILON));
        KRATOS_EXPECT_TRUE(p_geom->IsInside(PointInVertex, LocalCoords, EPSILON));
        KRATOS_EXPECT_TRUE(p_geom->IsInside(PointInEdge, LocalCoords, EPSILON));
    }

    /** Checks the point local coordinates for a given point respect to the
    * line. The baricentre of the triangle is selected due to its known
    * solution.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsDiagonalLine3D2();

        // Compute the global coordinates of the baricentre
        const Geometry<Point>::PointsArrayType geom_pts = p_geom->Points();
        Point centre = Point{geom_pts[0] + geom_pts[1]};
        centre *= 1.0/2.0;

        // Compute the centre local coordinates
        array_1d<double, 3> centre_local_coords;
        p_geom->PointLocalCoordinates(centre_local_coords, centre);

        KRATOS_EXPECT_NEAR(centre_local_coords(0), 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(centre_local_coords(1), 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(centre_local_coords(2), 0.0, TOLERANCE);
    }

        /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        const double ExpectedJacobian = 0.5;

        Vector JacobianDeterminants;
        p_geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_1 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i) {
            KRATOS_EXPECT_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianArray2, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        const double ExpectedJacobian = 0.5;

        Vector JacobianDeterminants;
        p_geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_2 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i) {
            KRATOS_EXPECT_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianArray3, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        const double ExpectedJacobian = 0.5;

        Vector JacobianDeterminants;
        p_geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_3 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i) {
            KRATOS_EXPECT_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianArray4, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
        const double ExpectedJacobian = 0.5;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_4 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i) {
            KRATOS_EXPECT_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianArray5, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        const double ExpectedJacobian = 0.5;

        Vector JacobianDeterminants;
        p_geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_5 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i) {
            KRATOS_EXPECT_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        const double ExpectedJacobian = 0.5;

        double JacobianDeterminant = p_geom->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_1 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianIndex2, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 0.5;

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_2 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_2 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianIndex3, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 0.5;

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_3 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_3 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 3, GeometryData::IntegrationMethod::GI_GAUSS_3 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianIndex4, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 0.5;

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_4 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_4 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 3, GeometryData::IntegrationMethod::GI_GAUSS_4 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 4, GeometryData::IntegrationMethod::GI_GAUSS_4 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianIndex5, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsUnitXDirectionLine3D2();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 0.5;

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_5 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_5 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 3, GeometryData::IntegrationMethod::GI_GAUSS_5 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 4, GeometryData::IntegrationMethod::GI_GAUSS_5 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = p_geom->DeterminantOfJacobian( 5, GeometryData::IntegrationMethod::GI_GAUSS_5 );
        KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Line3D2ShapeFunctionValue, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsDiagonalLine3D2();
        array_1d<double, 3> coord(3);
        coord[0] = 2.0/3.0;
        coord[1] = 2.0/3.0;
        coord[2] = 2.0/3.0;
        KRATOS_EXPECT_NEAR(p_geom->ShapeFunctionValue(0, coord), 1.0/6.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(p_geom->ShapeFunctionValue(1, coord), 5.0/6.0, TOLERANCE);
        auto& r_geom = *p_geom;
        Geometry<Node>::Pointer p_geom_nodes = Kratos::make_shared<Line3D2<Node>>(
        Kratos::make_intrusive<Node>(1, r_geom[0].X(), r_geom[0].Y(), r_geom[0].Z()),
        Kratos::make_intrusive<Node>(2, r_geom[1].X(), r_geom[1].Y(), r_geom[1].Z())
        );
        CrossCheckShapeFunctionsValues(*p_geom_nodes);
    }

    KRATOS_TEST_CASE_IN_SUITE(Line3D2ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsDiagonalLine3D2();
        array_1d<double, 3> coord(3);
        coord[0] = 2.0/3.0;
        coord[1] = 2.0/3.0;
        coord[2] = 2.0/3.0;
        Vector N_values;
        p_geom->ShapeFunctionsValues(N_values, coord);
        KRATOS_EXPECT_NEAR(N_values[0], 1.0/6.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values[1], 5.0/6.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Line3D2ShapeFunctionsValuesMatrix, KratosCoreGeometriesFastSuite) {

        Geometry<Point>::Pointer p_geom = GeneratePointsDiagonalLine3D2();
        Line3D2<Point>::Pointer p_line = GeneratePointsDiagonalLine3D2();

        const Matrix N_values_geom = p_geom->ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);
        const Matrix N_values_line = p_line->ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

        KRATOS_EXPECT_NEAR(N_values_geom(0, 0), 0.788675, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_geom(0, 1), 0.211325, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_geom(1, 0), 0.211325, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_geom(1, 1), 0.788675, TOLERANCE);

        KRATOS_EXPECT_NEAR(N_values_line(0, 0), 0.788675, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_line(0, 1), 0.211325, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_line(1, 0), 0.211325, TOLERANCE);
        KRATOS_EXPECT_NEAR(N_values_line(1, 1), 0.788675, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Line3D2ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        Geometry<Point>::Pointer p_geom = GeneratePointsDiagonalLine3D2();
        auto& r_geom = *p_geom;
        Geometry<Node>::Pointer p_geom_nodes = Kratos::make_shared<Line3D2<Node>>(
        Kratos::make_intrusive<Node>(1, r_geom[0].X(), r_geom[0].Y(), r_geom[0].Z()),
        Kratos::make_intrusive<Node>(2, r_geom[1].X(), r_geom[1].Y(), r_geom[1].Z())
        );
        TestAllShapeFunctionsLocalGradients(*p_geom_nodes);
    }

} // namespace Testing.
} // namespace Kratos.
