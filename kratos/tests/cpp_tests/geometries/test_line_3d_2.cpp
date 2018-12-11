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
    * Checks if the number of edges is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2EdgesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
        KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 2);
    }

    /** Checks if the number of faces is correct.
    * Checks if the number of faces is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2FacesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
        KRATOS_CHECK_EQUAL(geom->FacesNumber(), 2);
    }

    /** Checks if the length of the line is calculated correctly.
    * Checks if the length of the line is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(LengthLine3D2, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsDiagonalLine3D2();

        KRATOS_CHECK_NEAR(geom->Length(), std::sqrt(3.0), TOLERANCE);
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
        auto geom = GeneratePointsDiagonalLine3D2();

        Point PointInside(0.5, 0.5, 0.5);
        Point PointOutside(1.66, 0.66, 0.0);
        Point PointInVertex(0.0, 0.0, 0.0);
        Point PointInEdge(1.0, 1.0, 1.0);

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
    KRATOS_TEST_CASE_IN_SUITE(Line3D2PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsDiagonalLine3D2();

        // Compute the global coordinates of the baricentre
        const Geometry<Point>::PointsArrayType geom_pts = geom->Points();
        Point centre = geom_pts[0] + geom_pts[1];
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
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
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
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianArray2, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
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
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianArray3, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
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
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianArray4, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
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
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianArray5, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
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
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
        const double ExpectedJacobian = 0.5;

        double JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_1 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianIndex2, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
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
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianIndex3, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
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
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianIndex4, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
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
    KRATOS_TEST_CASE_IN_SUITE(Line3D2DeterminantOfJacobianIndex5, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsUnitXDirectionLine3D2();
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

    KRATOS_TEST_CASE_IN_SUITE(Line3D2ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsDiagonalLine3D2();
        array_1d<double, 3> coord(3);
        coord[0] = 2.0/3.0;
        coord[1] = 2.0/3.0;
        coord[2] = 2.0/3.0;
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), 1.0/6.0, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 5.0/6.0, TOLERANCE);
        auto& r_geom = *geom;
        auto p_geom_nodes = Kratos::make_shared<Line3D2<Node<3>>>(
        Kratos::make_shared<Node<3>>(1, r_geom[0].X(), r_geom[0].Y(), r_geom[0].Z()),
        Kratos::make_shared<Node<3>>(2, r_geom[1].X(), r_geom[1].Y(), r_geom[1].Z())
        );
        CrossCheckShapeFunctionsValues(*p_geom_nodes);
    }

    KRATOS_TEST_CASE_IN_SUITE(Line3D2ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsDiagonalLine3D2();
        auto& r_geom = *geom;
        auto p_geom_nodes = Kratos::make_shared<Line3D2<Node<3>>>(
        Kratos::make_shared<Node<3>>(1, r_geom[0].X(), r_geom[0].Y(), r_geom[0].Z()),
        Kratos::make_shared<Node<3>>(2, r_geom[1].X(), r_geom[1].Y(), r_geom[1].Z())
        );
        TestAllShapeFunctionsLocalGradients(*p_geom_nodes);
    }

} // namespace Testing.
} // namespace Kratos.
