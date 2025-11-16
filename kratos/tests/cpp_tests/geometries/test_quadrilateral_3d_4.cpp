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
#include "geometries/quadrilateral_3d_4.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos::Testing
{
    /// Factory functions

    /** Generates a sample Quadrilateral3D4.
    * Generates a quadrilateral defined by three random points in the space.
    * @return  Pointer to a Quadrilateral3D4
    */
    template<class TPointType>
    typename Quadrilateral3D4<TPointType>::Pointer GenerateQuadrilateral3D4(
        typename TPointType::Pointer PointA = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointB = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointC = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointD = GeneratePoint<TPointType>()
        )
    {
        return typename Quadrilateral3D4<TPointType>::Pointer(new Quadrilateral3D4<TPointType>(
            PointA,
            PointB,
            PointC,
            PointD
            ));
    }

    /** Generates a sample Quadrilateral3D4.
    * Generates a right quadrilateral with origin in the origin and leg size 1.
    * @return  Pointer to a Quadrilateral3D4
    */
    template<class TPointType>
    typename Quadrilateral3D4<TPointType>::Pointer GenerateRightQuadrilateral3D4()
    {
        return typename Quadrilateral3D4<TPointType>::Pointer(new Quadrilateral3D4<TPointType>(
        GeneratePoint<TPointType>(0.0, 0.0, 0.0),
        GeneratePoint<TPointType>(std::cos(Globals::Pi/4), 0.0, std::sin(Globals::Pi/4)),
        GeneratePoint<TPointType>(1.0, 1.0, 0.5),
        GeneratePoint<TPointType>(0.0, 1.0, 0.0)
        ));
    }

    /** Generates a sample Quadrilateral3D4.
    * Generates a right quadrilateral with origin in the origin and leg size 1.
    * @return  Pointer to a Quadrilateral3D4
    */
    template<class TPointType>
    typename Quadrilateral3D4<TPointType>::Pointer GenerateFlatQuadrilateral3D4()
    {
        return typename Quadrilateral3D4<TPointType>::Pointer(new Quadrilateral3D4<TPointType>(
        GeneratePoint<TPointType>( 0.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 0.0, 0.0),
        GeneratePoint<TPointType>( 1.0, 1.0, 0.0),
        GeneratePoint<TPointType>( 0.0, 1.0, 0.0)
        ));
    }

    /// Tests

    /** Checks if the number of edges is correct.
    * Checks if the number of edges is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4EdgesNumber, KratosCoreGeometriesFastSuite)
    {
        auto geom = GenerateRightQuadrilateral3D4<Node>();

        KRATOS_EXPECT_EQ(geom->EdgesNumber(), 4);
    }

    /** Checks if the number of edges is correct.
    * Checks if the number of edges is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4Edges, KratosCoreGeometriesFastSuite)
    {
        auto p_geom = GenerateRightQuadrilateral3D4<Node>();

        const auto& r_edges = p_geom->GenerateEdges();
        KRATOS_EXPECT_NEAR((r_edges[0])[0].X(), (p_geom->pGetPoint(0))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[0])[0].Y(), (p_geom->pGetPoint(0))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[0])[0].Z(), (p_geom->pGetPoint(0))->Z(), TOLERANCE);

        KRATOS_EXPECT_NEAR((r_edges[0])[1].X(), (p_geom->pGetPoint(1))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[0])[1].Y(), (p_geom->pGetPoint(1))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[0])[1].Z(), (p_geom->pGetPoint(1))->Z(), TOLERANCE);

        KRATOS_EXPECT_NEAR((r_edges[1])[0].X(), (p_geom->pGetPoint(1))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[1])[0].Y(), (p_geom->pGetPoint(1))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[1])[0].Z(), (p_geom->pGetPoint(1))->Z(), TOLERANCE);

        KRATOS_EXPECT_NEAR((r_edges[1])[1].X(), (p_geom->pGetPoint(2))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[1])[1].Y(), (p_geom->pGetPoint(2))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[1])[1].Z(), (p_geom->pGetPoint(2))->Z(), TOLERANCE);

        KRATOS_EXPECT_NEAR((r_edges[2])[0].X(), (p_geom->pGetPoint(2))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[2])[0].Y(), (p_geom->pGetPoint(2))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[2])[0].Z(), (p_geom->pGetPoint(2))->Z(), TOLERANCE);

        KRATOS_EXPECT_NEAR((r_edges[2])[1].X(), (p_geom->pGetPoint(3))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[2])[1].Y(), (p_geom->pGetPoint(3))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[2])[1].Z(), (p_geom->pGetPoint(3))->Z(), TOLERANCE);

        KRATOS_EXPECT_NEAR((r_edges[3])[0].X(), (p_geom->pGetPoint(3))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[3])[0].Y(), (p_geom->pGetPoint(3))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[3])[0].Z(), (p_geom->pGetPoint(3))->Z(), TOLERANCE);

        KRATOS_EXPECT_NEAR((r_edges[3])[1].X(), (p_geom->pGetPoint(0))->X(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[3])[1].Y(), (p_geom->pGetPoint(0))->Y(), TOLERANCE);
        KRATOS_EXPECT_NEAR((r_edges[3])[1].Z(), (p_geom->pGetPoint(0))->Z(), TOLERANCE);
    }

    /** Checks if the number of faces is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4FacesNumber, KratosCoreGeometriesFastSuite)
    {
        auto geom = GenerateRightQuadrilateral3D4<Node>();

        KRATOS_EXPECT_EQ(geom->FacesNumber(), 1);
    }

    /** Checks if the faces are correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4Faces, KratosCoreGeometriesFastSuite) {
        auto p_geom = GenerateRightQuadrilateral3D4<Node>();

        const auto& r_faces = p_geom->GenerateFaces();
        ASSERT_EQ(r_faces.size(), 1);
        for (std::size_t i = 0; i < r_faces.front().PointsNumber(); ++i) {
            KRATOS_EXPECT_NEAR(r_faces.front()[i].X(), (p_geom->pGetPoint(i))->X(), TOLERANCE);
            KRATOS_EXPECT_NEAR(r_faces.front()[i].Y(), (p_geom->pGetPoint(i))->Y(), TOLERANCE);
            KRATOS_EXPECT_NEAR(r_faces.front()[i].Z(), (p_geom->pGetPoint(i))->Z(), TOLERANCE);
        }
    }

    /** Checks if the area of the quadrilateral is calculated correctly.
    * Checks if the area of the quadrilateral is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4Area, KratosCoreGeometriesFastSuite)
    {
        auto geom = GenerateRightQuadrilateral3D4<Node>();
        KRATOS_EXPECT_NEAR(geom->Area(), 1.06948, 1.0e-5);
    }

    /** Tests the PointLocalCoordinates for Quadrilateral2D4.
     * Tests the PointLocalCoordinates for Quadrilateral2D4.
     */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateFlatQuadrilateral3D4<Node>();

        Point TestPointA(1.0, 1.0, 0.0);
        Point TestPointB(0.5, 0.5, 0.0);
        Point TestResultA(0.0, 0.0, 0.0);
        Point TestResultB(0.0, 0.0, 0.0);

        geom->PointLocalCoordinates(TestResultA, TestPointA);
        geom->PointLocalCoordinates(TestResultB, TestPointB);

        // Test transformation in the edge
        KRATOS_EXPECT_NEAR(TestResultA[0], 1.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(TestResultA[1], 1.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(TestResultA[2], 0.0, TOLERANCE);

        // Test transformation in the center
        KRATOS_EXPECT_NEAR(TestResultB[0], 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(TestResultB[1], 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(TestResultB[2], 0.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateFlatQuadrilateral3D4<Node>();
      array_1d<double, 3> coord(3);
      coord[0] = 1.0 / 2.0;
      coord[1] = 1.0 / 4.0;
      coord[2] = 0.0;
      KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(0, coord), 0.09375, TOLERANCE);
      KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(1, coord), 0.28125, TOLERANCE);
      KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(2, coord), 0.46875, TOLERANCE);
      KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(3, coord), 0.15625, TOLERANCE);
      CrossCheckShapeFunctionsValues(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateFlatQuadrilateral3D4<Node>();
        TestAllShapeFunctionsLocalGradients(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4CoplanarPointIntersection, KratosCoreGeometriesFastSuite) {
        Quadrilateral3D4<Point > quadrilateral_1(
            std::make_shared<Point>(0.0, 0.0, 0.0),
            std::make_shared<Point>(10., 0.0, 2.0),
            std::make_shared<Point>(0.0, 1.0, 0.0),
            std::make_shared<Point>(0.0, 1.0, 2.0)
            );
        Quadrilateral3D4<Point > quadrilateral_2(
            std::make_shared<Point>(0.00, 0.00, 0.0),
            std::make_shared<Point>(-10., 0.0, -2.0),
            std::make_shared<Point>(0.0, -1.0, 0.00),
            std::make_shared<Point>(0.0, -1.0, -2.00)
            );

        KRATOS_EXPECT_TRUE(quadrilateral_1.HasIntersection(quadrilateral_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4EdgeIntersection, KratosCoreGeometriesFastSuite) {
        Quadrilateral3D4<Point > quadrilateral_1(
            std::make_shared<Point>(0.0, 0.0, 0.0),
            std::make_shared<Point>(10., 0.0, 2.0),
            std::make_shared<Point>(0.0, 1.0, 0.0),
            std::make_shared<Point>(10.0, 1.0, 0.0)
            );
        Quadrilateral3D4<Point > quadrilateral_2(
            std::make_shared<Point>(0.00, 0.00, 0.0),
            std::make_shared<Point>(10., 0.0, 2.0),
            std::make_shared<Point>(0.0, -1.0, 0.00),
            std::make_shared<Point>(0.0, -1.0, 2.00)
            );

        KRATOS_EXPECT_TRUE(quadrilateral_1.HasIntersection(quadrilateral_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4InsideIntersection, KratosCoreGeometriesFastSuite) {
        Quadrilateral3D4<Point > quadrilateral_1(
            std::make_shared<Point>(0.0, 0.0, 0.0),
            std::make_shared<Point>(0.0, 0.0, 4.0),
            std::make_shared<Point>(0.0, 4.0, 0.0),
            std::make_shared<Point>(0.0, 4.0, 4.0)
            );
        Quadrilateral3D4<Point > quadrilateral_2(
            std::make_shared<Point>(0.0, 1.0, 1.0),
            std::make_shared<Point>(0.0, 1.0, 3.0),
            std::make_shared<Point>(0.0, 3.0, 1.0),
            std::make_shared<Point>(0.0, 3.0, 3.0)
            );

        KRATOS_EXPECT_TRUE(quadrilateral_1.HasIntersection(quadrilateral_2));
    }

    /** Checks the ProjectionPoint test for a given point respect to the quadrilateral
    * Checks the ProjectionPoint test for a given point respect to the quadrilateral
    */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4ProjectionPoint, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateFlatQuadrilateral3D4<NodeType>();

        Point point(0.25,0.35,0.2);

        Geometry<Point>::CoordinatesArrayType global_coords;
        Geometry<Point>::CoordinatesArrayType local_coords;

        geom->ProjectionPointGlobalToLocalSpace(point.Coordinates(), local_coords);
        geom->GlobalCoordinates(global_coords, local_coords);

        KRATOS_EXPECT_RELATIVE_NEAR(global_coords[0], 0.25, 1.0e-4);
        KRATOS_EXPECT_RELATIVE_NEAR(global_coords[1], 0.35, 1.0e-4);
        KRATOS_EXPECT_NEAR(global_coords[2], 0.0, 1.0e-4);

        KRATOS_EXPECT_RELATIVE_NEAR(local_coords[0], -0.5, 1.0e-4);
        KRATOS_EXPECT_RELATIVE_NEAR(local_coords[1], -0.3, 1.0e-4);
        KRATOS_EXPECT_NEAR(local_coords[2], 0.0, 1.0e-4);
    }

    /**
     * Checks the distance from a point to a quadrilateral
     */
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4CalculateDistance, KratosCoreGeometriesFastSuite)
    {
        auto geom = GenerateFlatQuadrilateral3D4<Node>();

        Point point1(0.0, 0.0, 0.0);
        KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point1), 0.0);

        Point point2(0.0, 0.0, 0.5);
        KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point2), 0.5);
    }

//     /** Checks if the volume of the quadrilateral is calculated correctly.
//     * Checks if the volume of the quadrilateral is calculated correctly.
//     * For quadrilateral 2D3 'volume()' call defaults to 'area()'
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4Volume, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//
//         KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Volume(), "Calling base class 'Volume' method instead of derived class one.");
//     }
//
//     /** Checks the inside test for a given point respect to the quadrilateral
//     * Checks the inside test for a given point respect to the quadrilateral
//     * It performs 4 tests:
//     * A Point inside the quadrilateral: Expected result TRUE
//     * A Point outside the quadrilateral: Expected result FALSE
//     * A Point over a vertex of the quadrilateral: Expected result TRUE
//     * A Point over an edge of the quadrilateral: Expected result TRUE
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4IsInside, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//
//         Point<3> PointInside(1.0/3.0, 2.0/3.0, 1.0/6.0);
//         Point<3> PointOutside(2.0/3.0, 2.0/3.0, 0.0);
//         Point<3> PointInVertex(0.0, 0.0, 0.0);
//         Point<3> PointInEdge(0.5, 0.5, 0.0);
//
//         Point<3> LocalCoords;
//
//         // It appears that the function checks whether the PROJECTION of the point is inside the geometry.
//         KRATOS_EXPECT_TRUE(geom->IsInside(PointInside, LocalCoords, EPSILON));
//         KRATOS_EXPECT_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
//         KRATOS_EXPECT_TRUE(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
//         KRATOS_EXPECT_TRUE(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
//     }
//
//     /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//         const double ExpectedJacobian = 1.0;
//
//         Vector JacobianDeterminants;
//         geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_1 );
//
//         for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
//         {
//             KRATOS_EXPECT_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
//         }
//     }
//
//     /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray2, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//         const double ExpectedJacobian = 1.0;
//
//         Vector JacobianDeterminants;
//         geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_2 );
//
//         for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
//         {
//             KRATOS_EXPECT_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
//         }
//     }
//
//     /** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray3, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//         const double ExpectedJacobian = 1.0;
//
//         Vector JacobianDeterminants;
//         geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_3 );
//
//         for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
//         {
//             KRATOS_EXPECT_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
//         }
//     }
//
//     /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray4, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//         const double ExpectedJacobian = 1.0;
//
//         Vector JacobianDeterminants;
//         geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_4 );
//
//         for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
//         {
//             KRATOS_EXPECT_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
//         }
//     }
//
//     /** Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianArray5, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//         const double ExpectedJacobian = 1.0;
//
//         Vector JacobianDeterminants;
//         geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_5 );
//
//         for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
//         {
//             KRATOS_EXPECT_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
//         }
//     }
//
//     /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//         const double ExpectedJacobian = 1.0;
//
//         double JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_1 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }
//
//     /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex2, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//         double JacobianDeterminant = 0.0;
//         const double ExpectedJacobian = 1.0;
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_2 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_2 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }
//
//     /** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex3, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//         double JacobianDeterminant = 0.0;
//         const double ExpectedJacobian = 1.0;
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_3 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_3 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 3, GeometryData::IntegrationMethod::GI_GAUSS_3 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }
//
//     /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex4, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//         double JacobianDeterminant = 0.0;
//         const double ExpectedJacobian = 1.0;
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_4 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_4 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 3, GeometryData::IntegrationMethod::GI_GAUSS_4 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 4, GeometryData::IntegrationMethod::GI_GAUSS_4 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }
//
//     /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
//     */
//     KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D4DeterminantOfJacobianIndex5, KratosCoreGeometriesFastSuite)
//     {
//         auto geom = GenerateRightQuadrilateral3D4<Node>();
//         double JacobianDeterminant = 0.0;
//         const double ExpectedJacobian = 1.0;
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::IntegrationMethod::GI_GAUSS_5 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::IntegrationMethod::GI_GAUSS_5 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 3, GeometryData::IntegrationMethod::GI_GAUSS_5 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 4, GeometryData::IntegrationMethod::GI_GAUSS_5 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//
//         JacobianDeterminant = geom->DeterminantOfJacobian( 5, GeometryData::IntegrationMethod::GI_GAUSS_5 );
//         KRATOS_EXPECT_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
//     }

} // namespace Kratos::Testing.
