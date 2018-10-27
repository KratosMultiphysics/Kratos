//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/triangle_3d_3.h"
#include "tests/geometries/test_geometry.h"
#include "tests/geometries/test_shape_function_derivatives.h"
#include "tests/geometries/cross_check_shape_functions_values.h"

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos 
{
namespace Testing 
{
    typedef Node<3> NodeType;

    /// Factory functions

    /** Generates a sample Triangle3D3.
    * Generates a triangle defined by three random points in the space.
    * @return  Pointer to a Triangle3D3
    */
    template<class TPointType>
    typename Triangle3D3<TPointType>::Pointer GenerateTriangle3D3(
        typename TPointType::Pointer PointA = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointB = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointC = GeneratePoint<TPointType>()) {
    return typename Triangle3D3<TPointType>::Pointer(new Triangle3D3<TPointType>(
        PointA,
        PointB,
        PointC
    ));
    }

    /** Generates a sample Triangle3D3.
    * Generates a right triangle with origin in the origin and leg size 1.
    * @return  Pointer to a Triangle3D3
    */
    template<class TPointType>
    typename Triangle3D3<TPointType>::Pointer GenerateRightTriangle3D3() {
    return typename Triangle3D3<TPointType>::Pointer(new Triangle3D3<TPointType>(
        GeneratePoint<TPointType>(0.0, 0.0, 0.0),
        GeneratePoint<TPointType>(std::cos(Globals::Pi/4), 0.0, std::sin(Globals::Pi/4)),
        GeneratePoint<TPointType>(0.0, 1.0, 0.0)
    ));
    }

    /** Generates a sample Triangle3D3.
    * Generates an equilateral triangle with vertices at each axis.
    * @return  Pointer to a Triangle3D3
    */
    template<class TPointType>
    typename Triangle3D3<TPointType>::Pointer GenerateEquilateralTriangle3D3() {
    return typename Triangle3D3<TPointType>::Pointer(new Triangle3D3<TPointType>(
        GeneratePoint<TPointType>(1.0, 0.0, 0.0),
        GeneratePoint<TPointType>(0.0, 1.0, 0.0),
        GeneratePoint<TPointType>(0.0, 0.0, 1.0)
    ));
    }

    /// Tests

    /** Checks if the number of edges is correct.
    * Checks if the number of edges is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3EdgesNumber, KratosCoreFastSuite) {
    auto geom = GenerateRightTriangle3D3<NodeType>();

    KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 3);
    }

    /** Checks if the number of faces is correct.
    * Checks if the number of faces is correct.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3FacesNumber, KratosCoreFastSuite) {
    auto geom = GenerateRightTriangle3D3<NodeType>();

    // Charlie: I will let this to 3 but probably 'FacesNumber' needs to be documented to state
    // that for planar geometries it also return the number of edges.
    KRATOS_CHECK_EQUAL(geom->FacesNumber(), 3);
    }

    /** Checks if the area of the triangle is calculated correctly.
    * Checks if the area of the triangle is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3Area, KratosCoreFastSuite) {
    auto geom = GenerateRightTriangle3D3<NodeType>();

    KRATOS_CHECK_NEAR(geom->Area(), 0.5, TOLERANCE);
    }

    /** Checks if the volume of the triangle is calculated correctly.
    * Checks if the volume of the triangle is calculated correctly.
    * For triangle 2D3 'volume()' call defaults to 'area()'
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3Volume, KratosCoreFastSuite) {
    auto geom = GenerateRightTriangle3D3<NodeType>();

    KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Volume(), "Calling base class 'Volume' method instead of derived class one.");
    }

    /** Checks if the minimum edge length is calculated correctly.
    * Checks if the minimum edge length is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3MinEdgeLength, KratosCoreFastSuite) {
    auto geom = GenerateRightTriangle3D3<NodeType>();

    KRATOS_CHECK_NEAR(geom->MinEdgeLength(), 1.0, TOLERANCE);
    }

    /** Checks if the maximum edge length is calculated correctly.
    * Checks if the maximum edge length is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3MaxEdgeLength, KratosCoreFastSuite) {
    auto geom = GenerateRightTriangle3D3<NodeType>();

    KRATOS_CHECK_NEAR(geom->MaxEdgeLength(), 1.414213, TOLERANCE);
    }

    /** Checks if the average edge length is calculated correctly.
    * Checks if the average edge length is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3AverageEdgeLength, KratosCoreFastSuite) {
    auto geom = GenerateRightTriangle3D3<NodeType>();

    KRATOS_CHECK_NEAR(geom->AverageEdgeLength(), 1.138071, TOLERANCE);
    }

    /** Checks if the circumradius is calculated correctly.
    * Checks if the circumradius is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3Circumradius, KratosCoreFastSuite) {
    auto geom = GenerateRightTriangle3D3<NodeType>();

    KRATOS_CHECK_NEAR(geom->Circumradius(), 0.707107, TOLERANCE);
    }

    /** Checks if the inradius is calculated correctly.
    * Checks if the inradius is calculated correctly.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3Inradius, KratosCoreFastSuite) {
    auto geom = GenerateRightTriangle3D3<NodeType>();

    KRATOS_CHECK_NEAR(geom->Inradius(), 0.292893, TOLERANCE);
    }

    /** Checks the inside test for a given point respect to the triangle
    * Checks the inside test for a given point respect to the triangle
    * It performs 4 tests:
    * A Point inside the triangle: Expected result TRUE
    * A Point outside the triangle: Expected result FALSE
    * A Point over a vertex of the triangle: Expected result TRUE
    * A Point over an edge of the triangle: Expected result TRUE
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3IsInside, KratosCoreFastSuite) {
    auto geom = GenerateRightTriangle3D3<NodeType>();

    Point PointInside(0.33, 0.33, 0.0);
    Point PointOutside(0.66, 0.66, 0.0);
    Point PointInVertex(0.0, 0.0, 0.0);
    Point PointInEdge(0.5, 0.5, 0.0);

    Point LocalCoords;

    // It appears that the function checks whether the PROJECTION of the point is inside the geometry.
    KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
    KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3DeterminantOfJacobianArray1, KratosCoreFastSuite) {
        auto geom = GenerateRightTriangle3D3<NodeType>();
        const double ExpectedJacobian = 1.0;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_1 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
        {
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3DeterminantOfJacobianArray2, KratosCoreFastSuite) {
        auto geom = GenerateRightTriangle3D3<NodeType>();
        const double ExpectedJacobian = 1.0;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_2 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
        {
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3DeterminantOfJacobianArray3, KratosCoreFastSuite) {
        auto geom = GenerateRightTriangle3D3<NodeType>();
        const double ExpectedJacobian = 1.0;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_3 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
        {
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_4' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3DeterminantOfJacobianArray4, KratosCoreFastSuite) {
        auto geom = GenerateRightTriangle3D3<NodeType>();
        const double ExpectedJacobian = 1.0;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_4 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
        {
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_5' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3DeterminantOfJacobianArray5, KratosCoreFastSuite) {
        auto geom = GenerateRightTriangle3D3<NodeType>();
        const double ExpectedJacobian = 1.0;

        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_5 );

        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
        {
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
        }
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
        * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
        */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3DeterminantOfJacobianIndex1, KratosCoreFastSuite) {
        auto geom = GenerateRightTriangle3D3<NodeType>();
        const double ExpectedJacobian = 1.0;

        double JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_1 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3DeterminantOfJacobianIndex2, KratosCoreFastSuite) {
        auto geom = GenerateRightTriangle3D3<NodeType>();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 1.0;

        JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_2 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

        JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::GI_GAUSS_2 );
        KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
    }

    /** Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    * Tests the Jacobian determinants using 'GI_GAUSS_3' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3DeterminantOfJacobianIndex3, KratosCoreFastSuite) {
        auto geom = GenerateRightTriangle3D3<NodeType>();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 1.0;

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
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3DeterminantOfJacobianIndex4, KratosCoreFastSuite) {
        auto geom = GenerateRightTriangle3D3<NodeType>();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 1.0;

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
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3DeterminantOfJacobianIndex5, KratosCoreFastSuite) {
        auto geom = GenerateRightTriangle3D3<NodeType>();
        double JacobianDeterminant = 0.0;
        const double ExpectedJacobian = 1.0;

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


    /** Tests two very near parallel triangles HasIntegration which should give false
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3ParallelNoIntersection, KratosCoreFastSuite) {
        Triangle3D3<Point > triangle_1(
            GeneratePoint<NodeType >(0.0, 0.0, 0.0),
            GeneratePoint<NodeType >(10., 0.0, 2.0),
            GeneratePoint<NodeType >(0.0, 1.0, 0.0)
            );
        Triangle3D3<Point > triangle_2(
            GeneratePoint<NodeType >(0.0, 0.0, 0.01),
            GeneratePoint<NodeType >(10., 0.0, 2.01),
            GeneratePoint<NodeType >(0.0, 1.0, 0.01)
            );

        KRATOS_CHECK_IS_FALSE(triangle_1.HasIntersection(triangle_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3ParallelNearIntersection, KratosCoreFastSuite) {
        Triangle3D3<Point > triangle_1(
            GeneratePoint<NodeType >(0.0, 0.0, 0.0),
            GeneratePoint<NodeType >(10., 0.0, 2.0),
            GeneratePoint<NodeType >(0.0, 1.0, 0.0)
            );
        Triangle3D3<Point > triangle_2(
            GeneratePoint<NodeType >(0.0, 0.0, 0.00000001),
            GeneratePoint<NodeType >(10., 0.0, 2.00000001),
            GeneratePoint<NodeType >(0.0, 1.0, 0.00000001)
            );

        KRATOS_CHECK_IS_FALSE(triangle_1.HasIntersection(triangle_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3CoplanarNoIntersection, KratosCoreFastSuite) {
        Triangle3D3<Point > triangle_1(
            GeneratePoint<NodeType >(0.0, 0.0, 0.0),
            GeneratePoint<NodeType >(10., 0.0, 2.0),
            GeneratePoint<NodeType >(0.0, 1.0, 0.0)
            );
        Triangle3D3<Point > triangle_2(
            GeneratePoint<NodeType >(0.00000001, 0.00000001, 0.00000001),
            GeneratePoint<NodeType >(-10., 0.0, -2.0),
            GeneratePoint<NodeType >(0.0, -1.0, 0.00)
            );

        KRATOS_CHECK_IS_FALSE(triangle_1.HasIntersection(triangle_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3CoplanarPointIntersection, KratosCoreFastSuite) {
        Triangle3D3<Point > triangle_1(
            GeneratePoint<NodeType >(0.0, 0.0, 0.0),
            GeneratePoint<NodeType >(10., 0.0, 2.0),
            GeneratePoint<NodeType >(0.0, 1.0, 0.0)
            );
        Triangle3D3<Point > triangle_2(
            GeneratePoint<NodeType >(0.00, 0.00, 0.0),
            GeneratePoint<NodeType >(-10., 0.0, -2.0),
            GeneratePoint<NodeType >(0.0, -1.0, 0.00)
            );

        KRATOS_CHECK(triangle_1.HasIntersection(triangle_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3EdgeIntersection, KratosCoreFastSuite) {
        Triangle3D3<Point > triangle_1(
            GeneratePoint<NodeType >(0.0, 0.0, 0.0),
            GeneratePoint<NodeType >(10., 0.0, 2.0),
            GeneratePoint<NodeType >(0.0, 1.0, 0.0)
            );
        Triangle3D3<Point > triangle_2(
            GeneratePoint<NodeType >(0.00, 0.00, 0.0),
            GeneratePoint<NodeType >(10., 0.0, 2.0),
            GeneratePoint<NodeType >(0.0, -1.0, 0.00)
            );

        KRATOS_CHECK(triangle_1.HasIntersection(triangle_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3InsideIntersection, KratosCoreFastSuite) {
        Triangle3D3<Point > triangle_1(
            GeneratePoint<NodeType >(0.0, 0.0, 0.0),
            GeneratePoint<NodeType >(0.0, 0.0, 4.0),
            GeneratePoint<NodeType >(0.0, 4.0, 0.0)
            );
        Triangle3D3<Point > triangle_2(
            GeneratePoint<NodeType >(0.0, 1.0, 1.0),
            GeneratePoint<NodeType >(0.0, 1.0, 3.0),
            GeneratePoint<NodeType >(0.0, 3.0, 1.0)
            );

        KRATOS_CHECK(triangle_1.HasIntersection(triangle_2));
    }

    /**
    * Test an overlaping box and triangle (intersects a triangle edge) HasIntersection
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3IntersectionBoxEdge, KratosCoreFastSuite) {
        auto geom = GenerateEquilateralTriangle3D3<NodeType>();
        Point point_1( 0.3, 0.3,-0.3);
        Point point_2( 1.0, 1.0, 1.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));

        Point point_3(-0.3, 0.3, 0.3);
        Point point_4( 1.0, 1.0, 1.0);
        KRATOS_CHECK(geom->HasIntersection(point_3, point_4));

        Point point_5( 0.3,-0.3, 0.3);
        Point point_6( 1.0, 1.0, 1.0);
        KRATOS_CHECK(geom->HasIntersection(point_5, point_6));
    }

    /**
    * Test an overlaping box and triangle (intersects a triangle node) HasIntersection
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3IntersectionBoxNode, KratosCoreFastSuite) {
        auto geom = GenerateEquilateralTriangle3D3<NodeType>();
        Point point_1(-0.5, 0.8,-0.3);
        Point point_2( 0.5, 1.2, 0.3);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));

        Point point_3(-0.3,-0.5, 0.8);
        Point point_4( 0.3, 0.5, 1.2);
        KRATOS_CHECK(geom->HasIntersection(point_3, point_4));

        Point point_5( 1.2, 0.3, 0.5);
        Point point_6( 0.8,-0.3,-0.5);
        KRATOS_CHECK(geom->HasIntersection(point_5, point_6));
    }

    /**
    * Test an overlaping box and triangle (intersects the triangle face) HasIntersection
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3IntersectionBoxPlane, KratosCoreFastSuite) {
        auto geom = GenerateEquilateralTriangle3D3<NodeType>();
        Point point_1( 0.0, 0.0, 0.0);
        Point point_2( 0.4, 0.5, 0.6);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

    /**
    * Test a non overlaping box and triangle HasIntersection
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3IntersectionBoxNoIntersect, KratosCoreFastSuite) {
        auto geom = GenerateEquilateralTriangle3D3<NodeType>();
        Point point_1( 0.4, 0.5, 0.6);
        Point point_2( 1.0, 1.0, 1.0);
        KRATOS_CHECK_IS_FALSE(geom->HasIntersection(point_1, point_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3ShapeFunctionsValues, KratosCoreFastSuite) {
        auto geom = GenerateEquilateralTriangle3D3<NodeType>();
        array_1d<double, 3> coord(3);
        coord[0] = 1.0 / 2.0;
        coord[1] = 1.0 / 8.0;
        coord[2] = 0.0;
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), 0.375, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(2, coord), 0.125, TOLERANCE);
        CrossCheckShapeFunctionsValues(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D3ShapeFunctionsLocalGradients, KratosCoreFastSuite) {
        auto geom = GenerateEquilateralTriangle3D3<NodeType>();
        TestAllShapeFunctionsLocalGradients(*geom);
    }

} // namespace Testing.
} // namespace Kratos.
