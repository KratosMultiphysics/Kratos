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
#include "geometries/triangle_2d_3.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos {
namespace Testing {

  /// Factory functions

  /** Generates a sample triangle2D3.
   * Generates a triangle defined by three random points in the space.
   * @return  Pointer to a triangle2D3
   */
  template<class TPointType>
  typename Triangle2D3<TPointType>::Pointer GenerateTriangle2D3(
      typename TPointType::Pointer PointA = GeneratePoint<TPointType>(),
      typename TPointType::Pointer PointB = GeneratePoint<TPointType>(),
      typename TPointType::Pointer PointC = GeneratePoint<TPointType>()) {
    return typename Triangle2D3<TPointType>::Pointer(new Triangle2D3<TPointType>(
      PointA,
      PointB,
      PointC
    ));
  }

  /** Generates a point type sample triangle2D3.
   * Generates a point type right triangle with origin in the origin and leg size 1.
   * @return  Pointer to a triangle2D3
   */
  Triangle2D3<Point>::Pointer GeneratePointsRightTriangle2D3() {
    return Triangle2D3<Point>::Pointer(new Triangle2D3<Point>(
      Point::Pointer(new Point(0.0, 0.0, 0.0)),
      Point::Pointer(new Point(1.0, 0.0, 0.0)),
      Point::Pointer(new Point(0.0, 1.0, 0.0))
    ));
  }

  /** Generates a node type sample triangle2D3.
   * Generates a point type right triangle with origin in the origin and leg size 1.
   * @return  Pointer to a triangle2D3
   */
  Triangle2D3<Node<3>>::Pointer GenerateNodesRightTriangle2D3() {
    return Triangle2D3<Node<3>>::Pointer(new Triangle2D3<Node<3>>(
      Node<3>::Pointer(new Node<3>(1, 0.0, 0.0, 0.0)),
      Node<3>::Pointer(new Node<3>(2, 1.0, 0.0, 0.0)),
      Node<3>::Pointer(new Node<3>(3, 0.0, 1.0, 0.0))
    ));
  }

  /** Generates a point type sample triangle2D3.
   * Generates a point type  irregular triangle.
   * @return  Pointer to a triangle2D3
   */
  Triangle2D3<Point>::Pointer GeneratePointsIrregularTriangle2D3() {
    return Triangle2D3<Point>::Pointer(new Triangle2D3<Point>(
      Point::Pointer(new Point(1.0, 1.0, 0.0)),
      Point::Pointer(new Point(3.0, 0.5, 0.0)),
      Point::Pointer(new Point(2.5, 2.0, 0.0))
    ));
  }

  /** Generates a sample triangle2D3.
   * Generates a node irregular triangle.
   * @return  Pointer to a triangle2D3
   */
  Triangle2D3<Node<3>>::Pointer GenerateNodesIrregularTriangle2D3() {
    return Triangle2D3<Node<3>>::Pointer(new Triangle2D3<Node<3>>(
      Node<3>::Pointer(new Node<3>(1, 1.0, 1.0, 0.0)),
      Node<3>::Pointer(new Node<3>(2, 3.0, 0.5, 0.0)),
      Node<3>::Pointer(new Node<3>(3, 2.5, 2.0, 0.0))
    ));
  }

  /// Tests

  /** Checks if the number of edges is correct.
   * Checks if the number of edges is correct.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3EdgesNumber, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();

    KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 3);
  }

  /** Checks if the number of faces is correct.
   * Checks if the number of faces is correct.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3FacesNumber, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();

    // Charlie: I will let this to 3 but probably 'FacesNumber' needs to be documented to state
    // that for planar geometries it also return the number of edges.
    KRATOS_CHECK_EQUAL(geom->FacesNumber(), 3);
  }

  /** Checks if the area of the triangle is calculated correctly.
   * Checks if the area of the triangle is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3Area, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();

    KRATOS_CHECK_NEAR(geom->Area(), 0.5, TOLERANCE);
  }

  /** Checks if the area of the triangle is calculated correctly.
   * Checks if the area of the triangle is calculated correctly using the Jaccobian matrix.
   * This test correctness is tied to the correctness of 'TestTriangle2D3Area'.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3AreaJaccobi, KratosCoreGeometriesFastSuite) {
		auto geom = GenerateNodesRightTriangle2D3();

    BoundedMatrix<double,3,2> DN_DX;
    array_1d<double,3> N;
    double Area;

    GeometryUtils::CalculateGeometryData(*geom, DN_DX, N, Area);

		KRATOS_CHECK_NEAR(Area, 0.5, TOLERANCE);
	}

  /** Checks if the volume of the triangle is calculated correctly.
   * Checks if the volume of the triangle is calculated correctly.
   * For triangle 2D3 'volume()' call defaults to 'area()'
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3Volume, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();

    KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Volume(), "Calling base class 'Volume' method instead of derived class one.");
	}

  /** Checks if the minimum edge length is calculated correctly.
   * Checks if the minimum edge length is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3MinEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();

    KRATOS_CHECK_NEAR(geom->MinEdgeLength(), 1.0, TOLERANCE);
  }

  /** Checks if the maximum edge length is calculated correctly.
   * Checks if the maximum edge length is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3MaxEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();

    KRATOS_CHECK_NEAR(geom->MaxEdgeLength(), 1.414213, TOLERANCE);
  }

  /** Checks if the average edge length is calculated correctly.
   * Checks if the average edge length is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3AverageEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();

    KRATOS_CHECK_NEAR(geom->AverageEdgeLength(), 1.138071, TOLERANCE);
  }

  /** Checks if the circumradius is calculated correctly.
   * Checks if the circumradius is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3Circumradius, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();

    KRATOS_CHECK_NEAR(geom->Circumradius(), 0.707107, TOLERANCE);
  }

  /** Checks if the inradius is calculated correctly.
   * Checks if the inradius is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3Inradius, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();

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
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3IsInside, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();

    Point PointInside(0.33, 0.33);
    Point PointOutside(0.66, 0.66);
    Point PointInVertex(0.0, 0.0);
    Point PointInEdge(0.5, 0.5);

    Point LocalCoords;

    KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
    KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
  }

  /** Checks the point local coordinates for a given point respect to the
   * triangle. The baricentre of the triangle is selected due to its known
   * solution.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsIrregularTriangle2D3();

    // Compute the global coordinates of the baricentre
    const Geometry<Point>::PointsArrayType geom_pts = geom->Points();
    auto baricentre = Point{geom_pts[0] + geom_pts[1] + geom_pts[2]};
    baricentre *= 1.0/3.0;

    // Compute the baricentre local coordinates
    array_1d<double, 3> baricentre_local_coords;
    geom->PointLocalCoordinates(baricentre_local_coords, baricentre);

    KRATOS_CHECK_NEAR(baricentre_local_coords(0), 1.0/3.0, TOLERANCE);
    KRATOS_CHECK_NEAR(baricentre_local_coords(1), 1.0/3.0, TOLERANCE);
    KRATOS_CHECK_NEAR(baricentre_local_coords(2), 0.0, TOLERANCE);
  }

  /** Tests the area using 'GI_GAUSS_1' integration method.
   * Tests the area using 'GI_GAUSS_1' integration method.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3GaussPoint1, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateNodesRightTriangle2D3();

    BoundedMatrix<double,3,2> DN_DX;
    array_1d<double,3> N;
    double ExpectedArea;

    GeometryUtils::CalculateGeometryData(*geom, DN_DX, N, ExpectedArea);

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_1), ExpectedArea, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_1);
  }

  /** Tests the area using 'GI_GAUSS_2' integration method.
   * Tests the area using 'GI_GAUSS_2' integration method.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3GaussPoint2, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateNodesRightTriangle2D3();

    BoundedMatrix<double,3,2> DN_DX;
    array_1d<double,3> N;
    double ExpectedArea;

    GeometryUtils::CalculateGeometryData(*geom, DN_DX, N, ExpectedArea);

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_2), ExpectedArea, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_2);
  }

  /** Tests the area using 'GI_GAUSS_3' integration method.
   * Tests the area using 'GI_GAUSS_3' integration method.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3GaussPoint3, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateNodesRightTriangle2D3();

    BoundedMatrix<double,3,2> DN_DX;
    array_1d<double,3> N;
    double ExpectedArea;

    GeometryUtils::CalculateGeometryData(*geom, DN_DX, N, ExpectedArea);

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_3), ExpectedArea, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_3);
  }

  /** Tests the area using 'GI_GAUSS_4' integration method.
   * Tests the area using 'GI_GAUSS_4' integration method.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3GaussPoint4, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateNodesRightTriangle2D3();

    BoundedMatrix<double,3,2> DN_DX;
    array_1d<double,3> N;
    double ExpectedArea;

    GeometryUtils::CalculateGeometryData(*geom, DN_DX, N, ExpectedArea);

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_4), ExpectedArea, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_4);
  }

  /** Tests the area using 'GI_GAUSS_5' integration method.
   * Tests the area using 'GI_GAUSS_5' integration method.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3GaussPoint5, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateNodesRightTriangle2D3();

    BoundedMatrix<double,3,2> DN_DX;
    array_1d<double,3> N;
    double ExpectedArea;

    GeometryUtils::CalculateGeometryData(*geom, DN_DX, N, ExpectedArea);

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_5), ExpectedArea, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_5);
  }

  /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
   * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();
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
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianArray2, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();
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
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianArray3, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();
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
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianArray4, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();
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
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianArray5, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();
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
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();
    const double ExpectedJacobian = 1.0;

    double JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_1 );
    KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
  }

  /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
   * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianIndex2, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();
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
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianIndex3, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();
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
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianIndex4, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();
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
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianIndex5, KratosCoreGeometriesFastSuite) {
    auto geom = GeneratePointsRightTriangle2D3();
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

    /**
     * Test an overlaping box and triangle (intersects a triangle edge) HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Triangle2D3IntersectionBoxEdge, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsRightTriangle2D3();
        Point point_1(-0.1, 0.1, 0.0);
        Point point_2( 0.1, 0.3, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));

        Point point_3( 0.1,-0.1, 0.0);
        Point point_4( 0.3, 0.1, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_3, point_4));

        Point point_5( 0.3, 0.2, 0.0);
        Point point_6( 1.0, 1.0, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_5, point_6));
    }

    /**
     * Test an overlaping box and triangle (intersects a triangle node) HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Triangle2D3IntersectionBoxNode, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsRightTriangle2D3();
        Point point_1(-0.5, 0.8, 0.0);
        Point point_2( 0.5, 1.2, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));

        Point point_3( 0.3,-0.5, 0.0);
        Point point_4( 1.2, 0.5, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_3, point_4));

        Point point_5( 0.2, 0.3, 0.0);
        Point point_6(-0.8,-0.3, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_5, point_6));
    }

    /**
     * Test a box inside a triangle HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Triangle2D3IntersectionBoxInside, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsRightTriangle2D3();
        Point point_1( 0.1, 0.1, 0.0);
        Point point_2( 0.3, 0.4, 0.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

    /**
     * Test a non overlaping box and triangle HasIntersection
     */
    KRATOS_TEST_CASE_IN_SUITE(Triangle2D3IntersectionBoxNoIntersect, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsRightTriangle2D3();
        Point point_1( 0.6, 0.5, 0.0);
        Point point_2( 1.0, 1.0, 0.0);
        KRATOS_CHECK_IS_FALSE(geom->HasIntersection(point_1, point_2));
    }

    /**
     * Test a box outside the triangle plane
     * HasIntersection should return true, because a 2D-space does not take care about Z-coordinates
     */
    KRATOS_TEST_CASE_IN_SUITE(Triangle2D3IntersectionBoxOutsidePlane, KratosCoreGeometriesFastSuite) {
        auto geom = GeneratePointsRightTriangle2D3();
        Point point_1( 0.2, 0.1, 0.1);
        Point point_2( 0.3, 0.5, 1.0);
        KRATOS_CHECK(geom->HasIntersection(point_1, point_2));
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D3ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateNodesRightTriangle2D3();
      array_1d<double, 3> coord(3);
      coord[0] = 1.0 / 2.0;
      coord[1] = 1.0 / 8.0;
      coord[2] = 0.0;
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), 0.375, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 0.5, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(2, coord), 0.125, TOLERANCE);
      CrossCheckShapeFunctionsValues(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D3ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateNodesRightTriangle2D3();
      TestAllShapeFunctionsLocalGradients(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D3Normal, KratosCoreGeometriesFastSuite) {

        auto geom = GeneratePointsRightTriangle2D3();
        Point::CoordinatesArrayType LocalCoord;
        LocalCoord.clear();
        auto normal = geom->Normal(LocalCoord);

        KRATOS_CHECK_NEAR(normal[0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(normal[1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(normal[2], 0.5, TOLERANCE);

        array_1d<double, 3> cross_norm;
        cross_norm[0] = 0.0;
        cross_norm[1] = 0.0;
        cross_norm[2] = 1.0;
        array_1d<double, 3> cross;
        MathUtils<double>::CrossProduct(cross, cross_norm, normal);

        KRATOS_CHECK_NEAR(cross[0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(cross[1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(cross[2], 0.0, TOLERANCE);

        normal /= norm_2(normal);

        auto unit_normal = geom->UnitNormal(LocalCoord);

        KRATOS_CHECK_NEAR(unit_normal[0], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(unit_normal[1], 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(unit_normal[2], 1.0, TOLERANCE);

        KRATOS_CHECK_NEAR(unit_normal[0], normal[0], TOLERANCE);
        KRATOS_CHECK_NEAR(unit_normal[1], normal[1], TOLERANCE);
        KRATOS_CHECK_NEAR(unit_normal[2], normal[2], TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.
