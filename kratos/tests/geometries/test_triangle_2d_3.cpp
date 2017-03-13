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
#include "tests/geometries/test_geometry.h"

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

  /** Generates a sample triangle2D3.
   * Generates a right triangle with origin in the origin and leg size 1.
   * @return  Pointer to a triangle2D3
   */
  template<class TPointType>
  typename Triangle2D3<TPointType>::Pointer GenerateRightTriangle2D3() {
    return typename Triangle2D3<TPointType>::Pointer(new Triangle2D3<TPointType>(
      GeneratePoint<TPointType>(0.0, 0.0, 0.0),
      GeneratePoint<TPointType>(1.0, 0.0, 0.0),
      GeneratePoint<TPointType>(0.0, 1.0, 0.0)
    ));
  }

  /// Tests

  /** Checks if the number of edges is correct.
   * Checks if the number of edges is correct.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3EdgesNumber, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 3);
  }

  /** Checks if the number of faces is correct.
   * Checks if the number of faces is correct.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3FacesNumber, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    // Charlie: I will let this to 3 but probably 'FacesNumber' needs to be documented to state
    // that for planar geometries it also return the number of edges.
    KRATOS_CHECK_EQUAL(geom->FacesNumber(), 3);
  }

  /** Checks if the area of the triangle is calculated correctly.
   * Checks if the area of the triangle is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3Area, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    KRATOS_CHECK_NEAR(geom->Area(), 0.5, TOLERANCE);
  }

  /** Checks if the area of the triangle is calculated correctly.
   * Checks if the area of the triangle is calculated correctly using the Jaccobian matrix.
   * This test correctness is tied to the correctness of 'TestTriangle2D3Area'.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3AreaJaccobi, KratosCoreGeometriesFastSuite) {
		auto geom = GenerateRightTriangle2D3<Node<3>>();

    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
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
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    KRATOS_CHECK_EXCEPTION_RAISED(geom->Volume(), Exception);
	}

  /** Checks if the minimum edge length is calculated correctly.
   * Checks if the minimum edge length is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3MinEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    KRATOS_CHECK_NEAR(geom->MinEdgeLength(), 1.0, TOLERANCE);
  }

  /** Checks if the maximum edge length is calculated correctly.
   * Checks if the maximum edge length is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3MaxEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    KRATOS_CHECK_NEAR(geom->MaxEdgeLength(), 1.414213, TOLERANCE);
  }

  /** Checks if the average edge length is calculated correctly.
   * Checks if the average edge length is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3AverageEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    KRATOS_CHECK_NEAR(geom->AverageEdgeLength(), 1.138071, TOLERANCE);
  }

  /** Checks if the circumradius is calculated correctly.
   * Checks if the circumradius is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3Circumradius, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    KRATOS_CHECK_NEAR(geom->Circumradius(), 0.707107, TOLERANCE);
  }

  /** Checks if the inradius is calculated correctly.
   * Checks if the inradius is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3Inradius, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();

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
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    Point<3> PointInside(0.33, 0.33);
    Point<3> PointOutside(0.66, 0.66);
    Point<3> PointInVertex(0.0, 0.0);
    Point<3> PointInEdge(0.5, 0.5);

    Point<3> LocalCoords;

    KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
    KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
  }

  /** Tests the area using 'GI_GAUSS_1' integration method.
   * Tests the area using 'GI_GAUSS_1' integration method.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3GaussPoint1, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
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
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
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
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
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
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
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
    auto geom = GenerateRightTriangle2D3<Node<3>>();

    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
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
    auto geom = GenerateRightTriangle2D3<Node<3>>();
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
    auto geom = GenerateRightTriangle2D3<Node<3>>();
    const double ExpectedJacobian = 1.0;

    Vector JacobianDeterminants;
    geom->DeterminantOfJacobian( JacobianDeterminants, GeometryData::GI_GAUSS_2 );

    for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
    {
        KRATOS_CHECK_NEAR(JacobianDeterminants[i], ExpectedJacobian, TOLERANCE);
    }
  }

  /** Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
   * Tests the Jacobian determinants using 'GI_GAUSS_1' integration method.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();
    const double ExpectedJacobian = 1.0;

    double JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_1 );
    KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
  }

  /** Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
   * Tests the Jacobian determinants using 'GI_GAUSS_2' integration method.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle2D3DeterminantOfJacobianIndex2, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle2D3<Node<3>>();
    double JacobianDeterminant = 0.0;
    const double ExpectedJacobian = 1.0;

    JacobianDeterminant = geom->DeterminantOfJacobian( 1, GeometryData::GI_GAUSS_2 );
    KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);

    JacobianDeterminant = geom->DeterminantOfJacobian( 2, GeometryData::GI_GAUSS_2 );
    KRATOS_CHECK_NEAR(JacobianDeterminant, ExpectedJacobian, TOLERANCE);
  }

} // namespace Testing.
} // namespace Kratos.
