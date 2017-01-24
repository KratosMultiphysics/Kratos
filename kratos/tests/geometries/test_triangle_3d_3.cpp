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

// Utility includes
#include "utilities/geometry_utilities.h"

namespace Kratos {
namespace Testing {

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
      GeneratePoint<TPointType>(std::cos(M_PI/4), 0.0, std::sin(M_PI/4)),
      GeneratePoint<TPointType>(0.0, 1.0, 0.0)
    ));
  }

  /// Tests

  /** Checks if the number of edges is correct.
   * Checks if the number of edges is correct.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle3D3EdgesNumber, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle3D3<Node<3>>();

    KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 3);
  }

  /** Checks if the number of faces is correct.
   * Checks if the number of faces is correct.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle3D3FacesNumber, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle3D3<Node<3>>();

    // Charlie: I will let this to 3 but probably 'FacesNumber' needs to be documented to state
    // that for planar geometries it also return the number of edges.
    KRATOS_CHECK_EQUAL(geom->FacesNumber(), 3);
  }

  /** Checks if the area of the triangle is calculated correctly.
   * Checks if the area of the triangle is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle3D3Area, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle3D3<Node<3>>();

    KRATOS_CHECK_NEAR(geom->Area(), 0.5, TOLERANCE);
  }

  /** Checks if the volume of the triangle is calculated correctly.
   * Checks if the volume of the triangle is calculated correctly.
   * For triangle 2D3 'volume()' call defaults to 'area()'
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle3D3Volume, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle3D3<Node<3>>();

    KRATOS_CHECK_EXCPETION_RAISED(geom->Volume(), Exception);
	}

  /** Checks if the minimum edge length is calculated correctly.
   * Checks if the minimum edge length is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle3D3MinEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle3D3<Node<3>>();

    KRATOS_CHECK_NEAR(geom->MinEdgeLength(), 1.0, TOLERANCE);
  }

  /** Checks if the maximum edge length is calculated correctly.
   * Checks if the maximum edge length is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle3D3MaxEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle3D3<Node<3>>();

    KRATOS_CHECK_NEAR(geom->MaxEdgeLength(), 1.414213, TOLERANCE);
  }

  /** Checks if the average edge length is calculated correctly.
   * Checks if the average edge length is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle3D3AverageEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle3D3<Node<3>>();

    KRATOS_CHECK_NEAR(geom->AverageEdgeLength(), 1.138071, TOLERANCE);
  }

  /** Checks if the circumradius is calculated correctly.
   * Checks if the circumradius is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle3D3Circumradius, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle3D3<Node<3>>();

    KRATOS_CHECK_NEAR(geom->Circumradius(), 0.707107, TOLERANCE);
  }

  /** Checks if the inradius is calculated correctly.
   * Checks if the inradius is calculated correctly.
   */
  KRATOS_TEST_CASE_IN_SUITE(Triangle3D3Inradius, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle3D3<Node<3>>();

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
  KRATOS_TEST_CASE_IN_SUITE(Triangle3D3IsInside, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRightTriangle3D3<Node<3>>();

    Point<3> PointInside(0.33, 0.33, 0.0);
    Point<3> PointOutside(0.66, 0.66, 0.0);
    Point<3> PointInVertex(0.0, 0.0, 0.0);
    Point<3> PointInEdge(0.5, 0.5, 0.0);

    Point<3> LocalCoords;

    // It appears that the function checks whether the PROJECTION of the point is inside the geometry.
    KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
    KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
  }

} // namespace Testing.
} // namespace Kratos.
