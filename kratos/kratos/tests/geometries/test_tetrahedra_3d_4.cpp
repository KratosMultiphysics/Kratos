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
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/tetrahedra_3d_4.h"
#include "tests/geometries/test_geometry.h"

namespace Kratos {
	namespace Testing {

    /** Generates a sample Tetrahedra3D4.
     * Generates a triangle defined by three random points in the space.
     * @return  Pointer to a Tetrahedra3D4
     */
    template<class TPointType>
    typename Tetrahedra3D4<TPointType>::Pointer GenerateTetrahedra3D4(
        typename TPointType::Pointer PointA = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointB = GeneratePoint<TPointType>(),
        typename TPointType::Pointer PointC = GeneratePoint<TPointType>()) {
      return typename Tetrahedra3D4<TPointType>::Pointer(new Tetrahedra3D4<TPointType>(
        PointA,
        PointB,
        PointC
      ));
    }

    /** Generates a sample Tetrahedra3D4.
     * Generates a right triangle with origin in the origin and leg size 1.
     * @return  Pointer to a Tetrahedra3D4
     */
    template<class TPointType>
    typename Tetrahedra3D4<TPointType>::Pointer GenerateTriRectangularTetrahedra3D4() {
      return typename Tetrahedra3D4<TPointType>::Pointer(new Tetrahedra3D4<TPointType>(
        GeneratePoint<TPointType>(0.0, 0.0, 0.0),
        GeneratePoint<TPointType>(1.0, 0.0, 0.0),
        GeneratePoint<TPointType>(0.0, 1.0, 0.0),
        GeneratePoint<TPointType>(0.0, 0.0, 1.0)
      ));
    }

    /** Checks if the number of edges is correct.
     * Checks if the number of edges is correct.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4EdgesNumber, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4<Node<3>>();

      KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 6);
    }

    /** Checks if the number of faces is correct.
     * Checks if the number of faces is correct.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4FacesNumber, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4<Node<3>>();

      // Charlie: I will let this to 3 but probably 'FacesNumber' needs to be documented to state
      // that for planar geometries it also return the number of edges.
      KRATOS_CHECK_EQUAL(geom->FacesNumber(), 4);
    }

    /** Checks if the area of the triangle is calculated correctly.
     * Checks if the area of the triangle is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Area, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4<Node<3>>();

      KRATOS_CHECK_NEAR(geom->Area(), 1.0/6.0, TOLERANCE);
    }

    /** Checks if the volume of the triangle is calculated correctly.
     * Checks if the volume of the triangle is calculated correctly.
     * For triangle 2D3 'volume()' call defaults to 'area()'
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Volume, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4<Node<3>>();

      KRATOS_CHECK_NEAR(geom->Volume(), 1.0/6.0, TOLERANCE);
  	}

    /** Checks if the minimum edge length is calculated correctly.
     * Checks if the minimum edge length is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4MinEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4<Node<3>>();

      KRATOS_CHECK_NEAR(geom->MinEdgeLength(), 1.0, TOLERANCE);
    }

    /** Checks if the maximum edge length is calculated correctly.
     * Checks if the maximum edge length is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4MaxEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4<Node<3>>();

      KRATOS_CHECK_NEAR(geom->MaxEdgeLength(), 1.414213, TOLERANCE);
    }

    /** Checks if the average edge length is calculated correctly.
     * Checks if the average edge length is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4AverageEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4<Node<3>>();

      KRATOS_CHECK_NEAR(geom->AverageEdgeLength(), 1.207106, TOLERANCE);
    }

    /** Checks if the circumradius is calculated correctly.
     * Checks if the circumradius is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Circumradius, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4<Node<3>>();

      KRATOS_CHECK_NEAR(geom->Circumradius(), 0.707107, TOLERANCE);
    }

    /** Checks if the inradius is calculated correctly.
     * Checks if the inradius is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Inradius, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4<Node<3>>();

      KRATOS_CHECK_NEAR(geom->Inradius(), 0.292893, TOLERANCE);
    }


	}
}  // namespace Kratos.
