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

    typedef Node<3>                   PointType;
    typedef Node<3>::Pointer          PointPtrType;
    typedef Tetrahedra3D4<PointType>  GeometryType;
    typedef GeometryType::Pointer     GeometryPtrType;

    /** Generates a sample Tetrahedra3D4.
     * Generates a triangle defined by three random points in the space.
     * @return  Pointer to a Tetrahedra3D4
     */
    GeometryPtrType GenerateTetrahedra3D4(
        PointPtrType PointA = GeneratePoint<PointType>(),
        PointPtrType PointB = GeneratePoint<PointType>(),
        PointPtrType PointC = GeneratePoint<PointType>(),
        PointPtrType PointD = GeneratePoint<PointType>()) {
      return GeometryPtrType(new GeometryType(PointA, PointB, PointC, PointD));
    }

    /** Generates a sample Tetrahedra3D4.
     * Generates a trirectangular tetrahedra on the origin with leg size 1.
     * @return  Pointer to a Tetrahedra3D4
     */
    GeometryPtrType GenerateTriRectangularTetrahedra3D4() {
      return GeometryPtrType(new GeometryType(
        GeneratePoint<PointType>(0.0, 0.0, 0.0),
        GeneratePoint<PointType>(1.0, 0.0, 0.0),
        GeneratePoint<PointType>(0.0, 1.0, 0.0),
        GeneratePoint<PointType>(0.0, 0.0, 1.0)
      ));
    }

    /** Generates a sample Tetrahedra3D4.
     * Generates a regular tetrahedra.
     * @return  Pointer to a Tetrahedra3D4
     */
    GeometryPtrType GenerateRegularTetrahedra3D4() {
      return GeometryPtrType(new GeometryType(
        GeneratePoint<PointType>(0.0, 0.0, 0.0),
        GeneratePoint<PointType>(1.0, 0.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(1.0, 1.0, 0.0)
      ));
    }

    /** Checks if the number of edges is correct.
     * Checks if the number of edges is correct.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4EdgesNumber, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 6);
    }

    /** Checks if the number of faces is correct.
     * Checks if the number of faces is correct.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4FacesNumber, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4();

      // Charlie: I will let this to 3 but probably 'FacesNumber' needs to be documented to state
      // that for planar geometries it also return the number of edges.
      KRATOS_CHECK_EQUAL(geom->FacesNumber(), 4);
    }

    /** Checks if the area of the triangle is calculated correctly.
     * Checks if the area of the triangle is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Area, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geom->Area(), 1.0/6.0, TOLERANCE);
    }

    /** Checks if the volume of the triangle is calculated correctly.
     * Checks if the volume of the triangle is calculated correctly.
     * For triangle 2D3 'volume()' call defaults to 'area()'
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Volume, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geom->Volume(), 1.0/6.0, TOLERANCE);
  	}

    /** Checks if the minimum edge length is calculated correctly.
     * Checks if the minimum edge length is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4MinEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geom->MinEdgeLength(), 1.0, TOLERANCE);
    }

    /** Checks if the maximum edge length is calculated correctly.
     * Checks if the maximum edge length is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4MaxEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geom->MaxEdgeLength(), 1.414213, TOLERANCE);
    }

    /** Checks if the average edge length is calculated correctly.
     * Checks if the average edge length is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4AverageEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geom->AverageEdgeLength(), 1.207106, TOLERANCE);
    }

    /** Checks if the circumradius is calculated correctly.
     * Checks if the circumradius is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Circumradius, KratosCoreGeometriesFastSuite) {
      auto geomRegular = GenerateRegularTetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      // Should be the same
      KRATOS_CHECK_NEAR(geomRegular->Circumradius(), 0.866025, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Circumradius(), 0.866025, TOLERANCE);
    }

    /** Checks if the inradius is calculated correctly.
     * Checks if the inradius is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Inradius, KratosCoreGeometriesFastSuite) {
      auto geomRegular = GenerateRegularTetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      // Should be the same
      KRATOS_CHECK_NEAR(geomRegular->Inradius(), 0.288675, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Inradius(), 0.211324, TOLERANCE);
    }

    /** Checks if the inradius to circumradius quality metric is correctly calculated.
     * Checks if the inradius to circumradius quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4InradiusToCircumradiusQuality, KratosCoreGeometriesFastSuite) {
      auto geomRegular = GenerateRegularTetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::INRADIUS_TO_CIRCUMRADIUS;

      KRATOS_CHECK_EQUAL(geomRegular->Quality(criteria), 1.0);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), -1.0, TOLERANCE);
    }

    /** Checks if the inradius to longest edge quality metric is correctly calculated.
     * Checks if the inradius to longest edge quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4InradiusToLongestEdgeQuality, KratosCoreGeometriesFastSuite) {
      auto geomRegular = GenerateRegularTetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::INRADIUS_TO_LONGEST_EDGE;

      KRATOS_CHECK_EQUAL(geomRegular->Quality(criteria), 1.0);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), -1.0, TOLERANCE);
    }

    /** Checks if the shortest to longest edge length quality metric is correctly calculated.
     * Checks if the shortest to longest edge length quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4ShortestToLongestEdgeQuality, KratosCoreGeometriesFastSuite) {
      auto geomRegular = GenerateRegularTetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::SHORTEST_TO_LONGEST_EDGE;

      KRATOS_CHECK_EQUAL(geomRegular->Quality(criteria), 1.0);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), 0.707106, TOLERANCE);
    }

    /** Checks if the regularity quality metric is correctly calculated.
     * Checks if the regularity quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4RegualrityQuiality, KratosCoreGeometriesFastSuite) {
      auto geomRegular = GenerateRegularTetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::REGULARITY;

      KRATOS_CHECK_EQUAL(geomRegular->Quality(criteria), 1.0);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), -1.0, TOLERANCE);
    }

    /** Checks if the volume to surface area quality metric is correctly calculated.
     * Checks if the volume to surface area quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4VolumeToSurfaceAreaQuality, KratosCoreGeometriesFastSuite) {
      auto geomRegular = GenerateRegularTetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::VOLUME_TO_SURFACE_AREA;

      KRATOS_CHECK_EQUAL(geomRegular->Quality(criteria), 1.0);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), -1.0, TOLERANCE);
    }

    /** Checks if the volume to edge length quality metric is correctly calculated.
     * Checks if the volume to edge length quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4VolumeToEdgeLengthQuality, KratosCoreGeometriesFastSuite) {
      auto geomRegular = GenerateRegularTetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::VOLUME_TO_EDGE_LENGTH;

      KRATOS_CHECK_EQUAL(geomRegular->Quality(criteria), 1.0);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), -1.0, TOLERANCE);
    }

    /** Checks if the volume to average edge length quality metric is correctly calculated.
     * Checks if the volume to average edge length quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4VolumeToAverageEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geomRegular = GenerateRegularTetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH;

      KRATOS_CHECK_EQUAL(geomRegular->Quality(criteria), 1.0);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), -1.0, TOLERANCE);
    }

    /** Checks if the volume to RMS edge length quality metric is correctly calculated.
     * Checks if the volume to RMS edge length quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4VolumeToRMSEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geomRegular = GenerateRegularTetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::VOLUME_TO_EDGE_LENGTH;

      KRATOS_CHECK_EQUAL(geomRegular->Quality(criteria), 1.0);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), -1.0, TOLERANCE);
    }



	}
}  // namespace Kratos.
