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
     * Generates a tetrahedra defined by three random points in the space.
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
     * Generates a trirectangular tetrahedra on the origin with positive volume and side 1.
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
     * Generates a regular tetrahedra with positive volume and side 1.
     * @return  Pointer to a Tetrahedra3D4
     */
    GeometryPtrType GenerateRegInvtLen1Tetrahedra3D4() {
      return GeometryPtrType(new GeometryType(
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(1.0, 0.0, 1.0),
        GeneratePoint<PointType>(1.0, 1.0, 0.0),
        GeneratePoint<PointType>(0.0, 0.0, 0.0)
      ));
    }

    /** Generates a sample Tetrahedra3D4.
     * Generates a regular tetrahedra with negative volume and side 1.
     * @return  Pointer to a Tetrahedra3D4
     */
    GeometryPtrType GenerateRegularLen1Tetrahedra3D4() {
      return GeometryPtrType(new GeometryType(
        GeneratePoint<PointType>(0.0, 0.0, 0.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(1.0, 0.0, 1.0),
        GeneratePoint<PointType>(1.0, 1.0, 0.0)
      ));
    }

    /** Generates a sample Tetrahedra3D4.
     * Generates a regular tetrahedra with positive volume with side 2.
     * @return  Pointer to a Tetrahedra3D4
     */
    GeometryPtrType GenerateRegularLen2Tetrahedra3D4() {
      return GeometryPtrType(new GeometryType(
        GeneratePoint<PointType>(0.0, 0.0, 0.0),
        GeneratePoint<PointType>(0.0, 2.0, 2.0),
        GeneratePoint<PointType>(2.0, 0.0, 2.0),
        GeneratePoint<PointType>(2.0, 2.0, 0.0)
      ));
    }

    /** Checks if the number of edges is correct.
     * Checks if the number of edges is correct.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4EdgesNumber, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_EQUAL(geomInvLen1->EdgesNumber(), 6);
      KRATOS_CHECK_EQUAL(geomRegLen1->EdgesNumber(), 6);
      KRATOS_CHECK_EQUAL(geomRegLen2->EdgesNumber(), 6);
      KRATOS_CHECK_EQUAL(geomTriRect->EdgesNumber(), 6);
    }

    /** Checks if the number of faces is correct.
     * Checks if the number of faces is correct.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4FacesNumber, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      // Charlie: I will let this to 3 but probably 'FacesNumber' needs to be documented to state
      // that for planar geometries it also return the number of edges.
      KRATOS_CHECK_EQUAL(geomInvLen1->FacesNumber(), 4);
      KRATOS_CHECK_EQUAL(geomRegLen1->FacesNumber(), 4);
      KRATOS_CHECK_EQUAL(geomRegLen2->FacesNumber(), 4);
      KRATOS_CHECK_EQUAL(geomTriRect->FacesNumber(), 4);
    }

    /** Checks if the characteristic length of the tetrahedra is calculated correctly.
     * Checks if the characteristic length of the tetrahedra is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Length, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geomInvLen1->Length(), 1.414213, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Length(), 1.414213, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Length(), 2.828427, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Length(), 1.122462, TOLERANCE);
    }

    /** Checks if the area of the tetrahedra is calculated correctly.
     * Checks if the area of the tetrahedra is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Area, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geomInvLen1->Area(), -1.0/3.0, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Area(),  1.0/3.0, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Area(),  8.0/3.0, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Area(),  1.0/6.0, TOLERANCE);
    }

    /** Checks if the volume of the tetrahedra is calculated correctly.
     * Checks if the volume of the tetrahedra is calculated correctly.
     * For tetrahedra 3D4 'volume()' call defaults to 'area()'
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Volume, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geomInvLen1->Volume(), -1.0/3.0, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Volume(),  1.0/3.0, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Volume(),  8.0/3.0, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Volume(),  1.0/6.0, TOLERANCE);
  	}

    /** Checks if the minimum edge length is calculated correctly.
     * Checks if the minimum edge length is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4MinEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geomInvLen1->MinEdgeLength(), 1.414213, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->MinEdgeLength(), 1.414213, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->MinEdgeLength(), 2.828427, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->MinEdgeLength(), 1.000000, TOLERANCE);
    }

    /** Checks if the maximum edge length is calculated correctly.
     * Checks if the maximum edge length is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4MaxEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geomInvLen1->MaxEdgeLength(), 1.414213, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->MaxEdgeLength(), 1.414213, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->MaxEdgeLength(), 2.828427, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->MaxEdgeLength(), 1.414213, TOLERANCE);
    }

    /** Checks if the average edge length is calculated correctly.
     * Checks if the average edge length is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4AverageEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geomInvLen1->AverageEdgeLength(), 1.414213, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->AverageEdgeLength(), 1.414213, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->AverageEdgeLength(), 2.828427, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->AverageEdgeLength(), 1.207106, TOLERANCE);
    }

    /** Checks if the circumradius is calculated correctly.
     * Checks if the circumradius is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Circumradius, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geomInvLen1->Circumradius(), 0.866025, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Circumradius(), 0.866025, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Circumradius(), 1.732050, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Circumradius(), 0.866025, TOLERANCE);
    }

    /** Checks if the inradius is calculated correctly.
     * Checks if the inradius is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4Inradius, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      KRATOS_CHECK_NEAR(geomInvLen1->Inradius(), 0.288675, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Inradius(), 0.288675, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Inradius(), 0.577350, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Inradius(), 0.211324, TOLERANCE);
    }

    /** Checks if the inradius to circumradius quality metric is correctly calculated.
     * Checks if the inradius to circumradius quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4InradiusToCircumradiusQuality, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::INRADIUS_TO_CIRCUMRADIUS;

      KRATOS_CHECK_NEAR(geomInvLen1->Quality(criteria), 1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Quality(criteria), 1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Quality(criteria), 1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), 0.732051, TOLERANCE);
    }

    /** Checks if the inradius to longest edge quality metric is correctly calculated.
     * Checks if the inradius to longest edge quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4InradiusToLongestEdgeQuality, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::INRADIUS_TO_LONGEST_EDGE;

      KRATOS_CHECK_NEAR(geomInvLen1->Quality(criteria), 1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Quality(criteria), 1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Quality(criteria), 1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), 0.732051, TOLERANCE);
    }

    /** Checks if the shortest to longest edge length quality metric is correctly calculated.
     * Checks if the shortest to longest edge length quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4ShortestToLongestEdgeQuality, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::SHORTEST_TO_LONGEST_EDGE;

      KRATOS_CHECK_NEAR(geomInvLen1->Quality(criteria), 1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Quality(criteria), 1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Quality(criteria), 1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), 0.707106, TOLERANCE);
    }

    /** Checks if the regularity quality metric is correctly calculated.
     * Checks if the regularity quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4RegularityQuality, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::REGULARITY;

      // KRATOS_CHECK_NEAR(geomRegLen1->Quality(criteria), 1.0, TOLERANCE);
      // KRATOS_CHECK_NEAR(geomRegLen2->Quality(criteria), 1.0, TOLERANCE);
      // KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), -1.0, TOLERANCE);

      KRATOS_CHECK_EXCEPTION_IS_THROWN(geomInvLen1->Quality(criteria), "Method 'RegularityQuality' is not yet implemented for Tetrahedra3D4");
      KRATOS_CHECK_EXCEPTION_IS_THROWN(geomRegLen1->Quality(criteria), "Method 'RegularityQuality' is not yet implemented for Tetrahedra3D4");
      KRATOS_CHECK_EXCEPTION_IS_THROWN(geomRegLen2->Quality(criteria), "Method 'RegularityQuality' is not yet implemented for Tetrahedra3D4");
      KRATOS_CHECK_EXCEPTION_IS_THROWN(geomTriRect->Quality(criteria), "Method 'RegularityQuality' is not yet implemented for Tetrahedra3D4");
    }

    /** Checks if the volume to surface area quality metric is correctly calculated.
     * Checks if the volume to surface area quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4VolumeToSurfaceAreaQuality, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::VOLUME_TO_SURFACE_AREA;

      // KRATOS_CHECK_NEAR(geomRegLen1->Quality(criteria), 1.0, TOLERANCE);
      // KRATOS_CHECK_NEAR(geomRegLen2->Quality(criteria), 1.0, TOLERANCE);
      // KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria), -1.0, TOLERANCE);

      KRATOS_CHECK_EXCEPTION_IS_THROWN(geomInvLen1->Quality(criteria), "Method 'VolumeToSurfaceAreaQuality' is not yet implemented for Tetrahedra3D4");
      KRATOS_CHECK_EXCEPTION_IS_THROWN(geomRegLen1->Quality(criteria), "Method 'VolumeToSurfaceAreaQuality' is not yet implemented for Tetrahedra3D4");
      KRATOS_CHECK_EXCEPTION_IS_THROWN(geomRegLen2->Quality(criteria), "Method 'VolumeToSurfaceAreaQuality' is not yet implemented for Tetrahedra3D4");
      KRATOS_CHECK_EXCEPTION_IS_THROWN(geomTriRect->Quality(criteria), "Method 'VolumeToSurfaceAreaQuality' is not yet implemented for Tetrahedra3D4");
    }

    /** Checks if the volume to edge length quality metric is correctly calculated.
     * Checks if the volume to edge length quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4VolumeToEdgeLengthQuality, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::VOLUME_TO_EDGE_LENGTH;

      KRATOS_CHECK_NEAR(geomInvLen1->Quality(criteria), -1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Quality(criteria),  1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Quality(criteria),  1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria),  0.839947, TOLERANCE);
    }

    /** Checks if the volume to average edge length quality metric is correctly calculated.
     * Checks if the volume to average edge length quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4VolumeToAverageEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH;

      KRATOS_CHECK_NEAR(geomInvLen1->Quality(criteria), -1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Quality(criteria),  1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Quality(criteria),  1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria),  0.804041, TOLERANCE);
    }

    /** Checks if the volume to RMS edge length quality metric is correctly calculated.
     * Checks if the volume to RMS edge length quality metric is correctly calculated for:
     * - Regular tetrahedra, which should return a perfect score.
     * - TriRectangular tetrahedra, which should return a sub-optimal score.
     */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4VolumeToRMSEdgeLength, KratosCoreGeometriesFastSuite) {
      auto geomInvLen1 = GenerateRegInvtLen1Tetrahedra3D4();
      auto geomRegLen1 = GenerateRegularLen1Tetrahedra3D4();
      auto geomRegLen2 = GenerateRegularLen2Tetrahedra3D4();
      auto geomTriRect = GenerateTriRectangularTetrahedra3D4();

      auto criteria = GeometryType::QualityCriteria::VOLUME_TO_RMS_EDGE_LENGTH;

      KRATOS_CHECK_NEAR(geomInvLen1->Quality(criteria), -1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Quality(criteria),  1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Quality(criteria),  1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria),  0.769800, TOLERANCE);
    }


    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4BoxIntersection, KratosCoreGeometriesFastSuite) {
      auto tetrahedron = GenerateTriRectangularTetrahedra3D4();

      //tetrahedron inside the box
      KRATOS_CHECK(tetrahedron->HasIntersection(Point(-.1,-.2,-.1), Point(1.1,1.1,1.2)));

      //tetrahedron contains the box
      KRATOS_CHECK(tetrahedron->HasIntersection(Point(.25,.25,.25), Point(.26,.26,.26)));

      //tetrahedron intersects the box
      KRATOS_CHECK(tetrahedron->HasIntersection(Point(.25,.25,.25), Point(1.1,1.1,1.2)));

      //tetrahedron not intersects the box
      KRATOS_CHECK_IS_FALSE(tetrahedron->HasIntersection(Point(.51,.51,.51), Point(1.1,1.1,1.2)));
    }


	}
}  // namespace Kratos.
