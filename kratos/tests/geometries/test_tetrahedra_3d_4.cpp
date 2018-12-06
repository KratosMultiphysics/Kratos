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
//                   Vicente Mataix Ferrandiz
//
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/tetrahedra_3d_4.h"
#include "tests/geometries/test_geometry.h"
#include "tests/geometries/test_shape_function_derivatives.h"
#include "tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos {
	namespace Testing {

    typedef Node<3>                   PointType;
    typedef Node<3>::Pointer          PointPtrType;
    typedef Tetrahedra3D4<PointType>  TetGeometryType;
    typedef TetGeometryType::Pointer  TetGeometryPtrType;

    /** Generates a sample Tetrahedra3D4.
     * Generates a tetrahedra defined by three random points in the space.
     * @return  Pointer to a Tetrahedra3D4
     */
    TetGeometryPtrType GenerateTetrahedra3D4(
        PointPtrType PointA = GeneratePoint<PointType>(),
        PointPtrType PointB = GeneratePoint<PointType>(),
        PointPtrType PointC = GeneratePoint<PointType>(),
        PointPtrType PointD = GeneratePoint<PointType>()) {
      return TetGeometryPtrType(new TetGeometryType(PointA, PointB, PointC, PointD));
    }

    /** Generates a sample Tetrahedra3D4.
     * Generates a trirectangular tetrahedra on the origin with positive volume and side 1.
     * @return  Pointer to a Tetrahedra3D4
     */
    TetGeometryPtrType GenerateTriRectangularTetrahedra3D4() {
      return TetGeometryPtrType(new TetGeometryType(
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
    TetGeometryPtrType GenerateRegInvtLen1Tetrahedra3D4() {
      return TetGeometryPtrType(new TetGeometryType(
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
    TetGeometryPtrType GenerateRegularLen1Tetrahedra3D4() {
      return TetGeometryPtrType(new TetGeometryType(
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
    TetGeometryPtrType GenerateRegularLen2Tetrahedra3D4() {
      return TetGeometryPtrType(new TetGeometryType(
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

      auto criteria = TetGeometryType::QualityCriteria::INRADIUS_TO_CIRCUMRADIUS;

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

      auto criteria = TetGeometryType::QualityCriteria::INRADIUS_TO_LONGEST_EDGE;

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

      auto criteria = TetGeometryType::QualityCriteria::SHORTEST_TO_LONGEST_EDGE;

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

      auto criteria = TetGeometryType::QualityCriteria::REGULARITY;

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

      auto criteria = TetGeometryType::QualityCriteria::VOLUME_TO_SURFACE_AREA;

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

      auto criteria = TetGeometryType::QualityCriteria::VOLUME_TO_EDGE_LENGTH;

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

      auto criteria = TetGeometryType::QualityCriteria::VOLUME_TO_AVERAGE_EDGE_LENGTH;

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

      auto criteria = TetGeometryType::QualityCriteria::VOLUME_TO_RMS_EDGE_LENGTH;

      KRATOS_CHECK_NEAR(geomInvLen1->Quality(criteria), -1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen1->Quality(criteria),  1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomRegLen2->Quality(criteria),  1.000000, TOLERANCE);
      KRATOS_CHECK_NEAR(geomTriRect->Quality(criteria),  0.769800, TOLERANCE);
    }

    /**
     * This test performs the check of the box intersection method
     */
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
    
    /** Checks the inside test for a given point respect to the tetrahedra
    * Checks the inside test for a given point respect to the tetrahedra
    * It performs 4 tests:
    * A Point inside the tetrahedra: Expected result TRUE
    * A Point outside the tetrahedra: Expected result FALSE
    * A Point over a vertex of the tetrahedra: Expected result TRUE
    * A Point over an edge of the tetrahedra: Expected result TRUE
    */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4IsInside, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateTriRectangularTetrahedra3D4();

        Point PointInside(0.1666, 0.1666, 0.1666);
        Point PointOutside(0.66, 0.66, 0.66);
        Point PointInVertex(0.0, 0.0, 0.0);
        Point PointInEdge(0.33, 0.33, 0.33);

        Point LocalCoords;
        
        KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
        KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
        KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
        KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
    }

    /** Checks the point local coordinates for a given point respect to the
    * tetrahedra. The baricentre of the tetrahedra is selected due to its known
    * solution.
    */
    KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateTriRectangularTetrahedra3D4();

        // Compute the global coordinates of the baricentre
        auto points = geom->Points();
        auto baricentre = Point{points[0] + points[1] + points[2] + points[3]};
        baricentre /= 3.0;

        // Compute the baricentre local coordinates
        array_1d<double, 3> baricentre_local_coords;
        geom->PointLocalCoordinates(baricentre_local_coords, baricentre);

        KRATOS_CHECK_NEAR(baricentre_local_coords(0), 1.0/3.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords(1), 1.0/3.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords(2), 1.0/3.0, TOLERANCE);
        
        Point baricentre_face_1;
        baricentre_face_1.Coordinates()[0] = 0.5;
        baricentre_face_1.Coordinates()[1] = 0.5;
        baricentre_face_1.Coordinates()[2] = 0.0;

        // Compute the baricentre local coordinates
        array_1d<double, 3> baricentre_local_coords_face_1;
        geom->PointLocalCoordinates(baricentre_local_coords_face_1, baricentre_face_1);
        
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_1(0), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_1(1), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_1(2), 0.0, TOLERANCE);
        
        Point baricentre_face_2;
        baricentre_face_2.Coordinates()[0] = 0.5;
        baricentre_face_2.Coordinates()[1] = 0.0;
        baricentre_face_2.Coordinates()[2] = 0.5;

        // Compute the baricentre local coordinates
        array_1d<double, 3> baricentre_local_coords_face_2;
        geom->PointLocalCoordinates(baricentre_local_coords_face_2, baricentre_face_2);
        
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_2(0), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_2(1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_2(2), 0.5, TOLERANCE);
        
        Point baricentre_face_3;
        baricentre_face_3.Coordinates()[0] = 0.0;
        baricentre_face_3.Coordinates()[1] = 0.5;
        baricentre_face_3.Coordinates()[2] = 0.5;

        // Compute the baricentre local coordinates
        array_1d<double, 3> baricentre_local_coords_face_3;
        geom->PointLocalCoordinates(baricentre_local_coords_face_3, baricentre_face_3);

        KRATOS_CHECK_NEAR(baricentre_local_coords_face_3(0), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_3(1), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(baricentre_local_coords_face_3(2), 0.5, TOLERANCE);

        Point outside_point;
        outside_point.Coordinates()[0] = 0.5;
        outside_point.Coordinates()[1] = 0.5;
        outside_point.Coordinates()[2] = 0.5;

        // Compute the baricentre local coordinates
        array_1d<double, 3> local_coords_outside_point;
        geom->PointLocalCoordinates(local_coords_outside_point, outside_point);

        KRATOS_CHECK_NEAR(local_coords_outside_point(0), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(local_coords_outside_point(1), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(local_coords_outside_point(2), 0.5, TOLERANCE);
    }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4();
      array_1d<double, 3> coord(3);
      coord[0] = 1.0 / 2.0;
      coord[1] = 1.0 / 4.0;
      coord[2] = 1.0 / 16.0;
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), 0.1875, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 0.5, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(2, coord), 0.25, TOLERANCE);
      KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(3, coord), 0.0625, TOLERANCE);
      CrossCheckShapeFunctionsValues(*geom);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D4ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
      auto geom = GenerateTriRectangularTetrahedra3D4();
      TestAllShapeFunctionsLocalGradients(*geom);
  }
	}
}  // namespace Kratos.
