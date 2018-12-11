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
//

// System includes
#include <set>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/prism_3d_6.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

    typedef Node<3>                                PointType;
    typedef Node<3>::Pointer                    PointPtrType;
    typedef Prism3D6<PointType>            PrismGeometryType;
    typedef PrismGeometryType::Pointer  PrismGeometryPtrType;

    /** Generates a sample Prism3D6.
     * Generates a prism defined by three random points in the space.
     * @return  Pointer to a Prism3D6
     */
    PrismGeometryPtrType GeneratePrism3D6(
        PointPtrType PointA = GeneratePoint<PointType>(),
        PointPtrType PointB = GeneratePoint<PointType>(),
        PointPtrType PointC = GeneratePoint<PointType>(),
        PointPtrType PointD = GeneratePoint<PointType>(),
        PointPtrType PointE = GeneratePoint<PointType>(),
        PointPtrType PointF = GeneratePoint<PointType>()) {
      return PrismGeometryPtrType(new PrismGeometryType(PointA, PointB, PointC, PointD, PointE, PointF));
    }

    /** Generates a sample Prism3D6.
     * Generates a trirectangular prism on the origin with positive volume and side 1.
     * @return  Pointer to a Prism3D6
     */
    PrismGeometryPtrType GenerateRegularPrism3D6() {
      return PrismGeometryPtrType(new PrismGeometryType(
        GeneratePoint<PointType>(0.0, 0.0, 0.0),
        GeneratePoint<PointType>(1.0, 0.0, 0.0),
        GeneratePoint<PointType>(0.0, 1.0, 0.0),
        GeneratePoint<PointType>(0.0, 0.0, 1.0),
        GeneratePoint<PointType>(1.0, 0.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0)
      ));
    }

    /** Checks if the number of edges is correct.
     * Checks if the number of edges is correct.
     */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6EdgesNumber, KratosCoreGeometriesFastSuite) {
        auto geomRegular = GenerateRegularPrism3D6();

        KRATOS_CHECK_EQUAL(geomRegular->EdgesNumber(), 9);
    }

    /** Checks if the number of faces is correct.
     * Checks if the number of faces is correct.
     */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6FacesNumber, KratosCoreGeometriesFastSuite) {
        auto geomRegular = GenerateRegularPrism3D6();

        KRATOS_CHECK_EQUAL(geomRegular->FacesNumber(), 5);
    }

    /** Checks if the characteristic length of the prism is calculated correctly.
     * Checks if the characteristic length of the prism is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6Length, KratosCoreGeometriesFastSuite) {
        auto geomRegular = GenerateRegularPrism3D6();

        KRATOS_CHECK_NEAR(geomRegular->Length(), 0.264567, TOLERANCE);
    }

    /** Checks if the area of the prism is calculated correctly.
     * Checks if the area of the prism is calculated correctly.
     */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6Area, KratosCoreGeometriesFastSuite) {
        auto geomRegular = GenerateRegularPrism3D6();

        KRATOS_CHECK_NEAR(geomRegular->Area(),  0.5, TOLERANCE);
    }

    /** Checks if the volume of the prism is calculated correctly.
     * Checks if the volume of the prism is calculated correctly.
     * For prism 3D6 'volume()' call defaults to 'area()'
     */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6Volume, KratosCoreGeometriesFastSuite) {
        auto geomRegular = GenerateRegularPrism3D6();

        KRATOS_CHECK_NEAR(geomRegular->Volume(),  0.5, TOLERANCE);
    }

    /** Checks the inside test for a given point respect to the prism
    * Checks the inside test for a given point respect to the prism
    * It performs 4 tests:
    * A Point inside the prism: Expected result TRUE
    * A Point outside the prism: Expected result FALSE
    * A Point over a vertex of the prism: Expected result TRUE
    * A Point over an edge of the prism: Expected result TRUE
    */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6IsInside, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRegularPrism3D6();

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
    * prism. The centre of the prism is selected due to its known
    * solution.
    */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRegularPrism3D6();

        // Compute the global coordinates of the centre
        auto points = geom->Points();
        Point centre = points[0] + points[1] + points[2] + points[3] + points[4] + points[5];
        centre /= 6.0;

        // Compute the centre local coordinates
        array_1d<double, 3> centre_local_coords;
        geom->PointLocalCoordinates(centre_local_coords, centre);

        KRATOS_CHECK_NEAR(centre_local_coords(0), 1.0/3.0, TOLERANCE);
        KRATOS_CHECK_NEAR(centre_local_coords(1), 1.0/3.0, TOLERANCE);
        KRATOS_CHECK_NEAR(centre_local_coords(2), 1.0/2.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Prism3D6ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRegularPrism3D6();
        array_1d<double, 3> coord(3);
        coord[0] = 1.0 / 2.0;
        coord[1] = 1.0 / 4.0;
        coord[2] = 1.0 / 16.0;
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), 0.234375, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 0.46875, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(2, coord), 0.234375, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(3, coord), 0.015625, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(4, coord), 0.03125, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(5, coord), 0.015625, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Prism3D6ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRegularPrism3D6();
        Matrix gradient;

        // Compute the global coordinates of the centre
        auto points = geom->Points();
        Point centre = points[0] + points[1] + points[2] + points[3] + points[4] + points[5];
        centre /= 6.0;

        // Compute the centre local coordinates
        array_1d<double, 3> centre_local_coords;
        geom->PointLocalCoordinates(centre_local_coords, centre);

        gradient = geom->ShapeFunctionsLocalGradients(gradient, centre_local_coords);

        KRATOS_CHECK_NEAR(gradient(0,0), -0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(0,1), -0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(0,2), -1.0/3.0, TOLERANCE);

        KRATOS_CHECK_NEAR(gradient(1,0), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(1,1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(1,2), -1.0/3.0, TOLERANCE);

        KRATOS_CHECK_NEAR(gradient(2,0), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(2,1), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(2,2), -1.0/3.0, TOLERANCE);

        KRATOS_CHECK_NEAR(gradient(3,0), -0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(3,1), -0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(3,2), 1.0/3.0, TOLERANCE);

        KRATOS_CHECK_NEAR(gradient(4,0), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(4,1), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(4,2), 1.0/3.0, TOLERANCE);

        KRATOS_CHECK_NEAR(gradient(5,0), 0.0, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(5,1), 0.5, TOLERANCE);
        KRATOS_CHECK_NEAR(gradient(5,2), 1.0/3.0, TOLERANCE);
    }

    /** Tests the area using 'GI_GAUSS_1' integration method.
    * Tests the area using 'GI_GAUSS_1' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6GaussPoint1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRegularPrism3D6();

        const double expected_vol = 0.5;

        KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_1), expected_vol, TOLERANCE);
        VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_1);
    }

    /** Tests the area using 'GI_GAUSS_2' integration method.
    * Tests the area using 'GI_GAUSS_2' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6GaussPoint2, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRegularPrism3D6();

        const double expected_vol = 0.5;

        KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_2), expected_vol, TOLERANCE);
        VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_2);
    }

    /** Tests the area using 'GI_GAUSS_3' integration method.
    * Tests the area using 'GI_GAUSS_3' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6GaussPoint3, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRegularPrism3D6();

        const double expected_vol = 0.5;

        KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_3), expected_vol, TOLERANCE);
        VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_3);
    }

    /** Tests the area using 'GI_GAUSS_4' integration method.
    * Tests the area using 'GI_GAUSS_4' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6GaussPoint4, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRegularPrism3D6();

        const double expected_vol = 0.5;

        KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_4), expected_vol, TOLERANCE);
        VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_4);
    }

    /** Tests the area using 'GI_GAUSS_5' integration method.
    * Tests the area using 'GI_GAUSS_5' integration method.
    */
    KRATOS_TEST_CASE_IN_SUITE(Prism3D6GaussPoint5, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateRegularPrism3D6();

        const double expected_vol = 0.5;

        KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_5), expected_vol, TOLERANCE);
        VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_5);
    }
}
}  // namespace Kratos.
