//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/pyramid_3d_13.h"
#include "tests/cpp_tests/geometries/test_geometry.h"


namespace Kratos {
namespace Testing {

typedef Node<3>                           PointType;
typedef Node<3>::Pointer                  PointPtrType;
typedef Geometry<PointType>               BaseGeometryType;
typedef BaseGeometryType::Pointer         BaseGeometryPtrType;
typedef Pyramid3D13<PointType>            Pyramid3D13GeometryType;
typedef Pyramid3D13GeometryType::Pointer  Pyramid3D13GeometryPtrType;

/** Generates a sample Pyramid3D5.
 * Generates a pyramid defined by three random points in the space.
 * @return Pointer to a Pyramid3D13
 */
BaseGeometryPtrType GeneratePyramid3D13(
    PointPtrType Point1 = GeneratePoint<PointType>(),
    PointPtrType Point2 = GeneratePoint<PointType>(),
    PointPtrType Point3 = GeneratePoint<PointType>(),
    PointPtrType Point4 = GeneratePoint<PointType>(),
    PointPtrType Point5 = GeneratePoint<PointType>(),
    PointPtrType Point6 = GeneratePoint<PointType>(),
    PointPtrType Point7 = GeneratePoint<PointType>(),
    PointPtrType Point8 = GeneratePoint<PointType>(),
    PointPtrType Point9 = GeneratePoint<PointType>(),
    PointPtrType Point10 = GeneratePoint<PointType>(),
    PointPtrType Point11 = GeneratePoint<PointType>(),
    PointPtrType Point12 = GeneratePoint<PointType>(),
    PointPtrType Point13 = GeneratePoint<PointType>()) {
    return BaseGeometryPtrType(new Pyramid3D13GeometryPtrType(
        Point1,
        Point2,
        Point3,
        Point4,
        Point5,
        Point6,
        Point7,
        Point8,
        Point9,
        Point10,
        Point11,
        Point12,
        Point13));
}

/** Generates a sample Pyramid3D5.
 * Generates a trirectangular pyramid on the origin with positive volume and side 1.
 * @return Pointer to a Pyramid3D13
 */
BaseGeometryPtrType GenerateRegularPyramid3D13() {
    return BaseGeometryPtrType(new Pyramid3D13GeometryPtrType(
        GeneratePoint<PointType>(0.0, 0.0, 0.0),
        GeneratePoint<PointType>(1.0, 0.0, 0.0),
        GeneratePoint<PointType>(0.0, 1.0, 0.0),
        GeneratePoint<PointType>(0.0, 0.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0),
        GeneratePoint<PointType>(0.0, 1.0, 1.0)
    ));
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13EdgesNumber, KratosCoreGeometriesFastSuite) {
    auto geomRegular = GenerateRegularPyramid3D13();

    KRATOS_CHECK_EQUAL(geomRegular->EdgesNumber(), 8);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13FacesNumber, KratosCoreGeometriesFastSuite) {
    auto geomRegular = GenerateRegularPyramid3D13();

    KRATOS_CHECK_EQUAL(geomRegular->FacesNumber(), 5);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13Volume, KratosCoreGeometriesFastSuite) {
    auto geomRegular = GenerateRegularPyramid3D13();

    KRATOS_CHECK_NEAR(geomRegular->Volume(),  0.5, TOLERANCE);
}

/** Checks the inside test for a given point respect to the pyramid
* Checks the inside test for a given point respect to the pyramid
* It performs 4 tests:
* A Point inside the pyramid: Expected result TRUE
* A Point outside the pyramid: Expected result FALSE
* A Point over a vertex of the pyramid: Expected result TRUE
* A Point over an edge of the pyramid: Expected result TRUE
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13IsInside, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPyramid3D13();

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
* pyramid. The centre of the pyramid is selected due to its known
* solution.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPyramid3D13();

    // Compute the global coordinates of the centre
    auto points = geom->Points();
    auto centre = Point{points[0] + points[1] + points[2] + points[3] + points[4] + points[5]};
    centre /= 6.0; // TODO use Center()

    // Compute the centre local coordinates
    array_1d<double, 3> centre_local_coords;
    geom->PointLocalCoordinates(centre_local_coords, centre);

    KRATOS_CHECK_NEAR(centre_local_coords(0), 1.0/3.0, TOLERANCE);
    KRATOS_CHECK_NEAR(centre_local_coords(1), 1.0/3.0, TOLERANCE);
    KRATOS_CHECK_NEAR(centre_local_coords(2), 1.0/2.0, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPyramid3D13();
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

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPyramid3D13();
    Matrix gradient;

    // Compute the global coordinates of the centre
    auto points = geom->Points();
    Point centre = Point{points[0] + points[1] + points[2] + points[3] + points[4] + points[5]};
    centre /= 6.0; // TODO use Center()

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
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GaussPoint1, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPyramid3D13();

    const double expected_vol = 0.5;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_1), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_1);
}

/** Tests the area using 'GI_GAUSS_2' integration method.
* Tests the area using 'GI_GAUSS_2' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GaussPoint2, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPyramid3D13();

    const double expected_vol = 0.5;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_2), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_2);
}

/** Tests the area using 'GI_GAUSS_3' integration method.
* Tests the area using 'GI_GAUSS_3' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GaussPoint3, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPyramid3D13();

    const double expected_vol = 0.5;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_3), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_3);
}

/** Tests the area using 'GI_GAUSS_4' integration method.
* Tests the area using 'GI_GAUSS_4' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GaussPoint4, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPyramid3D13();

    const double expected_vol = 0.5;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_4), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_4);
}

/** Tests the area using 'GI_GAUSS_5' integration method.
* Tests the area using 'GI_GAUSS_5' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GaussPoint5, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateRegularPyramid3D13();

    const double expected_vol = 0.5;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_5), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_5);
}

} // namespace Testing
} // namespace Kratos.
