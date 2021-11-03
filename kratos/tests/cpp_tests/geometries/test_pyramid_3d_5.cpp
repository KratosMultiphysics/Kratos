//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//	                 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Ashish Darekar
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/pyramid_3d_5.h"
#include "tests/cpp_tests/geometries/test_geometry.h"


namespace Kratos {
namespace Testing {

typedef GeometryType::Pointer            BaseGeometryPtrType;
typedef Pyramid3D5<NodeType>             Pyramid3D5GeometryType;

/** Generates a sample Pyramid3D5.
 * Generates a rectangular pyramid on the origin with positive volume and side 1.
 * @return Pointer to a Pyramid3D5
 */
BaseGeometryPtrType GenerateRegularPyramid3D5() {
    return BaseGeometryPtrType(new Pyramid3D5GeometryType(
        GeneratePoint<NodeType>(-1.0,  1.0, 0.0),
        GeneratePoint<NodeType>(-1.0, -1.0, 0.0),
        GeneratePoint<NodeType>( 1.0, -1.0, 0.0),
        GeneratePoint<NodeType>( 1.0,  1.0, 0.0),
        GeneratePoint<NodeType>( 0.0,  0.0, 1.5)
    ));
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5EdgesNumber, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D5();

    KRATOS_CHECK_EQUAL(geomRegular->EdgesNumber(), 8);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5FacesNumber, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D5();

    KRATOS_CHECK_EQUAL(geomRegular->FacesNumber(), 5);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5Volume, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D5();

    KRATOS_CHECK_NEAR(geomRegular->Volume(), 2.0, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5Center, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D5();

    array_1d<double, 3> center{0, 0, 0.3};

    KRATOS_CHECK_VECTOR_NEAR(geomRegular->Center(), center, TOLERANCE);
}

/** Checks the inside test for a given point respect to the pyramid
* Checks the inside test for a given point respect to the pyramid
* It performs 4 tests:
* A Point inside the pyramid: Expected result TRUE
* A Point outside the pyramid: Expected result FALSE
* A Point over a vertex of the pyramid: Expected result TRUE
* A Point over an edge of the pyramid: Expected result TRUE
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5IsInside, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    Point PointInside(0.0, 0.0, 0.3);
    Point PointOutside(0.0, 0.0, 1.6);
    Point PointOutside2(0.1, 0.1, 1.4);
    Point PointInVertex(-1.0, 1.0, 0.0);
    Point PointInEdge(-1.0, 0.0, 0.0);

    Point LocalCoords;

    KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
    KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
    KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside2, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
}

/** Checks the point local coordinates for a given point respect to the
* pyramid. The centre of the pyramid is selected due to its known
* solution.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5PointLocalCoordinates, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    // Compute the global coordinates of the centre
    auto points = geom->Points();
    auto centre = geom->Center();

    // Compute the centre local coordinates
    array_1d<double, 3> centre_local_coords;
    geom->PointLocalCoordinates(centre_local_coords, centre);

    KRATOS_CHECK_NEAR(centre_local_coords(0), 0.0, TOLERANCE);
    KRATOS_CHECK_NEAR(centre_local_coords(1), 0.0, TOLERANCE);
    KRATOS_CHECK_NEAR(centre_local_coords(2), -0.6, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5ShapeFunctionsValues, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();
    array_1d<double, 3> coord(3);
    coord[0] = 1.0 / 2.0;
    coord[1] = 1.0 / 4.0;
    coord[2] = 1.0 / 16.0;
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), 0.0439453, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 0.131836, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(2, coord), 0.219727, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(3, coord), 0.0732422, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(4, coord), 0.53125, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();
    Matrix gradient;

    // Compute the global coordinates of the centre
    auto points = geom->Points();
    auto centre = geom->Center();

    // Compute the centre local coordinates
    array_1d<double, 3> centre_local_coords;
    geom->PointLocalCoordinates(centre_local_coords, centre);

    gradient = geom->ShapeFunctionsLocalGradients(gradient, centre_local_coords);

    KRATOS_CHECK_NEAR(gradient(0,0), -0.2, TOLERANCE);
    KRATOS_CHECK_NEAR(gradient(0,1), -0.2, TOLERANCE);
    KRATOS_CHECK_NEAR(gradient(0,2), -0.125, TOLERANCE);

    KRATOS_CHECK_NEAR(gradient(1,0), +0.2, TOLERANCE);
    KRATOS_CHECK_NEAR(gradient(1,1), -0.2, TOLERANCE);
    KRATOS_CHECK_NEAR(gradient(1,2), -0.125, TOLERANCE);

    KRATOS_CHECK_NEAR(gradient(2,0), +0.2, TOLERANCE);
    KRATOS_CHECK_NEAR(gradient(2,1), +0.2, TOLERANCE);
    KRATOS_CHECK_NEAR(gradient(2,2), -0.125, TOLERANCE);

    KRATOS_CHECK_NEAR(gradient(3,0), -0.2, TOLERANCE);
    KRATOS_CHECK_NEAR(gradient(3,1), +0.2, TOLERANCE);
    KRATOS_CHECK_NEAR(gradient(3,2), -0.125, TOLERANCE);

    KRATOS_CHECK_NEAR(gradient(4,0), 0.00, TOLERANCE);
    KRATOS_CHECK_NEAR(gradient(4,1), 0.00, TOLERANCE);
    KRATOS_CHECK_NEAR(gradient(4,2), +0.5, TOLERANCE);
}

/** Tests the area using 'GI_GAUSS_1' integration method.
* Tests the area using 'GI_GAUSS_1' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GaussPoint1, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    const double expected_vol = 2.0;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_1), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_1);
}

/** Tests the area using 'GI_GAUSS_2' integration method.
* Tests the area using 'GI_GAUSS_2' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GaussPoint2, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    const double expected_vol = 2.0;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_2), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_2);
}

/** Tests the area using 'GI_GAUSS_3' integration method.
* Tests the area using 'GI_GAUSS_3' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GaussPoint3, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    const double expected_vol = 2.0;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_3), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_3);
}

/** Tests the area using 'GI_GAUSS_4' integration method.
* Tests the area using 'GI_GAUSS_4' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GaussPoint4, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    const double expected_vol = 2.0;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_4), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_4);
}

/** Tests the area using 'GI_GAUSS_5' integration method.
* Tests the area using 'GI_GAUSS_5' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GaussPoint5, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    const double expected_vol = 2.0;

    KRATOS_CHECK_NEAR(CalculateAreaByIntegration(*geom, GeometryData::GI_GAUSS_5), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::GI_GAUSS_5);
}

} // namespace Testing
} // namespace Kratos.
