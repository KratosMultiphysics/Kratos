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
#include "geometries/pyramid_3d_13.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos::Testing {

typedef GeometryType::Pointer             BaseGeometryPtrType;
typedef Pyramid3D13<NodeType>             Pyramid3D13GeometryType;

/** Generates a sample Pyramid3D13.
 * Generates a trirectangular pyramid on the origin with positive volume and side 1.
 * @return Pointer to a Pyramid3D13
 */
BaseGeometryPtrType GenerateRegularPyramid3D13() {
    return BaseGeometryPtrType(new Pyramid3D13GeometryType(
        GeneratePoint<NodeType>(-1.0,  1.0, 0.0),
        GeneratePoint<NodeType>(-1.0, -1.0, 0.0),
        GeneratePoint<NodeType>( 1.0, -1.0, 0.0),
        GeneratePoint<NodeType>( 1.0,  1.0, 0.0),
        GeneratePoint<NodeType>( 0.0,  0.0, 1.5),
        GeneratePoint<NodeType>(-1.0,  0.0, 0.0),
        GeneratePoint<NodeType>( 0.0, -1.0, 0.0),
        GeneratePoint<NodeType>( 1.0,  0.0, 0.0),
        GeneratePoint<NodeType>( 0.0,  1.0, 0.0),
        GeneratePoint<NodeType>(-0.5,  0.5, 0.75),
        GeneratePoint<NodeType>(-0.5, -0.5, 0.75),
        GeneratePoint<NodeType>( 0.5, -0.5, 0.75),
        GeneratePoint<NodeType>( 0.5,  0.5, 0.75)
    ));
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GetGeometryType, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D13();

    KRATOS_EXPECT_EQ(geomRegular->GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Pyramid);
    KRATOS_EXPECT_EQ(geomRegular->GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Pyramid3D13);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13EdgesNumber, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D13();

    KRATOS_EXPECT_EQ(geomRegular->EdgesNumber(), 8);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13FacesNumber, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D13();

    KRATOS_EXPECT_EQ(geomRegular->FacesNumber(), 5);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13Volume, KratosCoreGeometriesFastSuite)
{
    //KRATOS_SKIP_TEST << "NOT IMPLEMENTED!";
    auto geomRegular = GenerateRegularPyramid3D13();

    KRATOS_EXPECT_NEAR(geomRegular->Volume(), 2.0, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13Center, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D13();

    array_1d<double, 3> center{0, 0, 4.5/13};

    KRATOS_EXPECT_VECTOR_NEAR(geomRegular->Center(), center, TOLERANCE);
}

/** Checks the inside test for a given point respect to the pyramid
* Checks the inside test for a given point respect to the pyramid
* It performs 4 tests:
* A Point inside the pyramid: Expected result TRUE
* A Point outside the pyramid: Expected result FALSE
* A Point over a vertex of the pyramid: Expected result TRUE
* A Point over an edge of the pyramid: Expected result TRUE
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13IsInside, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D13();

    Point PointInside(0.0, 0.0, 0.3);
    Point PointOutside(0.0, 0.0, 1.6);
    Point PointOutside2(0.1, 0.1, 1.4);
    Point PointInVertex(-1.0, 1.0, 0.0);
    Point PointInEdge(-1.0, 0.0, 0.0);

    Point LocalCoords;

    KRATOS_EXPECT_TRUE(geom->IsInside(PointInside, LocalCoords, EPSILON));
    KRATOS_EXPECT_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
    KRATOS_EXPECT_FALSE(geom->IsInside(PointOutside2, LocalCoords, EPSILON));
    KRATOS_EXPECT_TRUE(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
    KRATOS_EXPECT_TRUE(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
}

/** Checks the point local coordinates for a given point respect to the
* pyramid. The centre of the pyramid is selected due to its known
* solution.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13PointLocalCoordinates, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D13();

    // Compute the global coordinates of the centre
    auto points = geom->Points();
    auto centre = geom->Center();

    // Compute the centre local coordinates
    array_1d<double, 3> centre_local_coords;
    geom->PointLocalCoordinates(centre_local_coords, centre);

    KRATOS_EXPECT_NEAR(centre_local_coords(0), 0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(centre_local_coords(1), 0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(centre_local_coords(2), -0.538461538, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13ShapeFunctionsValues, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D13();
    array_1d<double, 3> coord(3);
    coord[0] = 1.0 / 2.0;
    coord[1] = 1.0 / 4.0;
    coord[2] = 1.0 / 16.0;
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(0, coord), -0.146942, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(1, coord), -0.203934, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(2, coord), -0.230026, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(3, coord), -0.169373, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(4, coord), 0.0332031, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(5, coord), 0.149345, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(6, coord), 0.242043, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(7, coord), 0.190544, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(8, coord), 0.139046, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(9, coord), 0.0933838, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(10, coord), 0.280151, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(11, coord), 0.466919, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(12, coord), 0.15564, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D13();
    Matrix gradient;

    // Compute the global coordinates of the centre
    auto points = geom->Points();
    Point centre = geom->Center();

    // Compute the centre local coordinates
    array_1d<double, 3> centre_local_coords;
    geom->PointLocalCoordinates(centre_local_coords, centre);

    gradient = geom->ShapeFunctionsLocalGradients(gradient, centre_local_coords);

    KRATOS_EXPECT_NEAR(gradient(0,0), +0.0443787, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(0,1), +0.0443787, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(0,2), -0.00961538, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(1,0), -0.0443787, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(1,1), +0.0443787, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(1,2), -0.00961538, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(2,0), -0.0443787, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(2,1), -0.0443787, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(2,2), -0.00961538, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(3,0), +0.0443787, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(3,1), -0.0443787, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(3,2), -0.00961538, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(4,0), 0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(4,1), 0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(4,2), -0.0384615, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(5,0),  0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(5,1), -0.295858, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(5,2), -0.25, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(6,0),  0.295858, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(6,1),  0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(6,2), -0.25, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(7,0),  0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(7,1), +0.295858, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(7,2), -0.25, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(8,0), -0.295858, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(8,1),  0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(8,2), -0.25, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(9,0), -0.177515, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(9,1), -0.177515, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(9,2), +0.269231, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(10,0), +0.177515, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(10,1), -0.177515, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(10,2), +0.269231, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(11,0), +0.177515, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(11,1), +0.177515, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(11,2), +0.269231, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(12,0), -0.177515, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(12,1), +0.177515, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(12,2), +0.269231, TOLERANCE);

}

/** Tests the area using 'GI_GAUSS_1' integration method.
* Tests the area using 'GI_GAUSS_1' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GaussPoint1, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D13();

    const double expected_vol = 2.0;

    KRATOS_EXPECT_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_1), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_1);
}

/** Tests the area using 'GI_GAUSS_2' integration method.
* Tests the area using 'GI_GAUSS_2' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GaussPoint2, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D13();

    const double expected_vol = 2.0;

    KRATOS_EXPECT_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_2), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_2);
}

/** Tests the area using 'GI_GAUSS_3' integration method.
* Tests the area using 'GI_GAUSS_3' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GaussPoint3, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D13();

    const double expected_vol = 2.0;

    KRATOS_EXPECT_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_3), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_3);
}

/** Tests the area using 'GI_GAUSS_4' integration method.
* Tests the area using 'GI_GAUSS_4' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GaussPoint4, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D13();

    const double expected_vol = 2.0;

    KRATOS_EXPECT_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_4), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_4);
}

/** Tests the area using 'GI_GAUSS_5' integration method.
* Tests the area using 'GI_GAUSS_5' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D13GaussPoint5, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D13();

    const double expected_vol = 2.0;

    KRATOS_EXPECT_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_5), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_5);
}

} // namespace Kratos::Testing.
