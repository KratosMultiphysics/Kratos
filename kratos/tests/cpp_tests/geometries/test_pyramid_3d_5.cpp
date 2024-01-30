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


namespace Kratos::Testing {

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

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GetGeometryType, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D5();

    KRATOS_EXPECT_EQ(geomRegular->GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Pyramid);
    KRATOS_EXPECT_EQ(geomRegular->GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Pyramid3D5);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5EdgesNumber, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D5();

    KRATOS_EXPECT_EQ(geomRegular->EdgesNumber(), 8);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GenerateEdges, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();
    auto edges = geom->GenerateEdges();

    // Edge 1
    KRATOS_EXPECT_EQ(edges[0].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(edges[0].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Line3D2);
    KRATOS_EXPECT_VECTOR_EQ(edges[0].GetPoint( 0 ).Coordinates(), geom->GetPoint( 0 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(edges[0].GetPoint( 1 ).Coordinates(), geom->GetPoint( 1 ).Coordinates());

    // Edge 2
    KRATOS_EXPECT_EQ(edges[1].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(edges[1].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Line3D2);
    KRATOS_EXPECT_VECTOR_EQ(edges[1].GetPoint( 0 ).Coordinates(), geom->GetPoint( 1 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(edges[1].GetPoint( 1 ).Coordinates(), geom->GetPoint( 2 ).Coordinates());

    // Edge 3
    KRATOS_EXPECT_EQ(edges[2].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(edges[2].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Line3D2);
    KRATOS_EXPECT_VECTOR_EQ(edges[2].GetPoint( 0 ).Coordinates(), geom->GetPoint( 2 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(edges[2].GetPoint( 1 ).Coordinates(), geom->GetPoint( 3 ).Coordinates());

    // Edge 4
    KRATOS_EXPECT_EQ(edges[3].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(edges[3].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Line3D2);
    KRATOS_EXPECT_VECTOR_EQ(edges[3].GetPoint( 0 ).Coordinates(), geom->GetPoint( 3 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(edges[3].GetPoint( 1 ).Coordinates(), geom->GetPoint( 0 ).Coordinates());

    // Edge 5
    KRATOS_EXPECT_EQ(edges[4].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(edges[4].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Line3D2);
    KRATOS_EXPECT_VECTOR_EQ(edges[4].GetPoint( 0 ).Coordinates(), geom->GetPoint( 0 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(edges[4].GetPoint( 1 ).Coordinates(), geom->GetPoint( 4 ).Coordinates());

    // Edge 6
    KRATOS_EXPECT_EQ(edges[5].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(edges[5].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Line3D2);
    KRATOS_EXPECT_VECTOR_EQ(edges[5].GetPoint( 0 ).Coordinates(), geom->GetPoint( 1 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(edges[5].GetPoint( 1 ).Coordinates(), geom->GetPoint( 4 ).Coordinates());

    // Edge 7
    KRATOS_EXPECT_EQ(edges[6].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(edges[6].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Line3D2);
    KRATOS_EXPECT_VECTOR_EQ(edges[6].GetPoint( 0 ).Coordinates(), geom->GetPoint( 2 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(edges[6].GetPoint( 1 ).Coordinates(), geom->GetPoint( 4 ).Coordinates());

    // Edge 8
    KRATOS_EXPECT_EQ(edges[7].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Linear);
    KRATOS_EXPECT_EQ(edges[7].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Line3D2);
    KRATOS_EXPECT_VECTOR_EQ(edges[7].GetPoint( 0 ).Coordinates(), geom->GetPoint( 3 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(edges[7].GetPoint( 1 ).Coordinates(), geom->GetPoint( 4 ).Coordinates());
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5FacesNumber, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D5();

    KRATOS_EXPECT_EQ(geomRegular->FacesNumber(), 5);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GenerateFaces, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();
    auto faces = geom->GenerateFaces();

    // Face 1
    KRATOS_EXPECT_EQ(faces[0].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Triangle);
    KRATOS_EXPECT_EQ(faces[0].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Triangle3D3);
    KRATOS_EXPECT_VECTOR_EQ(faces[0].GetPoint( 0 ).Coordinates(), geom->GetPoint( 0 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[0].GetPoint( 1 ).Coordinates(), geom->GetPoint( 1 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[0].GetPoint( 2 ).Coordinates(), geom->GetPoint( 4 ).Coordinates());

    // Face 2
    KRATOS_EXPECT_EQ(faces[1].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Triangle);
    KRATOS_EXPECT_EQ(faces[1].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Triangle3D3);
    KRATOS_EXPECT_VECTOR_EQ(faces[1].GetPoint( 0 ).Coordinates(), geom->GetPoint( 1 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[1].GetPoint( 1 ).Coordinates(), geom->GetPoint( 2 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[1].GetPoint( 2 ).Coordinates(), geom->GetPoint( 4 ).Coordinates());

    // Face 3
    KRATOS_EXPECT_EQ(faces[2].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Quadrilateral);
    KRATOS_EXPECT_EQ(faces[2].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4);
    KRATOS_EXPECT_VECTOR_EQ(faces[2].GetPoint( 0 ).Coordinates(), geom->GetPoint( 0 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[2].GetPoint( 1 ).Coordinates(), geom->GetPoint( 1 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[2].GetPoint( 2 ).Coordinates(), geom->GetPoint( 2 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[2].GetPoint( 3 ).Coordinates(), geom->GetPoint( 3 ).Coordinates());

    // Face 4
    KRATOS_EXPECT_EQ(faces[3].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Triangle);
    KRATOS_EXPECT_EQ(faces[3].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Triangle3D3);
    KRATOS_EXPECT_VECTOR_EQ(faces[3].GetPoint( 0 ).Coordinates(), geom->GetPoint( 2 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[3].GetPoint( 1 ).Coordinates(), geom->GetPoint( 3 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[3].GetPoint( 2 ).Coordinates(), geom->GetPoint( 4 ).Coordinates());

    // Face 5
    KRATOS_EXPECT_EQ(faces[4].GetGeometryFamily(), GeometryData::KratosGeometryFamily::Kratos_Triangle);
    KRATOS_EXPECT_EQ(faces[4].GetGeometryType(), GeometryData::KratosGeometryType::Kratos_Triangle3D3);
    KRATOS_EXPECT_VECTOR_EQ(faces[4].GetPoint( 0 ).Coordinates(), geom->GetPoint( 3 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[4].GetPoint( 1 ).Coordinates(), geom->GetPoint( 0 ).Coordinates());
    KRATOS_EXPECT_VECTOR_EQ(faces[4].GetPoint( 2 ).Coordinates(), geom->GetPoint( 4 ).Coordinates());
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5Volume, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D5();

    KRATOS_EXPECT_NEAR(geomRegular->Volume(), 2.0, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5Center, KratosCoreGeometriesFastSuite)
{
    auto geomRegular = GenerateRegularPyramid3D5();

    array_1d<double, 3> center{0, 0, 0.3};

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
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5IsInside, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

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
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5PointLocalCoordinates, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    // Compute the global coordinates of the centre
    auto points = geom->Points();
    auto centre = geom->Center();

    // Compute the centre local coordinates
    array_1d<double, 3> centre_local_coords;
    geom->PointLocalCoordinates(centre_local_coords, centre);

    KRATOS_EXPECT_NEAR(centre_local_coords(0), 0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(centre_local_coords(1), 0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(centre_local_coords(2), -0.6, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5ShapeFunctionsValues, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();
    array_1d<double, 3> coord(3);
    coord[0] = 1.0 / 2.0;
    coord[1] = 1.0 / 4.0;
    coord[2] = 1.0 / 16.0;
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(0, coord), 0.0439453, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(1, coord), 0.131836, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(2, coord), 0.219727, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(3, coord), 0.0732422, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(4, coord), 0.53125, TOLERANCE);
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

    KRATOS_EXPECT_NEAR(gradient(0,0), -0.2, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(0,1), -0.2, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(0,2), -0.125, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(1,0), +0.2, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(1,1), -0.2, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(1,2), -0.125, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(2,0), +0.2, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(2,1), +0.2, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(2,2), -0.125, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(3,0), -0.2, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(3,1), +0.2, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(3,2), -0.125, TOLERANCE);

    KRATOS_EXPECT_NEAR(gradient(4,0), 0.00, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(4,1), 0.00, TOLERANCE);
    KRATOS_EXPECT_NEAR(gradient(4,2), +0.5, TOLERANCE);
}

/** Tests the area using 'GI_GAUSS_1' integration method.
* Tests the area using 'GI_GAUSS_1' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GaussPoint1, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    const double expected_vol = 2.0;

    KRATOS_EXPECT_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_1), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_1);
}

/** Tests the area using 'GI_GAUSS_2' integration method.
* Tests the area using 'GI_GAUSS_2' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GaussPoint2, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    const double expected_vol = 2.0;

    KRATOS_EXPECT_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_2), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_2);
}

/** Tests the area using 'GI_GAUSS_3' integration method.
* Tests the area using 'GI_GAUSS_3' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GaussPoint3, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    const double expected_vol = 2.0;

    KRATOS_EXPECT_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_3), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_3);
}

/** Tests the area using 'GI_GAUSS_4' integration method.
* Tests the area using 'GI_GAUSS_4' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GaussPoint4, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    const double expected_vol = 2.0;

    KRATOS_EXPECT_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_4), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_4);
}

/** Tests the area using 'GI_GAUSS_5' integration method.
* Tests the area using 'GI_GAUSS_5' integration method.
*/
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5GaussPoint5, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    const double expected_vol = 2.0;

    KRATOS_EXPECT_NEAR(CalculateAreaByIntegration(*geom, GeometryData::IntegrationMethod::GI_GAUSS_5), expected_vol, TOLERANCE);
    VerifyStrainExactness(*geom, GeometryData::IntegrationMethod::GI_GAUSS_5);
}

/** Checks if CalculateDistance is correct.
 * Checks if CalculateDistance is correct.
 */
KRATOS_TEST_CASE_IN_SUITE(Pyramid3D5CalculateDistance, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateRegularPyramid3D5();

    Point point1(0.0, 0.0, 0.5);
    KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point1), 0.0);
    KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point1), GeometryUtils::CalculateDistanceFrom3DGeometry(*geom, point1));

    Point point2(0.0, 0.0, -0.5);
    KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point2), 0.5);
    KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point2), GeometryUtils::CalculateDistanceFrom3DGeometry(*geom, point2));
}

} // namespace Kratos::Testing.
