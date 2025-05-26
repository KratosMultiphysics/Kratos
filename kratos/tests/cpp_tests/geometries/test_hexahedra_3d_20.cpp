//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/hexahedra_3d_20.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos::Testing
{
using PointType = Node;
using PointPtrType = Node::Pointer;
using Hexa20GeometryType = Hexahedra3D20<PointType>;
using Hexa20GeometryPtrType = Hexa20GeometryType::Pointer;

/** Generates a sample Hexahedra3D20.
 * Generates a hexahedra defined by eight random points in the space.
 * @return  Pointer to a Hexahedra3D20
 */
Hexa20GeometryPtrType GenerateHexahedra3D20(
    PointPtrType PointA = GeneratePoint<PointType>(),
    PointPtrType PointB = GeneratePoint<PointType>(),
    PointPtrType PointC = GeneratePoint<PointType>(),
    PointPtrType PointD = GeneratePoint<PointType>(),
    PointPtrType PointE = GeneratePoint<PointType>(),
    PointPtrType PointF = GeneratePoint<PointType>(),
    PointPtrType PointG = GeneratePoint<PointType>(),
    PointPtrType PointH = GeneratePoint<PointType>(),
    PointPtrType PointI = GeneratePoint<PointType>(),
    PointPtrType PointJ = GeneratePoint<PointType>(),
    PointPtrType PointK = GeneratePoint<PointType>(),
    PointPtrType PointM = GeneratePoint<PointType>(),
    PointPtrType PointN = GeneratePoint<PointType>(),
    PointPtrType PointO = GeneratePoint<PointType>(),
    PointPtrType PointP = GeneratePoint<PointType>(),
    PointPtrType PointQ = GeneratePoint<PointType>(),
    PointPtrType PointR = GeneratePoint<PointType>(),
    PointPtrType PointS = GeneratePoint<PointType>(),
    PointPtrType PointT = GeneratePoint<PointType>(),
    PointPtrType PointU = GeneratePoint<PointType>()) {
  return Hexa20GeometryPtrType(new Hexa20GeometryType(PointA, PointB, PointC, PointD, PointE, PointF, PointG, PointH, PointI, PointJ, PointK, PointM, PointN, PointO, PointP, PointQ, PointR, PointS, PointT, PointU));
}

/** Generates a sample Hexahedra3D20.
 * Generates a hexahedra with center on the origin with positive volume and side 1.
 * @return  Pointer to a Hexahedra3D20
 */
Hexa20GeometryPtrType GenerateCanonicalHexahedra3D20()
{
  auto p1 = GeneratePoint<PointType>(0.0,1.0,1.0);
  auto p2 = GeneratePoint<PointType>(0.0,1.0,0.5);
  auto p3 = GeneratePoint<PointType>(0.5,1.0,1.0);
  auto p4 = GeneratePoint<PointType>(0.0,0.5,1.0);
  auto p5 = GeneratePoint<PointType>(0.0,0.0,1.0);
  auto p6 = GeneratePoint<PointType>(1.0,1.0,1.0);
  auto p7 = GeneratePoint<PointType>(0.0,1.0,0.0);
  auto p8 = GeneratePoint<PointType>(1.0,0.5,1.0);
  auto p9 = GeneratePoint<PointType>(1.0,1.0,0.5);
  auto p10 = GeneratePoint<PointType>(0.0,0.0,0.5);
  auto p11 = GeneratePoint<PointType>(0.5,1.0,0.0);
  auto p12 = GeneratePoint<PointType>(0.0,0.5,0.0);
  auto p13 = GeneratePoint<PointType>(0.5,0.0,1.0);
  auto p14 = GeneratePoint<PointType>(1.0,1.0,0.0);
  auto p15 = GeneratePoint<PointType>(1.0,0.0,1.0);
  auto p16 = GeneratePoint<PointType>(0.0,0.0,0.0);
  auto p17 = GeneratePoint<PointType>(0.5,0.0,0.0);
  auto p18 = GeneratePoint<PointType>(1.0,0.5,0.0);
  auto p19 = GeneratePoint<PointType>(1.0,0.0,0.5);
  auto p20 = GeneratePoint<PointType>(1.0,0.0,0.0);
  return Hexa20GeometryPtrType(new Hexa20GeometryType(p6, p1, p7, p14, p15, p5, p16, p20, p3, p2, p11, p9, p8, p4, p12, p18, p13, p10, p17, p19));
}

/** Checks if the number of edges is correct.
 * Checks if the number of edges is correct.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20EdgesNumber, KratosCoreGeometriesFastSuite)
{
  auto geom = GenerateCanonicalHexahedra3D20();

  KRATOS_EXPECT_EQ(geom->EdgesNumber(), 12);
}

/** Checks if the number of faces is correct.
 * Checks if the number of faces is correct.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20FacesNumber, KratosCoreGeometriesFastSuite)
{
  auto geom = GenerateCanonicalHexahedra3D20();

  KRATOS_EXPECT_EQ(geom->FacesNumber(), 6);
}

/** Checks if the characteristic length of the hexahedra is calculated correctly.
 * Checks if the characteristic length of the hexahedra is calculated correctly.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20Length, KratosCoreGeometriesFastSuite)
{
  auto geom = GenerateCanonicalHexahedra3D20();

  KRATOS_EXPECT_NEAR(geom->Length(), 0.353553, TOLERANCE);
}

/** Checks if the area of the hexahedra is calculated correctly.
 * Checks if the area of the hexahedra is calculated correctly.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20Area, KratosCoreGeometriesFastSuite)
{
  auto geom = GenerateCanonicalHexahedra3D20();

  KRATOS_EXPECT_NEAR(geom->Area(), 1.0, TOLERANCE);
}

/** Checks if the volume of the hexahedra is calculated correctly.
 * Checks if the volume of the hexahedra is calculated correctly.
 * For hexahedra 3D8 'volume()' call defaults to 'area()'
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20Volume, KratosCoreGeometriesFastSuite)
{
  auto geom = GenerateCanonicalHexahedra3D20();

  KRATOS_EXPECT_NEAR(geom->Volume(), 1.0, TOLERANCE);
}

/** Checks the values of the shape functions at a local coordinate */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20ShapeFunctionsValues, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D20();
    array_1d<double, 3> coord(3);
    coord[0] = 1.0 / 2.0;
    coord[1] = 1.0 / 4.0;
    coord[2] = 1.0 / 16.0;

    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(0, coord), -0.12359619140625, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(1, coord), -0.23895263671875, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(2, coord), -0.28839111328125, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(3, coord), -0.16937255859375, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(4, coord), -0.13385009765625, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(5, coord), -0.25213623046875, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(6, coord), -0.29571533203125, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(7, coord), -0.18157958984375, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(8, coord), 0.1318359375, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(9, coord), 0.32958984375, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(10, coord), 0.2197265625, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(11, coord), 0.10986328125, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(12, coord), 0.0933837890625, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(13, coord), 0.2801513671875, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(14, coord), 0.4669189453125, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(15, coord), 0.1556396484375, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(16, coord), 0.1494140625, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(17, coord), 0.37353515625, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(18, coord), 0.2490234375, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(19, coord), 0.12451171875, TOLERANCE);

    CrossCheckShapeFunctionsValues(*geom);
}

KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D20();
    TestAllShapeFunctionsLocalGradients(*geom);
}

/** Checks if the coordinates of the generate faces are correct
*/
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20GenerateFaces, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D20();

    auto faces = geom->GenerateFaces();

    // create array of expected coordinates per face
    std::vector<std::vector<array_1d<double, 3>>> all_expected_face_coordinates;
    all_expected_face_coordinates.push_back(
        {array_1d<double, 3>{1.0, 1.0, 0.0}, array_1d<double, 3>{0.0, 1.0, 0.0},
         array_1d<double, 3>{0.0, 1.0, 1.0}, array_1d<double, 3>{1.0, 1.0, 1.0},
         array_1d<double, 3>{0.5, 1.0, 0.0}, array_1d<double, 3>{0.0, 1.0, 0.5},
         array_1d<double, 3>{0.5, 1.0, 1.0}, array_1d<double, 3>{1.0, 1.0, 0.5}});
    all_expected_face_coordinates.push_back(
        {array_1d<double, 3>{1.0, 1.0, 1.0}, array_1d<double, 3>{0.0, 1.0, 1.0},
         array_1d<double, 3>{0.0, 0.0, 1.0}, array_1d<double, 3>{1.0, 0.0, 1.0},
         array_1d<double, 3>{0.5, 1.0, 1.0}, array_1d<double, 3>{0.0, 0.5, 1.0},
         array_1d<double, 3>{0.5, 0.0, 1.0}, array_1d<double, 3>{1.0, 0.5, 1.0}});

    all_expected_face_coordinates.push_back(
        {array_1d<double, 3>{0.0, 1.0, 0.0}, array_1d<double, 3>{0.0, 0.0, 0.0},
         array_1d<double, 3>{0.0, 0.0, 1.0}, array_1d<double, 3>{0.0, 1.0, 1.0},
         array_1d<double, 3>{0.0, 0.5, 0.0}, array_1d<double, 3>{0.0, 0.0, 0.5},
         array_1d<double, 3>{0.0, 0.5, 1.0}, array_1d<double, 3>{0.0, 1.0, 0.5}});

    all_expected_face_coordinates.push_back(
        {array_1d<double, 3>{1.0, 0.0, 0.0}, array_1d<double, 3>{0.0, 0.0, 0.0},
         array_1d<double, 3>{0.0, 1.0, 0.0}, array_1d<double, 3>{1.0, 1.0, 0.0},
         array_1d<double, 3>{0.5, 0.0, 0.0}, array_1d<double, 3>{0.0, 0.5, 0.0},
         array_1d<double, 3>{0.5, 1.0, 0.0}, array_1d<double, 3>{1.0, 0.5, 0.0}});

    all_expected_face_coordinates.push_back(
        {array_1d<double, 3>{1.0, 0.0, 0.0}, array_1d<double, 3>{1.0, 1.0, 0.0},
         array_1d<double, 3>{1.0, 1.0, 1.0}, array_1d<double, 3>{1.0, 0.0, 1.0},
         array_1d<double, 3>{1.0, 0.5, 0.0}, array_1d<double, 3>{1.0, 1.0, 0.5},
         array_1d<double, 3>{1.0, 0.5, 1.0}, array_1d<double, 3>{1.0, 0.0, 0.5}});

    all_expected_face_coordinates.push_back(
        {array_1d<double, 3>{1.0, 0.0, 1.0}, array_1d<double, 3>{0.0, 0.0, 1.0},
         array_1d<double, 3>{0.0, 0.0, 0.0}, array_1d<double, 3>{1.0, 0.0, 0.0},
         array_1d<double, 3>{0.5, 0.0, 1.0}, array_1d<double, 3>{0.0, 0.0, 0.5},
         array_1d<double, 3>{0.5, 0.0, 0.0}, array_1d<double, 3>{1.0, 0.0, 0.5}});

    // check every coordinate of every generate face
    for (unsigned int i = 0; i < faces.size(); i++) {
        // loop over points of face
        std::vector<array_1d<double, 3>> face_coordinates;
        for (unsigned int j = 0; j < faces[i].Points().size(); j++) {
            KRATOS_EXPECT_VECTOR_NEAR(faces[i].Points()[j].Coordinates(),
                                      all_expected_face_coordinates[i][j], TOLERANCE);
        }
    }
}

} // namespace Kratos::Testing.
