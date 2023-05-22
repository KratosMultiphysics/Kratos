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

  KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 12);
}

/** Checks if the number of faces is correct.
 * Checks if the number of faces is correct.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20FacesNumber, KratosCoreGeometriesFastSuite)
{
  auto geom = GenerateCanonicalHexahedra3D20();

  KRATOS_CHECK_EQUAL(geom->FacesNumber(), 6);
}

/** Checks if the characteristic length of the hexahedra is calculated correctly.
 * Checks if the characteristic length of the hexahedra is calculated correctly.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20Length, KratosCoreGeometriesFastSuite)
{
  auto geom = GenerateCanonicalHexahedra3D20();

  KRATOS_CHECK_NEAR(geom->Length(), 0.353553, TOLERANCE);
}

/** Checks if the area of the hexahedra is calculated correctly.
 * Checks if the area of the hexahedra is calculated correctly.
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20Area, KratosCoreGeometriesFastSuite)
{
  auto geom = GenerateCanonicalHexahedra3D20();

  KRATOS_CHECK_NEAR(geom->Area(), 1.0, TOLERANCE);
}

/** Checks if the volume of the hexahedra is calculated correctly.
 * Checks if the volume of the hexahedra is calculated correctly.
 * For hexahedra 3D8 'volume()' call defaults to 'area()'
 */
KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20Volume, KratosCoreGeometriesFastSuite)
{
  auto geom = GenerateCanonicalHexahedra3D20();

  KRATOS_CHECK_NEAR(geom->Volume(), 1.0, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(Hexahedra3D20ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite)
{
    auto geom = GenerateCanonicalHexahedra3D20();
    TestAllShapeFunctionsLocalGradients(*geom);
}

}  // namespace Kratos::Testing.
