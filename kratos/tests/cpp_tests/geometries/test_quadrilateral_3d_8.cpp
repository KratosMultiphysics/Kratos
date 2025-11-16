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

#include "testing/testing.h"
#include "geometries/quadrilateral_3d_8.h"
#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos::Testing {

namespace {
/// Test utility functions

/** Generates a sample Quadrilateral3D8.
* Generates a right quadrilateral with origin in the origin and leg size 1.
* @return  Pointer to a Quadrilateral3D8
*/
template<class TPointType>
typename Quadrilateral3D8<TPointType>::Pointer GenerateFlatQuadrilateral3D8()
{
  return Kratos::make_shared<Quadrilateral3D8<TPointType>>(
  GeneratePoint<TPointType>( 0.0, 0.0, 0.0),
  GeneratePoint<TPointType>( 1.0, 0.0, 0.0),
  GeneratePoint<TPointType>( 1.0, 1.0, 0.0),
  GeneratePoint<TPointType>( 0.0, 1.0, 0.0),
  GeneratePoint<TPointType>( 0.5, 0.0, 0.0),
  GeneratePoint<TPointType>( 1.0, 0.5, 0.0),
  GeneratePoint<TPointType>( 0.5, 1.0, 0.0),
  GeneratePoint<TPointType>( 0.0, 0.5, 0.0)
  );
}

}

/** Checks if the number of faces is correct.
*/
KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D8FacesNumber, KratosCoreGeometriesFastSuite)
{
  auto p_geom = GenerateFlatQuadrilateral3D8<Node>();

  KRATOS_EXPECT_EQ(p_geom->FacesNumber(), 1);
}

/** Checks if the faces are correct.
*/
KRATOS_TEST_CASE_IN_SUITE(Quadrilateral3D8Faces, KratosCoreGeometriesFastSuite) {
  auto p_geom = GenerateFlatQuadrilateral3D8<Node>();

  const auto& r_faces = p_geom->GenerateFaces();
  ASSERT_EQ(r_faces.size(), 1);
  for (std::size_t i = 0; i < r_faces.front().PointsNumber(); ++i) {
    KRATOS_EXPECT_NEAR(r_faces.front()[i].X(), (p_geom->pGetPoint(i))->X(), TOLERANCE);
    KRATOS_EXPECT_NEAR(r_faces.front()[i].Y(), (p_geom->pGetPoint(i))->Y(), TOLERANCE);
    KRATOS_EXPECT_NEAR(r_faces.front()[i].Z(), (p_geom->pGetPoint(i))->Z(), TOLERANCE);
  }
}

}
