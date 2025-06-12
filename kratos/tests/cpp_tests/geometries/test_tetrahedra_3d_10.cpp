//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//                   Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/tetrahedra_3d_10.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos::Testing {

  // /// Factory functions

namespace{
  typedef Node NodeType;

  Geometry<NodeType>::PointsArrayType GenerateReferenceNodes3D10()
  {
      Geometry<NodeType>::PointsArrayType points;
      points.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(4, 0.0, 0.0, 1.0)));
      points.push_back(NodeType::Pointer(new NodeType(5, 0.5, 0.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(6, 0.5, 0.5, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(7, 0.0, 0.5, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(8, 0.0, 0.0, 0.5)));
      points.push_back(NodeType::Pointer(new NodeType(9, 0.5, 0.0, 0.5)));
      points.push_back(NodeType::Pointer(new NodeType(10, 0.0, 0.5, 0.5)));
      return points;
  }

  Geometry<NodeType>::Pointer GenerateReferenceTetrahedra3D10() {
      return Geometry<NodeType>::Pointer(
          new Tetrahedra3D10<NodeType>(GenerateReferenceNodes3D10()));
  }

  Geometry<NodeType>::Pointer GenerateCurvedTetrahedra3D10() {
      auto nodes = GenerateReferenceNodes3D10();
      nodes[5].X0() = nodes[5].X() = 0.4;
      nodes[5].Y0() = nodes[5].Y() = 0.4;
      nodes[6].X0() = nodes[6].X() = 0.1;
      nodes[6].Z0() = nodes[6].Z() = 0.1;
      nodes[7].X0() = nodes[7].X() = 0.1;
      nodes[7].Y0() = nodes[7].Y() = 0.1;
      nodes[8].X0() = nodes[8].X() = 0.4;
      nodes[8].Z0() = nodes[8].Z() = 0.4;
      nodes[9].Y0() = nodes[9].Y() = 0.4;
      nodes[9].Z0() = nodes[9].Z() = 0.4;
      return Geometry<NodeType>::Pointer(new Tetrahedra3D10<NodeType>(nodes));
  }
}

  // /// Tests

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10EdgesNumber, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_EXPECT_EQ(geom->EdgesNumber(), 6);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10FacesNumber, KratosCoreGeometriesFastSuite) {
    Geometry<NodeType>::Pointer geom = GenerateReferenceTetrahedra3D10();
    KRATOS_EXPECT_EQ(geom->FacesNumber(), 4);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10Volume, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_EXPECT_NEAR(geom->Volume(), 1.0 / 6.0, TOLERANCE);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10MinEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->MinEdgeLength(), "Calling base class 'MinEdgeLength' method instead of derived class one.");
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10MaxEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->MaxEdgeLength(), "Calling base class 'MaxEdgeLength' method instead of derived class one.");
  }


  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10Circumradius, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Circumradius(), "Calling base class 'Circumradius' method instead of derived class one.");
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10Inradius, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Inradius(), "Calling base class 'Inradius' method instead of derived class one.");
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10IsInside, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateCurvedTetrahedra3D10();
    Point PointInside(0.1, 0.1, 0.1);
    Point PointOutside(0.5, 0.5, 0.0);
    Point PointInVertex(0.0, 0.0, 0.0);
    Point PointInEdge(0.6, 0.0, 0.0);
    Point LocalCoords;
    KRATOS_EXPECT_TRUE(geom->IsInside(PointInside, LocalCoords, EPSILON));
    KRATOS_EXPECT_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
    KRATOS_EXPECT_TRUE(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
    geom->IsInside(PointInEdge, LocalCoords, EPSILON);
    KRATOS_WATCH(LocalCoords);
    KRATOS_EXPECT_TRUE(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    Vector JacobianDeterminants;
    geom->DeterminantOfJacobian(JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_2);
    for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
        KRATOS_EXPECT_NEAR(JacobianDeterminants[i], 1.0, TOLERANCE);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10IntegrationPointsNumber, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_EXPECT_EQ(geom->IntegrationPointsNumber(), 4);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    for (unsigned g = 0; g < geom->IntegrationPointsNumber(); ++g)
      KRATOS_EXPECT_NEAR(geom->DeterminantOfJacobian(g), 1.0, TOLERANCE);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    array_1d<double, 3> coord(3);
    coord[0] = 1.0 / 2.0;
    coord[1] = 1.0 / 4.0;
    coord[2] = 1.0 / 16.0;
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(0, coord), -0.1171875, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(1, coord), 0.0, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(2, coord), -0.125, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(3, coord), -0.0546875, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(4, coord), 0.375, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(5, coord), 0.5, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(6, coord), 0.1875, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(7, coord), 0.046875, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(8, coord), 0.125, TOLERANCE);
    KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(9, coord), 0.0625, TOLERANCE);
    CrossCheckShapeFunctionsValues(*geom);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    TestAllShapeFunctionsLocalGradients(*geom);
  }

KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10ShapeFunctionsSecondDerivatives, KratosCoreGeometriesFastSuite)
{
    // Set a second order isoparametric tetrahedron so we avoid the Jacobian transformation
    auto p_node_0 = Kratos::make_intrusive<NodeType>(1, 0.0, 0.0, 0.0);
    auto p_node_1 = Kratos::make_intrusive<NodeType>(2, 1.0, 0.0, 0.0);
    auto p_node_2 = Kratos::make_intrusive<NodeType>(3, 0.0, 1.0, 0.0);
    auto p_node_3 = Kratos::make_intrusive<NodeType>(4, 0.0, 0.0, 1.0);
    auto p_node_4 = Kratos::make_intrusive<NodeType>(5, 0.5, 0.0, 0.0);
    auto p_node_5 = Kratos::make_intrusive<NodeType>(6, 0.5, 0.5, 0.0);
    auto p_node_6 = Kratos::make_intrusive<NodeType>(7, 0.0, 0.5, 0.0);
    auto p_node_7 = Kratos::make_intrusive<NodeType>(8, 0.0, 0.0, 0.5);
    auto p_node_8 = Kratos::make_intrusive<NodeType>(9, 0.5, 0.0, 0.5);
    auto p_node_9 = Kratos::make_intrusive<NodeType>(10, 0.0, 0.5, 0.5);
    auto p_geom = Kratos::make_unique<Tetrahedra3D10<NodeType>>(
        p_node_0, p_node_1, p_node_2, p_node_3, p_node_4, p_node_5, p_node_6, p_node_7, p_node_8, p_node_9);

    // Compute the shape functions second derivatives in the barycencer
    GeometryType::ShapeFunctionsSecondDerivativesType DDN_DDX;
    array_1d<double,3> barycenter_coords({1.0/3.0, 1.0/3.0, 1.0/3.0});
    p_geom->ShapeFunctionsSecondDerivatives(DDN_DDX, barycenter_coords);

    // Evaluate the second derivatives of the auxiliary field xi^2 + eta^2 + zeta^2 at the barycenter
    BoundedMatrix<double, 3, 3> expected_result = ZeroMatrix(3,3);
    expected_result(0,0) = 2.0;
    expected_result(1,1) = 2.0;
    expected_result(2,2) = 2.0;

    const auto field_func = [](array_1d<double,3>& rCoords){return std::pow(rCoords[0],2) + std::pow(rCoords[1],2) + std::pow(rCoords[2],2);};
    BoundedMatrix<double, 3, 3> result = ZeroMatrix(3,3);
    for (IndexType i = 0; i < 10; ++i) {
        const auto& r_DDN_DDX_i = DDN_DDX[i];
        const double field_i = field_func((*p_geom)[i]);
        for (IndexType d1 = 0; d1 < 3; ++d1) {
            for (IndexType d2 = 0; d2 < 3; ++d2) {
                result(d1,d2) += r_DDN_DDX_i(d1,d2) * field_i;
            }
        }
    }

    KRATOS_CHECK_MATRIX_NEAR(expected_result, result, 1.0e-12)
}

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10AverageEdgeLength, KratosCoreGeometriesFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_EXPECT_NEAR(geom->AverageEdgeLength(), 1.20710678119, 1e-7);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10HasIntersection, KratosCoreGeometriesFastSuite) {
    const Point LowPoint(0.1,0.1,-0.1);
    const Point HighPoint(0.1,0.1,1.1);

    const Point OutLowPoint(1.1,0.1,-0.1);
    const Point OutHighPoint(1.1,0.1,1.1);

    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_EXPECT_TRUE(geom->HasIntersection(LowPoint, HighPoint));
    KRATOS_EXPECT_FALSE(geom->HasIntersection(OutLowPoint, OutHighPoint));

    auto curved_geom = GenerateCurvedTetrahedra3D10();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
      curved_geom->HasIntersection(LowPoint, HighPoint),
      "\"HasIntersection\" is not implemented for non-planar 10 noded tetrahedra.");
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10CalculateDistance, KratosCoreGeometriesFastSuite)
  {
      auto geom = GenerateReferenceTetrahedra3D10();

      Point point1(0.25, 0.25, 0.25);
      KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point1), 0.0);
      KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point1), GeometryUtils::CalculateDistanceFrom3DGeometry(*geom, point1));

      Point point2(1.5, 0.0, 0.0);
      KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point2), 0.5);
      KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point2), GeometryUtils::CalculateDistanceFrom3DGeometry(*geom, point2));
  }

} // namespace Kratos::Testing.
