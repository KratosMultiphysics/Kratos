//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/tetrahedra_3d_10.h"
#include "tests/geometries/test_geometry.h"
#include "tests/geometries/test_shape_function_derivatives.h"
#include "tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos {
namespace Testing {

  // /// Factory functions

namespace{
  typedef Node<3> NodeType;

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

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10EdgesNumber, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 6);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10FacesNumber, KratosCoreFastSuite) {
    Geometry<NodeType>::Pointer geom = GenerateReferenceTetrahedra3D10();
    KRATOS_CHECK_EQUAL(geom->FacesNumber(), 4);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10Volume, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_CHECK_NEAR(geom->Volume(), 1.0 / 6.0, TOLERANCE);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10MinEdgeLength, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->MinEdgeLength(), "Calling base class 'MinEdgeLength' method instead of derived class one.");
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10MaxEdgeLength, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->MaxEdgeLength(), "Calling base class 'MaxEdgeLength' method instead of derived class one.");
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10AverageEdgeLength, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->AverageEdgeLength(), "Calling base class 'AverageEdgeLength' method instead of derived class one.");
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10Circumradius, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Circumradius(), "Calling base class 'Circumradius' method instead of derived class one.");
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10Inradius, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Inradius(), "Calling base class 'Inradius' method instead of derived class one.");
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10IsInside, KratosCoreFastSuite) {
    auto geom = GenerateCurvedTetrahedra3D10();
    Point PointInside(0.1, 0.1, 0.1);
    Point PointOutside(0.5, 0.5, 0.0);
    Point PointInVertex(0.0, 0.0, 0.0);
    Point PointInEdge(0.6, 0.0, 0.0);
    Point LocalCoords;
    KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
    KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
    KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
    geom->IsInside(PointInEdge, LocalCoords, EPSILON);
    KRATOS_WATCH(LocalCoords);
    KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10DeterminantOfJacobianArray1, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    Vector JacobianDeterminants;
    geom->DeterminantOfJacobian(JacobianDeterminants, GeometryData::GI_GAUSS_2);
    for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
        KRATOS_CHECK_NEAR(JacobianDeterminants[i], 1.0, TOLERANCE);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10IntegrationPointsNumber, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    KRATOS_CHECK_EQUAL(geom->IntegrationPointsNumber(), 4);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10DeterminantOfJacobianIndex1, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    for (unsigned g = 0; g < geom->IntegrationPointsNumber(); ++g)
      KRATOS_CHECK_NEAR(geom->DeterminantOfJacobian(g), 1.0, TOLERANCE);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10ShapeFunctionsValues, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    array_1d<double, 3> coord(3);
    coord[0] = 1.0 / 2.0;
    coord[1] = 1.0 / 4.0;
    coord[2] = 1.0 / 16.0;
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), -0.1171875, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), 0.0, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(2, coord), -0.125, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(3, coord), -0.0546875, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(4, coord), 0.375, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(5, coord), 0.5, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(6, coord), 0.1875, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(7, coord), 0.046875, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(8, coord), 0.125, TOLERANCE);
    KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(9, coord), 0.0625, TOLERANCE);
    CrossCheckShapeFunctionsValues(*geom);
  }

  KRATOS_TEST_CASE_IN_SUITE(Tetrahedra3D10ShapeFunctionsLocalGradients, KratosCoreFastSuite) {
    auto geom = GenerateReferenceTetrahedra3D10();
    TestAllShapeFunctionsLocalGradients(*geom);
  }

} // namespace Testing.
} // namespace Kratos.
