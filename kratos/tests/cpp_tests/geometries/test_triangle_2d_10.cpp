//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohamed Nabi
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/triangle_2d_10.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos {
namespace Testing {
typedef Node NodeType;
const double oneThird = 1.0 / 3.0;
const double twoThird = 2.0 / 3.0;

// /// Factory functions
namespace {
  Geometry<NodeType>::PointsArrayType GenerateReferenceNodes2D10()
  {
      Geometry<NodeType>::PointsArrayType points;
      points.push_back(NodeType::Pointer(new NodeType( 1, 0.0     , 0.0     , 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 2, 1.0     , 0.0     , 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 3, 0.0     , 1.0     , 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 4, oneThird, 0.0     , 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 5, twoThird, 0.0     , 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 6, twoThird, oneThird, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 7, oneThird, twoThird, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 8, 0.0     , twoThird, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 9, 0.0     , oneThird, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(10, oneThird, oneThird, 0.0)));
      return points;
  }

  Geometry<NodeType>::Pointer GenerateReferenceTriangle2D10() {
      return Geometry<NodeType>::Pointer(
          new Triangle2D10<NodeType>(GenerateReferenceNodes2D10()));
  }

  Geometry<NodeType>::Pointer GenerateCurvedTriangle2D10() {
      auto nodes = GenerateReferenceNodes2D10();
      nodes[3].Y0() = nodes[3].Y() = 0.10;
      nodes[4].Y0() = nodes[4].Y() = 0.10;
      nodes[5].X0() = nodes[5].X() = 0.56;
      nodes[5].Y0() = nodes[5].Y() = 0.23;
      nodes[6].X0() = nodes[6].X() = 0.23;
      nodes[6].Y0() = nodes[6].Y() = 0.56;
      nodes[7].X0() = nodes[7].X() = 0.10;
      nodes[8].X0() = nodes[8].X() = 0.10;
      return Geometry<NodeType>::Pointer(new Triangle2D10<NodeType>(nodes));
  }
}

    ///// Tests

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10EdgesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        KRATOS_EXPECT_EQ(geom->EdgesNumber(), 3);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10FacesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        KRATOS_EXPECT_EQ(geom->FacesNumber(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10Area, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        KRATOS_EXPECT_NEAR(geom->Area(), 0.5, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10Volume, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Volume(), "Calling base class 'Volume' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10MinEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->MinEdgeLength(), "Calling base class 'MinEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10MaxEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->MaxEdgeLength(), "Calling base class 'MaxEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10AverageEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->AverageEdgeLength(), "Calling base class 'AverageEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10Circumradius, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Circumradius(), "Calling base class 'Circumradius' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10Inradius, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Inradius(), "Calling base class 'Inradius' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10IsInside, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateCurvedTriangle2D10();
        Point PointInside(0.1, 0.1);
        Point PointOutside(0.5, 0.5);
        Point PointInVertex(0.0, 0.0);
        Point PointInEdge(oneThird, 0.1);
        Point LocalCoords;
        KRATOS_EXPECT_TRUE(geom->IsInside(PointInside, LocalCoords, EPSILON));
        KRATOS_EXPECT_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
        KRATOS_EXPECT_TRUE(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
        KRATOS_EXPECT_TRUE(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian(JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_2);
        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
            KRATOS_EXPECT_NEAR(JacobianDeterminants[i], 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10IntegrationPointsNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        KRATOS_EXPECT_EQ(geom->IntegrationPointsNumber(), 6);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        for (unsigned g = 0; g < geom->IntegrationPointsNumber(); ++g)
        KRATOS_EXPECT_NEAR(geom->DeterminantOfJacobian(g), 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        array_1d<double, 3> coord(3);
        coord[0] = 1.0 / 2.0;
        coord[1] = 1.0 / 8.0;
        coord[2] = 0.0;
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(0, coord), -0.0205078125, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(1, coord), -0.0625      , TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(2, coord),  0.0634765625, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(3, coord),  0.10546875  , TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(4, coord),  0.421875    , TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(5, coord),  0.140625    , TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(6, coord), -0.17578125  , TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(7, coord), -0.1318359375, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(8, coord),  0.0263671875, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(9, coord),  0.6328125   , TOLERANCE);
        CrossCheckShapeFunctionsValues(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();
        TestAllShapeFunctionsLocalGradients(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D10LumpingFactorsRegularShape, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D10();

        Vector lumping_factors(10);
        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::ROW_SUM);

        KRATOS_EXPECT_NEAR(lumping_factors[0], 1.0/30.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[1], 1.0/30.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[2], 1.0/30.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[3], 3.0/40.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[4], 3.0/40.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[5], 3.0/40.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[6], 3.0/40.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[7], 3.0/40.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[8], 3.0/40.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[9], 9.0/20.0, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);

        KRATOS_EXPECT_NEAR(lumping_factors[0], 0.014042867701404, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[1], 0.014042867701404, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[2], 0.014042867701404, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[3], 0.099778270509978, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[4], 0.099778270509978, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[5], 0.099778270509978, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[6], 0.099778270509978, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[7], 0.099778270509978, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[8], 0.099778270509978, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[9], 0.359201773835920, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::QUADRATURE_ON_NODES);

        KRATOS_EXPECT_NEAR(lumping_factors[0], 0.1, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[1], 0.1, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[2], 0.1, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[3], 0.1, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[4], 0.1, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[5], 0.1, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[6], 0.1, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[7], 0.1, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[8], 0.1, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[9], 0.1, TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.
