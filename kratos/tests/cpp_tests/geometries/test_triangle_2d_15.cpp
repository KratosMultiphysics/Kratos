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
#include "geometries/triangle_2d_15.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos {
namespace Testing {
typedef Node NodeType;

// /// Factory functions
namespace {
  Geometry<NodeType>::PointsArrayType GenerateReferenceNodes2D15()
  {
      Geometry<NodeType>::PointsArrayType points;
      points.push_back(NodeType::Pointer(new NodeType( 1, 0.00, 0.00, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 2, 1.00, 0.00, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 3, 0.00, 1.00, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 4, 0.25, 0.00, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 5, 0.50, 0.00, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 6, 0.75, 0.00, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 7, 0.75, 0.25, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 8, 0.50, 0.50, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType( 9, 0.25, 0.75, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(10, 0.00, 0.75, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(11, 0.00, 0.50, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(12, 0.00, 0.25, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(13, 0.25, 0.25, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(14, 0.50, 0.25, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(15, 0.25, 0.50, 0.0)));
      return points;
  }

  Geometry<NodeType>::Pointer GenerateReferenceTriangle2D15() {
      return Geometry<NodeType>::Pointer(
          new Triangle2D15<NodeType>(GenerateReferenceNodes2D15()));
  }

  Geometry<NodeType>::Pointer GenerateCurvedTriangle2D15() {
      auto nodes = GenerateReferenceNodes2D15();
      nodes[ 3].Y0() = nodes[ 3].Y() = 0.10;
      nodes[ 4].Y0() = nodes[ 4].Y() = 0.15;
      nodes[ 5].Y0() = nodes[ 5].Y() = 0.10;
      nodes[ 6].X0() = nodes[ 6].X() = 0.65;
      nodes[ 6].Y0() = nodes[ 6].Y() = 0.15;
      nodes[ 7].X0() = nodes[ 7].X() = 0.35;
      nodes[ 7].Y0() = nodes[ 7].Y() = 0.35;
      nodes[ 6].X0() = nodes[ 6].X() = 0.15;
      nodes[ 6].Y0() = nodes[ 6].Y() = 0.65;
      nodes[ 9].X0() = nodes[ 9].X() = 0.10;
      nodes[10].X0() = nodes[10].X() = 0.15;
      nodes[11].X0() = nodes[11].X() = 0.10;
      return Geometry<NodeType>::Pointer(new Triangle2D15<NodeType>(nodes));
  }
}

    ///// Tests

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15EdgesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 3);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15FacesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        KRATOS_CHECK_EQUAL(geom->FacesNumber(), 3);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15Area, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        KRATOS_CHECK_NEAR(geom->Area(), 0.5, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15Volume, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Volume(), "Calling base class 'Volume' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15MinEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->MinEdgeLength(), "Calling base class 'MinEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15MaxEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->MaxEdgeLength(), "Calling base class 'MaxEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15AverageEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->AverageEdgeLength(), "Calling base class 'AverageEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15Circumradius, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Circumradius(), "Calling base class 'Circumradius' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15Inradius, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Inradius(), "Calling base class 'Inradius' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15IsInside, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateCurvedTriangle2D15();
        Point PointInside(0.1, 0.1);
        Point PointOutside(0.5, 0.5);
        Point PointInVertex(0.0, 0.0);
        Point PointInEdge(0.5, 0.15);
        Point LocalCoords;
        KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
        KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
        KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
        KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian(JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_2);
        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15IntegrationPointsNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        KRATOS_CHECK_EQUAL(geom->IntegrationPointsNumber(), 12);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        for (unsigned g = 0; g < geom->IntegrationPointsNumber(); ++g)
        KRATOS_CHECK_NEAR(geom->DeterminantOfJacobian(g), 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        array_1d<double, 3> coord(3);
        coord[0] = 1.0 / 2.0;
        coord[1] = 1.0 / 8.0;
        coord[2] = 0.0;
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue( 0, coord),  0.0234375, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue( 1, coord),  0.0      , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue( 2, coord), -0.0390625, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue( 3, coord), -0.125    , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue( 4, coord),  0.375    , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue( 5, coord),  0.0      , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue( 6, coord),  0.0      , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue( 7, coord), -0.125    , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue( 8, coord),  0.125    , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue( 9, coord),  0.09375  , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(10, coord), -0.046875 , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(11, coord), -0.03125  , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(12, coord),  0.375    , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(13, coord),  0.75     , TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(14, coord), -0.375    , TOLERANCE);
        CrossCheckShapeFunctionsValues(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();
        TestAllShapeFunctionsLocalGradients(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D15LumpingFactorsRegularShape, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D15();

        Vector lumping_factors(15);
        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::ROW_SUM);

        KRATOS_CHECK_NEAR(lumping_factors[ 0],  0.0     , TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 1],  0.0     , TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 2],  0.0     , TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 3],  4.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 4], -1.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 5],  4.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 6],  4.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 7], -1.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 8],  4.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 9],  4.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[10], -1.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[11],  4.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[12],  8.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[13],  8.0/45.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[14],  8.0/45.0, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);

        KRATOS_CHECK_NEAR(lumping_factors[ 0], 0.0049479514208965, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 1], 0.0049479514208965, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 2], 0.0049479514208965, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 3], 0.0415105750486979, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 4], 0.0411365548334816, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 5], 0.0415105750486979, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 6], 0.0415105750486979, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 7], 0.0411365548334816, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 8], 0.0415105750486979, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 9], 0.0415105750486979, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[10], 0.0411365548334816, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[11], 0.0415105750486979, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[12], 0.2042276769815594, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[13], 0.2042276769815594, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[14], 0.2042276769815594, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::QUADRATURE_ON_NODES);

        KRATOS_CHECK_NEAR(lumping_factors[ 0], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 1], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 2], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 3], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 4], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 5], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 6], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 7], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 8], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[ 9], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[10], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[11], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[12], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[13], 1.0/15.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[14], 1.0/15.0, TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.
