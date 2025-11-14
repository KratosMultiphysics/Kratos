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
#include "geometries/triangle_2d_6.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos {
namespace Testing {
typedef Node NodeType;

// /// Factory functions
namespace {
  Geometry<NodeType>::PointsArrayType GenerateReferenceNodes2D6()
  {
      Geometry<NodeType>::PointsArrayType points;
      points.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(3, 0.0, 1.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(4, 0.5, 0.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(5, 0.5, 0.5, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(6, 0.0, 0.5, 0.0)));
      return points;
  }

  Geometry<NodeType>::Pointer GenerateReferenceTriangle2D6() {
      return Geometry<NodeType>::Pointer(
          new Triangle2D6<NodeType>(GenerateReferenceNodes2D6()));
  }

  Geometry<NodeType>::Pointer GenerateCurvedTriangle2D6() {
      auto nodes = GenerateReferenceNodes2D6();
      nodes[3].Y0() = nodes[3].Y() = 0.1;
      nodes[4].X0() = nodes[4].X() = 0.4;
      nodes[4].Y0() = nodes[4].Y() = 0.4;
      nodes[5].X0() = nodes[5].X() = 0.1;
      return Geometry<NodeType>::Pointer(new Triangle2D6<NodeType>(nodes));
  }
}

    ///// Tests

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6EdgesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        KRATOS_EXPECT_EQ(geom->EdgesNumber(), 3);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6FacesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        KRATOS_EXPECT_EQ(geom->FacesNumber(), 1);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6Faces, KratosCoreGeometriesFastSuite) {
        auto p_geom = GenerateReferenceTriangle2D6();

        const auto& r_faces = p_geom->GenerateFaces();
        ASSERT_EQ(r_faces.size(), 1);
        for (std::size_t i = 0; i < r_faces.front().PointsNumber(); ++i) {
            KRATOS_EXPECT_NEAR(r_faces.front()[i].X(), (p_geom->pGetPoint(i))->X(), TOLERANCE);
            KRATOS_EXPECT_NEAR(r_faces.front()[i].Y(), (p_geom->pGetPoint(i))->Y(), TOLERANCE);
        }
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6Area, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        KRATOS_EXPECT_NEAR(geom->Area(), 0.5, ZERO_TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6Volume, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        // TODO: Remove code in June 2023
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Volume(), "Calling base class 'Volume' method instead of derived class one.");
        // TODO: Activate code in June 2023
        //KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Volume(), "Triangle2D6:: Method not well defined. Replace with DomainSize() instead.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6MinEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->MinEdgeLength(), "Calling base class 'MinEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6MaxEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->MaxEdgeLength(), "Calling base class 'MaxEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6AverageEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->AverageEdgeLength(), "Calling base class 'AverageEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6Circumradius, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Circumradius(), "Calling base class 'Circumradius' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6Inradius, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Inradius(), "Calling base class 'Inradius' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6IsInside, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateCurvedTriangle2D6();
        Point PointInside(0.1, 0.1);
        Point PointOutside(0.5, 0.5);
        Point PointInVertex(0.0, 0.0);
        Point PointInEdge(0.5, 0.1);
        Point LocalCoords;
        KRATOS_EXPECT_TRUE(geom->IsInside(PointInside, LocalCoords, EPSILON));
        KRATOS_EXPECT_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
        KRATOS_EXPECT_TRUE(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
        KRATOS_EXPECT_TRUE(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian(JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_2);
        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
            KRATOS_EXPECT_NEAR(JacobianDeterminants[i], 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6IntegrationPointsNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        KRATOS_EXPECT_EQ(geom->IntegrationPointsNumber(), 3);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        for (unsigned g = 0; g < geom->IntegrationPointsNumber(); ++g)
        KRATOS_EXPECT_NEAR(geom->DeterminantOfJacobian(g), 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        array_1d<double, 3> coord(3);
        coord[0] = 1.0 / 2.0;
        coord[1] = 1.0 / 8.0;
        coord[2] = 0.0;
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(0, coord), -0.09375, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(1, coord), 0.0,      TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(2, coord), -0.09375, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(3, coord), 0.75,     TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(4, coord), 0.25,     TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(5, coord), 0.1875,   TOLERANCE);
        CrossCheckShapeFunctionsValues(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();
        TestAllShapeFunctionsLocalGradients(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle2D6LumpingFactorsRegularShape, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle2D6();

        Vector lumping_factors(6);
        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::ROW_SUM);

        KRATOS_EXPECT_NEAR(lumping_factors[0], 0.0,        TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[1], 0.0,        TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[2], 0.0,        TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[3], 1.0/3.0,    TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[4], 1.0/3.0,    TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[5], 1.0/3.0,    TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);

        KRATOS_EXPECT_NEAR(lumping_factors[0], 0.0328638,  TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[1], 0.0328638,  TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[2], 0.0328638,  TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[3], 0.300469,   TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[4], 0.300469,   TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[5], 0.300469,   TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::QUADRATURE_ON_NODES);

        KRATOS_EXPECT_NEAR(lumping_factors[0], 1.0/6.0,   TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[1], 1.0/6.0,   TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[2], 1.0/6.0,   TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[3], 1.0/6.0,   TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[4], 1.0/6.0,   TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[5], 1.0/6.0,   TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.
