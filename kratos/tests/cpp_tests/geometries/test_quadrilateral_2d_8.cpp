//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/quadrilateral_2d_8.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos {
namespace Testing {
typedef Node NodeType;

// /// Factory functions
namespace {
  Geometry<NodeType>::PointsArrayType GenerateReferenceNodes2D8()
  {
      Geometry<NodeType>::PointsArrayType points;
      points.push_back(NodeType::Pointer(new NodeType(1, 0.0, 0.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(2, 1.0, 0.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(3, 1.0, 1.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(4, 0.0, 1.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(5, 0.5, 0.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(6, 1.0, 0.5, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(7, 0.5, 1.0, 0.0)));
      points.push_back(NodeType::Pointer(new NodeType(8, 0.0, 0.5, 0.0)));
      return points;
  }

  Geometry<NodeType>::Pointer GenerateReferenceQuadrilateral2D8() {
      return Geometry<NodeType>::Pointer(
          new Quadrilateral2D8<NodeType>(GenerateReferenceNodes2D8()));
  }
}

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8EdgesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        KRATOS_CHECK_EQUAL(geom->EdgesNumber(), 4);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8Area, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        KRATOS_CHECK_NEAR(geom->Area(), 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8Volume, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        // TODO: Remove code in June 2023
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Volume(), "Calling base class 'Volume' method instead of derived class one.");
        // TODO: Activate code in June 2023
        //KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Volume(), "Quadrilateral2D8:: Method not well defined. Replace with DomainSize() instead.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8MinEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->MinEdgeLength(), "Calling base class 'MinEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8MaxEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->MaxEdgeLength(), "Calling base class 'MaxEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8AverageEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->AverageEdgeLength(), "Calling base class 'AverageEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8Circumradius, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Circumradius(), "Calling base class 'Circumradius' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8Inradius, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        KRATOS_CHECK_EXCEPTION_IS_THROWN(geom->Inradius(), "Calling base class 'Inradius' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8IsInside, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        Point PointInside(0.1, 0.1);
        Point PointOutside(1.5, 0.5);
        Point PointInVertex(0.0, 0.0);
        Point PointInEdge(0.5, 0.0);
        Point LocalCoords;
        KRATOS_CHECK(geom->IsInside(PointInside, LocalCoords, EPSILON));
        KRATOS_CHECK_IS_FALSE(geom->IsInside(PointOutside, LocalCoords, EPSILON));
        KRATOS_CHECK(geom->IsInside(PointInVertex, LocalCoords, EPSILON));
        KRATOS_CHECK(geom->IsInside(PointInEdge, LocalCoords, EPSILON));
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian(JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_2);
        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
            KRATOS_CHECK_NEAR(JacobianDeterminants[i], 0.25, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8IntegrationPointsNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        KRATOS_CHECK_EQUAL(geom->IntegrationPointsNumber(), 9);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        for (unsigned g = 0; g < geom->IntegrationPointsNumber(); ++g)
        KRATOS_CHECK_NEAR(geom->DeterminantOfJacobian(g), 0.25, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        array_1d<double, 3> coord(3);
        coord[0] = 1.0 / 2.0;
        coord[1] = 1.0 / 8.0;
        coord[2] = 0.0;

        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(0, coord), -0.177734, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(1, coord), -0.205078, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(2, coord), -0.158203, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(3, coord), -0.193359, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(4, coord), 0.328125, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(5, coord), 0.738281, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(6, coord), 0.421875, TOLERANCE);
        KRATOS_CHECK_NEAR(geom->ShapeFunctionValue(7, coord), 0.246094, TOLERANCE);
        CrossCheckShapeFunctionsValues(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();
        TestAllShapeFunctionsLocalGradients(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Quadrilateral2D8LumpingFactorsRegularShape, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceQuadrilateral2D8();

        Vector lumping_factors(8);
        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::ROW_SUM);

        KRATOS_CHECK_NEAR(lumping_factors[0], -1/12.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[1], -1/12.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[2], -1/12.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[3], -1/12.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[4], 1.0/3.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[5], 1.0/3.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[6], 1.0/3.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[7], 1.0/3.0, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);

        KRATOS_CHECK_NEAR(lumping_factors[0], 0.0394737, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[1], 0.0394737, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[2], 0.0394737, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[3], 0.0394737, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[4], 0.210526, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[5], 0.210526, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[6], 0.210526, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[7], 0.210526, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::QUADRATURE_ON_NODES);

        KRATOS_CHECK_NEAR(lumping_factors[0], 1.0/8.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[1], 1.0/8.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[2], 1.0/8.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[3], 1.0/8.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[4], 1.0/8.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[5], 1.0/8.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[6], 1.0/8.0, TOLERANCE);
        KRATOS_CHECK_NEAR(lumping_factors[7], 1.0/8.0, TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.
