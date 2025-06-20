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
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/triangle_3d_6.h"
#include "tests/cpp_tests/geometries/test_geometry.h"
#include "tests/cpp_tests/geometries/test_shape_function_derivatives.h"
#include "tests/cpp_tests/geometries/cross_check_shape_functions_values.h"

namespace Kratos::Testing {
typedef Node NodeType;

// /// Factory functions
namespace {
  Geometry<NodeType>::PointsArrayType GenerateReferenceNodes3D6()
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

  Geometry<NodeType>::Pointer GenerateReferenceTriangle3D6() {
      return Geometry<NodeType>::Pointer(
          new Triangle3D6<NodeType>(GenerateReferenceNodes3D6()));
  }

}

    ///// Tests

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6EdgesNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        KRATOS_EXPECT_EQ(geom->EdgesNumber(), 3);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6Area, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        KRATOS_EXPECT_NEAR(geom->Area(), 0.5, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6Volume, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        // TODO: Remove code in June 2023
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Volume(), "Calling base class 'Volume' method instead of derived class one.");
        // TODO: Activate code in June 2023
        //KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Volume(), "Triangle3D6:: Method not well defined. Replace with DomainSize() instead.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6MinEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->MinEdgeLength(), "Calling base class 'MinEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6MaxEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->MaxEdgeLength(), "Calling base class 'MaxEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6AverageEdgeLength, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->AverageEdgeLength(), "Calling base class 'AverageEdgeLength' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6Circumradius, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Circumradius(), "Calling base class 'Circumradius' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6Inradius, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(geom->Inradius(), "Calling base class 'Inradius' method instead of derived class one.");
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6DeterminantOfJacobianArray1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        Vector JacobianDeterminants;
        geom->DeterminantOfJacobian(JacobianDeterminants, GeometryData::IntegrationMethod::GI_GAUSS_2);
        for (unsigned int i=0; i<JacobianDeterminants.size(); ++i)
            KRATOS_EXPECT_NEAR(JacobianDeterminants[i], 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6IntegrationPointsNumber, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        KRATOS_EXPECT_EQ(geom->IntegrationPointsNumber(), 3);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6DeterminantOfJacobianIndex1, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        for (unsigned g = 0; g < geom->IntegrationPointsNumber(); ++g)
        KRATOS_EXPECT_NEAR(geom->DeterminantOfJacobian(g), 1.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6ShapeFunctionsLocalGradients, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        TestAllShapeFunctionsLocalGradients(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6ShapeFunctionsValues, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();
        array_1d<double, 3> coord(3);
        coord[0] = 0.5;
        coord[1] = 0.125;
        coord[2] = 0.0;
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(0, coord), -0.09375, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(1, coord), 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(2, coord), -0.09375, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(3, coord), 0.75, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(4, coord), 0.25, TOLERANCE);
        KRATOS_EXPECT_NEAR(geom->ShapeFunctionValue(5, coord), 0.1875, TOLERANCE);
        CrossCheckShapeFunctionsValues(*geom);
    }

    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6LumpingFactorsRegularShape, KratosCoreGeometriesFastSuite) {
        auto geom = GenerateReferenceTriangle3D6();

        Vector lumping_factors(6);
        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::ROW_SUM);

        KRATOS_EXPECT_NEAR(lumping_factors[0], 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[1], 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[2], 0.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[3], 1.0/3.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[4], 1.0/3.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[5], 1.0/3.0, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::DIAGONAL_SCALING);

        KRATOS_EXPECT_NEAR(lumping_factors[0], 0.0328638, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[1], 0.0328638, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[2], 0.0328638, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[3], 0.300469, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[4], 0.300469, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[5], 0.300469, TOLERANCE);

        geom->LumpingFactors(lumping_factors, Geometry<Node>::LumpingMethods::QUADRATURE_ON_NODES);

        KRATOS_EXPECT_NEAR(lumping_factors[0], 1.0/6.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[1], 1.0/6.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[2], 1.0/6.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[3], 1.0/6.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[4], 1.0/6.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(lumping_factors[5], 1.0/6.0, TOLERANCE);
    }

    /**
     * Checks the distance from a point to a triangle
     */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6CalculateDistance, KratosCoreGeometriesFastSuite)
    {
        auto geom = GenerateReferenceTriangle3D6();

        Point point1(0.0, 0.0, 0.0);
        KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point1), 0.0);

        Point point2(0.0, 0.0, 0.5);
        KRATOS_EXPECT_DOUBLE_EQ(geom->CalculateDistance(point2), 0.5);
    }

    /*
    * Computes point local coordinates from a given point.
    * his triangle has straight edges so we can use the same
    * implementation as the Triangle3D3
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6PointLocalCoordinates, KratosCoreGeometriesFastSuite) {
        Triangle3D6<Point> geom(
            Kratos::make_shared<Point>(0.0, 0.0, 0.0),
            Kratos::make_shared<Point>(1.0, 0.0, 0.0),
            Kratos::make_shared<Point>(0.0, 1.0, 0.0),
            Kratos::make_shared<Point>(0.5, 0.0, 0.0),
            Kratos::make_shared<Point>(0.5, 0.5, 0.0),
            Kratos::make_shared<Point>(0.0, 0.5, 0.0)
        );

        // Compute the global coordinates of the baricentre
        array_1d<double, 3> baricentre;
        baricentre[0] = 1.0/3.0; baricentre[1] = 1.0/3.0; baricentre[2] = 0.0;

        // Compute the baricentre local coordinates
        array_1d<double, 3> baricentre_local_coords;
        geom.PointLocalCoordinates(baricentre_local_coords, baricentre);

        KRATOS_EXPECT_NEAR(baricentre_local_coords(0), 1.0/3.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(baricentre_local_coords(1), 1.0/3.0, TOLERANCE);
        KRATOS_EXPECT_NEAR(baricentre_local_coords(2), 0.0, TOLERANCE);
    }

    /*
    * Computes point local coordinates from a given point.
    * This triangle does not ghave straight edges so this will
    * throw an exception
    */
    KRATOS_TEST_CASE_IN_SUITE(Triangle3D6PointLocalCoordinatesError, KratosCoreGeometriesFastSuite) {
        Triangle3D6<Point> geom(
            Kratos::make_shared<Point>(1.0, 0.0, 0.0),
            Kratos::make_shared<Point>(0.0, 1.0, 0.0),
            Kratos::make_shared<Point>(0.0, 0.0, 1.0),
            Kratos::make_shared<Point>(0.6, 0.5, 0.0),
            Kratos::make_shared<Point>(0.0, 0.5, 0.5),
            Kratos::make_shared<Point>(0.5, 0.0, 0.5)
        );

        // Compute the global coordinates of the baricentre
        array_1d<double, 3> baricentre;
        baricentre[0] = 1.0/3.0; baricentre[1] = 1.0/3.0; baricentre[2] = 1.0/3.0;

        // Compute the baricentre local coordinates
        array_1d<double, 3> baricentre_local_coords;
        KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        geom.PointLocalCoordinates(baricentre_local_coords, baricentre),
            "ERROR:: Attention, the Point Local Coordinates must be specialized for the current geometry");
    }

} // namespace Kratos::Testing.
