//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Thomas Oberbichler
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/pointer_vector.h"
#include "geometries/nurbs_curve_geometry.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

typedef Node<3> NodeType;

// /// Factory functions
//namespace {
    NurbsCurveGeometry<2, PointerVector<Point>> GenerateReferenceCurve2d()
    {
        PointerVector<Point> points;

        points.push_back(Point::Pointer(new Point(0, 0, 0)));
        points.push_back(Point::Pointer(new Point(3.3333333333333335, 1.6666666666666667, 0)));
        points.push_back(Point::Pointer(new Point(6.6666666666666661, 3.333333333333333, 0)));
        points.push_back(Point::Pointer(new Point(10, 5, 0)));

        Vector knot_vector = ZeroVector(6);
        knot_vector[0] = 0.0;
        knot_vector[1] = 0.0;
        knot_vector[2] = 0.0;
        knot_vector[3] = 11.180339887498949;
        knot_vector[4] = 11.180339887498949;
        knot_vector[5] = 11.180339887498949;

        int p = 3;

        auto curve = NurbsCurveGeometry<2, PointerVector<Point>>(points, p, knot_vector);

        return curve;
    }

    NurbsCurveGeometry<3, PointerVector<Point>> GenerateCircle()
    {
        PointerVector<Point> points(9);

        points(0) = Kratos::make_shared<Point>(9.4868329805051381,  5, 3.1622776601683795);
        points(1) = Kratos::make_shared<Point>(9.4868329805051381, 10, 3.16227766016838);
        points(2) = Kratos::make_shared<Point>(4.7434164902525691, 10, 1.5811388300841898);
        points(3) = Kratos::make_shared<Point>(0,                  10, 0);
        points(4) = Kratos::make_shared<Point>(0,                   5, 0);
        points(5) = Kratos::make_shared<Point>(0,                   0, 0);
        points(6) = Kratos::make_shared<Point>(4.7434164902525682,  0, 1.5811388300841895);
        points(7) = Kratos::make_shared<Point>(9.4868329805051381,  0, 3.1622776601683791);
        points(8) = Kratos::make_shared<Point>(9.4868329805051381,  5, 3.1622776601683795);

        Vector knot_vector = ZeroVector(10);
        knot_vector[0] = 0.0;
        knot_vector[1] = 0.0;
        knot_vector[2] = 7.8539816339744828;
        knot_vector[3] = 7.8539816339744828;
        knot_vector[4] = 15.707963267948966;
        knot_vector[5] = 15.707963267948966;
        knot_vector[6] = 23.561944901923447;
        knot_vector[7] = 23.561944901923447;
        knot_vector[8] = 31.415926535897931;
        knot_vector[9] = 31.415926535897931;

        Vector weights = ZeroVector(9);
        weights[0] = 1.0;
        weights[1] = 0.70710678118654757;
        weights[2] = 1.0;
        weights[3] = 0.70710678118654757;
        weights[4] = 1.0;
        weights[5] = 0.70710678118654757;
        weights[6] = 1.0;
        weights[7] = 0.70710678118654757;
        weights[8] = 1.0;

        int p = 2;

        auto curve = NurbsCurveGeometry<3, PointerVector<Point>>(points, p, knot_vector, weights);

        return curve;
    }

    NurbsCurveGeometry<2, PointerVector<NodeType>> GenerateReferenceCurve2dNodes()
    {
        PointerVector<NodeType> points;

        points.push_back(NodeType::Pointer(new NodeType(1, 0, 0, 0)));
        points.push_back(NodeType::Pointer(new NodeType(2, 3.3333333333333335, 1.6666666666666667, 0)));
        points.push_back(NodeType::Pointer(new NodeType(3, 6.6666666666666661, 3.333333333333333, 0)));
        points.push_back(NodeType::Pointer(new NodeType(4, 10, 5, 0)));

        Vector knot_vector = ZeroVector(6);
        knot_vector[0] = 0.0;
        knot_vector[1] = 0.0;
        knot_vector[2] = 0.0;
        knot_vector[3] = 11.180339887498949;
        knot_vector[4] = 11.180339887498949;
        knot_vector[5] = 11.180339887498949;

        int p = 3;

        auto curve = NurbsCurveGeometry<2, PointerVector<NodeType>>(points, p, knot_vector);

        return curve;
    }

    NurbsCurveGeometry<3, PointerVector<NodeType>> GenerateReferenceCurve3d()
    {
        PointerVector<NodeType> points;

        points.push_back(NodeType::Pointer(new NodeType(1, 0, -25, -5)));
        points.push_back(NodeType::Pointer(new NodeType(2, -15, -15, 0)));
        points.push_back(NodeType::Pointer(new NodeType(3, 5, -5, -3)));
        points.push_back(NodeType::Pointer(new NodeType(4, 15, -15, 3)));
        points.push_back(NodeType::Pointer(new NodeType(5, 25, 0, 6)));
        points.push_back(NodeType::Pointer(new NodeType(6, 15, 15, 6)));
        points.push_back(NodeType::Pointer(new NodeType(7, -5, -5, -3)));
        points.push_back(NodeType::Pointer(new NodeType(8, -25, 15, 4)));

        Vector knot_vector = ZeroVector(11);
        knot_vector[0] = 0.0;
        knot_vector[1] = 0.0;
        knot_vector[2] = 0.0;
        knot_vector[3] = 0.0;
        knot_vector[4] = 32.9731425998736;
        knot_vector[5] = 65.9462851997473;
        knot_vector[6] = 98.9194277996209;
        knot_vector[7] = 131.892570399495;
        knot_vector[8] = 131.892570399495;
        knot_vector[9] = 131.892570399495;
        knot_vector[10] = 131.892570399495;

        int p = 4;

        Vector weights = ZeroVector(8);
        weights[0] = 1.0;
        weights[1] = 3.0;
        weights[2] = 1.0;
        weights[3] = 2.5;
        weights[4] = 1.0;
        weights[5] = 0.5;
        weights[6] = 1.0;
        weights[7] = 2.0;

        auto curve = NurbsCurveGeometry<3, PointerVector<NodeType>>(points, p, knot_vector, weights);

        return curve;
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurve2d, KratosCoreNurbsGeometriesFastSuite) {
        auto curve = GenerateReferenceCurve2d();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(curve.Dimension(), 1);
        KRATOS_CHECK_EQUAL(curve.WorkingSpaceDimension(), 2);
        KRATOS_CHECK_EQUAL(curve.LocalSpaceDimension(), 1);
        KRATOS_CHECK_EQUAL(curve.IsRational(), false);

        KRATOS_CHECK_EQUAL(curve.PolynomialDegree(0), 3);
        KRATOS_CHECK_EQUAL(curve.NumberOfKnots(), 6);
        KRATOS_CHECK_EQUAL(curve.PointsNumber(), 4);

        KRATOS_CHECK_EQUAL(curve.DomainInterval().GetT0(), 0);
        KRATOS_CHECK_EQUAL(curve.DomainInterval().GetT1(), 11.180339887498949);

        // Check length of 2D curve
        KRATOS_CHECK_NEAR(curve.Length(), 11.180339887498949, TOLERANCE);

        array_1d<double, 3> parameter(0.0);
        parameter[0] = 1.0;
        array_1d<double, 3> result(0.0);

        curve.GlobalCoordinates(result, parameter);

        KRATOS_CHECK_NEAR(result[0], 0.894427, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], 0.447214, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], 0.0, TOLERANCE);

        std::vector<array_1d<double, 3>> derivatives;
        curve.GlobalSpaceDerivatives(derivatives, parameter, 5);

        KRATOS_CHECK_NEAR(derivatives[4][1], 0.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsCurve3dCircleLength, KratosCoreNurbsGeometriesFastSuite) {
        auto curve = GenerateCircle();
        KRATOS_CHECK_NEAR(curve.Length(), 31.415926535897931, 1e-2);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsCurve3d, KratosCoreNurbsGeometriesFastSuite) {

        auto curve = GenerateReferenceCurve3d();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(curve.Dimension(), 1);
        KRATOS_CHECK_EQUAL(curve.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(curve.LocalSpaceDimension(), 1);
        KRATOS_CHECK_EQUAL(curve.IsRational(), true);

        KRATOS_CHECK_EQUAL(curve.PolynomialDegree(0), 4);
        KRATOS_CHECK_EQUAL(curve.NumberOfKnots(), 11);
        KRATOS_CHECK_EQUAL(curve.PointsNumber(), 8);

        KRATOS_CHECK_EQUAL(curve.DomainInterval().GetT0(), 0);
        KRATOS_CHECK_EQUAL(curve.DomainInterval().GetT1(), 131.892570399495);

        // Check the curve length
        KRATOS_CHECK_NEAR(curve.Length(), 105.464152102341, TOLERANCE);

        // Check location information

        // check point at t = 0
        {
            array_1d<double, 3> parameter(0.0);
            parameter[0] = 0.0;

            array_1d<double, 3> result(0.0);
            curve.GlobalCoordinates(result, parameter);

            KRATOS_CHECK_NEAR(result[0], 0, TOLERANCE);
            KRATOS_CHECK_NEAR(result[1], -25, TOLERANCE);
            KRATOS_CHECK_NEAR(result[2], -5, TOLERANCE);
        }

        // check derivatives at t = 0
        {
            array_1d<double, 3> parameter(0.0);
            parameter[0] = 0.0;

            std::vector<array_1d<double, 3>> derivatives;
            curve.GlobalSpaceDerivatives(derivatives, parameter, 3);

            KRATOS_CHECK_NEAR(derivatives[0][0], 0, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[0][1], -25, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[0][2], -5, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[1][0], -5.458988, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[1][1], 3.639326, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[1][2], 1.819663, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[2][0], 3.421545, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[2][1], -2.152262, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[2][2], -1.12028, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[3][0], -3.084298, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[3][1], 1.953733, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[3][2], 1.014747, TOLERANCE);
        }

        // check point at t = 65.9462851997
        {
            array_1d<double, 3> parameter(0.0);
            parameter[0] = 65.9462851997;

            array_1d<double, 3> result(0.0);
            curve.GlobalCoordinates(result, parameter);

            KRATOS_CHECK_NEAR(result[0], 17.372881, TOLERANCE);
            KRATOS_CHECK_NEAR(result[1], -10.084746, TOLERANCE);
            KRATOS_CHECK_NEAR(result[2], 3.661017, TOLERANCE);
        }

        // check derivatives at t = 65.9462851997
        {
            array_1d<double, 3> parameter(0.0);
            parameter[0] = 65.9462851997;

            std::vector<array_1d<double, 3>> derivatives;
            curve.GlobalSpaceDerivatives(derivatives, parameter, 3);

            KRATOS_CHECK_NEAR(derivatives[0][0], 17.372881, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[0][1], -10.084746, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[0][2], 3.661017, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[1][0], 0.157519, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[1][1], 0.214672, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[1][2], 0.065029, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[2][0], -0.001173, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[2][1], 0.013599, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[2][2], -0.00044, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[3][0], -0.000212, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[3][1], -0.000031, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[3][2], 0.000078, TOLERANCE);
        }

        // check point at t = 125
        {
            array_1d<double, 3> parameter(0.0);
            parameter[0] = 125;

            array_1d<double, 3> result(0.0);
            curve.GlobalCoordinates(result, parameter);

            KRATOS_CHECK_NEAR(result[0], -15.801248, TOLERANCE);
            KRATOS_CHECK_NEAR(result[1], 7.432826, TOLERANCE);
            KRATOS_CHECK_NEAR(result[2], 1.456648, TOLERANCE);
        }

        // check derivatives at t = 125
        {
            array_1d<double, 3> parameter(0.0);
            parameter[0] = 125;

            std::vector<array_1d<double, 3>> derivatives;
            curve.GlobalSpaceDerivatives(derivatives, parameter, 3);

            KRATOS_CHECK_NEAR(derivatives[0][0], -15.801248, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[0][1], 7.432826, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[0][2], 1.456648, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[1][0], -1.44436, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[1][1], 0.927174, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[1][2], 0.287317, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[2][0], 0.026303, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[2][1], 0.065941, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[2][2], 0.031612, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[3][0], 0.00346, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[3][1], -0.006864, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[3][2], -0.003416, TOLERANCE);
        }

        // check point at t = 131.892570399495
        {
            array_1d<double, 3> parameter(0.0);
            parameter[0] = 131.892570399495;

            array_1d<double, 3> result(0.0);
            curve.GlobalCoordinates(result, parameter);

            KRATOS_CHECK_NEAR(result[0], -25, TOLERANCE);
            KRATOS_CHECK_NEAR(result[1], 15, TOLERANCE);
            KRATOS_CHECK_NEAR(result[2], 4, TOLERANCE);
        }

        // check derivatives at t = 131.892570399495
        {
            array_1d<double, 3> parameter(0.0);
            parameter[0] = 131.892570399495;

            std::vector<array_1d<double, 3>> derivatives;
            curve.GlobalSpaceDerivatives(derivatives, parameter, 3);

            KRATOS_CHECK_NEAR(derivatives[0][0], -25, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[0][1], 15, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[0][2], 4, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[1][0], -1.213109, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[1][1], 1.213109, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[1][2], 0.424588, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[2][0], 0.036791, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[2][1], 0.018395, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[2][2], 0.009198, TOLERANCE);

            KRATOS_CHECK_NEAR(derivatives[3][0], 0, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[3][1], -0.005858, TOLERANCE);
            KRATOS_CHECK_NEAR(derivatives[3][2], -0.00265, TOLERANCE);
        }
    }

    ///// Test integration points of nurbs curve
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurve2dCreateIntegrationPoints, KratosCoreNurbsGeometriesFastSuite) {
        auto curve = GenerateReferenceCurve2d();

        // Check general information, input to ouput
        typename Geometry<Node<3>>::IntegrationPointsArrayType integration_points;
        curve.CreateIntegrationPoints(integration_points);

        KRATOS_CHECK_EQUAL(integration_points.size(), 4);
        double length = 0;
        for (IndexType i = 0; i < integration_points.size(); ++i) {
            length += integration_points[i].Weight();
        }
        KRATOS_CHECK_NEAR(length, 11.180339887498949, TOLERANCE);
    }


    // test quadrature points of curve on surface
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurve2dCreateQuadraturePoints, KratosCoreNurbsGeometriesFastSuite)
    {
        // Nurbs curve on a Nurbs surface
        auto curve = GenerateReferenceCurve3d();

        typename Geometry<Node<3>>::IntegrationPointsArrayType integration_points;
        curve.CreateIntegrationPoints(integration_points);

        typename Geometry<Node<3>>::GeometriesArrayType quadrature_points;
        curve.CreateQuadraturePointGeometries(quadrature_points, 3, integration_points);

        KRATOS_CHECK_EQUAL(quadrature_points.size(), 20);
        double length = 0;
        for (IndexType i = 0; i < quadrature_points.size(); ++i) {
            for (IndexType j = 0; j < quadrature_points[i].IntegrationPointsNumber(); ++j) {
                length += quadrature_points[i].IntegrationPoints()[j].Weight();
            }
        }
        KRATOS_CHECK_NEAR(length, 131.892570399495, TOLERANCE);

        auto element = Element(0, quadrature_points(2));

        // Check shape functions. This is to guarantee that the information does not get lost.
        KRATOS_CHECK_MATRIX_NEAR(
            element.pGetGeometry()->ShapeFunctionsValues(),
            quadrature_points(2)->ShapeFunctionsValues(),
            TOLERANCE);

        // Check first derivatives
        KRATOS_CHECK_MATRIX_NEAR(
            element.GetGeometry().ShapeFunctionDerivatives(1, 0),
            quadrature_points(2)->ShapeFunctionLocalGradient(0),
            TOLERANCE);

        // Check second derivatives
        KRATOS_CHECK_MATRIX_NEAR(
            element.GetGeometry().ShapeFunctionDerivatives(2, 0),
            quadrature_points(2)->ShapeFunctionDerivatives(2, 0),
            TOLERANCE);

        array_1d<double, 3> global_coords;
        array_1d<double, 3> local_coords;
        local_coords[0] = integration_points[10][0];
        local_coords[1] = integration_points[10][1];
        curve.GlobalCoordinates(global_coords, local_coords);

        KRATOS_CHECK_VECTOR_NEAR(quadrature_points[10].Center(), global_coords, TOLERANCE);
    }
} // namespace Testing.
} // namespace Kratos.
