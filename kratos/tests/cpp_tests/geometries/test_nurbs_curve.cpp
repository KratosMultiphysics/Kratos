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

        KRATOS_CHECK_EQUAL(curve.PolynomialDegree(), 3);
        KRATOS_CHECK_EQUAL(curve.NumberOfKnots(), 6);
        KRATOS_CHECK_EQUAL(curve.PointsNumber(), 4);

        KRATOS_CHECK_EQUAL(curve.DomainInterval().GetT0(), 0);
        KRATOS_CHECK_EQUAL(curve.DomainInterval().GetT1(), 11.180339887498949);

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

    KRATOS_TEST_CASE_IN_SUITE(NurbsCurve3d, KratosCoreNurbsGeometriesFastSuite) {

        auto curve = GenerateReferenceCurve3d();

        // Check general information, input to ouput
        KRATOS_CHECK_EQUAL(curve.Dimension(), 1);
        KRATOS_CHECK_EQUAL(curve.WorkingSpaceDimension(), 3);
        KRATOS_CHECK_EQUAL(curve.LocalSpaceDimension(), 1);
        KRATOS_CHECK_EQUAL(curve.IsRational(), true);

        KRATOS_CHECK_EQUAL(curve.PolynomialDegree(), 4);
        KRATOS_CHECK_EQUAL(curve.NumberOfKnots(), 11);
        KRATOS_CHECK_EQUAL(curve.PointsNumber(), 8);

        KRATOS_CHECK_EQUAL(curve.DomainInterval().GetT0(), 0);
        KRATOS_CHECK_EQUAL(curve.DomainInterval().GetT1(), 131.892570399495);

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

} // namespace Testing.
} // namespace Kratos.
