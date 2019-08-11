//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "geometries/nurbs_curve_geometry.h"

#include "tests/cpp_tests/geometries/test_geometry.h"

namespace Kratos {
namespace Testing {

typedef Node<3> NodeType;

// /// Factory functions
//namespace {
    NurbsCurveGeometry<2, Point> GenerateReferenceCurve2d()
    {
        NurbsCurveGeometry<2, Point>::PointsArrayType points;

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

        auto curve = NurbsCurveGeometry<2, Point>(points, p, knot_vector);

        return curve;
    }

    NurbsCurveGeometry<3, NodeType> GenerateReferenceCurve3d()
    {
        NurbsCurveGeometry<3, NodeType>::PointsArrayType points;

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
        weights[1] = 1.0;
        weights[2] = 1.0;
        weights[3] = 2.5;
        weights[4] = 1.0;
        weights[5] = 1.0;
        weights[6] = 1.0;
        weights[7] = 1.0;

        auto curve = NurbsCurveGeometry<3, NodeType>(points, p, knot_vector, weights);

        return curve;
    }

    ///// Tests
    KRATOS_TEST_CASE_IN_SUITE(NurbsCurve2d, KratosCoreGeometriesFastSuite) {
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

        auto derivatives = curve.GlobalDerivatives(parameter, 5);

        KRATOS_CHECK_NEAR(derivatives[4][1], 0.0, TOLERANCE);
    }

    KRATOS_TEST_CASE_IN_SUITE(NurbsCurve3d, KratosCoreGeometriesFastSuite) {

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

        array_1d<double, 3> result(0.0);

        // Check location information

        // parameter t=0.0
        array_1d<double, 3> parameter(0.0);
        parameter[0] = 0.0;
        curve.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], 0, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], -25, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], -5, TOLERANCE);

        const auto derivatives_1 = curve.GlobalDerivatives(parameter, 3);
        KRATOS_CHECK_NEAR(derivatives_1[0][0], 0, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_1[0][1], -25, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_1[0][2], -5, TOLERANCE);

        KRATOS_CHECK_NEAR(derivatives_1[1][0], -1.81966277, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_1[1][1], 1.2131085134, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_1[1][2], 0.6065542567, TOLERANCE);

        KRATOS_CHECK_NEAR(derivatives_1[2][0], 0.2759310497, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_1[2][1], -0.0551862099, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_1[2][2], -0.0717420729, TOLERANCE);

        KRATOS_CHECK_NEAR(derivatives_1[3][0], -0.0189682773, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_1[3][1], 0.0005578905, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_1[3][2], 0.005523116, TOLERANCE);


        // parameter t=65.9462851997
        parameter[0] = 65.9462851997;
        result = ZeroVector(3);

        curve.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], 21.333333, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], -3.6666667, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], 4.9, TOLERANCE);

        const auto derivatives_2 = curve.GlobalDerivatives(parameter, 3);
        KRATOS_CHECK_NEAR(derivatives_2[0][0], 21.33333333, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[0][1], -3.66666667, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[0][2], 4.9, TOLERANCE);

        KRATOS_CHECK_NEAR(derivatives_2[1][0], 0.20218475, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[1][1], 0.33697459, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[1][2], 0.10109238, TOLERANCE);

        KRATOS_CHECK_NEAR(derivatives_2[2][0], -0.0122636, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[2][1], 0.0153295, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[2][2], -0.00367908, TOLERANCE);

        KRATOS_CHECK_NEAR(derivatives_2[3][0], -5.57890509e-04, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[3][1], -6.50872261e-04, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_2[3][2], 5.57890509e-05, TOLERANCE);

        // parameter t=131.892570399495
        parameter[0] = 131.892570399495;
        result = ZeroVector(3);

        curve.GlobalCoordinates(result, parameter);
        KRATOS_CHECK_NEAR(result[0], -25, TOLERANCE);
        KRATOS_CHECK_NEAR(result[1], 15, TOLERANCE);
        KRATOS_CHECK_NEAR(result[2], 4, TOLERANCE);

        const auto derivatives_3 = curve.GlobalDerivatives(parameter, 3);
        KRATOS_CHECK_NEAR(derivatives_3[0][0], -25, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_3[0][1], 15, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_3[0][2], 4, TOLERANCE);

        KRATOS_CHECK_NEAR(derivatives_3[1][0], -2.4262170267, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_3[1][1], 2.4262170267, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_3[1][2], 0.8491759593, TOLERANCE);

        KRATOS_CHECK_NEAR(derivatives_3[2][0], -0.1103724199, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_3[2][1], 0.3311172597, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_3[2][2], 0.1269282829, TOLERANCE);

        KRATOS_CHECK_NEAR(derivatives_3[3][0], -0.0044631241, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_3[3][1], 0.0251050729, TOLERANCE);
        KRATOS_CHECK_NEAR(derivatives_3[3][2], 0.0092051934, TOLERANCE);
    }

} // namespace Testing.
} // namespace Kratos.
