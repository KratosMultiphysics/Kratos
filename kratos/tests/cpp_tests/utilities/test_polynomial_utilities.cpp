//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela Dalmau
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "testing/testing.h"
#include "utilities/polynomial_utilities.h"

namespace Kratos::Testing::PolynomialUtilitiesTests {

namespace {
    using Polynomial = PolynomialUtilities::PolynomialType;
    using Interval = PolynomialUtilities::IntervalType;
    using RootIntervals = std::vector<PolynomialUtilities::IntervalType>;
    constexpr double TOLERANCE = 1e-9;

    std::size_t CountIntervalsContaining(const RootIntervals& rIntervals, double Coordinate) {
        std::size_t containing = 0;
        for (const auto& range: rIntervals) {
            if (range[0] <= Coordinate && range[1] >= Coordinate) {
                ++containing;
            }
        }
        return containing;
    }
}

KRATOS_TEST_CASE_IN_SUITE(PolynomialUtilitiesDegree, KratosCoreFastSuite) {
    Polynomial p{1.0, 2.0, 3.0, 4.0}; // x^3 + 2x^2 + 3x + 4
    KRATOS_CHECK_EQUAL(PolynomialUtilities::Degree(p), 3);
}

KRATOS_TEST_CASE_IN_SUITE(PolynomialUtilitiesEvaluate, KratosCoreFastSuite) {
    Polynomial p{1.0, 2.0, 3.0, 4.0}; // x^3 + 2x^2 + 3x + 4

    KRATOS_CHECK_NEAR(PolynomialUtilities::Evaluate(p, 0.0), 4.0, TOLERANCE);
    KRATOS_CHECK_NEAR(PolynomialUtilities::Evaluate(p, 0.5), 6.125, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(PolynomialUtilitiesDifferentiate, KratosCoreFastSuite) {
    Polynomial p{1.0, 2.0, 3.0, 4.0}; // x^3 + 2x^2 + 3x + 4
    Polynomial d{3.0, 4.0, 3.0}; // 3x^2 + 4x + 3

    KRATOS_CHECK_VECTOR_NEAR(PolynomialUtilities::Differentiate(p), d, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(PolynomialUtilitiesMultiply, KratosCoreFastSuite) {
    Polynomial a{1.0, 2.0, 3.0}; // x^2 + 2x + 3
    Polynomial b{4.0, 5.0}; // 4x + 5
    Polynomial c{4.0, 13.0, 22.0, 15.0}; // 4x^3 + 13x^2 + 22x + 15

    KRATOS_CHECK_VECTOR_NEAR(PolynomialUtilities::Multiply(a, b), c, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(PolynomialUtilitiesDivide, KratosCoreFastSuite) {
    Polynomial a{1.0, 2.0, 3.0}; // x^2 + 2x + 3
    Polynomial b{4.0, 5.0}; // 4x + 5
    Polynomial q, r;

    PolynomialUtilities::Divide(q, r, a, b);

    Polynomial expected_quotient{1.0/4.0, 3.0/16.0};
    Polynomial expected_remainder{33.0/16.0};

    KRATOS_CHECK_VECTOR_NEAR(q, expected_quotient, TOLERANCE);
    KRATOS_CHECK_VECTOR_NEAR(r, expected_remainder, TOLERANCE);

    Polynomial c{4.0, 13.0, 22.0, 15.0}; // a*b
    PolynomialUtilities::Divide(q, r, c, a);
    KRATOS_CHECK_VECTOR_NEAR(q, b, TOLERANCE);
    KRATOS_CHECK_VECTOR_NEAR(r, Polynomial{0.0}, TOLERANCE);
}

KRATOS_TEST_CASE_IN_SUITE(PolynomialUtilitiesIsolateRoots, KratosCoreFastSuite) {
    Polynomial a{1.06, 0.75, -0.26, -0.175}; // Has three roots in [-1, 1]
    RootIntervals intervals;

    PolynomialUtilities::IsolateRoots(intervals, a, Interval{-1, 1});
    KRATOS_CHECK_EQUAL(intervals.size(), 3);
    KRATOS_CHECK_EQUAL(CountIntervalsContaining(intervals, -0.7360625845831237), 1);
    KRATOS_CHECK_EQUAL(CountIntervalsContaining(intervals, -0.45955361708183773), 1);
    KRATOS_CHECK_EQUAL(CountIntervalsContaining(intervals,  0.4880690318514098), 1);

    PolynomialUtilities::IsolateRoots(intervals, a, Interval{0, 1});
    KRATOS_CHECK_EQUAL(intervals.size(), 1);
    KRATOS_CHECK_EQUAL(CountIntervalsContaining(intervals, 0.4880690318514098), 1);

    PolynomialUtilities::IsolateRoots(intervals, a, Interval{-2, -1});
    KRATOS_CHECK_EQUAL(intervals.size(), 0);
}

KRATOS_TEST_CASE_IN_SUITE(PolynomialUtilitiesFindRoot, KratosCoreFastSuite) {
    Polynomial a{1.06, 0.75, -0.26, -0.175}; // Has three roots in [-1, 1]

    KRATOS_CHECK_NEAR(PolynomialUtilities::FindRoot(a, Interval{-1, -0.5}), -0.7360625845831237, TOLERANCE);
    KRATOS_CHECK_NEAR(PolynomialUtilities::FindRoot(a, Interval{-0.5, 0}), -0.45955361708183773, TOLERANCE);
    KRATOS_CHECK_NEAR(PolynomialUtilities::FindRoot(a, Interval{0, 1}), 0.4880690318514098, TOLERANCE);
}

}