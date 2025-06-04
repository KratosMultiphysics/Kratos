//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

// Project includes
#include "testing/testing.h"
#include "utilities/comparison.h" // Comparison


namespace Kratos::Testing {


KRATOS_TEST_CASE_IN_SUITE(IntegerComparison, KratosCoreFastSuite)
{
    Comparison<int>::Equal equality_comparison;
    Comparison<int>::Less ordering;

    for (const auto& [left, right, equality_reference, ordering_reference] :
        std::vector<std::tuple<int,int,bool,bool>> {{-4, -4,  true, false},
                                                    {-2, -4, false, false},
                                                    {-1, -4, false, false},
                                                    { 0, -4, false, false},
                                                    { 1, -4, false, false},
                                                    { 2, -4, false, false},
                                                    { 4, -4, false, false},

                                                    {-4, -2, false,  true},
                                                    {-2, -2,  true, false},
                                                    {-1, -2, false, false},
                                                    { 0, -2, false, false},
                                                    { 1, -2, false, false},
                                                    { 2, -2, false, false},
                                                    { 4, -2, false, false},

                                                    {-4, -1, false,  true},
                                                    {-2, -1, false,  true},
                                                    {-1, -1,  true, false},
                                                    { 0, -1, false, false},
                                                    { 1, -1, false, false},
                                                    { 2, -1, false, false},
                                                    { 4, -1, false, false},

                                                    {-4,  0, false,  true},
                                                    {-2,  0, false,  true},
                                                    {-1,  0, false,  true},
                                                    { 0,  0,  true, false},
                                                    { 1,  0, false, false},
                                                    { 2,  0, false, false},
                                                    { 4,  0, false, false},

                                                    {-4,  1, false,  true},
                                                    {-2,  1, false,  true},
                                                    {-1,  1, false,  true},
                                                    { 0,  1, false,  true},
                                                    { 1,  1,  true, false},
                                                    { 2,  1, false, false},
                                                    { 4,  1, false, false},

                                                    {-4,  2, false,  true},
                                                    {-2,  2, false,  true},
                                                    {-1,  2, false,  true},
                                                    { 0,  2, false,  true},
                                                    { 1,  2, false,  true},
                                                    { 2,  2,  true, false},
                                                    { 4,  2, false, false},

                                                    {-4,  4, false,  true},
                                                    {-2,  4, false,  true},
                                                    {-1,  4, false,  true},
                                                    { 0,  4, false,  true},
                                                    { 1,  4, false,  true},
                                                    { 2,  4, false,  true},
                                                    { 4,  4,  true, false}}) {
        KRATOS_EXPECT_TRUE(equality_comparison(left, right) == equality_reference);
        KRATOS_EXPECT_TRUE(ordering(left, right) == ordering_reference);
    }
}


KRATOS_TEST_CASE_IN_SUITE(FloatComparison, KratosCoreFastSuite)
{
    Comparison<double>::Equal equality_comparison(1.0, 0.5);
    Comparison<double>::Less ordering(1.0, 0.5);

    for (const auto& [left, right, equality_reference, ordering_reference] :
        std::vector<std::tuple<double,double,bool,bool>> {{-4.0, -4.0,  true, false},
                                                          {-2.0, -4.0,  true, false},
                                                          {-1.0, -4.0, false, false},
                                                          { 0.0, -4.0, false, false},
                                                          { 1.0, -4.0, false, false},
                                                          { 2.0, -4.0, false, false},
                                                          { 4.0, -4.0, false, false},

                                                          {-4.0, -2.0,  true, false},
                                                          {-2.0, -2.0,  true, false},
                                                          {-1.0, -2.0,  true, false},
                                                          { 0.0, -2.0, false, false},
                                                          { 1.0, -2.0, false, false},
                                                          { 2.0, -2.0, false, false},
                                                          { 4.0, -2.0, false, false},

                                                          {-4.0, -1.0, false,  true},
                                                          {-2.0, -1.0,  true, false},
                                                          {-1.0, -1.0,  true, false},
                                                          { 0.0, -1.0, false, false},
                                                          { 1.0, -1.0, false, false},
                                                          { 2.0, -1.0, false, false},
                                                          { 4.0, -1.0, false, false},

                                                          {-4.0,  0.0, false,  true},
                                                          {-2.0,  0.0, false,  true},
                                                          {-1.0,  0.0, false,  true},
                                                          { 0.0,  0.0,  true, false},
                                                          { 1.0,  0.0, false, false},
                                                          { 2.0,  0.0, false, false},
                                                          { 4.0,  0.0, false, false},

                                                          {-4.0,  1.0, false,  true},
                                                          {-2.0,  1.0, false,  true},
                                                          {-1.0,  1.0, false,  true},
                                                          { 0.0,  1.0, false,  true},
                                                          { 1.0,  1.0,  true, false},
                                                          { 2.0,  1.0,  true, false},
                                                          { 4.0,  1.0, false, false},

                                                          {-4.0,  2.0, false,  true},
                                                          {-2.0,  2.0, false,  true},
                                                          {-1.0,  2.0, false,  true},
                                                          { 0.0,  2.0, false,  true},
                                                          { 1.0,  2.0,  true, false},
                                                          { 2.0,  2.0,  true, false},
                                                          { 4.0,  2.0,  true, false},

                                                          {-4.0,  4.0, false,  true},
                                                          {-2.0,  4.0, false,  true},
                                                          {-1.0,  4.0, false,  true},
                                                          { 0.0,  4.0, false,  true},
                                                          { 1.0,  4.0, false,  true},
                                                          { 2.0,  4.0,  true, false},
                                                          { 4.0,  4.0,  true, false}}) {
        KRATOS_EXPECT_TRUE(equality_comparison(left, right) == equality_reference);
        KRATOS_EXPECT_TRUE(ordering(left, right) == ordering_reference);
    }
}


KRATOS_TEST_CASE_IN_SUITE(FloatComparisonConsistency, KratosCoreFastSuite)
{
    const float absolute_tolerance = std::numeric_limits<float>::min();
    const float relative_tolerance = 1e-4f;

    Comparison<float>::Equal equality_comparison(absolute_tolerance, relative_tolerance);

    const std::array<float,3> test_values {4.9303807e-32f,  //< reference value
                                           4.9303810e-32f,
                                           4.9309825e-32f};

    // Example of an inconsistent comparison.
    {
        const auto inconsistent_comparison = [absolute_tolerance, relative_tolerance](float left, float right) {
            const float diff = std::abs(left - right);
            if (left == 0.0f || right == 0.0f || diff < absolute_tolerance) {
                return diff < relative_tolerance * absolute_tolerance;
            } else {
                return diff < relative_tolerance * (std::abs(left) + std::abs(right));
            }
        };

        // Even though these values are very close to each other, they are already
        // too far apart for the absolute tolerance, but not far enough for
        // the relative comparison to kick in yet.
        KRATOS_EXPECT_FALSE(inconsistent_comparison(test_values[0], test_values[1]));

        // Even though these values are farther apart than the previous two, they're
        // picked up by the relative comparison.
        KRATOS_EXPECT_TRUE(inconsistent_comparison(test_values[0], test_values[2]));
    }

    // Make sure the implemented comparison deals with the case
    // the inconsistent one fails at.
    KRATOS_EXPECT_TRUE(equality_comparison(test_values[0], test_values[1]));
    KRATOS_EXPECT_TRUE(equality_comparison(test_values[0], test_values[2]));

    // Classic example of finite precision inaccuracies.
    KRATOS_EXPECT_TRUE(3.0f * 0.3f != 1.0f - 0.1f
                    && equality_comparison(3.0f * 0.3f, 1.0f - 0.1f));
}


} // namespace Kratos::Testing
