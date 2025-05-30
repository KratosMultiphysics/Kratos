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


} // namespace Kratos::Testing
