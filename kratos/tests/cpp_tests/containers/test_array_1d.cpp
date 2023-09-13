//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//                   Vicente Mataix Ferrandiz
//
//

// System includes
#include <functional>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "testing/testing.h"

namespace Kratos::Testing 
{

/**
 * Test case for initializing a 1D array with a specific value.
 */
KRATOS_TEST_CASE_IN_SUITE(Array1DInitializationValue, KratosCoreFastSuite) 
{
    const array_1d<double, 3> arr(2, 2.2);
    Vector ref(3);
    ref[0] = 2.2;
    ref[1] = 2.2;
    KRATOS_EXPECT_DOUBLE_EQ(arr[0], ref[0]);
    KRATOS_EXPECT_DOUBLE_EQ(arr[1], ref[1]);
}

/**
 * A test case to check if two 1D arrays are equal or not.
 */
KRATOS_TEST_CASE_IN_SUITE(Array1DTest, KratosCoreFastSuite) 
{
    const array_1d<double, 3> arr1{1.0, 2.0, 3.0};
    const array_1d<double, 3> arr2{1.0, 2.0, 3.0};
    const array_1d<double, 3> arr3{1.0, 2.0, 4.0};
    KRATOS_EXPECT_TRUE(arr1 == arr2);
    KRATOS_EXPECT_FALSE(arr1 == arr3);
}

/**
 * Tests that the size function returns the expected value for a given array.
 */
KRATOS_TEST_CASE_IN_SUITE(SizeFunction, KratosCoreFastSuite)
{
    const array_1d<double, 3> arr{1.0, 2.0, 3.0};
    KRATOS_EXPECT_EQ(arr.size(), 3);
}

/**
 * Test case for the index operator.
 */
KRATOS_TEST_CASE_IN_SUITE(IndexOperator, KratosCoreFastSuite)
{
    const array_1d<double, 3> arr{1.0, 2.0, 3.0};
    KRATOS_EXPECT_EQ(arr[0], 1.0);
    KRATOS_EXPECT_EQ(arr[1], 2.0);
    KRATOS_EXPECT_EQ(arr[2], 3.0);
}

/**
 * Tests the initialization of a 1D array with given values.
 */
KRATOS_TEST_CASE_IN_SUITE(Array1DInitializerList, KratosCoreFastSuite) 
{
    const array_1d<double, 3> arr {1.1, -2.3, 3.4};
    Vector ref(3);
    ref[0] = 1.1;
    ref[1] = -2.3;
    ref[2] = 3.4;
    KRATOS_EXPECT_VECTOR_EQ(arr, ref);
}

/**
 * Tests the addition operator between two arrays of 3 doubles.
 */
KRATOS_TEST_CASE_IN_SUITE(AdditionOperator, KratosCoreFastSuite) 
{
    const array_1d<double, 3> arr1{1.0, 2.0, 3.0};
    const array_1d<double, 3> arr2{4.0, 5.0, 6.0};
    const array_1d<double, 3> sum = arr1 + arr2;
    const array_1d<double, 3> expected_sum{5.0, 7.0, 9.0};
    KRATOS_EXPECT_VECTOR_EQ(sum, expected_sum);
}

/**
 * Tests the substraction operator between two arrays of 3 doubles.
 */
KRATOS_TEST_CASE_IN_SUITE(SubtractionOperator, KratosCoreFastSuite)
{
    const array_1d<double, 3> arr1{1.0, 2.0, 3.0};
    const array_1d<double, 3> arr2{4.0, 5.0, 6.0};
    const array_1d<double, 3> diff = arr2 - arr1;
    const array_1d<double, 3> expected_diff{3.0, 3.0, 3.0};
    KRATOS_EXPECT_VECTOR_EQ(diff, expected_diff);
}

/**
 * Test case to check if a Kratos array of doubles is hashed correctly.
 */
KRATOS_TEST_CASE_IN_SUITE(HashesCorrectly, KratosCoreFastSuite) 
{
    const array_1d<double, 3> arr1{1.0, 2.0, 3.0};
    const array_1d<double, 3> arr2{1.0, 2.0, 3.0};
    const array_1d<double, 3> arr3{1.0, 2.0, 4.0};
    std::hash<array_1d<double, 3>> hasher;
    const std::size_t hash_value = hasher(arr1);
    const std::size_t expected_hash_value = hasher(arr2);
    const std::size_t not_expected_hash_value = hasher(arr3);
    KRATOS_EXPECT_EQ(hash_value, expected_hash_value);
    KRATOS_EXPECT_NE(hash_value, not_expected_hash_value);
}

} // namespace Kratos::Testing.
