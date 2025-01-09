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
//

// System includes
#include <vector>

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"


namespace Kratos::Testing {

/// @details Loop over an index range larger than an array, overwriting
///          its components with identical values from all threads.
KRATOS_TEST_CASE(AtomicSet)
{
    constexpr std::size_t array_size = 5;
    constexpr std::size_t iteration_count = 1e5;
    std::vector<std::int8_t> array(array_size, 0);

    IndexPartition<std::size_t>(iteration_count).for_each([&array](const std::size_t i){
        AtomicSet(array[i % array_size], i % array_size);
    });

    for (unsigned i=0u; i<array_size; ++i) {
        KRATOS_EXPECT_EQ(array[i], i);
    } 
}

KRATOS_TEST_CASE(AtomicAdd)
{
    constexpr std::size_t size = 12345;
    double sum = 0;

    IndexPartition(size).for_each(
        [&sum](auto i){
            AtomicAdd(sum, 1.0);
            }
        );

    KRATOS_EXPECT_DOUBLE_EQ(static_cast<double>(size), sum);
}

KRATOS_TEST_CASE(AtomicSub)
{
    constexpr std::size_t size = 12345;
    double sum = 0;

    IndexPartition(size).for_each(
        [&sum](auto i){
            AtomicSub(sum, 1.0);
            }
        );

    KRATOS_EXPECT_DOUBLE_EQ(static_cast<double>(size), -sum);
}

KRATOS_TEST_CASE(AtomicMult)
{
    constexpr std::size_t size = 12345;
    const double exp = 1.001;
    double sum = 5;

    IndexPartition(size).for_each(
        [&sum, exp](auto i){
            AtomicMult(sum, exp);
            }
        );

    KRATOS_EXPECT_NEAR(5 * std::pow(exp, size), sum, 1e-3);
}

}  // namespace Kratos::Testing
