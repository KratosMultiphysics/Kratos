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
#include <ranges>
#include <execution>
#include <algorithm>

// External includes

// Project includes
#include "testing/testing.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"


namespace Kratos::Testing {

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

#ifdef KRATOS_USE_TBB
    auto range = std::ranges::views::iota(0uz, size);
    std::for_each(std::execution::par_unseq, range.begin(), range.end(), [&sum, exp](int x) {
        AtomicMult(sum, exp);
    });
#else
    IndexPartition(size).for_each(
        [&sum, exp](auto i) {
            AtomicMult(sum, exp);
        }
    );
#endif

    KRATOS_EXPECT_NEAR(5 * std::pow(exp, size), sum, 1e-3);
}

}  // namespace Kratos::Testing
