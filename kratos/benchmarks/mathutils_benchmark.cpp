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

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "utilities/math_utils.h"
#include "includes/ublas_interface.h"  // Provides ZeroVector and Vector

namespace Kratos
{

// Template benchmark for the Dot product with different vector sizes
template <std::size_t Size>
static void BM_MathUtilsDot(benchmark::State& state) {
    // Initialize two vectors of the given size
    Vector a = ZeroVector(Size);
    Vector b = ZeroVector(Size);

    // Set some values to avoid zero vector dot product
    if (Size > 1) {
        a[1] = 1.0;  // Only the second component is non-zero
        b[0] = 1.0;  // Only the first component is non-zero
    }

    double result = 0.0;
    for (auto _ : state) {
        // Compute the dot product; result should be 0.0 for these specific vectors
        result = MathUtils<double>::Dot(a, b);
        // Prevent the compiler from optimizing away the computation
        benchmark::DoNotOptimize(result);
    }
}

// Register the function as benchmarks for different sizes
BENCHMARK_TEMPLATE(BM_MathUtilsDot, 3);
BENCHMARK_TEMPLATE(BM_MathUtilsDot, 10);
BENCHMARK_TEMPLATE(BM_MathUtilsDot, 100);
BENCHMARK_TEMPLATE(BM_MathUtilsDot, 1000);

}  // namespace Kratos

BENCHMARK_MAIN();
