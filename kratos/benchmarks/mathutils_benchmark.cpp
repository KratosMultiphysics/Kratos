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

// Benchmark for the Dot product
static void BM_MathUtilsDot(benchmark::State& state) {
    // Initialize two 3D vectors
    Vector a = ZeroVector(3);
    a[1] = 1.0;  // Only the second component is non-zero
    Vector b = ZeroVector(3);
    b[0] = 1.0;  // Only the first component is non-zero

    double result = 0.0;
    for (auto _ : state) {
        // Compute the dot product; result should be 0.0
        result = MathUtils<double>::Dot(a, b);
        // Prevent the compiler from optimizing away the computation
        benchmark::DoNotOptimize(result);
    }
}

// Register the function as a benchmark
BENCHMARK(BM_MathUtilsDot);

}  // namespace Kratos

BENCHMARK_MAIN();
