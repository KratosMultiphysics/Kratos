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
#include <utility>
#include <numeric>
#include <iostream>
#include <unordered_map>

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

namespace Kratos
{
// Template class for testing
template<std::size_t TSize>
class RHSElement {
public:
    explicit RHSElement(const double Val) : mRHSVal(Val) {}
    void CalculateRHS(std::vector<double>& rVector) {
        if (rVector.size() != TSize) { rVector.resize(TSize); }
        std::fill(rVector.begin(), rVector.end(), mRHSVal);
    }
    double GetAccumRHSValue() { return mAccumRHSValue; }
    void SetAccumRHSValue(double Value) { mAccumRHSValue = Value; }

private:
    double mRHSVal;
    double mAccumRHSValue = 0.0;
};

// Benchmark for power operation on a vector
static void BM_VectorPower(benchmark::State& state) {
    int nsize = state.range(0);
    std::vector<double> data_vector(nsize, 5.0);

    for (auto _ : state) {
        block_for_each(data_vector, [](double& item) {
            item = std::pow(item, 0.1);
        });
    }
}

// Benchmark for reduction
static void BM_VectorReduction(benchmark::State& state) {
    int nsize = state.range(0);
    std::vector<double> data_vector(nsize, 5.0);

    for (auto _ : state) {
        auto final_sum = BlockPartition<std::vector<double>::iterator>(
            data_vector.begin(), data_vector.end()
        ).for_each<SumReduction<double>>([](double& item) {
            return item;
        });

        benchmark::DoNotOptimize(final_sum); // <-- Fixes [[nodiscard]] warning
    }
}

// Benchmark for element-wise operations with thread-local storage
static void BM_ThreadLocalStorage(benchmark::State& state) {
    constexpr std::size_t vec_size = 6;
    std::size_t n_elems = state.range(0);

    using RHSElementType = RHSElement<vec_size>;

    std::vector<double> rhs_vals(n_elems);
    for (std::size_t i = 0; i < n_elems; ++i) {
        rhs_vals[i] = (i % 12) * 1.889;
    }

    std::vector<RHSElementType> elements;
    for (std::size_t i = 0; i < rhs_vals.size(); ++i) {
        elements.push_back(RHSElementType(rhs_vals[i]));
    }

    auto tls_lambda_manual_reduction = [](RHSElementType& rElem, std::vector<double>& rTLS)
    {
        rElem.CalculateRHS(rTLS);
        double rhs_sum = std::accumulate(rTLS.begin(), rTLS.end(), 0.0);
        rElem.SetAccumRHSValue(rhs_sum);
    };

    for (auto _ : state) {
        BlockPartition<std::vector<RHSElementType>::iterator>(
            elements.begin(), elements.end()
        ).for_each(std::vector<double>(), tls_lambda_manual_reduction);

        const double sum_elem_rhs_vals = std::accumulate(
            elements.begin(), elements.end(), 0.0,
            [](double acc, RHSElementType& rElem) {
                return acc + rElem.GetAccumRHSValue();
            }
        );

        benchmark::DoNotOptimize(sum_elem_rhs_vals);
    }
}

// Register benchmarks and provide input size as a command-line option
BENCHMARK(BM_VectorPower)->Arg(1e3)->Arg(1e5)->Arg(1e6);
BENCHMARK(BM_VectorReduction)->Arg(1e3)->Arg(1e5)->Arg(1e6);
BENCHMARK(BM_ThreadLocalStorage)->Arg(1e3)->Arg(1e5)->Arg(1e6);

}  // namespace Kratos

BENCHMARK_MAIN();