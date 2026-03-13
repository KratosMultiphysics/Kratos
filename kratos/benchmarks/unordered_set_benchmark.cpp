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
#include <random>

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "containers/unordered_set.h"

namespace Kratos
{

// Function to generate random integers
std::vector<int> generate_random_integers(size_t size) {
    std::vector<int> data(size);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(1, 1000000);
    for (size_t i = 0; i < size; ++i) {
        data[i] = dis(gen);
    }
    return data;
}

// Benchmark for insertion
static void BM_Unordered_Set_Insert(benchmark::State& state) {
    std::vector<int> data = generate_random_integers(state.range(0));
    for (auto _ : state) {
        Kratos::unordered_set<int> uset;
        for (const auto& val : data) {
            uset.insert(val);
        }
    }
}
BENCHMARK(BM_Unordered_Set_Insert)->Range(1 << 10, 1 << 20);

// Benchmark for lookup
static void BM_Unordered_Set_Lookup(benchmark::State& state) {
    std::vector<int> data = generate_random_integers(state.range(0));
    Kratos::unordered_set<int> uset(data.begin(), data.end());
    for (auto _ : state) {
        for (const auto& val : data) {
            benchmark::DoNotOptimize(uset.find(val));
        }
    }
}
BENCHMARK(BM_Unordered_Set_Lookup)->Range(1 << 10, 1 << 20);

// Benchmark for deletion
static void BM_Unordered_Set_Delete(benchmark::State& state) {
    std::vector<int> data = generate_random_integers(state.range(0));
    for (auto _ : state) {
        Kratos::unordered_set<int> uset(data.begin(), data.end());
        for (const auto& val : data) {
            uset.erase(val);
        }
    }
}
BENCHMARK(BM_Unordered_Set_Delete)->Range(1 << 10, 1 << 20);


// Benchmark for iteration
static void BM_Unordered_Set_Iteration(benchmark::State& state) {
    std::vector<int> data = generate_random_integers(state.range(0));
    Kratos::unordered_set<int> uset(data.begin(), data.end());
    for (auto _ : state) {
        for (const auto& val : uset) {
            benchmark::DoNotOptimize(val);
        }
    }
}
BENCHMARK(BM_Unordered_Set_Iteration)->Range(1 << 10, 1 << 20);

}  // namespace Kratos

BENCHMARK_MAIN();