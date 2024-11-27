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
#include <vector>

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "containers/unordered_map.h"

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
static void BM_Unordered_Map_Insert(benchmark::State& state) {
    std::vector<int> data = generate_random_integers(state.range(0));
    for (auto _ : state) {
        Kratos::unordered_map<int, int> umap;
        for (const auto& val : data) {
            umap[val] = val;
        }
    }
}
BENCHMARK(BM_Unordered_Map_Insert)->Range(1 << 10, 1 << 20);

// Benchmark for lookup
static void BM_Unordered_Map_Lookup(benchmark::State& state) {
    std::vector<int> data = generate_random_integers(state.range(0));
    Kratos::unordered_map<int, int> umap;
    for (const auto& val : data) {
        umap[val] = val;
    }
    for (auto _ : state) {
        for (const auto& val : data) {
            benchmark::DoNotOptimize(umap.find(val));
        }
    }
}
BENCHMARK(BM_Unordered_Map_Lookup)->Range(1 << 10, 1 << 20);

// Benchmark for deletion
static void BM_Unordered_Map_Delete(benchmark::State& state) {
    std::vector<int> data = generate_random_integers(state.range(0));
    for (auto _ : state) {
        Kratos::unordered_map<int, int> umap;
        for (const auto& val : data) {
            umap[val] = val;
        }
        for (const auto& val : data) {
            umap.erase(val);
        }
    }
}
BENCHMARK(BM_Unordered_Map_Delete)->Range(1 << 10, 1 << 20);

// Benchmark for rehashing
static void BM_Unordered_Map_Rehash(benchmark::State& state) {
    std::vector<int> data = generate_random_integers(state.range(0));
    Kratos::unordered_map<int, int> umap;
    for (auto _ : state) {
        umap.insert(data.begin(), data.end());
        umap.rehash(umap.size() * 2);
    }
}
BENCHMARK(BM_Unordered_Map_Rehash)->Range(1 << 10, 1 << 20);

// Benchmark for iteration
static void BM_Unordered_Map_Iteration(benchmark::State& state) {
    std::vector<int> data = generate_random_integers(state.range(0));
    Kratos::unordered_map<int, int> umap;
    for (const auto& val : data) {
        umap[val] = val;
    }
    for (auto _ : state) {
        for (const auto& pair : umap) {
            benchmark::DoNotOptimize(pair);
        }
    }
}
BENCHMARK(BM_Unordered_Map_Iteration)->Range(1 << 10, 1 << 20);

}  // namespace Kratos

BENCHMARK_MAIN();