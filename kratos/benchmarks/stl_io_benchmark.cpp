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
#include "containers/model.h"
#include "input_output/stl_io.h"

namespace Kratos
{

static void BM_StlIO(benchmark::State& state) {
    // 1. Prepare the raw data string outside the loop so it's not measured
    const std::string stl_data = R"input(
        solid CUBE
            facet normal 0 0 1
            outer loop
                vertex -35 60 20
                vertex -55 60 20
                vertex -35 40 20
            endloop
            endfacet
            facet normal 0 0 1
            outer loop
                vertex -35 40 20
                vertex -55 60 20
                vertex -55 40 20
            endloop
            endfacet
            facet normal 0 0 -1
            outer loop
                vertex -35 40 0
                vertex -55 40 0
                vertex -35 60 0
            endloop
            endfacet
            facet normal -0 0 -1
            outer loop
                vertex -35 60 0
                vertex -55 40 0
                vertex -55 60 0
            endloop
            endfacet
            facet normal 0 -1 0
            outer loop
                vertex -55 40 20
                vertex -55 40 0
                vertex -35 40 20
            endloop
            endfacet
            facet normal 0 -1 -0
            outer loop
                vertex -35 40 20
                vertex -55 40 0
                vertex -35 40 0
            endloop
            endfacet
            facet normal -1 -0 -0
            outer loop
                vertex -55 60 20
                vertex -55 60 0
                vertex -55 40 20
            endloop
            endfacet
            facet normal -1 0 0
            outer loop
                vertex -55 40 20
                vertex -55 60 0
                vertex -55 40 0
            endloop
            endfacet
            facet normal 0 1 0
            outer loop
                vertex -35 60 20
                vertex -35 60 0
                vertex -55 60 20
            endloop
            endfacet
            facet normal 0 1 0
            outer loop
                vertex -55 60 20
                vertex -35 60 0
                vertex -55 60 0
            endloop
            endfacet
            facet normal 1 -0 0
            outer loop
                vertex -35 40 20
                vertex -35 40 0
                vertex -35 60 20
            endloop
            endfacet
            facet normal 1 0 0
            outer loop
                vertex -35 60 20
                vertex -35 40 0
                vertex -35 60 0
            endloop
            endfacet
        endsolid CUBE
        )input";

    Model current_model;

    for (auto _ : state) {
        // 2. Setup that we don't want to profile as "IO Time"
        state.PauseTiming();
        // We create a new ModelPart to ensure we start from a clean slate
        // and aren't measuring the cost of searching through existing nodes.
        ModelPart& r_model_part = current_model.CreateModelPart("Main_" + std::to_string(state.iterations()));
        auto p_input = Kratos::make_shared<std::stringstream>(stl_data);
        StlIO stl_io(p_input);
        state.ResumeTiming();

        // 3. The actual operation being benchmarked
        stl_io.ReadModelPart(r_model_part);

        // Cleanup to prevent memory bloat during high-iteration benchmarks
        state.PauseTiming();
        current_model.DeleteModelPart("Main_" + std::to_string(state.iterations()));
        state.ResumeTiming();
    }
}

BENCHMARK(BM_StlIO);

}  // namespace Kratos

BENCHMARK_MAIN();