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
#include <fstream>
#include <stdexcept>

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "containers/model.h"
#include "input_output/stl_io.h"

namespace Kratos
{

const std::string file_name = "file.stl";

static void BM_StlIO(benchmark::State& state) {
    // Reading the file once to ensure it exists and is accessible, and to get past any initial overhead of file access.
    Kratos::shared_ptr<std::iostream> p_input_file =
        Kratos::make_shared<std::fstream>(file_name, std::ios::in);

    if (!p_input_file || !(*p_input_file)) {
        state.SkipWithError("Could not open file.stl in current folder.");
        return;
    }

    std::dynamic_pointer_cast<std::fstream>(p_input_file)->close();

    // Settings to read the STL as elements instead of conditions, which is a common use case and more expensive than reading as conditions.
    Parameters settings(R"({
        "new_entity_type" : "element"
    })");

    for (auto _ : state) {
        // Setup that we don't want to profile as "IO Time"
        // state.PauseTiming();
        // We create a new ModelPart to ensure we start from a clean slate
        // and aren't measuring the cost of searching through existing nodes.
        Model current_model;
        ModelPart& r_model_part = current_model.CreateModelPart("Main");

        // state.ResumeTiming();

        // The actual operation being benchmarked
        StlIO stl_io(p_input_file, settings);
        stl_io.ReadModelPart(r_model_part);
    }
}

BENCHMARK(BM_StlIO);

}  // namespace Kratos

BENCHMARK_MAIN();