//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
//

// System includes
#include <filesystem>
#include <fstream>
#include <string>

// External includes
#include <benchmark/benchmark.h>

// Project includes
#include "containers/model.h"
#include "geometries/quadrilateral_2d_4.h"
#include "includes/kratos_application.h"
#include "includes/model_part_io.h"
#include "processes/structured_mesh_generator_process.h"

namespace Kratos {

namespace {

// Register KratosCore components exactly once per process.
void RegisterKratos()
{
    static KratosApplication application("KratosApplication");
    static bool registered = false;
    if (!registered) {
        application.Register();
        registered = true;
    }
}

// Build a Cartesian 2D mesh of triangles on the unit square [0,1]².
// NumberOfDivisions n produces:
//   - 2·n²     Element2D3N  elements
//   - (n+1)²   nodes
void GenerateCartesianMesh(ModelPart& rModelPart, int NumberOfDivisions)
{
    auto p1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto p2 = Kratos::make_intrusive<Node>(2, 1.0, 0.0, 0.0);
    auto p3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    auto p4 = Kratos::make_intrusive<Node>(4, 0.0, 1.0, 0.0);
    Quadrilateral2D4<Node> domain(p1, p2, p3, p4);

    const std::string params_str =
        "{ \"number_of_divisions\": " + std::to_string(NumberOfDivisions) + ","
        "  \"element_name\": \"Element2D3N\","
        "  \"condition_name\": \"LineCondition2D2N\","
        "  \"create_skin_sub_model_part\": false }";

    StructuredMeshGeneratorProcess(domain, rModelPart, Parameters(params_str)).Execute();
}

std::string TempFileStem()
{
    return (std::filesystem::temp_directory_path() / "kratos_io_benchmark").string();
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// BM_ModelPartIOWrite
// Setup  (not timed): generate a Cartesian mesh in memory.
// Measured: serialise the full model part to an .mdpa file on disk.
// ---------------------------------------------------------------------------
static void BM_ModelPartIOWrite(benchmark::State& rState)
{
    RegisterKratos();

    const int n_div = static_cast<int>(rState.range(0));
    const std::string stem = TempFileStem();

    Model model;
    auto& r_model_part = model.CreateModelPart("Mesh");
    GenerateCartesianMesh(r_model_part, n_div);

    const std::size_t n_nodes    = r_model_part.NumberOfNodes();
    const std::size_t n_elements = r_model_part.NumberOfElements();

    rState.SetLabel("nodes=" + std::to_string(n_nodes) +
                    " elems=" + std::to_string(n_elements));

    for (auto _ : rState) {
        ModelPartIO io(stem, IO::WRITE);
        io.WriteModelPart(r_model_part);
    }

    rState.counters["nodes"]    = benchmark::Counter(static_cast<double>(n_nodes));
    rState.counters["elements"] = benchmark::Counter(static_cast<double>(n_elements));
    rState.counters["nodes/s"]  = benchmark::Counter(
        static_cast<double>(n_nodes), benchmark::Counter::kIsRate);
    rState.counters["elems/s"]  = benchmark::Counter(
        static_cast<double>(n_elements), benchmark::Counter::kIsRate);

    std::filesystem::remove(stem + ".mdpa");
}

// ---------------------------------------------------------------------------
// BM_ModelPartIORead
// Setup  (not timed): generate a Cartesian mesh and write it to disk.
// Measured: deserialise the model part from the .mdpa file.
// ---------------------------------------------------------------------------
static void BM_ModelPartIORead(benchmark::State& rState)
{
    RegisterKratos();

    const int n_div = static_cast<int>(rState.range(0));
    const std::string stem = TempFileStem();

    std::size_t n_nodes    = 0;
    std::size_t n_elements = 0;
    {
        Model setup_model;
        auto& r_mp = setup_model.CreateModelPart("Setup");
        GenerateCartesianMesh(r_mp, n_div);
        n_nodes    = r_mp.NumberOfNodes();
        n_elements = r_mp.NumberOfElements();

        ModelPartIO io_write(stem, IO::WRITE);
        io_write.WriteModelPart(r_mp);
    }

    rState.SetLabel("nodes=" + std::to_string(n_nodes) +
                    " elems=" + std::to_string(n_elements));

    for (auto _ : rState) {
        rState.PauseTiming();
        Model model;
        auto& r_model_part = model.CreateModelPart("Mesh");
        rState.ResumeTiming();

        ModelPartIO io(stem);
        io.ReadModelPart(r_model_part);

        benchmark::DoNotOptimize(r_model_part.NumberOfNodes());
    }

    rState.counters["nodes"]    = benchmark::Counter(static_cast<double>(n_nodes));
    rState.counters["elements"] = benchmark::Counter(static_cast<double>(n_elements));
    rState.counters["nodes/s"]  = benchmark::Counter(
        static_cast<double>(n_nodes), benchmark::Counter::kIsRate);
    rState.counters["elems/s"]  = benchmark::Counter(
        static_cast<double>(n_elements), benchmark::Counter::kIsRate);

    std::filesystem::remove(stem + ".mdpa");
}

// Mesh sizes by number of Cartesian divisions (2D unit square → triangles):
//   n_div=100  →    20 000 elements /   10 201 nodes
//   n_div=300  →   180 000 elements /   90 601 nodes
//   n_div=700  →   980 000 elements /  491 401 nodes
//   n_div=1000 → 2 000 000 elements / 1 002 001 nodes
BENCHMARK(BM_ModelPartIOWrite)->Arg(100)->Arg(300)->Arg(700)->Arg(1000);
BENCHMARK(BM_ModelPartIORead) ->Arg(100)->Arg(300)->Arg(700)->Arg(1000);

} // namespace Kratos

BENCHMARK_MAIN();
