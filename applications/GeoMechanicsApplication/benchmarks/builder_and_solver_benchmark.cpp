// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "containers/model.h"
#include "geo_mechanics_application.h"
#include "geometries/quadrilateral_2d_4.h"
#include "linear_solvers/linear_solver.h"
#include "processes/structured_mesh_generator_process.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

#include <benchmark/benchmark.h>

namespace Kratos
{
using SparseSpaceType  = UblasSpace<double, CompressedMatrix, Vector>;
using DenseSpaceType   = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, DenseSpaceType>;

void BenchmarkBuildRHS(benchmark::State& rState)
{
    KratosApplication application("KratosApplication");
    application.Register();

    const auto number_of_cores = rState.range();

    ParallelUtilities::SetNumThreads(number_of_cores);

    // Create the test model part
    Model test_model;
    auto& r_test_model_part = test_model.CreateModelPart("TestModelPart");

    // Set up the test model part mesh
    auto                   p_point_1 = Kratos::make_intrusive<Node>(1, 0.0, 0.0, 0.0);
    auto                   p_point_2 = Kratos::make_intrusive<Node>(2, 0.0, 1.0, 0.0);
    auto                   p_point_3 = Kratos::make_intrusive<Node>(3, 1.0, 1.0, 0.0);
    auto                   p_point_4 = Kratos::make_intrusive<Node>(4, 1.0, 0.0, 0.0);
    Quadrilateral2D4<Node> domain_geometry(p_point_1, p_point_2, p_point_3, p_point_4);

    Parameters mesher_parameters(R"({
        "number_of_divisions": 10,
        "element_name": "Element2D3N",
        "condition_name": "LineCondition",
        "create_skin_sub_model_part": true
    })");
    StructuredMeshGeneratorProcess(domain_geometry, r_test_model_part, mesher_parameters).Execute();
    std::cout << "Number of elements = " << r_test_model_part.Elements().size() << std::endl;


    ResidualBasedBlockBuilderAndSolver<SparseSpaceType, DenseSpaceType, LinearSolverType> builder_and_solver{};
    Scheme<SparseSpaceType, DenseSpaceType>::Pointer p_scheme = std::make_shared<Scheme<SparseSpaceType, DenseSpaceType>>();


    ResidualBasedBlockBuilderAndSolver<SparseSpaceType, DenseSpaceType, LinearSolverType>::TSystemVectorPointerType b;
    ResidualBasedBlockBuilderAndSolver<SparseSpaceType, DenseSpaceType, LinearSolverType>::TSystemVectorPointerType Dx;
    ResidualBasedBlockBuilderAndSolver<SparseSpaceType, DenseSpaceType, LinearSolverType>::TSystemMatrixPointerType A;
    builder_and_solver.InitializeSolutionStep(r_test_model_part, *A, *Dx, *b);

    for (auto _ : rState) {
        builder_and_solver.BuildRHS(p_scheme, r_test_model_part, *b);
    }
}

BENCHMARK(BenchmarkBuildRHS)->Arg(1)->Arg(2)->Arg(4)->Arg(8)->Arg(16);

} // namespace Kratos

BENCHMARK_MAIN();