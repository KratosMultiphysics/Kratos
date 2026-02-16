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

#include "linear_solvers/linear_solver.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

#include <benchmark/benchmark.h>


namespace Kratos
{
using SparseSpaceType   = UblasSpace<double, CompressedMatrix, Vector>;
using DenseSpaceType    = UblasSpace<double, Matrix, Vector>;
using LinearSolverType = LinearSolver<SparseSpaceType, DenseSpaceType>;

void BenchmarkBuildRHS(benchmark::State& rState)
{
    ResidualBasedBlockBuilderAndSolver<SparseSpaceType, DenseSpaceType,LinearSolverType> builder_and_solver();


    for (auto _ : rState) {
        double test = 0.0;
        test *= 2;
    }
}

BENCHMARK(BenchmarkBuildRHS);

} // namespace Kratos

BENCHMARK_MAIN();