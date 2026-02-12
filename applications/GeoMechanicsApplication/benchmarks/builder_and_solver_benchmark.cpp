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

#include <benchmark/benchmark.h>


namespace Kratos
{


void BenchmarkBuildRHS(benchmark::State& rState)
{
    for (auto _ : rState) {
        double test = 0.0;
        test *= 2;
    }
}

BENCHMARK(BenchmarkBuildRHS);

} // namespace Kratos

BENCHMARK_MAIN();