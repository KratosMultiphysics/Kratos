// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include <string>

#include "containers/model.h"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy_with_seepage.hpp"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "spaces/ublas_space.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite.h"

using namespace Kratos;
using namespace std::string_literals;

namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(GeoMechanicsNewtonRaphsonStrategyWithSeepage_ThrowsWhenNumberOfCyclesIsGreaterThanOne,
                          KratosGeoMechanicsFastSuiteWithoutKernel)
{
    // At present, the fixity of degrees of freedom is not properly restored when a solution step
    // does not converge and another cycle is attempted. This test ensures that we don't
    // accidentally run a seepage analysis with multiple cycles, which would lead to incorrect
    // results. See also method `_RevertStateToStartOfStep` in `geomechanics_analysis.py`, which
    // omits restoring the fixity of degrees of freedom.

    // The settings below are supposed to be contained by the "solver_settings" object
    const auto solver_settings = Parameters{R"(
    {
        "number_cycles": 2
    })"s};

    auto  model        = Model{};
    auto& r_model_part = model.CreateModelPart("Main"s);

    using SparseSpaceType  = UblasSpace<double, CompressedMatrix, Vector>;
    using LocalSpaceType   = UblasSpace<double, Matrix, Vector>;
    using LinearSolverType = LinearSolver<SparseSpaceType, LocalSpaceType>;
    auto p_convergence_criterion = std::make_shared<ConvergenceCriteria<SparseSpaceType, LocalSpaceType>>();
    auto p_builder_and_solver =
        std::make_shared<BuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>();
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(
        (GeoMechanicsNewtonRaphsonStrategyWithSeepage<SparseSpaceType, LocalSpaceType, LinearSolverType>{r_model_part, nullptr, p_convergence_criterion, p_builder_and_solver, solver_settings}), "GeoMechanicsNewtonRaphsonStrategyWithSeepage does not support multiple cycles. This is because the fixity of degrees of freedom is not properly restored when a solution step does not converge and another cycle is attempted.")
}

} // namespace Kratos::Testing