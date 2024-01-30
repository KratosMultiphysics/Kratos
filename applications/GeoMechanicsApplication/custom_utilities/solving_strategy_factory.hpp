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

#pragma once

#include <memory>
#include <string>

#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/line_search_strategy.h"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy.hpp"
#include "factories/standard_linear_solver_factory.h"
#include "scheme_factory.hpp"
#include "convergence_criteria_factory.hpp"
#include "builder_and_solver_factory.hpp"
#include "parameters_utilities.h"

namespace Kratos {

template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class SolvingStrategyFactory {
public:
    [[nodiscard]] static std::unique_ptr<SolvingStrategy<TSparseSpace, TDenseSpace>>
    Create(const Parameters& rSolverSettings, ModelPart& rModelPart) {
        const std::string strategy_type = "strategy_type";

        KRATOS_ERROR_IF_NOT(rSolverSettings.Has(strategy_type))
                    << "The parameter strategy_type is undefined, aborting.";

        const auto echo_level = rSolverSettings["echo_level"].GetInt();

        auto solver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(
                rSolverSettings["linear_solver_settings"]);
        KRATOS_ERROR_IF_NOT(solver) << "Failed to create a linear solver" << std::endl;

        auto scheme = SchemeFactory<TSparseSpace, TDenseSpace>::Create(rSolverSettings);
        KRATOS_ERROR_IF_NOT(scheme) << "Failed to create a scheme" << std::endl;

        auto builder_and_solver = BuilderAndSolverFactory<TSparseSpace, TDenseSpace, TLinearSolver>::Create(
                rSolverSettings, solver);
        KRATOS_ERROR_IF_NOT(builder_and_solver) << "Failed to create a builder-and-solver" << std::endl;
        builder_and_solver->SetEchoLevel(echo_level);

        auto criteria = ConvergenceCriteriaFactory<TSparseSpace, TDenseSpace>::Create(rSolverSettings);
        KRATOS_ERROR_IF_NOT(criteria) << "Failed to create convergence criteria" << std::endl;
        criteria->SetEchoLevel(echo_level);

        if (rSolverSettings[strategy_type].GetString() == "newton_raphson") {
            const auto max_iterations = rSolverSettings["max_iterations"].GetInt();
            const auto calculate_reactions = rSolverSettings["calculate_reactions"].GetBool();
            const auto reform_dof_set_at_each_step = rSolverSettings["reform_dofs_at_each_step"].GetBool();
            const auto move_mesh_flag = rSolverSettings["move_mesh_flag"].GetBool();

            const std::vector<std::string> strategy_entries = {"loads_sub_model_part_list",
                                                               "loads_variable_list"};
            auto strategy_parameters = ParametersUtilities::CopyOptionalParameters(
                rSolverSettings, strategy_entries);
            auto result = std::make_unique<GeoMechanicsNewtonRaphsonStrategy<TSparseSpace,
                    TDenseSpace,
                    TLinearSolver>>(rModelPart,
                                    scheme,
                                    solver,
                                    criteria,
                                    builder_and_solver,
                                    strategy_parameters,
                                    max_iterations,
                                    calculate_reactions,
                                    reform_dof_set_at_each_step,
                                    move_mesh_flag);
            result->SetEchoLevel(echo_level);
            return result;
        } else if (rSolverSettings[strategy_type].GetString() == "line_search") {
            const std::vector<std::string> strategy_entries = {"max_iteration",
                                                               "compute_reactions",
                                                               "max_line_search_iterations",
                                                               "first_alpha_value",
                                                               "second_alpha_value",
                                                               "min_alpha",
                                                               "max_alpha",
                                                               "line_search_tolerance",
                                                               "move_mesh_flag",
                                                               "reform_dofs_at_each_step",
                                                               "echo_level"};

            auto strategy_parameters = ParametersUtilities::CopyOptionalParameters(
                rSolverSettings, strategy_entries);
            auto result = std::make_unique<LineSearchStrategy<TSparseSpace,
                    TDenseSpace,
                    TLinearSolver>>(rModelPart,
                                    scheme,
                                    solver,
                                    criteria,
                                    strategy_parameters);
            return result;
        }

        return nullptr;
    }
};

}