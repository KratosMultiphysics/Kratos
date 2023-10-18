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
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy.hpp"
#include "factories/standard_linear_solver_factory.h"
#include "scheme_factory.hpp"
#include "convergence_criteria_factory.hpp"
#include "builder_and_solver_factory.hpp"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class SolvingStrategyFactory
{
public:
    [[nodiscard]] static std::unique_ptr<SolvingStrategy<TSparseSpace, TDenseSpace>> Create(const Parameters& rSolverSettings, ModelPart& rModelPart)
    {
        const std::string strategy_type = "strategy_type";

        KRATOS_ERROR_IF_NOT(rSolverSettings.Has(strategy_type))
        << "The parameter strategy_type is undefined, aborting.";

        auto solver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(rSolverSettings["linear_solver_settings"]);
        auto scheme = SchemeFactory<TSparseSpace, TDenseSpace>::Create(rSolverSettings);
        auto builder_and_solver = BuilderAndSolverFactory<TSparseSpace, TDenseSpace, TLinearSolver>::Create(rSolverSettings, solver);
        auto criteria = ConvergenceCriteriaFactory<TSparseSpace, TDenseSpace>::Create(rSolverSettings);

        if (rSolverSettings[strategy_type].GetString() ==
            GeoMechanicsNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>::Name())
        {
            const int max_iterations = rSolverSettings["max_iterations"].GetInt();
            const bool calculate_reactions = rSolverSettings["calculate_reactions"].GetBool();
            const bool reform_dof_set_at_each_step = rSolverSettings["reform_dofs_at_each_step"].GetBool();
            const bool move_mesh_flag = rSolverSettings["move_mesh_flag"].GetBool();

            auto strategy_parameters = ExtractStrategyParameters(rSolverSettings);
            return std::make_unique<GeoMechanicsNewtonRaphsonStrategy<TSparseSpace,
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
        }

        return nullptr;
    }

private:
    static Parameters ExtractStrategyParameters(const Parameters &rSolverSettings)
    {
        Parameters result;
        std::vector<std::string> strategy_entries = {"min_iteration",
                                                     "number_cycles",
                                                     "increase_factor",
                                                     "reduction_factor",
                                                     "max_piping_iterations",
                                                     "desired_iterations",
                                                     "max_radius_factor",
                                                     "min_radius_factor",
                                                     "search_neighbours_step",
                                                     "body_domain_sub_model_part_list",
                                                     "loads_sub_model_part_list",
                                                     "loads_variable_list",
                                                     "rebuild_level"};
        for (const std::string& entry : strategy_entries)
        {
            if (rSolverSettings.Has(entry))
            {
                result.AddValue(entry, rSolverSettings[entry]);
            }
        }

        return result;
    }
};

}
