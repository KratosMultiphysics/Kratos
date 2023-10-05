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
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/convergencecriterias/mixed_generic_criteria.h"
#include "solving_strategy_factory.hpp"
#include "factories/standard_linear_solver_factory.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy.hpp"
#include "scheme_factory.hpp"
#include "convergence_criteria_factory.h"
#include "builder_and_solver_factory.hpp"

namespace Kratos
{

template<class TSparseSpace, class TDenseSpace, class TLinearSolver>
class SolvingStrategyFactory
{
public:
    [[nodiscard]] static std::unique_ptr<SolvingStrategy<TSparseSpace, TDenseSpace>> Create(Parameters& rSolverSettings, ModelPart& rModelPart)
    {
        if (rSolverSettings["strategy_type"].GetString() == ResidualBasedLinearStrategy<TSparseSpace,
                TDenseSpace,
                TLinearSolver>::Name())
        {
            return std::make_unique<ResidualBasedLinearStrategy<TSparseSpace,
                    TDenseSpace,
                    TLinearSolver>>();
        }
        else if (rSolverSettings["strategy_type"].GetString() == GeoMechanicsNewtonRaphsonStrategy<TSparseSpace,
                TDenseSpace,
                TLinearSolver>::Name())
        {
            auto solver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(rSolverSettings["linear_solver_settings"]);
            auto scheme = SchemeFactory<TSparseSpace, TDenseSpace>::Create(rSolverSettings);
            auto builder_and_solver = BuilderAndSolverFactory<TSparseSpace, TDenseSpace, TLinearSolver>::Create(rSolverSettings, solver);
            auto criteria = ConvergenceCriteriaFactory<TSparseSpace, TDenseSpace>::Create(rSolverSettings);

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

            Parameters parameters;
            for (const std::string& entry : strategy_entries)
            {
                if (rSolverSettings.Has(entry))
                {
                    parameters.AddValue(entry, rSolverSettings[entry]);
                }
            }

            return std::make_unique<GeoMechanicsNewtonRaphsonStrategy<TSparseSpace,
                    TDenseSpace,
                    TLinearSolver>>(rModelPart,
                    scheme,
                    solver,
                    criteria,
                    builder_and_solver,
                    parameters);
        }

        return nullptr;
    }

};

}
