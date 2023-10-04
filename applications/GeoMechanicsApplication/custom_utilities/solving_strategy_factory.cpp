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
#include "solving_strategy_factory.h"
#include "factories/standard_linear_solver_factory.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "custom_strategies/strategies/geo_mechanics_newton_raphson_strategy.hpp"
#include "scheme_factory.h"
#include "convergence_criteria_factory.h"

namespace Kratos
{

unique_ptr<SolvingStrategy<SolvingStrategyFactory::SparseSpaceType, SolvingStrategyFactory::LocalSpaceType>>
SolvingStrategyFactory::Create(Parameters& rSolverSettings, ModelPart& rModelPart)
{
    if (rSolverSettings["strategy_type"].GetString() == ResidualBasedLinearStrategy<SparseSpaceType,
            LocalSpaceType,
            LinearSolverType>::Name())
    {
        return std::make_unique<ResidualBasedLinearStrategy<SparseSpaceType,
                LocalSpaceType,
                LinearSolverType>>();
    }
    else if (rSolverSettings["strategy_type"].GetString() == GeoMechanicsNewtonRaphsonStrategy<SparseSpaceType,
            LocalSpaceType,
            LinearSolverType>::Name())
    {
        auto solver = LinearSolverFactory<SparseSpaceType, LocalSpaceType>().Create(rSolverSettings["linear_solver_settings"]);
        auto scheme = SchemeFactory<SparseSpaceType, LocalSpaceType>::Create(rSolverSettings);
        auto builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>(solver);
        auto criteria = ConvergenceCriteriaFactory<SparseSpaceType, LocalSpaceType>::Create(rSolverSettings);

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

        return std::make_unique<GeoMechanicsNewtonRaphsonStrategy<SparseSpaceType,
                LocalSpaceType,
                LinearSolverType>>(rModelPart,
                        scheme,
                        solver,
                        criteria,
                        builder_and_solver,
                                   parameters);
    }

    return nullptr;
}

}
