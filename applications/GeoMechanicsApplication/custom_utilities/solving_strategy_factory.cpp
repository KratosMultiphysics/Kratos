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
#include "custom_strategies/schemes/backward_euler_quasistatic_Pw_scheme.hpp"

namespace Kratos
{

unique_ptr<SolvingStrategy<SolvingStrategyFactory::SparseSpaceType, SolvingStrategyFactory::LocalSpaceType>>
SolvingStrategyFactory::Create(Parameters& rSolverSettings, ModelPart& rModelPart) const
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
        KRATOS_ERROR_IF(solver == nullptr);

        Scheme<SparseSpaceType, LocalSpaceType>::Pointer scheme = Kratos::make_shared<BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>>();
        auto builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>(solver);

        const double rel_tol = 1.0e-4;
        const double abs_tol = 1.0e-9;
        VariableData *p_water_pres = &WATER_PRESSURE;
        ConvergenceVariableListType convergence_settings{std::make_tuple(p_water_pres, rel_tol, abs_tol)};
        auto criteria = std::make_shared<ConvergenceCriteriaType>(MixedGenericCriteriaType(convergence_settings));
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
