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
const Parameters default_linear_solver_settings(
        R"(
    {
            "solver_type":     "amgcl",
            "smoother_type":   "ilu0",
            "krylov_type":     "gmres",
            "coarsening_type": "aggregation",
            "max_iteration":   100,
            "verbosity":       0,
            "tolerance":       1.0e-6,
            "scaling":         false
        }
    )");

unique_ptr<SolvingStrategy<SolvingStrategyFactory::SparseSpaceType, SolvingStrategyFactory::LocalSpaceType>>
SolvingStrategyFactory::Create(Parameters& rSolverSettings, ModelPart& rModelPart) const
{
    auto solver = LinearSolverFactory<SparseSpaceType, LocalSpaceType>().Create(default_linear_solver_settings);
    KRATOS_ERROR_IF(solver == nullptr);


    Scheme<SparseSpaceType, LocalSpaceType>::Pointer scheme = Kratos::make_shared<BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>>();
    auto builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, LinearSolverType>>(solver);

    const double rel_tol = 1.0e-4;
    const double abs_tol = 1.0e-9;
    VariableData *p_water_pres = &WATER_PRESSURE;
    ConvergenceVariableListType convergence_settings{std::make_tuple(p_water_pres, rel_tol, abs_tol)};
    auto criteria = std::make_shared<ConvergenceCriteriaType>(MixedGenericCriteriaType(convergence_settings));


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
        Parameters p_parameters(R"(
    {
        "min_iteration":    6,
        "number_cycles":    100,
        "increase_factor":  2.0,
        "reduction_factor": 0.5,
		"max_piping_iterations": 500,
        "desired_iterations": 4,
        "max_radius_factor": 10.0,
        "min_radius_factor": 0.1,
        "search_neighbours_step": false,
        "body_domain_sub_model_part_list": [],
        "loads_sub_model_part_list": [],
        "loads_variable_list" : []
    }  )");
        return std::make_unique<GeoMechanicsNewtonRaphsonStrategy<SparseSpaceType,
                LocalSpaceType,
                LinearSolverType>>(rModelPart,
                        scheme,
                        solver,
                        criteria,
                        builder_and_solver,
                                                           p_parameters);
    }




    return nullptr;
}

}
