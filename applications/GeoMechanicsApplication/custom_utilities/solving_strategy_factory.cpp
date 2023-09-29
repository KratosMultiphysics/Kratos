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
SolvingStrategyFactory::Create(const Parameters& rSolverSettings) const
{
    auto solver = LinearSolverFactory<SolvingStrategyFactory::SparseSpaceType, SolvingStrategyFactory::LocalSpaceType>().Create(default_linear_solver_settings);
    KRATOS_ERROR_IF(solver == nullptr);

    if (rSolverSettings["strategy_type"].GetString() == ResidualBasedLinearStrategy<SolvingStrategyFactory::SparseSpaceType,
            SolvingStrategyFactory::LocalSpaceType,
            SolvingStrategyFactory::LinearSolverType>::Name())
    {
        return std::make_unique<ResidualBasedLinearStrategy<SolvingStrategyFactory::SparseSpaceType,
                SolvingStrategyFactory::LocalSpaceType,
                SolvingStrategyFactory::LinearSolverType>>();
    }

    return nullptr;
}

//SolvingStrategyFactory::LinearSolverType::Pointer SolvingStrategyFactory::setup_solver_dgeoflow()
//{
//    return Kratos::make_shared<SkylineLUFactorizationSolverType>();
//}
//
//SolvingStrategyFactory::ResidualBasedLinearStrategy::Pointer SolvingStrategyFactory::setup_strategy_dgeoflow(ModelPart &model_part)
//{
//    // Create the linear strategy
//    auto p_solver = setup_solver_dgeoflow();
//
//    Scheme<SparseSpaceType, LocalSpaceType>::Pointer p_scheme = Kratos::make_shared<BackwardEulerQuasistaticPwScheme<SparseSpaceType, LocalSpaceType>>();
//
//    auto p_builder_and_solver = Kratos::make_shared<ResidualBasedBlockBuilderAndSolver<SparseSpaceType, LocalSpaceType, KratosExecute::LinearSolverType>>(p_solver);
//    p_builder_and_solver->SetEchoLevel(0);
//
//    auto p_criteria = setup_criteria_dgeoflow();
//    p_criteria->SetEchoLevel(0);
//
//    Parameters p_parameters(R"(
//    {
//        "min_iteration":    6,
//        "number_cycles":    100,
//        "increase_factor":  2.0,
//        "reduction_factor": 0.5,
//		"max_piping_iterations": 500,
//        "desired_iterations": 4,
//        "max_radius_factor": 10.0,
//        "min_radius_factor": 0.1,
//        "search_neighbours_step": false,
//        "body_domain_sub_model_part_list": [],
//        "loads_sub_model_part_list": [],
//        "loads_variable_list" : []
//    }  )");
//
//    int MaxIterations = 15;
//    bool CalculateReactions = true;
//    bool ReformDofSetAtEachStep = false;
//    bool MoveMeshFlag = false;
//
//    auto p_solving_strategy = Kratos::make_unique<GeoMechanicsNewtonRaphsonErosionProcessStrategy<SparseSpaceType, LocalSpaceType, KratosExecute::LinearSolverType>>(
//            model_part,
//                    p_scheme,
//                    p_solver,
//                    p_criteria,
//                    p_builder_and_solver,
//                    p_parameters,
//                    MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag);
//
//    p_solving_strategy->Check();
//    return p_solving_strategy;
//}



}
