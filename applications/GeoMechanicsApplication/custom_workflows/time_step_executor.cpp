// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#include "time_step_executor.h"

namespace Kratos
{

void TimeStepExecutor::SetSolverStrategy(std::shared_ptr<StrategyWrapper> SolverStrategy)
{
    mStrategyWrapper = std::move(SolverStrategy);
}

void TimeStepExecutor::SetProcessObservables(const std::vector<std::weak_ptr<Process>>& rProcessObservables)
{
    mProcessObservables = rProcessObservables;
}

TimeStepEndState TimeStepExecutor::Run(double Time)
{
    KRATOS_INFO("TimeStepExecutor") << "Running time step at time " << Time << std::endl;

    mStrategyWrapper->InitializeSolutionStep();

    for (const auto& process_observable : mProcessObservables) {
        auto process = process_observable.lock();
        if (process) process->ExecuteInitializeSolutionStep();
    }

    mStrategyWrapper->Predict();

    // Here we directly use the SolveSolutionStep() method of the strategy wrapper
    // to get the convergence state, instead of the IsConverged method of the strategy.
    // This is because the IsConverged method of the strategy makes unforeseen changes
    // to the results.
    TimeStepEndState result;
    result.time              = Time;
    result.convergence_state = mStrategyWrapper->SolveSolutionStep();
    result.num_of_iterations = mStrategyWrapper->GetNumberOfIterations();

    for (const auto& process_observable : mProcessObservables) {
        auto process = process_observable.lock();
        if (process) process->ExecuteFinalizeSolutionStep();
    }

    return result;
}

} // namespace Kratos