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

void TimeStepExecutor::SetProcessReferences(ProcessRefVec ProcessRefs)
{
    mProcessRefs = std::move(ProcessRefs);
}

TimeStepEndState TimeStepExecutor::Run(double Time)
{
    mStrategyWrapper->Initialize();
    mStrategyWrapper->InitializeSolutionStep();

    for (const auto& process : mProcessRefs)
    {
        process.get().ExecuteInitializeSolutionStep();
    }

    mStrategyWrapper->Predict();
    mStrategyWrapper->SolveSolutionStep();

    for (const auto& process : mProcessRefs)
    {
         process.get().ExecuteFinalizeSolutionStep();
    }

    mStrategyWrapper->FinalizeSolutionStep();

    TimeStepEndState result;
    result.time              = Time;
    result.convergence_state = mStrategyWrapper->GetConvergenceState();
    result.num_of_iterations = mStrategyWrapper->GetNumberOfIterations();
    return result;
}

}