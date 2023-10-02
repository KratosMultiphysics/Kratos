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

#pragma once

#include "processes/process.h"
#include "time_step_end_state.hpp"
#include "strategy_wrapper.hpp"
#include "geo_mechanics_application_variables.h"

#include <functional>
#include <memory>


namespace Kratos
{

class TimeStepExecutor
{
public:
    using ProcessRef    = std::reference_wrapper<Process>;
    using ProcessRefVec = std::vector<ProcessRef>;

    void SetSolverStrategy(std::shared_ptr<StrategyWrapper> SolverStrategy)
    {
        mStrategyWrapper = std::move(SolverStrategy);
    }

    void SetProcessReferences(ProcessRefVec ProcessRefs)
    {
        mProcessRefs = std::move(ProcessRefs);
    }

    TimeStepEndState Run(double Time)
    {
        mStrategyWrapper->Initialize();
        mStrategyWrapper->InitializeSolutionStep();

        for (auto& process : mProcessRefs)
        {
            process.get().ExecuteInitializeSolutionStep();
        }

        mStrategyWrapper->Predict();
        mStrategyWrapper->SolveSolutionStep();

        for (auto& process : mProcessRefs)
        {
            process.get().ExecuteFinalizeSolutionStep();
        }

        mStrategyWrapper->FinalizeSolutionStep();

        TimeStepEndState result;
        result.time = Time;
        result.convergence_state = mStrategyWrapper->IsConverged() ? TimeStepEndState::ConvergenceState::converged :
                                   TimeStepEndState::ConvergenceState::non_converged;
        //result.number_of_nonlinear_iterations = mStrategyWrapper->GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER];
        result.number_of_nonlinear_iterations = mStrategyWrapper->GetNumberOfIterations();
        return result;
    }

private:
    std::shared_ptr<StrategyWrapper> mStrategyWrapper;
    ProcessRefVec                 mProcessRefs;
};

}
