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

#include <functional>
#include <memory>


namespace Kratos
{

template <typename StrategyType>
class TimeStepExecutor
{
public:
    using ProcessRef    = std::reference_wrapper<Process>;
    using ProcessRefVec = std::vector<ProcessRef>;

    void SetSolverStrategy(std::shared_ptr<StrategyType> SolverStrategy)
    {
        mStrategy = std::move(SolverStrategy);
    }

    void SetProcessReferences(ProcessRefVec ProcessRefs)
    {
        mProcessRefs = std::move(ProcessRefs);
    }

    TimeStepEndState Run(double Time)
    {
        mStrategy->Initialize();
        mStrategy->InitializeSolutionStep();

        for (auto& process : mProcessRefs)
        {
            process.get().ExecuteInitializeSolutionStep();
        }

        mStrategy->Predict();
        mStrategy->SolveSolutionStep();

        for (auto& process : mProcessRefs)
        {
            process.get().ExecuteFinalizeSolutionStep();
        }

        mStrategy->FinalizeSolutionStep();

        TimeStepEndState result;
        result.time = Time;
        result.convergence_state = mStrategy->IsConverged() ? TimeStepEndState::ConvergenceState::converged :
                                                              TimeStepEndState::ConvergenceState::non_converged;
        return result;
    }

private:
    std::shared_ptr<StrategyType> mStrategy;
    ProcessRefVec                 mProcessRefs;
};

}
