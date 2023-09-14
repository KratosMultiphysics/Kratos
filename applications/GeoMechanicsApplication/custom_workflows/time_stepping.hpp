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

#include <functional>
#include <memory>


namespace Kratos
{

template <typename StrategyType>
class TimeStepExecuter
{
public:
    enum class ConvergenceState {converged, non_converged};

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

    ConvergenceState Run()
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

        return mStrategy->IsConverged() ? ConvergenceState::converged : ConvergenceState::non_converged;
    }

private:
    std::shared_ptr<StrategyType> mStrategy;
    ProcessRefVec                 mProcessRefs;
};

}
