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
//                   Anne van de Graaf
//                   Wijtze Pieter Kikstra
//

#pragma once

#include <vector>
#include <memory>

#include "time_step_end_state.hpp"
#include "time_incrementor.h."
#include "processes/process.h"
#include "strategy_wrapper.hpp"

namespace Kratos
{

class Process;

class TimeLoopExecutor{
public :
    virtual ~TimeLoopExecutor() = default;
    virtual void SetProcessReferences(const std::vector<std::weak_ptr<Process>>& rProcessRefs) {}

    void SetTimeIncrementor(std::unique_ptr<TimeIncrementor> ATimeIncrementor)
    {
        mTimeIncrementor = std::move(ATimeIncrementor);
    }

    void SetSolverStrategyTimeIncrementor(std::unique_ptr<StrategyWrapper> AStrategyWrapper)
    {
        mSolverStrategy = std::move(AStrategyWrapper);
    }

    std::vector<TimeStepEndState> Run(TimeStepEndState end_state)
    {
        std::vector<TimeStepEndState> result;
        while (mTimeIncrementor->WantNextStep(end_state)) {
            end_state.convergence_state = mSolverStrategy->GetConvergenceState();
            end_state.time              += mTimeIncrementor->GetIncrement();
            result.emplace_back(end_state);
            mTimeIncrementor->PostTimeStepExecution(end_state);
        }
        return result;
    }

private:
    std::unique_ptr<TimeIncrementor> mTimeIncrementor;
    std::unique_ptr<StrategyWrapper> mSolverStrategy;
};

}