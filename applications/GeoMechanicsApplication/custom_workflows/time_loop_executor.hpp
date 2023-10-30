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

#include "time_loop_executor_interface.h"
#include "time_step_end_state.hpp"
#include "time_incrementor.h"
#include "time_step_executor.h"
#include "strategy_wrapper.hpp"

namespace Kratos
{

class Process;

class TimeLoopExecutor : public TimeLoopExecutorInterface {
public :
    void SetCancelDelegate(const std::function<bool()>& rCancelDelegate) override
    {
        mCancelDelegateAvailable = true;
        mCancelDelegate = rCancelDelegate;
    }

    void SetProcessObservables(const std::vector<std::weak_ptr<Process>>& rProcessObservables) override
    {
        mTimeStepExecutor->SetProcessObservables(rProcessObservables);
    }

    void SetTimeIncrementor(std::unique_ptr<TimeIncrementor> pTimeIncrementor) override
    {
        mTimeIncrementor = std::move(pTimeIncrementor);
    }

    void SetSolverStrategyWrapper(std::shared_ptr<StrategyWrapper> pStrategyWrapper) override
    {
        mStrategyWrapper = std::move(pStrategyWrapper);
        mTimeStepExecutor->SetSolverStrategy(mStrategyWrapper);
    }

    std::vector<TimeStepEndState> Run(const TimeStepEndState& EndState) override
    {
        mStrategyWrapper->Initialize();
        mStrategyWrapper->SaveTotalDisplacementFieldAtStartOfTimeLoop();
        std::vector<TimeStepEndState> result;
        TimeStepEndState NewEndState = EndState;
        while (mTimeIncrementor->WantNextStep(NewEndState) && mCancelDelegateAvailable && !mCancelDelegate()) {
            mStrategyWrapper->IncrementStepNumber();
            // clone without end time, the end time is overwritten anyway
            mStrategyWrapper->CloneTimeStep();
            NewEndState = RunCycleLoop(NewEndState);
            mStrategyWrapper->AccumulateTotalDisplacementField();
            mStrategyWrapper->FinalizeSolutionStep();
            mStrategyWrapper->OutputProcess();
            result.emplace_back(NewEndState);
        }

        mStrategyWrapper->FinalizeOutput();

        return result;
    }

private:
    TimeStepEndState RunCycle(double PreviousTime)
    {
        // Setting the time and time increment may be needed for the processes
        const auto time_increment = mTimeIncrementor->GetIncrement();
        mStrategyWrapper->SetTimeIncrement(time_increment);
        const auto end_time = PreviousTime + time_increment;
        mStrategyWrapper->SetEndTime(end_time);

        auto end_state = mTimeStepExecutor->Run(end_time);
        mTimeIncrementor->PostTimeStepExecution(end_state);
        return end_state;
    }

    TimeStepEndState RunCycleLoop(const TimeStepEndState& previous_state)
    {
        auto cycle_number = 0;
        auto end_state = previous_state;
        while (mTimeIncrementor->WantRetryStep(cycle_number, end_state) && !mCancelDelegate()) {
            if (cycle_number > 0) mStrategyWrapper->RestorePositionsAndDOFVectorToStartOfStep();
            end_state = RunCycle(previous_state.time);
            ++cycle_number;
        }

        end_state.num_of_cycles = cycle_number;
        return end_state;
    }

    std::unique_ptr<TimeIncrementor>  mTimeIncrementor;
    bool                              mCancelDelegateAvailable = false;
    std::function<bool()>             mCancelDelegate;
    std::unique_ptr<TimeStepExecutor> mTimeStepExecutor = std::make_unique<TimeStepExecutor>();
    std::shared_ptr<StrategyWrapper>  mStrategyWrapper;
};

}