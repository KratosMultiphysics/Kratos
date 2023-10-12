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

class TimeLoopExecutor : public TimeLoopExecutorInterface{
public :
    TimeLoopExecutor() : mTimeStepExecutor{std::make_unique<TimeStepExecutor>()} {}

    virtual void SetProcessObservables(const std::vector<std::weak_ptr<Process>>& rProcessObservables) override
    {
        mTimeStepExecutor->SetProcessObservables(rProcessObservables);
    }

    void SetTimeIncrementor(std::unique_ptr<TimeIncrementor> pTimeIncrementor) override
    {
        mTimeIncrementor = std::move(pTimeIncrementor);
    }

    void SetSolverStrategyTimeStepExecutor(std::shared_ptr<StrategyWrapper> pStrategyWrapper) override
    {
        mStrategyWrapper = std::move(pStrategyWrapper);
        mTimeStepExecutor->SetSolverStrategy(mStrategyWrapper);
    }

    std::vector<TimeStepEndState> Run(TimeStepEndState EndState) override
    {
        mStrategyWrapper->SaveTotalDisplacementFieldAtStartOfStage();
        std::vector<TimeStepEndState> result;
        while (mTimeIncrementor->WantNextStep(EndState)) {
            mStrategyWrapper->IncrementStepNumber();
            // clone without end time, the end time is overwritten anyway
            mStrategyWrapper->CloneTimeStep();
            EndState = RunCycleLoop(EndState);
            mStrategyWrapper->AccumulateTotalDisplacementField();
            mStrategyWrapper->FinalizeSolutionStep();
            result.emplace_back(EndState);
            mStrategyWrapper->OutputProcess();
        }
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
        while (mTimeIncrementor->WantRetryStep(cycle_number, end_state)) {
            if (cycle_number > 0) mStrategyWrapper->RestorePositionsAndDOFVectorToStartOfStep();
            end_state = RunCycle(previous_state.time);
            ++cycle_number;
        }

        end_state.num_of_cycles = cycle_number;
        return end_state;
    }

    std::unique_ptr<TimeIncrementor>  mTimeIncrementor;
    std::unique_ptr<TimeStepExecutor> mTimeStepExecutor;
    std::shared_ptr<StrategyWrapper>  mStrategyWrapper;
};

}