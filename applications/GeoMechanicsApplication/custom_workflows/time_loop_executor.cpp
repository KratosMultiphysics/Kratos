// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//                   Wijtze Pieter Kikstra
//                   Richard Faasse
//

#include "time_loop_executor.h"
#include "processes/process.h"
#include "scoped_output_file_access.h"
#include "strategy_wrapper.hpp"
#include "time_incrementor.h"

namespace Kratos
{

void TimeLoopExecutor::SetCancelDelegate(const std::function<bool()>& rCancelDelegate)
{
    mCancelDelegate = rCancelDelegate;
}

void TimeLoopExecutor::SetProgressDelegate(const std::function<void(double)>& rProgressDelegate)
{
    mProgressDelegate = rProgressDelegate;
}

void TimeLoopExecutor::SetProcessObservables(const std::vector<std::weak_ptr<Process>>& rProcessObservables)
{
    mProcessObservables = rProcessObservables;
    mTimeStepExecutor->SetProcessObservables(rProcessObservables);
}

void TimeLoopExecutor::SetTimeIncrementor(std::unique_ptr<TimeIncrementor> pTimeIncrementor)
{
    mTimeIncrementor = std::move(pTimeIncrementor);
}

void TimeLoopExecutor::SetSolverStrategyWrapper(std::shared_ptr<StrategyWrapper> pStrategyWrapper)
{
    mStrategyWrapper = std::move(pStrategyWrapper);
    mTimeStepExecutor->SetSolverStrategy(mStrategyWrapper);
}

std::vector<TimeStepEndState> TimeLoopExecutor::Run(const TimeStepEndState& EndState)
{
    mStrategyWrapper->Initialize();
    CallExecuteBeforeSolutionLoopOnProcesses();

    std::vector<TimeStepEndState> result;
    TimeStepEndState              NewEndState = EndState;
    ScopedOutputFileAccess        limit_output_file_access_to_this_scope{*mStrategyWrapper};
    while (mTimeIncrementor->WantNextStep(NewEndState) && !IsCancelled()) {
        mStrategyWrapper->IncrementStepNumber();
        // clone without end time, the end time is overwritten anyway
        mStrategyWrapper->CloneTimeStep();
        NewEndState = RunCycleLoop(NewEndState);
        KRATOS_ERROR_IF_NOT(NewEndState.Converged())
            << "The calculation exited without converging.\n";
        mStrategyWrapper->ComputeIncrementalDisplacementField();
        mStrategyWrapper->AccumulateTotalDisplacementField();
        mStrategyWrapper->FinalizeSolutionStep();
        mStrategyWrapper->OutputProcess();
        result.emplace_back(NewEndState);
    }

    return result;
}

void TimeLoopExecutor::CallExecuteBeforeSolutionLoopOnProcesses() const
{
    for (const auto& r_process_observable : mProcessObservables) {
        auto p_process = r_process_observable.lock();
        if (p_process) p_process->ExecuteBeforeSolutionLoop();
    }
}

TimeStepEndState TimeLoopExecutor::RunCycle(double PreviousTime)
{
    // Setting the time and time increment may be needed for the processes
    const auto time_increment = mTimeIncrementor->GetIncrement();
    mStrategyWrapper->SetTimeIncrement(time_increment);
    const auto end_time = PreviousTime + time_increment;
    mStrategyWrapper->SetEndTime(end_time);

    auto end_state = mTimeStepExecutor->Run(end_time);
    mTimeIncrementor->PostTimeStepExecution(end_state);
    UpdateProgress(end_time);
    return end_state;
}

TimeStepEndState TimeLoopExecutor::RunCycleLoop(const TimeStepEndState& previous_state)
{
    auto cycle_number = 0;
    auto end_state    = previous_state;
    while (mTimeIncrementor->WantRetryStep(cycle_number, end_state) && !IsCancelled()) {
        if (cycle_number > 0) mStrategyWrapper->RestorePositionsAndDOFVectorToStartOfStep();
        end_state = RunCycle(previous_state.time);
        ++cycle_number;
    }

    end_state.num_of_cycles = cycle_number;
    return end_state;
}

bool TimeLoopExecutor::IsCancelled() const { return mCancelDelegate && mCancelDelegate(); }

void TimeLoopExecutor::UpdateProgress(double Time) const
{
    if (mProgressDelegate) {
        mProgressDelegate(Time);
    }
}

// This default destructor is added in the cpp to be able to forward declaremember variables in
// a unique_ptr
TimeLoopExecutor::~TimeLoopExecutor() = default;

// This default constructor is added in the cpp to make sure the destructor of the forward declared
// member variables is known when moving a constructed object
TimeLoopExecutor::TimeLoopExecutor() = default;

} // namespace Kratos