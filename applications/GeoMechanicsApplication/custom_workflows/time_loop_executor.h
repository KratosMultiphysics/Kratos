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

#include <memory>
#include <vector>

#include "time_loop_executor_interface.h"
#include "time_step_end_state.h"
#include "time_step_executor.h"

namespace Kratos
{

class Process;
class StrategyWrapper;
class TimeIncrementor;

class KRATOS_API(GEO_MECHANICS_APPLICATION) TimeLoopExecutor : public TimeLoopExecutorInterface
{
public:
    TimeLoopExecutor();
    ~TimeLoopExecutor() override;

    void SetCancelDelegate(const std::function<bool()>& rCancelDelegate) override;
    void SetProgressDelegate(const std::function<void(double)>& rProgressDelegate) override;
    void SetProcessObservables(const std::vector<std::weak_ptr<Process>>& rProcessObservables) override;
    void SetTimeIncrementor(std::unique_ptr<TimeIncrementor> pTimeIncrementor) override;
    void SetSolverStrategyWrapper(std::shared_ptr<StrategyWrapper> pStrategyWrapper) override;
    std::vector<TimeStepEndState> Run(const TimeStepEndState& EndState) override;

private:
    void             CallExecuteBeforeSolutionLoopOnProcesses() const;
    TimeStepEndState RunCycle(double PreviousTime);
    TimeStepEndState RunCycleLoop(const TimeStepEndState& previous_state);
    bool             IsCancelled() const;
    void             UpdateProgress(double Time) const;

    std::unique_ptr<TimeIncrementor>    mTimeIncrementor;
    std::function<bool()>               mCancelDelegate;
    std::function<void(double)>         mProgressDelegate;
    std::unique_ptr<TimeStepExecutor>   mTimeStepExecutor = std::make_unique<TimeStepExecutor>();
    std::shared_ptr<StrategyWrapper>    mStrategyWrapper;
    std::vector<std::weak_ptr<Process>> mProcessObservables;
};

} // namespace Kratos