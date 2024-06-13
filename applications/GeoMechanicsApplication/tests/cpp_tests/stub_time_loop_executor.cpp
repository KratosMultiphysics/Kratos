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
//
#include "stub_time_loop_executor.h"
#include "processes/process.h"
#include "testing/testing.h"

namespace Kratos {

StubTimeLoopExecutor::StubTimeLoopExecutor(std::size_t NumberOfExpectedProcesses) :
        mNumberOfExpectedProcesses{NumberOfExpectedProcesses}
{
}

void StubTimeLoopExecutor::SetProcessObservables(const std::vector<std::weak_ptr<Process>>& rProcessObservables)
{
    KRATOS_EXPECT_EQ(rProcessObservables.size(), mNumberOfExpectedProcesses);
    for (const auto& process_observable : rProcessObservables)
    {
        KRATOS_EXPECT_FALSE(process_observable.expired())

        std::shared_ptr<Process> process = process_observable.lock();
        KRATOS_EXPECT_EQ(process->Check(), 0);
    }
}

void StubTimeLoopExecutor::SetCancelDelegate(const std::function<bool()>& rCancelDelegate) {
    // intentionally empty
}

void StubTimeLoopExecutor::SetProgressDelegate(const std::function<void(double)>& rProgressDelegate) {
    // intentionally empty
}

void StubTimeLoopExecutor::SetTimeIncrementor(std::unique_ptr<TimeIncrementor> pTimeIncrementor) {
    // intentionally empty
}

void StubTimeLoopExecutor::SetSolverStrategyWrapper(std::shared_ptr<StrategyWrapper> pStrategyWrapper) {
    // intentionally empty
}

std::vector<TimeStepEndState> StubTimeLoopExecutor::Run(const TimeStepEndState& EndState) {
    return {};
}

}
