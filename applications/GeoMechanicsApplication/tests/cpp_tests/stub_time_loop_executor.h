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

#pragma once

#include "custom_workflows/time_loop_executor_interface.h"

namespace Kratos {

class StubTimeLoopExecutor : public TimeLoopExecutorInterface
{
public:
    explicit StubTimeLoopExecutor(size_t NumberOfExpectedProcesses = 0);

    void SetCancelDelegate(const std::function<bool()>& rCancelDelegate) override;

    void SetProgressDelegate(const std::function<void(double)>& rProgressDelegate) override;
    
    void SetProcessObservables(const std::vector<std::weak_ptr<Process>>& rProcessObservables) override;

    void SetTimeIncrementor(std::unique_ptr<TimeIncrementor> pTimeIncrementor) override;

    void SetSolverStrategyWrapper(std::shared_ptr<StrategyWrapper> pStrategyWrapper) override;

    std::vector<TimeStepEndState> Run(const TimeStepEndState& EndState) override;

private:
    std::size_t mNumberOfExpectedProcesses;
};

}
