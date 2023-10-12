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

#include "custom_workflows/time_loop_executor.hpp"

namespace Kratos {

class StubTimeLoopExecutor : public TimeLoopExecutor
{
public:
    explicit StubTimeLoopExecutor(size_t NumberOfExpectedProcesses = 0);
    void SetProcessObservables(const std::vector<std::weak_ptr<Process>>& rProcessObservables) override;

private:
    std::size_t mNumberOfExpectedProcesses;
};

}
