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

#include "custom_workflows/time_loop_executor.h"

namespace Kratos {

class StubTimeLoopExecutor : public TimeLoopExecutor
{
public:
    explicit StubTimeLoopExecutor(int NumberOfExpectedProcesses = 0);
    void SetProcessReferences(const std::vector<std::reference_wrapper<Process>>& rProcessRefs) override;

private:
    int mNumberOfExpectedProcesses;
};

}