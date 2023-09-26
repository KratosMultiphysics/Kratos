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

void StubTimeLoopExecutor::SetProcessReferences(const std::vector<std::weak_ptr<Process>>& rProcessRefs)
{
    KRATOS_EXPECT_EQ(rProcessRefs.size(), mNumberOfExpectedProcesses);
    for (const auto& process_ref : rProcessRefs)
    {
        KRATOS_EXPECT_FALSE(process_ref.expired())

        std::shared_ptr<Process> process = process_ref.lock();
        KRATOS_EXPECT_EQ(process->Check(), 0);
    }
}

}
