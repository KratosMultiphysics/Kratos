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

void StubTimeLoopExecutor::SetProcessReferences(std::vector<std::reference_wrapper<Process>> ProcessRefs)
{
    KRATOS_EXPECT_EQ(ProcessRefs.size(), mNumberOfExpectedProcesses);
    for (const auto process_ref : ProcessRefs)
    {
        KRATOS_EXPECT_EQ(process_ref.get().Check(), 0);
    }
}

StubTimeLoopExecutor::StubTimeLoopExecutor(const int NumberOfExpectedProcesses) :
    mNumberOfExpectedProcesses{NumberOfExpectedProcesses}
{
}

}
