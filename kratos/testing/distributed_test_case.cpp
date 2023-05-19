//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// System includes
#include <string>

// External includes

// Project includes
#include "includes/kernel.h"
#include "includes/parallel_environment.h"
#include "testing/distributed_test_case.h"

namespace Kratos::Testing {

DistributedTestCase::DistributedTestCase(std::string const& Name):
    TestCase(Name)
{}

DistributedTestCase::~DistributedTestCase() {}

void DistributedTestCase::Run()
{
    TestCase::Run();
    CheckRemoteFailure();
}

void DistributedTestCase::Profile()
{
    TestCase::Profile();
    CheckRemoteFailure();
}

bool DistributedTestCase::IsEnabled() const
{
    return TestCase::IsEnabled() && Kernel::IsDistributedRun();
}

bool DistributedTestCase::IsDisabled() const
{
    return !IsEnabled();
}

bool DistributedTestCase::IsDistributedTest() const
{
    return true;
}

std::string DistributedTestCase::Info() const
{
    return "Distributed test case " + Name();
}

void DistributedTestCase::CheckRemoteFailure()
{
    // The data communicator
    const DataCommunicator& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();

    // First we check if tests are run on all ranks
    const bool run_on_this_rank = GetResult().IsRun();
    const bool global_run = r_comm.AndReduceAll(run_on_this_rank);         // It is run on all ranks
    if (!global_run) {
        if (run_on_this_rank) {
            TestCaseResult remote_run(GetResult());
            remote_run.Reset();
            remote_run.SetErrorMessage("Test was reported as run on this rank, but did not run on a different rank.");
            SetResult(remote_run);
        } else {
            TestCaseResult rank_result(GetResult());
            rank_result.Reset();
            rank_result.SetErrorMessage("Test was not run on this rank.");
            SetResult(rank_result);
        }

        return;
    }

    // Test skipped on at least one rank.
    const bool skipped_on_this_rank = GetResult().IsSkipped();
    const bool global_skipped = r_comm.OrReduceAll(skipped_on_this_rank);  // It is skipped on some rank
    if (global_skipped) {
        if (!skipped_on_this_rank) { // This rank skipped, but another rank skipped.
            TestCaseResult remote_skip(GetResult());
            remote_skip.SetToSkipped();
            remote_skip.SetErrorMessage("Test was reported as not skipped on this rank, but skipped on a different rank.");
            SetResult(remote_skip);
        } else { // Test failed at least on this rank
            TestCaseResult rank_result(GetResult());
            rank_result.SetToSkipped();
            rank_result.SetErrorMessage("Test skipped at least on this rank.");
            SetResult(rank_result);
        }

        return;
    }

    // Test failed on at least one rank.
    const bool success_on_this_rank = GetResult().IsSucceed();
    const bool global_success = r_comm.AndReduceAll(success_on_this_rank); // It succeeds on all ranks
    if (!global_success) {
        if (success_on_this_rank) { // This rank succeeded, but another rank failed.
            TestCaseResult remote_failure(GetResult());
            remote_failure.SetToFailed();
            remote_failure.SetErrorMessage("Test was reported as successful on this rank, but failed on a different rank.");
            SetResult(remote_failure);
        } else { // Test failed at least on this rank
            TestCaseResult rank_result(GetResult());
            rank_result.SetToFailed();
            rank_result.SetErrorMessage("Test failed at least on this rank.");
            SetResult(rank_result);
        }
    }
}

///@}

} // namespace Kratos::Testing.
