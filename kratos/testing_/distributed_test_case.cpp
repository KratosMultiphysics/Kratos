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
#include "testing/testing.h"
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
    const bool success_on_this_rank = GetResult().IsSucceed();
    const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();
    const bool global_success = r_comm.AndReduceAll(success_on_this_rank);
    if (success_on_this_rank && !global_success) {
        TestCaseResult remote_failure(GetResult());
        remote_failure.SetToFailed();
        remote_failure.SetErrorMessage("Test was reported as successful on this rank, but failed on a different rank.");
        SetResult(remote_failure);
    }
}

///@}

} // namespace Kratos::Testing.
