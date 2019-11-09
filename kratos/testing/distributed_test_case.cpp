//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#include <string>

// External includes

// Project includes
#include "includes/kernel.h"
#include "testing/distributed_test_case.h"

namespace Kratos {
namespace Testing {

DistributedTestCase::DistributedTestCase(std::string const& Name):
    TestCase(Name)
{}

DistributedTestCase::~DistributedTestCase() {}

bool DistributedTestCase::IsEnabled() const
{
    return TestCase::IsEnabled() && Kernel::IsDistributedRun();
}

bool DistributedTestCase::IsDisabled() const
{
    return !IsEnabled();
}

///@}
///@name Input and output
///@{

/// Turn back information as a string.
std::string DistributedTestCase::Info() const
{
    return "Distributed test case " + Name();
}

///@}

}
}
