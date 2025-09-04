//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

// System includes

// External includes

// Project includes
#include "includes/parallel_environment.h"

#include "tests/test_utilities/test_environment.h"

namespace Kratos::Testing 
{

void KratosTestEnv::SetUp() {
    // This is intentionally left empty.
    // Add here your code if you need to set up some stuff for your tests or debuging
}

void KratosTestEnv::TearDown() {
    // This is intentionally left empty.
    // Add here your code if you need to tear down some stuff for your tests or debuging
}

KratosTestEnv::KratosTestEnv() {
    // This is intentionally left empty.
    // Add here your code if you need to do any operation in initialization or debuging
}

DataCommunicator& GetDefaultDataCommunicator()
{
    return ParallelEnvironment::GetDefaultDataCommunicator();
}

} // namespace Kratos::Testing