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
#include "testing/testing.h"
#include "includes/parallel_environment.h"

namespace Kratos::Testing 
{

void KratosTestEnv::SetUp() {
    std::cout << "KratosTestEnv::SetUp" << std::endl;
}

void KratosTestEnv::TearDown() {
    std::cout << "KratosTestEnv::TearDown" << std::endl;
}

KratosTestEnv::KratosTestEnv() {
    std::cout << "KratosTestEnv::KratosTestEnv" << std::endl;
}

void KratosCoreFastSuite::SetUp() { 
    mCoutBuffer = std::cout.rdbuf();
    mCerrBuffer = std::cerr.rdbuf();

    std::cout.rdbuf(mStream.rdbuf());
    std::cerr.rdbuf(mStream.rdbuf());   
}

void KratosCoreFastSuite::TearDown() { 
    // std::cout.rdbuf(mCoutBuffer);
    // std::cerr.rdbuf(mCerrBuffer); 
}

DataCommunicator& GetDefaultDataCommunicator()
{
    return ParallelEnvironment::GetDefaultDataCommunicator();
}

} // namespace Kratos::Testing