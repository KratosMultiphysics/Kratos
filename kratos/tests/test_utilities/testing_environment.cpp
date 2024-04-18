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

KratosTestEnv::KratosTestEnv() {
    std::cout << "KratosTestEnv::KratosTestEnv" << std::endl;
}

void KratosCoreFastSuite::SetUp() { 
    std::cout.rdbuf(mStreamBuffer.rdbuf());
    std::cerr.rdbuf(mStreamBuffer.rdbuf());   
}

void KratosCoreFastSuite::TearDown() { 

}

void KratosTestEnv::SetUp() {
    std::cout << "KratosTestEnv::SetUp" << std::endl;
}

void KratosTestEnv::TearDown() {
    std::cout << "KratosTestEnv::TearDown" << std::endl;
}

DataCommunicator& GetDefaultDataCommunicator()
{
    return ParallelEnvironment::GetDefaultDataCommunicator();
}

} // namespace Kratos::Testing