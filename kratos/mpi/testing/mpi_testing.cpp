//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//
//

// System includes

// External includes

// Project includes
#include "mpi/testing/mpi_testing.h"

// This is tipically neot needed, but for mpi we need to register the new testing environment
int main(int argc, char* argv[]) 
{
    ::testing::InitGoogleTest(&argc, argv);
    ::testing::AddGlobalTestEnvironment(new Kratos::Testing::KratosMpiTestEnv);
    return RUN_ALL_TESTS();
}
