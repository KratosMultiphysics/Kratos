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
#include "mpi.h"
#include "mpi/testing/mpi_testing.h"

// Create a custom main with the MPI environment and custom listeners for the test output
int main(int argc, char* argv[]) 
{
    return Kratos::Testing::MPIGTestMain::InitializeMPITesting(argc, argv);
}
