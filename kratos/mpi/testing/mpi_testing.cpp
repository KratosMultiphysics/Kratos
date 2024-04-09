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
    // Initialize MPI
    // int err = MPI_Init(&argc, &argv);
    MPI_Init(&argc, &argv);

    // Initialize the tests
    ::testing::InitGoogleTest(&argc, argv);

    // Get the size and rank 
    // TODO: we should ask this from the KratosMpiTestEnv, not from MPI directly, but for now it will have to suffice.
    int rank = 0;
    int size = 0;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Remove the default listener
    testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
    auto default_printer = listeners.Release(listeners.default_result_printer());

    // Create a configurable listener
    Kratos::Testing::ConfigurableEventListener *listener = new Kratos::Testing::ConfigurableEventListener(default_printer);

    // Set the listener configuration (by default false to all if rank != 0)
    if (rank != 0) {
        listener->showStart = false;
        listener->showIterations = false;
        listener->showEnvironment = false;
        listener->showTestCases = false;
        listener->showTestNames = false;
        listener->showSuccesses = false;
        listener->showInlineFailures = false;
        listener->showResult = false;
        listener->showEnd = false;
    }

    std::cout << "Initializing GTEST MPI environment with " << size << " ranks (" << rank << ")" << std::endl;
    
    // Add the MPI environment to the test 
    ::testing::AddGlobalTestEnvironment(new Kratos::Testing::KratosMpiTestEnv);

    // Add our listener
    listeners.Append(listener);

    // Run the tests
    return RUN_ALL_TESTS();

    // Finalize MPI
    // err = MPI_Finalize();
    MPI_Finalize();
}
