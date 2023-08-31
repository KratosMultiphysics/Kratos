//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

// External includes
#include "mpi.h"

// Project includes
#include "includes/parallel_environment.h"
#include "mpi/testing/mpi_testing.h"

namespace Kratos::Testing {

KRATOS_TEST_CASE_IN_SUITE(MPIEnvironmentSetUp, KratosMPICoreFastSuite)
{
    KRATOS_EXPECT_EQ(ParallelEnvironment::HasDataCommunicator("World"), true);
    KRATOS_EXPECT_EQ(ParallelEnvironment::HasDataCommunicator("Serial"), true);
    KRATOS_EXPECT_EQ(ParallelEnvironment::HasDataCommunicator("NotReallyACommunicator"), false);
}

KRATOS_TEST_CASE_IN_SUITE(MPIEnvironmentDefaultComms, KratosMPICoreFastSuite)
{
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // r_default_comm should be "World" in MPI runs
    DataCommunicator& r_default_comm = ParallelEnvironment::GetDefaultDataCommunicator();
    DataCommunicator& r_world = ParallelEnvironment::GetDataCommunicator("World");
    DataCommunicator& r_serial = ParallelEnvironment::GetDataCommunicator("Serial");

    KRATOS_EXPECT_EQ(r_default_comm.IsDistributed(), true);
    KRATOS_EXPECT_EQ(r_world.IsDistributed(), true);
    KRATOS_EXPECT_EQ(r_serial.IsDistributed(), false);

    KRATOS_EXPECT_EQ(r_default_comm.Rank(), rank);
    KRATOS_EXPECT_EQ(r_world.Rank(), rank);
    KRATOS_EXPECT_EQ(r_serial.Rank(), 0);

    KRATOS_EXPECT_EQ(r_default_comm.Size(), size);
    KRATOS_EXPECT_EQ(r_world.Size(), size);
    KRATOS_EXPECT_EQ(r_serial.Size(), 1);
}

}