//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//
//

#include "mpi.h"
#include "mpi/mpi_environment.h"
#include "includes/parallel_environment.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(MPIEnvironmentSetUp, KratosMPICoreFastSuite)
{
    ParallelEnvironment& parallel_environment = ParallelEnvironment::GetInstance();
    KRATOS_CHECK_EQUAL(parallel_environment.HasDataCommunicator("World"), true);
    KRATOS_CHECK_EQUAL(parallel_environment.HasDataCommunicator("Serial"), true);
    KRATOS_CHECK_EQUAL(parallel_environment.HasDataCommunicator("NotReallyACommunicator"), false);
}

KRATOS_TEST_CASE_IN_SUITE(MPIEnvironmentDefaultComms, KratosMPICoreFastSuite)
{
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    ParallelEnvironment& parallel_environment = ParallelEnvironment::GetInstance();

    // r_default_comm should be "World" in MPI runs
    DataCommunicator& r_default_comm = parallel_environment.GetDefaultDataCommunicator();
    DataCommunicator& r_world = parallel_environment.GetDataCommunicator("World");
    DataCommunicator& r_serial = parallel_environment.GetDataCommunicator("Serial");

    KRATOS_CHECK_EQUAL(r_default_comm.IsDistributed(), true);
    KRATOS_CHECK_EQUAL(r_world.IsDistributed(), true);
    KRATOS_CHECK_EQUAL(r_serial.IsDistributed(), false);

    KRATOS_CHECK_EQUAL(r_default_comm.Rank(), rank);
    KRATOS_CHECK_EQUAL(r_world.Rank(), rank);
    KRATOS_CHECK_EQUAL(r_serial.Rank(), 0);

    KRATOS_CHECK_EQUAL(r_default_comm.Size(), size);
    KRATOS_CHECK_EQUAL(r_world.Size(), size);
    KRATOS_CHECK_EQUAL(r_serial.Size(), 1);
}

}
}