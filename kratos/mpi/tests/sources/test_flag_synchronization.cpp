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

#include "includes/kratos_flags.h"
#include "mpi/includes/mpi_data_communicator.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsAndAll, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int rank = mpi_world_communicator.Rank();
    const int size = mpi_world_communicator.Size();

    Kratos::Flags test_flag;
    test_flag.Set(STRUCTURE, rank == 0);

    Kratos::Flags synchronized_flag = mpi_world_communicator.AndAll(test_flag, STRUCTURE);

    if (size > 1)
    {
        KRATOS_CHECK_EQUAL(synchronized_flag.Is(STRUCTURE), false);
    }
    else // if there is only one rank, result should be true
    {
        KRATOS_CHECK_EQUAL(synchronized_flag.Is(STRUCTURE), true);
    }
    KRATOS_CHECK_EQUAL(synchronized_flag.IsDefined(PERIODIC), false);
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorFlagsOrAll, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_SELF);
    const int rank = mpi_world_communicator.Rank();
    //const int size = mpi_world_communicator.Size();

    Kratos::Flags test_flag;
    test_flag.Set(STRUCTURE, rank == 0);

    Kratos::Flags synchronized_flag = mpi_world_communicator.OrAll(test_flag, STRUCTURE);

    KRATOS_CHECK_EQUAL(synchronized_flag.Is(STRUCTURE), true);
    KRATOS_CHECK_EQUAL(synchronized_flag.IsDefined(PERIODIC), false);
}

}

}