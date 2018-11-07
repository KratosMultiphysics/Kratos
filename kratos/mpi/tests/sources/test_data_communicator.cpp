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

#include "includes/data_communicator.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/mpi_environment.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorRankAndSize, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();

    KRATOS_CHECK_EQUAL(serial_communicator.Rank(), 0);
    KRATOS_CHECK_EQUAL(serial_communicator.Size(), 1);

    MPIDataCommunicator mpi_self_communicator = MPIDataCommunicator(MPI_COMM_SELF);

    KRATOS_CHECK_EQUAL(mpi_self_communicator.Rank(), 0);
    KRATOS_CHECK_EQUAL(mpi_self_communicator.Size(), 1);

    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);

    KRATOS_CHECK_EQUAL(mpi_world_communicator.Rank(), world_rank);
    KRATOS_CHECK_EQUAL(mpi_world_communicator.Size(), world_size);
}

KRATOS_TEST_CASE_IN_SUITE(MPICommRetrieval, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_self_communicator = MPIDataCommunicator(MPI_COMM_SELF);
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    KRATOS_CHECK_EQUAL(MPIEnvironment::GetMPICommunicator(serial_communicator), MPI_COMM_SELF);
    KRATOS_CHECK_EQUAL(MPIEnvironment::GetMPICommunicator(mpi_self_communicator), MPI_COMM_SELF);
    KRATOS_CHECK_EQUAL(MPIEnvironment::GetMPICommunicator(mpi_world_communicator), MPI_COMM_WORLD);
    KRATOS_CHECK_NOT_EQUAL(MPIEnvironment::GetMPICommunicator(mpi_world_communicator), MPI_COMM_SELF);
}

}
}