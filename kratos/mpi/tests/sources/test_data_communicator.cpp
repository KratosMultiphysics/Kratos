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

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSumAll, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    int int_sum = 1;
    double double_sum = 2.0f;
    array_1d<double,3> array_sum;
    array_sum[0] = -1.0f;
    array_sum[1] =  0.0f;
    array_sum[2] =  1.0f;

    // local version: do nothing
    serial_communicator.SumAll(int_sum);
    KRATOS_CHECK_EQUAL(int_sum, 1);

    serial_communicator.SumAll(double_sum);
    KRATOS_CHECK_EQUAL(double_sum, 2.0f);

    serial_communicator.SumAll(array_sum);
    KRATOS_CHECK_EQUAL(array_sum[0], -1.0f);
    KRATOS_CHECK_EQUAL(array_sum[1],  0.0f);
    KRATOS_CHECK_EQUAL(array_sum[2],  1.0f);

    // MPI version
    int world_size = mpi_world_communicator.Size();
    mpi_world_communicator.SumAll(int_sum);
    KRATOS_CHECK_EQUAL(int_sum, world_size);

    mpi_world_communicator.SumAll(double_sum);
    KRATOS_CHECK_EQUAL(double_sum, 2.*world_size);

   mpi_world_communicator.SumAll(array_sum);
   KRATOS_CHECK_EQUAL(array_sum[0], -1.0f*world_size);
   KRATOS_CHECK_EQUAL(array_sum[1],  0.0f);
   KRATOS_CHECK_EQUAL(array_sum[2],  1.0f*world_size);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMinAll, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    int world_rank = mpi_world_communicator.Rank();

    int int_min = world_rank;
    double double_min = 2.0f*world_rank;
    array_1d<double,3> array_min;
    array_min[0] = -1.0f*world_rank;
    array_min[1] =  0.0f;
    array_min[2] =  1.0f*world_rank;

    // local version: do nothing
    serial_communicator.MinAll(int_min);
    KRATOS_CHECK_EQUAL(int_min, world_rank);

    serial_communicator.MinAll(double_min);
    KRATOS_CHECK_EQUAL(double_min, 2.0f*world_rank);

    serial_communicator.MinAll(array_min);
    KRATOS_CHECK_EQUAL(array_min[0], -1.0f*world_rank);
    KRATOS_CHECK_EQUAL(array_min[1],  0.0f);
    KRATOS_CHECK_EQUAL(array_min[2],  1.0f*world_rank);

    // MPI version
    mpi_world_communicator.MinAll(int_min);
    KRATOS_CHECK_EQUAL(int_min, 0.0f);

    mpi_world_communicator.MinAll(double_min);
    KRATOS_CHECK_EQUAL(double_min, 0.0f);

    mpi_world_communicator.MinAll(array_min);
    int world_size = mpi_world_communicator.Size();
    KRATOS_CHECK_EQUAL(array_min[0], -1.0f*(world_size-1));
    KRATOS_CHECK_EQUAL(array_min[1],  0.0f);
    KRATOS_CHECK_EQUAL(array_min[2],  0.0f);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMaxAll, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    int world_rank = mpi_world_communicator.Rank();

    int int_max = world_rank;
    double double_max = 2.0f*world_rank;
    array_1d<double,3> array_max;
    array_max[0] = -1.0f*world_rank;
    array_max[1] =  0.0f;
    array_max[2] =  1.0f*world_rank;

    // local version: do nothing
    serial_communicator.MaxAll(int_max);
    KRATOS_CHECK_EQUAL(int_max, world_rank);

    serial_communicator.MaxAll(double_max);
    KRATOS_CHECK_EQUAL(double_max, 2.0f*world_rank);

    serial_communicator.MaxAll(array_max);
    KRATOS_CHECK_EQUAL(array_max[0], -1.0f*world_rank);
    KRATOS_CHECK_EQUAL(array_max[1],  0.0f);
    KRATOS_CHECK_EQUAL(array_max[2],  1.0f*world_rank);

    // MPI version
    int world_size = mpi_world_communicator.Size();
    mpi_world_communicator.MaxAll(int_max);
    KRATOS_CHECK_EQUAL(int_max, world_size-1);

    mpi_world_communicator.MaxAll(double_max);
    KRATOS_CHECK_EQUAL(double_max, 2.0f*(world_size-1));

    mpi_world_communicator.MaxAll(array_max);
    KRATOS_CHECK_EQUAL(array_max[0], 0.0f);
    KRATOS_CHECK_EQUAL(array_max[1], 0.0f);
    KRATOS_CHECK_EQUAL(array_max[2], 1.0f*(world_size-1));
}

}
}