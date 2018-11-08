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

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSum, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);
    constexpr int root = 0;

    int int_sum = 1;
    double double_sum = 2.0;
    array_1d<double,3> array_sum;
    array_sum[0] = -1.0;
    array_sum[1] =  0.0;
    array_sum[2] =  1.0;

    // local version: do nothing
    serial_communicator.Sum(int_sum, root);
    KRATOS_CHECK_EQUAL(int_sum, 1);

    serial_communicator.Sum(double_sum, root);
    KRATOS_CHECK_EQUAL(double_sum, 2.0);

    serial_communicator.Sum(array_sum, root);
    KRATOS_CHECK_EQUAL(array_sum[0], -1.0);
    KRATOS_CHECK_EQUAL(array_sum[1],  0.0);
    KRATOS_CHECK_EQUAL(array_sum[2],  1.0);

    // MPI version
    mpi_world_communicator.Sum(double_sum, root);
    mpi_world_communicator.Sum(int_sum, root);
    mpi_world_communicator.Sum(array_sum, root);

    int world_size = mpi_world_communicator.Size();
    int world_rank = mpi_world_communicator.Rank();
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(int_sum, world_size);

        KRATOS_CHECK_EQUAL(double_sum, 2.*world_size);

        KRATOS_CHECK_EQUAL(array_sum[0], -1.0*world_size);
        KRATOS_CHECK_EQUAL(array_sum[1],  0.0);
        KRATOS_CHECK_EQUAL(array_sum[2],  1.0*world_size);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMin, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);
    constexpr int root = 0;

    int world_rank = mpi_world_communicator.Rank();

    int int_min = world_rank;
    double double_min = 2.0*world_rank;
    array_1d<double,3> array_min;
    array_min[0] = -1.0*world_rank;
    array_min[1] =  0.0;
    array_min[2] =  1.0*world_rank;

    // local version: do nothing
    serial_communicator.Min(int_min, root);
    KRATOS_CHECK_EQUAL(int_min, world_rank);

    serial_communicator.Min(double_min, root);
    KRATOS_CHECK_EQUAL(double_min, 2.0*world_rank);

    serial_communicator.Min(array_min, root);
    KRATOS_CHECK_EQUAL(array_min[0], -1.0*world_rank);
    KRATOS_CHECK_EQUAL(array_min[1],  0.0);
    KRATOS_CHECK_EQUAL(array_min[2],  1.0*world_rank);

    // MPI version
    mpi_world_communicator.Min(int_min, root);
    mpi_world_communicator.Min(double_min, root);
    mpi_world_communicator.Min(array_min, root);

    if (mpi_world_communicator.Rank() == root)
    {
        KRATOS_CHECK_EQUAL(int_min, 0.0);

        KRATOS_CHECK_EQUAL(double_min, 0.0);

        int world_size = mpi_world_communicator.Size();
        KRATOS_CHECK_EQUAL(array_min[0], -1.0*(world_size-1));
        KRATOS_CHECK_EQUAL(array_min[1],  0.0);
        KRATOS_CHECK_EQUAL(array_min[2],  0.0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMax, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);
    constexpr int root = 0;

    int world_rank = mpi_world_communicator.Rank();

    int int_max = world_rank;
    double double_max = 2.0*world_rank;
    array_1d<double,3> array_max;
    array_max[0] = -1.0*world_rank;
    array_max[1] =  0.0;
    array_max[2] =  1.0*world_rank;

    // local version: do nothing
    serial_communicator.Max(int_max, root);
    KRATOS_CHECK_EQUAL(int_max, world_rank);

    serial_communicator.Max(double_max, root);
    KRATOS_CHECK_EQUAL(double_max, 2.0*world_rank);

    serial_communicator.Max(array_max, root);
    KRATOS_CHECK_EQUAL(array_max[0], -1.0*world_rank);
    KRATOS_CHECK_EQUAL(array_max[1],  0.0);
    KRATOS_CHECK_EQUAL(array_max[2],  1.0*world_rank);

    // MPI version
    mpi_world_communicator.Max(int_max, root);
    mpi_world_communicator.Max(double_max, root);
    mpi_world_communicator.Max(array_max, root);

    if (world_rank == root)
    {
        int world_size = mpi_world_communicator.Size();
        KRATOS_CHECK_EQUAL(int_max, world_size-1);

        KRATOS_CHECK_EQUAL(double_max, 2.0*(world_size-1));

        KRATOS_CHECK_EQUAL(array_max[0], 0.0);
        KRATOS_CHECK_EQUAL(array_max[1], 0.0);
        KRATOS_CHECK_EQUAL(array_max[2], 1.0*(world_size-1));
    }
}


KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSumAll, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    int int_sum = 1;
    double double_sum = 2.0;
    array_1d<double,3> array_sum;
    array_sum[0] = -1.0;
    array_sum[1] =  0.0;
    array_sum[2] =  1.0;

    // local version: do nothing
    serial_communicator.SumAll(int_sum);
    KRATOS_CHECK_EQUAL(int_sum, 1);

    serial_communicator.SumAll(double_sum);
    KRATOS_CHECK_EQUAL(double_sum, 2.0);

    serial_communicator.SumAll(array_sum);
    KRATOS_CHECK_EQUAL(array_sum[0], -1.0);
    KRATOS_CHECK_EQUAL(array_sum[1],  0.0);
    KRATOS_CHECK_EQUAL(array_sum[2],  1.0);

    // MPI version
    int world_size = mpi_world_communicator.Size();
    mpi_world_communicator.SumAll(int_sum);
    KRATOS_CHECK_EQUAL(int_sum, world_size);

    mpi_world_communicator.SumAll(double_sum);
    KRATOS_CHECK_EQUAL(double_sum, 2.*world_size);

   mpi_world_communicator.SumAll(array_sum);
   KRATOS_CHECK_EQUAL(array_sum[0], -1.0*world_size);
   KRATOS_CHECK_EQUAL(array_sum[1],  0.0);
   KRATOS_CHECK_EQUAL(array_sum[2],  1.0*world_size);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMinAll, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    int world_rank = mpi_world_communicator.Rank();

    int int_min = world_rank;
    double double_min = 2.0*world_rank;
    array_1d<double,3> array_min;
    array_min[0] = -1.0*world_rank;
    array_min[1] =  0.0;
    array_min[2] =  1.0*world_rank;

    // local version: do nothing
    serial_communicator.MinAll(int_min);
    KRATOS_CHECK_EQUAL(int_min, world_rank);

    serial_communicator.MinAll(double_min);
    KRATOS_CHECK_EQUAL(double_min, 2.0*world_rank);

    serial_communicator.MinAll(array_min);
    KRATOS_CHECK_EQUAL(array_min[0], -1.0*world_rank);
    KRATOS_CHECK_EQUAL(array_min[1],  0.0);
    KRATOS_CHECK_EQUAL(array_min[2],  1.0*world_rank);

    // MPI version
    mpi_world_communicator.MinAll(int_min);
    KRATOS_CHECK_EQUAL(int_min, 0.0);

    mpi_world_communicator.MinAll(double_min);
    KRATOS_CHECK_EQUAL(double_min, 0.0);

    mpi_world_communicator.MinAll(array_min);
    int world_size = mpi_world_communicator.Size();
    KRATOS_CHECK_EQUAL(array_min[0], -1.0*(world_size-1));
    KRATOS_CHECK_EQUAL(array_min[1],  0.0);
    KRATOS_CHECK_EQUAL(array_min[2],  0.0);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMaxAll, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    int world_rank = mpi_world_communicator.Rank();

    int int_max = world_rank;
    double double_max = 2.0*world_rank;
    array_1d<double,3> array_max;
    array_max[0] = -1.0*world_rank;
    array_max[1] =  0.0;
    array_max[2] =  1.0*world_rank;

    // local version: do nothing
    serial_communicator.MaxAll(int_max);
    KRATOS_CHECK_EQUAL(int_max, world_rank);

    serial_communicator.MaxAll(double_max);
    KRATOS_CHECK_EQUAL(double_max, 2.0*world_rank);

    serial_communicator.MaxAll(array_max);
    KRATOS_CHECK_EQUAL(array_max[0], -1.0*world_rank);
    KRATOS_CHECK_EQUAL(array_max[1],  0.0);
    KRATOS_CHECK_EQUAL(array_max[2],  1.0*world_rank);

    // MPI version
    int world_size = mpi_world_communicator.Size();
    mpi_world_communicator.MaxAll(int_max);
    KRATOS_CHECK_EQUAL(int_max, world_size-1);

    mpi_world_communicator.MaxAll(double_max);
    KRATOS_CHECK_EQUAL(double_max, 2.0*(world_size-1));

    mpi_world_communicator.MaxAll(array_max);
    KRATOS_CHECK_EQUAL(array_max[0], 0.0);
    KRATOS_CHECK_EQUAL(array_max[1], 0.0);
    KRATOS_CHECK_EQUAL(array_max[2], 1.0*(world_size-1));
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScanSum, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    int local_total_int = 1;
    double local_total_double = 2.0;

    // local version: do nothing
    int partial_sum_int = serial_communicator.ScanSum(local_total_int);
    double partial_sum_double = serial_communicator.ScanSum(local_total_double);
    KRATOS_CHECK_EQUAL(partial_sum_int, local_total_int);
    KRATOS_CHECK_EQUAL(partial_sum_double, local_total_double);

    // MPI version
    partial_sum_int = mpi_world_communicator.ScanSum(local_total_int);
    partial_sum_double = mpi_world_communicator.ScanSum(local_total_double);
    int mpi_world_rank = mpi_world_communicator.Rank();
    KRATOS_CHECK_EQUAL(partial_sum_int, mpi_world_rank + 1);
    KRATOS_CHECK_EQUAL(partial_sum_double, 2.0*(mpi_world_rank + 1) );
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSendRecv, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_rank + 1 == world_size ? 0 : world_rank + 1;
    const int recv_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;

    std::vector<int> send_buffer_int = {world_rank, world_rank};
    std::vector<int> recv_buffer_int = {0, 0};
    std::vector<double> send_buffer_double = {2.0*world_rank, 2.0*world_rank};
    std::vector<double> recv_buffer_double = {0.0, 0.0};

    serial_communicator.SendRecv(send_buffer_int, send_rank, recv_buffer_int, recv_rank);
    serial_communicator.SendRecv(send_buffer_double, send_rank, recv_buffer_double, recv_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer_int[i], 0);
        KRATOS_CHECK_EQUAL(recv_buffer_double[i], 0.0);
    }

    if (world_size > 1)
    {
        mpi_world_communicator.SendRecv(send_buffer_int, send_rank, recv_buffer_int, recv_rank);
        mpi_world_communicator.SendRecv(send_buffer_double, send_rank, recv_buffer_double, recv_rank);

        const int expected_recv_int = world_rank > 0 ? world_rank - 1 : world_size - 1;
        const double expected_recv_double = 2.0*expected_recv_int;

        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(recv_buffer_int[i], expected_recv_int);
            KRATOS_CHECK_EQUAL(recv_buffer_double[i], expected_recv_double);
        }
    }
}


KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScatter, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator = DataCommunicator();
    MPIDataCommunicator mpi_world_communicator = MPIDataCommunicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = 0;

    std::vector<int> send_buffer_int(0);
    std::vector<double> send_buffer_double(0);
    std::vector<int> recv_buffer_int = {0, 0};
    std::vector<double> recv_buffer_double = {0.0, 0.0};

    if (world_rank == send_rank)
    {
        send_buffer_int.resize(2*world_size);
        send_buffer_double.resize(2*world_size);
        for (unsigned int i = 0; i < 2*world_size; i++)
        {
            send_buffer_int[i] = 1;
            send_buffer_double[i] = 2.0;
        }
    }

    serial_communicator.Scatter(send_buffer_int, recv_buffer_int, send_rank);
    serial_communicator.Scatter(send_buffer_double, recv_buffer_double, send_rank);

    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer_int[i], 0);
        KRATOS_CHECK_EQUAL(recv_buffer_double[i], 0.0);
    }

    mpi_world_communicator.Scatter(send_buffer_int, recv_buffer_int, send_rank);
    mpi_world_communicator.Scatter(send_buffer_double, recv_buffer_double, send_rank);

    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer_int[i], 1);
        KRATOS_CHECK_EQUAL(recv_buffer_double[i], 2.0);
    }
}

}
}