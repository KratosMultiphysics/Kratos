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
#include "includes/kratos_components.h"
#include "mpi/mpi_environment.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "includes/parallel_environment.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorRankAndSize, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    KRATOS_CHECK_EQUAL(serial_communicator.Rank(), 0);
    KRATOS_CHECK_EQUAL(serial_communicator.Size(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorFromKratosComponents, KratosMPICoreFastSuite)
{
    KRATOS_CHECK_EQUAL(KratosComponents<DataCommunicator>::Has("Serial"), true);
    const DataCommunicator& r_serial = KratosComponents<DataCommunicator>::Get("Serial");
    KRATOS_CHECK_EQUAL(r_serial.IsDistributed(), false);
}

// Sum ////////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorSumInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    int local = 1;
    int result = serial_communicator.Sum(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, local);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorSumDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    double local = 2.0;
    double result = serial_communicator.Sum(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, local);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorSumArray1d, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    array_1d<double,3> local;
    local[0] = -1.0;
    local[1] =  0.0;
    local[2] =  1.0;
    array_1d<double,3> result = serial_communicator.Sum(local, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 3; i++)
        {
            KRATOS_CHECK_EQUAL(result[i], local[i]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSumIntVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    std::vector<int> local{1, 1};
    std::vector<int> output{-1, -1};

    // two-buffer version
    serial_communicator.Sum(local, output, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(output[i], -1);
        }
    }

    // return buffer version
    std::vector<int> returned_result = serial_communicator.Sum(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSumDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    std::vector<double> local{2.0, 2.0};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    serial_communicator.Sum(local, output, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(output[i], -1.0);
        }
    }

    // return buffer version
    std::vector<double> returned_result = serial_communicator.Sum(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
        }
    }
}

// Min ////////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMinInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    int local = 1;
    int result = serial_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, local);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMinDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    double local = 2.0;
    double result = serial_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, local);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMinArray1d, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    array_1d<double,3> local;
    local[0] = -1.0;
    local[1] =  0.0;
    local[2] =  1.0;
    array_1d<double,3> result = serial_communicator.Min(local, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 3; i++)
        {
            KRATOS_CHECK_EQUAL(result[i], local[i]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMinIntVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    std::vector<int> local{1, 1};
    std::vector<int> output{-1, -1};

    // two-buffer version
    serial_communicator.Min(local, output, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(output[i], -1);
        }
    }

    // return buffer version
    std::vector<int> returned_result = serial_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMinDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    std::vector<double> local{2.0, 2.0};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    serial_communicator.Min(local, output, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(output[i], -1.0);
        }
    }

    // return buffer version
    std::vector<double> returned_result = serial_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
        }
    }
}

// Max ////////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMaxInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    int local = 1;
    int result = serial_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, local);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMaxDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    double local = 2.0;
    double result = serial_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, local);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMaxArray1d, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    array_1d<double,3> local;
    local[0] = -1.0;
    local[1] =  0.0;
    local[2] =  1.0;
    array_1d<double,3> result = serial_communicator.Max(local, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 3; i++)
        {
            KRATOS_CHECK_EQUAL(result[i], local[i]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMaxIntVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    std::vector<int> local{1, 1};
    std::vector<int> output{-1, -1};

    // two-buffer version
    serial_communicator.Max(local, output, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(output[i], -1);
        }
    }

    // return buffer version
    std::vector<int> returned_result = serial_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
        }
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMaxDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    constexpr int root = 0;

    const int world_rank = r_world.Rank();

    std::vector<double> local{2.0, 2.0};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    serial_communicator.Max(local, output, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(output[i], -1.0);
        }
    }

    // return buffer version
    std::vector<double> returned_result = serial_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
        }
    }
}

// SumAll /////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorSumAllInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    int local = 1;
    int result = serial_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(result, local);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorSumAllDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    double local = 2.0;
    double result = serial_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(result, local);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorSumAllArray1d, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    array_1d<double,3> local;
    local[0] = -1.0;
    local[1] =  0.0;
    local[2] =  1.0;
    array_1d<double,3> result = serial_communicator.SumAll(local);
    for (int i = 0; i < 3; i++)
    {
        KRATOS_CHECK_EQUAL(result[i], local[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSumAllIntVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    std::vector<int> local{1, 1};
    std::vector<int> output{-1, -1};

    // two-buffer version
    serial_communicator.SumAll(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], -1);
    }

    // return buffer version
    std::vector<int> returned_result = serial_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSumAllDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    std::vector<double> local{2.0, 2.0};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    serial_communicator.SumAll(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], -1.0);
    }

    // return buffer version
    std::vector<double> returned_result = serial_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }
}

// MinAll /////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMinAllInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    int local = 1;
    int result = serial_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(result, local);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMinAllDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    double local = 2.0;
    double result = serial_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(result, local);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMinAllArray1d, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    array_1d<double,3> local;
    local[0] = -1.0;
    local[1] =  0.0;
    local[2] =  1.0;
    array_1d<double,3> result = serial_communicator.MinAll(local);
    for (int i = 0; i < 3; i++)
    {
        KRATOS_CHECK_EQUAL(result[i], local[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMinAllIntVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    std::vector<int> local{1, 1};
    std::vector<int> output{-1, -1};

    // two-buffer version
    serial_communicator.MinAll(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], -1);
    }

    // return buffer version
    std::vector<int> returned_result = serial_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMinAllDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    std::vector<double> local{2.0, 2.0};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    serial_communicator.MinAll(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], -1.0);
    }

    // return buffer version
    std::vector<double> returned_result = serial_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }
}

// MaxAll /////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMaxAllInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    int local = 1;
    int result = serial_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(result, local);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMaxAllDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    double local = 2.0;
    double result = serial_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(result, local);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMaxAllArray1d, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    array_1d<double,3> local;
    local[0] = -1.0;
    local[1] =  0.0;
    local[2] =  1.0;
    array_1d<double,3> result = serial_communicator.MaxAll(local);
    for (int i = 0; i < 3; i++)
    {
        KRATOS_CHECK_EQUAL(result[i], local[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMaxAllIntVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    std::vector<int> local{1, 1};
    std::vector<int> output{-1, -1};

    // two-buffer version
    serial_communicator.MaxAll(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], -1);
    }

    // return buffer version
    std::vector<int> returned_result = serial_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMaxAllDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    std::vector<double> local{2.0, 2.0};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    serial_communicator.MaxAll(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], -1.0);
    }

    // return buffer version
    std::vector<double> returned_result = serial_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }
}

// ScanSum ////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorScanSumInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    int local = 1;
    int result = serial_communicator.ScanSum(local);
    KRATOS_CHECK_EQUAL(result, local);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorScanSumDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    double local = 2.0;
    double result = serial_communicator.ScanSum(local);
    KRATOS_CHECK_EQUAL(result, local);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScanSumIntVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    std::vector<int> local{1, 1};
    std::vector<int> output{-1, -1};

    // two-buffer version
    serial_communicator.ScanSum(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], -1);
    }

    // return buffer version
    std::vector<int> returned_result = serial_communicator.ScanSum(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScanSumDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    std::vector<double> local{2.0, 2.0};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    serial_communicator.ScanSum(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], -1.0);
    }

    // return buffer version
    std::vector<double> returned_result = serial_communicator.ScanSum(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }
}

// SendRecv ///////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSendRecvInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();

    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();
    const int send_rank = world_rank + 1 == world_size ? 0 : world_rank + 1;
    const int recv_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;

    std::vector<int> send_buffer = {world_rank, world_rank};
    std::vector<int> recv_buffer = {-1, -1};

    // two-buffer version
    serial_communicator.SendRecv(send_buffer, send_rank, recv_buffer, recv_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], -1);
    }

    // return version
    std::vector<int> return_buffer = serial_communicator.SendRecv(send_buffer, send_rank, recv_rank);

    if (send_rank == world_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer[i]);
        }
    }
    else
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), 0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSendRecvDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();

    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();
    const int send_rank = world_rank + 1 == world_size ? 0 : world_rank + 1;
    const int recv_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;

    std::vector<double> send_buffer{2.0*world_rank, 2.0*world_rank};
    std::vector<double> recv_buffer{-1.0, -1.0};

    // two-buffer version
    serial_communicator.SendRecv(send_buffer, send_rank, recv_buffer, recv_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], -1.0);
    }

    // return version
    std::vector<double> return_buffer = serial_communicator.SendRecv(send_buffer, send_rank, recv_rank);

    if (send_rank == world_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer[i]);
        }
    }
    else
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), 0);
    }
}


KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorBroadcast, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_size-1;

    int send_buffer_int = 0;
    double send_buffer_double = 0.0;

    if (world_rank == send_rank)
    {
        send_buffer_int = 1;
        send_buffer_double = 2.0;
    }
    int send_buffer_int_reference = send_buffer_int;
    double send_buffer_double_reference = send_buffer_double;

    serial_communicator.Broadcast(send_buffer_int,send_rank);
    serial_communicator.Broadcast(send_buffer_double,send_rank);

    // Check that serial_communicator does nothing
    KRATOS_CHECK_EQUAL(send_buffer_int, send_buffer_int_reference);
    KRATOS_CHECK_EQUAL(send_buffer_double, send_buffer_double_reference);

    mpi_world_communicator.Broadcast(send_buffer_int,send_rank);
    mpi_world_communicator.Broadcast(send_buffer_double,send_rank);

    KRATOS_CHECK_EQUAL(send_buffer_int, 1);
    KRATOS_CHECK_EQUAL(send_buffer_double, 2.0);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorBroadcastVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_size-1;

    std::vector<int> send_buffer_int{0, 0};
    std::vector<double> send_buffer_double{0.0, 0.0};

    if (world_rank == send_rank)
    {
        for (int i = 0; i < 2; i++)
        {
            send_buffer_int[i] = 1;
            send_buffer_double[i] = 2.0;
        }
    }
    std::vector<int> send_buffer_int_reference(send_buffer_int);
    std::vector<double> send_buffer_double_reference(send_buffer_double);

    serial_communicator.Broadcast(send_buffer_int,send_rank);
    serial_communicator.Broadcast(send_buffer_double,send_rank);

    // Check that serial_communicator does nothing
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(send_buffer_int[i], send_buffer_int_reference[i]);
        KRATOS_CHECK_EQUAL(send_buffer_double[i], send_buffer_double_reference[i]);
    }

    mpi_world_communicator.Broadcast(send_buffer_int,send_rank);
    mpi_world_communicator.Broadcast(send_buffer_double,send_rank);

    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(send_buffer_int[i], 1);
        KRATOS_CHECK_EQUAL(send_buffer_double[i], 2.0);
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        if (world_rank == 0)
        {
            send_buffer_int.resize(3);
            send_buffer_int = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Broadcast(send_buffer_int, send_rank),"Input error in call to MPI_Bcast");
    }
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScatter, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

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
        for (int i = 0; i < 2*world_size; i++)
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

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        if (world_rank == 0)
        {
            send_buffer_int.resize(3);
            send_buffer_int = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Scatter(send_buffer_int, recv_buffer_int, send_rank),"Error");
    }
    // send rank has wrong size
    send_buffer_double.push_back(0.0);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Scatter(send_buffer_double, recv_buffer_double, send_rank),"Error");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScatterv, KratosMPICoreFastSuite)
{
    /* send message for ints is {0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, ...} (max 6 values per rank)
     * read only first <rank> values of the message per rank (up to 5 values per rank)
     * message containing doubles is double of message containing ints
     */
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_size-1;

    std::vector<int> send_buffer_int(0);
    std::vector<double> send_buffer_double(0);
    std::vector<int> send_sizes(0);
    std::vector<int> send_offsets(0);

    auto make_message_size = [](int rank) { return rank < 5 ? rank : 5; };
    auto make_message_distance = [](int rank, int padding) {
        return rank < 5 ? ((rank-1)*rank)/2 + rank*padding : rank*(5+padding) - 15;
    };

    const int recv_size = make_message_size(world_rank);
    std::vector<int> recv_buffer_int(recv_size, -1);
    std::vector<double> recv_buffer_double(recv_size, -1.0);

    const int message_padding = 1;
    if (world_rank == send_rank)
    {
        const int send_size = make_message_distance(world_size, message_padding);
        send_buffer_int.resize(send_size);
        send_buffer_double.resize(send_size);
        send_sizes.resize(world_size);
        send_offsets.resize(world_size);
        int counter = 0;
        for (int rank = 0; rank < world_size; rank++)
        {
            send_sizes[rank] = make_message_size(rank);
            send_offsets[rank] = make_message_distance(rank, message_padding);
            for (int i = 0; i < send_sizes[rank] + message_padding; i++, counter++)
            {
                send_buffer_int[counter] = rank;
                send_buffer_double[counter] = 2.0*rank;
            }
        }
    }

    serial_communicator.Scatterv(send_buffer_int, send_sizes, send_offsets, recv_buffer_int, send_rank);
    serial_communicator.Scatterv(send_buffer_double, send_sizes, send_offsets, recv_buffer_double, send_rank);

    for (int i = 0; i < recv_size; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer_int[i], -1);
        KRATOS_CHECK_EQUAL(recv_buffer_double[i], -1.0);
    }

    mpi_world_communicator.Scatterv(send_buffer_int, send_sizes, send_offsets, recv_buffer_int, send_rank);
    mpi_world_communicator.Scatterv(send_buffer_double, send_sizes, send_offsets, recv_buffer_double, send_rank);

    for (int i = 0; i < recv_size; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer_int[i], world_rank);
        KRATOS_CHECK_EQUAL(recv_buffer_double[i], 2.0*world_rank);
    }

    #ifdef KRATOS_DEBUG
    // send sizes do not match
    std::vector<int> wrong_send_sizes = send_sizes;
    if (world_rank == send_rank) {
        wrong_send_sizes[0] += 1;
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Scatterv(send_buffer_int, wrong_send_sizes, send_offsets, recv_buffer_int, send_rank),
        "Error");

    // sent message is too large
    std::vector<int> wrong_recv_message;
    if (world_rank == send_rank)
    {
        wrong_recv_message.resize(recv_buffer_int.size()-1);
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Scatterv(send_buffer_int, send_sizes, send_offsets, wrong_recv_message, send_rank),
        "Error");

    // sent offsets overflow
    std::vector<int> wrong_send_offsets = send_offsets;
    if (world_rank == send_rank)
    {
        wrong_send_offsets[world_rank - 1] += 5;
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Scatterv(send_buffer_int, send_sizes, wrong_send_offsets, recv_buffer_int, send_rank),
        "Error");

    #endif
}


KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorGather, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int recv_rank = 0;

    std::vector<int> send_buffer_int{world_rank, world_rank};
    std::vector<double> send_buffer_double{2.0*world_rank, 2.0*world_rank};

    std::vector<int> recv_buffer_int;
    std::vector<double> recv_buffer_double;

    if (world_rank == recv_rank)
    {
        recv_buffer_int = std::vector<int>(2*world_size, -1);
        recv_buffer_double = std::vector<double>(2*world_size, -1.0);
    }

    serial_communicator.Gather(send_buffer_int, recv_buffer_int, recv_rank);
    serial_communicator.Gather(send_buffer_double, recv_buffer_double, recv_rank);

    if (world_rank == recv_rank)
    {
        for (int i = 0; i < 2*world_size; i++)
        {
            KRATOS_CHECK_EQUAL(recv_buffer_int[i], -1);
            KRATOS_CHECK_EQUAL(recv_buffer_double[i], -1.0);
        }
    }

    mpi_world_communicator.Gather(send_buffer_int, recv_buffer_int, recv_rank);
    mpi_world_communicator.Gather(send_buffer_double, recv_buffer_double, recv_rank);

    if (world_rank == recv_rank)
    {
        for (int rank = 0; rank < world_size; rank++)
        {
            for (int j = 2*rank; j < 2*rank+2; j++)
            {
                KRATOS_CHECK_EQUAL(recv_buffer_int[j], rank);
                KRATOS_CHECK_EQUAL(recv_buffer_double[j], 2.0*rank);
            }
        }
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        if (world_rank == 0)
        {
            send_buffer_int.resize(3);
            send_buffer_int = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Gather(send_buffer_int, recv_buffer_int, recv_rank),"Error");
    }
    // recv rank has wrong size
    recv_buffer_double.push_back(0.0);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Gather(send_buffer_double, recv_buffer_double, recv_rank),"Error");
    #endif
}


KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorGatherv, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int recv_rank = world_size-1;

    auto make_message_size = [](int rank) { return rank < 5 ? rank : 5; };
    auto make_message_distance = [](int rank, int padding) {
        return rank < 5 ? ((rank-1)*rank)/2 + rank*padding : rank*(5+padding) - 15;
    };

    const int message_padding = 1;
    const int send_size = make_message_size(world_rank);
    const int recv_size = make_message_distance(world_size, message_padding);

    std::vector<int> send_buffer_int(send_size, world_rank);
    std::vector<double> send_buffer_double(send_size, 2.0*world_rank);
    std::vector<int> recv_buffer_int(0);
    std::vector<double> recv_buffer_double(0);
    std::vector<int> recv_sizes(0);
    std::vector<int> recv_offsets(0);

    if (world_rank == recv_rank)
    {
        recv_buffer_int.resize(recv_size, -1);
        recv_buffer_double.resize(recv_size, -1.0);
        recv_sizes.resize(world_size);
        recv_offsets.resize(world_size);
        for (int rank = 0; rank < world_size; rank++)
        {
            recv_sizes[rank] = make_message_size(rank);
            recv_offsets[rank] = make_message_distance(rank, message_padding);
        }
    }

    serial_communicator.Gatherv(send_buffer_int, recv_buffer_int, recv_sizes, recv_offsets, recv_rank);
    serial_communicator.Gatherv(send_buffer_double, recv_buffer_double, recv_sizes, recv_offsets, recv_rank);

    if (world_rank == recv_rank)
    {
        for (int i = 0; i < recv_size; i++)
        {
            KRATOS_CHECK_EQUAL(recv_buffer_int[i], -1);
            KRATOS_CHECK_EQUAL(recv_buffer_double[i], -1.0);
        }
    }

    mpi_world_communicator.Gatherv(send_buffer_int, recv_buffer_int, recv_sizes, recv_offsets, recv_rank);
    mpi_world_communicator.Gatherv(send_buffer_double, recv_buffer_double, recv_sizes, recv_offsets, recv_rank);

    /* send message is {rank,} repeated <rank> times (up to 5) for ints and {2.*rank,} for doubles.
     * read message assumes 1 extra position per rank, so that
     * there are some uninitialized padding values on the recv message.
     * This is essentially the inverse of the test DataCommunicatorScatterv
     */

    if (world_rank == recv_rank)
    {
        for (int rank = 0; rank < world_size; rank++)
        {
            int recv_size = make_message_size(rank);
            int recv_offset = make_message_distance(rank, message_padding);
            // the message from this rank...
            for (int i = recv_offset; i < recv_offset + recv_size; i++)
            {
                KRATOS_CHECK_EQUAL(recv_buffer_int[i], rank);
                KRATOS_CHECK_EQUAL(recv_buffer_double[i], 2.0*rank);
            }
            // ...followed by the expected padding.
            for (int i = recv_offset + recv_size; i < recv_offset + recv_size + message_padding; i++)
            {
                KRATOS_CHECK_EQUAL(recv_buffer_int[i], -1);
                KRATOS_CHECK_EQUAL(recv_buffer_double[i], -1.0);
            }

        }
    }

    #ifdef KRATOS_DEBUG
    // recv sizes do not match
    std::vector<int> wrong_recv_sizes = recv_sizes;
    if (world_rank == recv_rank) {
        wrong_recv_sizes[0] += 1;
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Gatherv(
            send_buffer_int, recv_buffer_int, wrong_recv_sizes, recv_offsets, recv_rank),
            "Error");

    // recv message is too small
    std::vector<int> wrong_recv_message;
    if (world_rank == recv_size)
    {
        wrong_recv_message.resize(recv_buffer_int.size()-1);
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Gatherv(send_buffer_int, wrong_recv_message, recv_sizes, recv_offsets, recv_rank),
        "Error");

    // sent offsets overflow
    std::vector<int> wrong_recv_offsets = recv_offsets;
    if (world_rank == recv_rank)
    {
        wrong_recv_offsets[world_rank - 1] += 5;
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Gatherv(send_buffer_int, recv_buffer_int, recv_sizes, wrong_recv_offsets, recv_rank),
        "Error");

    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorAllGather, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();

    std::vector<int> send_buffer_int{world_rank, world_rank};
    std::vector<double> send_buffer_double{2.0*world_rank, 2.0*world_rank};

    std::vector<int> recv_buffer_int(2*world_size, -1);
    std::vector<double> recv_buffer_double(2*world_size, -1.0);

    serial_communicator.AllGather(send_buffer_int, recv_buffer_int);
    serial_communicator.AllGather(send_buffer_double, recv_buffer_double);

    for (int i = 0; i < 2*world_size; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer_int[i], -1);
        KRATOS_CHECK_EQUAL(recv_buffer_double[i], -1.0);
    }

    mpi_world_communicator.AllGather(send_buffer_int, recv_buffer_int);
    mpi_world_communicator.AllGather(send_buffer_double, recv_buffer_double);

    for (int rank = 0; rank < world_size; rank++)
    {
        for (int j = 2*rank; j < 2*rank+2; j++)
        {
            KRATOS_CHECK_EQUAL(recv_buffer_int[j], rank);
            KRATOS_CHECK_EQUAL(recv_buffer_double[j], 2.0*rank);
        }
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        send_buffer_int.resize(3);
        send_buffer_int = {1,2,3};
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            mpi_world_communicator.AllGather(send_buffer_int, recv_buffer_int),
            "Input error in call to MPI_Allgather");
    }
    // recv rank has wrong size
    recv_buffer_double.push_back(0.0);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.AllGather(send_buffer_double, recv_buffer_double),
        "Input error in call to MPI_Allgather");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorErrorBroadcasting, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // The serial communicator does not throw,
    // since there are no "other ranks" to broadcast the error to.
    // All these functions need to do is to pass along the bool condition.
    KRATOS_CHECK_EQUAL(serial_communicator.BroadcastErrorIfTrue(true, 0), true);
    KRATOS_CHECK_EQUAL(serial_communicator.BroadcastErrorIfTrue(false, 0), false);
    KRATOS_CHECK_EQUAL(serial_communicator.BroadcastErrorIfFalse(true, 0), true);
    KRATOS_CHECK_EQUAL(serial_communicator.BroadcastErrorIfFalse(false, 0), false);

    KRATOS_CHECK_EQUAL(serial_communicator.ErrorIfTrueOnAnyRank(true), true);
    KRATOS_CHECK_EQUAL(serial_communicator.ErrorIfTrueOnAnyRank(false), false);
    KRATOS_CHECK_EQUAL(serial_communicator.ErrorIfFalseOnAnyRank(true), true);
    KRATOS_CHECK_EQUAL(serial_communicator.ErrorIfFalseOnAnyRank(false), false);
}

}
}