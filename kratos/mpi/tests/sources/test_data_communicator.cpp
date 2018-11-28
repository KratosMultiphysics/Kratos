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

// Broadcast //////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorBroadcastInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_rank = r_world.Rank();
    const int send_rank = 0;

    int send = 1 + world_rank;
    serial_communicator.Broadcast(send,send_rank);
    KRATOS_CHECK_EQUAL(send, 1 + world_rank);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorBroadcastDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_rank = r_world.Rank();
    const int send_rank = 0;

    double send = 1.0 + world_rank;
    serial_communicator.Broadcast(send,send_rank);
    KRATOS_CHECK_EQUAL(send, 1.0 + world_rank);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorBroadcastIntVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_rank = r_world.Rank();
    const int send_rank = 0;

    std::vector<int> send = {world_rank, 1 + world_rank};
    serial_communicator.Broadcast(send,send_rank);
    KRATOS_CHECK_EQUAL(send[0], world_rank);
    KRATOS_CHECK_EQUAL(send[1], 1 + world_rank);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorBroadcastDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_rank = r_world.Rank();
    const int send_rank = 0;

    std::vector<double> send = {1.0*world_rank, 1.0 + world_rank};
    serial_communicator.Broadcast(send,send_rank);
    KRATOS_CHECK_EQUAL(send[0], world_rank);
    KRATOS_CHECK_EQUAL(send[1], 1 + world_rank);
}

// Scatter ////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScatterIntVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();
    const int send_rank = 0;

    std::vector<int> send_buffer(0);
    std::vector<int> recv_buffer = {-1, -1};

    if (world_rank == send_rank)
    {
        send_buffer.resize(2*world_size);
        for (int i = 0; i < 2*world_size; i++)
        {
            send_buffer[i] = 1;
        }
    }

    // two-buffer version
    serial_communicator.Scatter(send_buffer, recv_buffer, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], -1);
    }

    // return version
    std::vector<int> return_buffer = serial_communicator.Scatter(send_buffer, send_rank);
    if (world_rank == send_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
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

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScatterDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();
    const int send_rank = 0;

    std::vector<double> send_buffer(0);
    std::vector<double> recv_buffer{-1.0, -1.0};

    if (world_rank == send_rank)
    {
        send_buffer.resize(2*world_size);
        for (int i = 0; i < 2*world_size; i++)
        {
            send_buffer[i] = 2.0;
        }
    }

    // two-buffer version
    serial_communicator.Scatter(send_buffer, recv_buffer, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], -1.0);
    }

    // return version
    std::vector<double> return_buffer = serial_communicator.Scatter(send_buffer, send_rank);
    if (world_rank == send_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
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

// Scatterv ///////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScattervInt, KratosMPICoreFastSuite)
{
    /* send message for ints is {0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, ...} (max 6 values per rank)
     * read only first <rank> values of the message per rank (up to 5 values per rank)
     * message containing doubles is double of message containing ints
     */
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();
    const int send_rank = world_size-1;

    std::vector<int> send_buffer(0);
    std::vector<int> send_sizes(0);
    std::vector<int> send_offsets(0);

    auto make_message_size = [](int rank) { return rank < 5 ? rank : 5; };
    auto make_message_distance = [](int rank, int padding) {
        return rank < 5 ? ((rank-1)*rank)/2 + rank*padding : rank*(5+padding) - 15;
    };

    const int recv_size = make_message_size(world_rank);
    std::vector<int> recv_buffer(recv_size, -1);

    const int message_padding = 1;
    if (world_rank == send_rank)
    {
        const int send_size = make_message_distance(world_size, message_padding);
        send_buffer.resize(send_size);
        send_sizes.resize(world_size);
        send_offsets.resize(world_size);
        int counter = 0;
        for (int rank = 0; rank < world_size; rank++)
        {
            send_sizes[rank] = make_message_size(rank);
            send_offsets[rank] = make_message_distance(rank, message_padding);
            for (int i = 0; i < send_sizes[rank] + message_padding; i++, counter++)
            {
                send_buffer[counter] = rank;
            }
        }
    }

    // two-buffer version
    serial_communicator.Scatterv(send_buffer, send_sizes, send_offsets, recv_buffer, send_rank);

    for (int i = 0; i < recv_size; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], -1);
    }

    // return buffer version

    // pre-process: prepare input vector-of-vectors
    std::vector<std::vector<int>> scatterv_message;
    if (serial_communicator.Rank() == send_rank)
    {
        scatterv_message.resize(world_size);
        for (int rank = 0; rank < world_size; rank++)
        {
            // note: single-buffer version does not support padding, entire message is sent.
            scatterv_message[rank].resize(make_message_size(rank));
            for (int i = 0; i < send_sizes[rank]; i++)
            {
                scatterv_message[rank][i] = rank;
            }
        }
    }
    std::vector<int> result = serial_communicator.Scatterv(scatterv_message, send_rank);
    if (serial_communicator.Rank() == send_rank)
    {
        unsigned int expected_size = make_message_size(send_rank);
        KRATOS_CHECK_EQUAL(result.size(), expected_size);
        for (unsigned int i = 0; i < expected_size; i++)
        {
            KRATOS_CHECK_EQUAL(result[i], send_rank);
        }
    }
    else
    {
        KRATOS_CHECK_EQUAL(result.size(), 0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScattervDouble, KratosMPICoreFastSuite)
{
    /* send message for ints is {0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, ...} (max 6 values per rank)
     * read only first <rank> values of the message per rank (up to 5 values per rank)
     * message containing doubles is double of message containing ints
     */
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();
    const int send_rank = world_size-1;

    std::vector<double> send_buffer(0);
    std::vector<int> send_sizes(0);
    std::vector<int> send_offsets(0);

    auto make_message_size = [](int rank) { return rank < 5 ? rank : 5; };
    auto make_message_distance = [](int rank, int padding) {
        return rank < 5 ? ((rank-1)*rank)/2 + rank*padding : rank*(5+padding) - 15;
    };

    const int recv_size = make_message_size(world_rank);
    std::vector<double> recv_buffer(recv_size, -1.0);

    const int message_padding = 1;
    if (world_rank == send_rank)
    {
        const int send_size = make_message_distance(world_size, message_padding);
        send_buffer.resize(send_size);
        send_sizes.resize(world_size);
        send_offsets.resize(world_size);
        int counter = 0;
        for (int rank = 0; rank < world_size; rank++)
        {
            send_sizes[rank] = make_message_size(rank);
            send_offsets[rank] = make_message_distance(rank, message_padding);
            for (int i = 0; i < send_sizes[rank] + message_padding; i++, counter++)
            {
                send_buffer[counter] = 2.0*rank;
            }
        }
    }

    // two-buffer version
    serial_communicator.Scatterv(send_buffer, send_sizes, send_offsets, recv_buffer, send_rank);

    for (int i = 0; i < recv_size; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], -1.0);
    }

    // return buffer version

    // pre-process: prepare input vector-of-vectors
    std::vector<std::vector<double>> scatterv_message;
    if (serial_communicator.Rank() == send_rank)
    {
        scatterv_message.resize(world_size);
        for (int rank = 0; rank < world_size; rank++)
        {
            // note: single-buffer version does not support padding, entire message is sent.
            scatterv_message[rank].resize(make_message_size(rank));
            for (int i = 0; i < send_sizes[rank]; i++)
            {
                scatterv_message[rank][i] = 2.0*rank;
            }
        }
    }
    std::vector<double> result = serial_communicator.Scatterv(scatterv_message, send_rank);
    if (serial_communicator.Rank() == send_rank)
    {
        unsigned int expected_size = make_message_size(send_rank);
        KRATOS_CHECK_EQUAL(result.size(), expected_size);
        for (unsigned int i = 0; i < expected_size; i++)
        {
            KRATOS_CHECK_EQUAL(result[i], 2.0*send_rank);
        }
    }
    else
    {
        KRATOS_CHECK_EQUAL(result.size(), 0);
    }
}

// Gather /////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorGatherInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();
    const int recv_rank = 0;

    std::vector<int> send_buffer{world_rank, world_rank};
    std::vector<int> recv_buffer;

    if (world_rank == recv_rank)
    {
        recv_buffer = std::vector<int>(2*world_size, -1);
    }

    // two-buffer version
    serial_communicator.Gather(send_buffer, recv_buffer, recv_rank);

    if (world_rank == recv_rank)
    {
        for (int rank = 0; rank < world_size; rank++)
        {
            for (int j = 2*rank; j < 2*rank+2; j++)
            {
                KRATOS_CHECK_EQUAL(recv_buffer[j], -1);
            }
        }
    }

    // return buffer version
    std::vector<int> return_buffer = serial_communicator.Gather(send_buffer, recv_rank);
    if (serial_communicator.Rank() == recv_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
        for (unsigned int j = 0; j < return_buffer.size(); j++)
        {
            KRATOS_CHECK_EQUAL(return_buffer[j], send_buffer[j]);
        }
    }
    else
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), 0);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorGatherDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();
    const int recv_rank = 0;

    std::vector<double> send_buffer{2.0*world_rank, 2.0*world_rank};
    std::vector<double> recv_buffer;

    if (world_rank == recv_rank)
    {
        recv_buffer = std::vector<double>(2*world_size, -1.0);
    }

    // two-buffer version
    serial_communicator.Gather(send_buffer, recv_buffer, recv_rank);

    if (world_rank == recv_rank)
    {
        for (int rank = 0; rank < world_size; rank++)
        {
            for (int j = 2*rank; j < 2*rank+2; j++)
            {
                KRATOS_CHECK_EQUAL(recv_buffer[j], -1.0);
            }
        }
    }

    // return buffer version
    std::vector<double> return_buffer = serial_communicator.Gather(send_buffer, recv_rank);
    if (serial_communicator.Rank() == recv_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
        for (unsigned int j = 0; j < return_buffer.size(); j++)
        {
            KRATOS_CHECK_EQUAL(return_buffer[j], send_buffer[j]);
        }
    }
    else
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), 0);
    }
}

// Gatherv ////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorGathervInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();
    const int recv_rank = world_size-1;

    auto make_message_size = [](int rank) { return rank < 5 ? rank : 5; };
    auto make_message_distance = [](int rank, int padding) {
        return rank < 5 ? ((rank-1)*rank)/2 + rank*padding : rank*(5+padding) - 15;
    };

    const int send_size = make_message_size(world_rank);
    std::vector<int> send_buffer(send_size, world_rank);

    // two-buffer version
    const int message_padding = 1;
    const int recv_size = make_message_distance(world_size, message_padding);
    std::vector<int> recv_buffer(0);
    std::vector<int> recv_sizes(0);
    std::vector<int> recv_offsets(0);

    if (world_rank == recv_rank)
    {
        recv_buffer.resize(recv_size, -1);
        recv_sizes.resize(world_size);
        recv_offsets.resize(world_size);
        for (int rank = 0; rank < world_size; rank++)
        {
            recv_sizes[rank] = make_message_size(rank);
            recv_offsets[rank] = make_message_distance(rank, message_padding);
        }
    }

    serial_communicator.Gatherv(send_buffer, recv_buffer, recv_sizes, recv_offsets, recv_rank);

    if (world_rank == recv_rank)
    {
        for (int i = 0; i < recv_size; i++)
        {
            KRATOS_CHECK_EQUAL(recv_buffer[i], -1);
        }
    }
    else
    {
        KRATOS_CHECK_EQUAL(recv_buffer.size(), 0);
    }

    // return buffer version
    std::vector<std::vector<int>> return_buffer = serial_communicator.Gatherv(send_buffer, recv_rank);

    if (serial_communicator.Rank() == recv_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), 1);
        unsigned int expected_size = make_message_size(world_rank);
        KRATOS_CHECK_EQUAL(return_buffer[0].size(), expected_size);
        for (unsigned int i = 0; i < expected_size; i++)
        {
            KRATOS_CHECK_EQUAL(return_buffer[0][i], send_buffer[i]);
        }
        // no padding in return version
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorGathervDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();
    const int recv_rank = world_size-1;

    auto make_message_size = [](int rank) { return rank < 5 ? rank : 5; };
    auto make_message_distance = [](int rank, int padding) {
        return rank < 5 ? ((rank-1)*rank)/2 + rank*padding : rank*(5+padding) - 15;
    };

    const int send_size = make_message_size(world_rank);
    std::vector<double> send_buffer(send_size, 2.0*world_rank);

    // two-buffer version
    const int message_padding = 1;
    const int recv_size = make_message_distance(world_size, message_padding);
    std::vector<double> recv_buffer(0);
    std::vector<int> recv_sizes(0);
    std::vector<int> recv_offsets(0);

    if (world_rank == recv_rank)
    {
        recv_buffer.resize(recv_size, -1.0);
        recv_sizes.resize(world_size);
        recv_offsets.resize(world_size);
        for (int rank = 0; rank < world_size; rank++)
        {
            recv_sizes[rank] = make_message_size(rank);
            recv_offsets[rank] = make_message_distance(rank, message_padding);
        }
    }

    serial_communicator.Gatherv(send_buffer, recv_buffer, recv_sizes, recv_offsets, recv_rank);

    if (world_rank == recv_rank)
    {
        for (int i = 0; i < recv_size; i++)
        {
            KRATOS_CHECK_EQUAL(recv_buffer[i], -1.0);
        }
    }
    else
    {
        KRATOS_CHECK_EQUAL(recv_buffer.size(), 0);
    }

    // return buffer version
    std::vector<std::vector<double>> return_buffer = serial_communicator.Gatherv(send_buffer, recv_rank);

    if (serial_communicator.Rank() == recv_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), 1);
        unsigned int expected_size = make_message_size(world_rank);
        KRATOS_CHECK_EQUAL(return_buffer[0].size(), expected_size);
        for (unsigned int i = 0; i < expected_size; i++)
        {
            KRATOS_CHECK_EQUAL(return_buffer[0][i], send_buffer[i]);
        }
        // no padding in return version
    }
}

// AllGather //////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorAllGatherInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();

    std::vector<int> send_buffer{world_rank, world_rank};

    // two-buffer version
    std::vector<int> recv_buffer(2*world_size, -1);
    serial_communicator.AllGather(send_buffer, recv_buffer);

    for (int i = 0; i < 2*world_size; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], -1);
    }

    // return buffer version
    std::vector<int> return_buffer = serial_communicator.AllGather(send_buffer);

    KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
    for (unsigned int i = 0; i < return_buffer.size(); i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorAllGatherDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = ParallelEnvironment::GetDefaultDataCommunicator();
    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();

    std::vector<double> send_buffer{2.0*world_rank, 2.0*world_rank};

    // two-buffer version
    std::vector<double> recv_buffer(2*world_size, -1.0);
    serial_communicator.AllGather(send_buffer, recv_buffer);

    for (int i = 0; i < 2*world_size; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], -1.0);
    }

    // return buffer version
    std::vector<double> return_buffer = serial_communicator.AllGather(send_buffer);

    KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
    for (unsigned int i = 0; i < return_buffer.size(); i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }
}

// Error broadcasting methods /////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorErrorBroadcasting, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // The serial communicator does not throw,
    // since it does not know about "other ranks" to broadcast the error to.
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