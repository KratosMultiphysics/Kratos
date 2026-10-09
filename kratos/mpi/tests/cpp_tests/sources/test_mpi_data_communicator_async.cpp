//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <string>
#include <vector>

// External includes
#include "mpi.h"

// Project includes
#include "includes/data_communicator_request.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/testing/mpi_testing.h"

namespace Kratos::Testing {

// Non-blocking communication is tested in a ring: each rank sends to the next one and receives from the previous one.
// Contrary to the blocking tests, this also works for a single rank (the message is sent to itself).

namespace {

int NextRank(const DataCommunicator& rComm)
{
    return (rComm.Rank() + 1) % rComm.Size();
}

int PreviousRank(const DataCommunicator& rComm)
{
    return (rComm.Rank() + rComm.Size() - 1) % rComm.Size();
}

template<class T>
void MPIDataCommunicatorISendIRecvScalarTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int send_rank = NextRank(mpi_world_communicator);
    const int recv_rank = PreviousRank(mpi_world_communicator);
    const T expected_recv = static_cast<T>(recv_rank + 1);

    // value version
    const T send_value = static_cast<T>(mpi_world_communicator.Rank() + 1);
    T recv_value = 0;
    auto recv_request = mpi_world_communicator.IRecv(recv_value, recv_rank, 0);
    auto send_request = mpi_world_communicator.ISend(send_value, send_rank, 0);
    KRATOS_EXPECT_FALSE(recv_request.IsCompleted());
    recv_request.Wait();
    send_request.Wait();
    KRATOS_EXPECT_TRUE(recv_request.IsCompleted());
    KRATOS_EXPECT_TRUE(send_request.IsCompleted());
    KRATOS_EXPECT_EQ(recv_value, expected_recv);

    // vector version
    const std::vector<T> send_values{send_value, static_cast<T>(2 * send_value)};
    std::vector<T> recv_values(2, 0);
    std::vector<DataCommunicatorRequest> requests;
    requests.push_back(mpi_world_communicator.IRecv(recv_values, recv_rank, 1));
    requests.push_back(mpi_world_communicator.ISend(send_values, send_rank, 1));
    DataCommunicatorRequest::WaitAll(requests);
    KRATOS_EXPECT_EQ(recv_values[0], expected_recv);
    KRATOS_EXPECT_EQ(recv_values[1], static_cast<T>(2 * expected_recv));
}

}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorISendIRecvInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorISendIRecvScalarTypeTest<int>();
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorISendIRecvUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorISendIRecvScalarTypeTest<unsigned int>();
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorISendIRecvLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorISendIRecvScalarTypeTest<long unsigned int>();
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorISendIRecvDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorISendIRecvScalarTypeTest<double>();
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorISendIRecvArray1d, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int send_rank = NextRank(mpi_world_communicator);
    const int recv_rank = PreviousRank(mpi_world_communicator);
    const double send_factor = mpi_world_communicator.Rank() + 1.0;
    const double recv_factor = recv_rank + 1.0;

    const array_1d<double, 3> send_value{2.0 * send_factor, 3.0 * send_factor, 5.0 * send_factor};
    const array_1d<double, 3> expected_recv{2.0 * recv_factor, 3.0 * recv_factor, 5.0 * recv_factor};
    const std::vector<array_1d<double, 3>> send_values{send_value, 2.0 * send_value};

    array_1d<double, 3> recv_value(3, -1.0);
    std::vector<array_1d<double, 3>> recv_values(2, array_1d<double, 3>(3, -1.0));

    std::vector<DataCommunicatorRequest> requests;
    requests.push_back(mpi_world_communicator.IRecv(recv_value, recv_rank, 0));
    requests.push_back(mpi_world_communicator.IRecv(recv_values, recv_rank, 1));
    requests.push_back(mpi_world_communicator.ISend(send_value, send_rank, 0));
    requests.push_back(mpi_world_communicator.ISend(send_values, send_rank, 1));
    DataCommunicatorRequest::WaitAll(requests);

    KRATOS_EXPECT_VECTOR_EQ(recv_value, expected_recv);
    KRATOS_EXPECT_VECTOR_EQ(recv_values[0], expected_recv);
    KRATOS_EXPECT_VECTOR_EQ(recv_values[1], (2.0 * expected_recv));
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorISendIRecvVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int send_rank = NextRank(mpi_world_communicator);
    const int recv_rank = PreviousRank(mpi_world_communicator);

    Vector send_value(4), expected_recv(4);
    for (std::size_t i = 0; i < 4; ++i) {
        send_value[i] = (i + 1.0) * (mpi_world_communicator.Rank() + 1.0);
        expected_recv[i] = (i + 1.0) * (recv_rank + 1.0);
    }
    // std::vector<Vector> is not contiguous in memory: this checks the data is copied back on completion.
    const std::vector<Vector> send_values{send_value, 2.0 * send_value};

    Vector recv_value(4, -1.0);
    std::vector<Vector> recv_values(2, Vector(4, -1.0));

    std::vector<DataCommunicatorRequest> requests;
    requests.push_back(mpi_world_communicator.IRecv(recv_value, recv_rank, 0));
    requests.push_back(mpi_world_communicator.IRecv(recv_values, recv_rank, 1));
    requests.push_back(mpi_world_communicator.ISend(send_value, send_rank, 0));
    requests.push_back(mpi_world_communicator.ISend(send_values, send_rank, 1));
    DataCommunicatorRequest::WaitAll(requests);

    KRATOS_EXPECT_VECTOR_EQ(recv_value, expected_recv);
    KRATOS_EXPECT_VECTOR_EQ(recv_values[0], expected_recv);
    KRATOS_EXPECT_VECTOR_EQ(recv_values[1], (2.0 * expected_recv));
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorISendIRecvMatrix, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int send_rank = NextRank(mpi_world_communicator);
    const int recv_rank = PreviousRank(mpi_world_communicator);

    Matrix send_value(2, 3), expected_recv(2, 3);
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            send_value(i, j) = (3.0 * i + j) * (mpi_world_communicator.Rank() + 1.0);
            expected_recv(i, j) = (3.0 * i + j) * (recv_rank + 1.0);
        }
    }
    const std::vector<Matrix> send_values{send_value, 2.0 * send_value};

    Matrix recv_value(2, 3, -1.0);
    std::vector<Matrix> recv_values(2, Matrix(2, 3, -1.0));

    std::vector<DataCommunicatorRequest> requests;
    requests.push_back(mpi_world_communicator.IRecv(recv_value, recv_rank, 0));
    requests.push_back(mpi_world_communicator.IRecv(recv_values, recv_rank, 1));
    requests.push_back(mpi_world_communicator.ISend(send_value, send_rank, 0));
    requests.push_back(mpi_world_communicator.ISend(send_values, send_rank, 1));
    DataCommunicatorRequest::WaitAll(requests);

    KRATOS_EXPECT_MATRIX_EQ(recv_value, expected_recv);
    KRATOS_EXPECT_MATRIX_EQ(recv_values[0], expected_recv);
    KRATOS_EXPECT_MATRIX_EQ(recv_values[1], (2.0 * expected_recv));
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorISendIRecvString, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int send_rank = NextRank(mpi_world_communicator);
    const int recv_rank = PreviousRank(mpi_world_communicator);

    // All messages have the same length, so that the receiving buffer can be sized in advance.
    const std::string send_string = "rank_" + std::to_string(mpi_world_communicator.Rank() % 10);
    const std::string expected_recv = "rank_" + std::to_string(recv_rank % 10);
    std::string recv_string(expected_recv.size(), ' ');

    auto recv_request = mpi_world_communicator.IRecv(recv_string, recv_rank, 0);
    auto send_request = mpi_world_communicator.ISend(send_string, send_rank, 0);
    recv_request.Wait();
    send_request.Wait();

    KRATOS_EXPECT_EQ(recv_string, expected_recv);
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorIRecvTest, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int send_rank = NextRank(mpi_world_communicator);
    const int recv_rank = PreviousRank(mpi_world_communicator);

    // non-contiguous type, to check that Test also copies the received data back
    const std::vector<Vector> send_values{Vector(3, mpi_world_communicator.Rank() + 1.0)};
    std::vector<Vector> recv_values(1, Vector(3, -1.0));

    auto recv_request = mpi_world_communicator.IRecv(recv_values, recv_rank, 0);
    auto send_request = mpi_world_communicator.ISend(send_values, send_rank, 0);
    while (!recv_request.Test()) {}
    while (!send_request.Test()) {}

    KRATOS_EXPECT_TRUE(recv_request.IsCompleted());
    KRATOS_EXPECT_TRUE(send_request.IsCompleted());
    KRATOS_EXPECT_VECTOR_EQ(recv_values[0], Vector(3, recv_rank + 1.0));
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorRequestDestructorCompletes, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int send_rank = NextRank(mpi_world_communicator);
    const int recv_rank = PreviousRank(mpi_world_communicator);

    const std::vector<Vector> send_values{Vector(2, mpi_world_communicator.Rank() + 1.0)};
    std::vector<Vector> recv_values(1, Vector(2, -1.0));
    {
        // requests going out of scope wait for the pending operations
        auto recv_request = mpi_world_communicator.IRecv(recv_values, recv_rank, 0);
        auto send_request = mpi_world_communicator.ISend(send_values, send_rank, 0);
    }

    KRATOS_EXPECT_VECTOR_EQ(recv_values[0], Vector(2, recv_rank + 1.0));
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorIBarrier, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    auto request = mpi_world_communicator.IBarrier();
    request.Wait();
    KRATOS_EXPECT_TRUE(request.IsCompleted());
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorIBroadcast, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int source_rank = mpi_world_communicator.Size() - 1;
    const bool is_source = world_rank == source_rank;

    int int_value = is_source ? 7 : -1;
    std::vector<double> double_values = is_source ? std::vector<double>{1.5, 2.5} : std::vector<double>{-1.0, -1.0};
    Vector vector_value = is_source ? Vector(3, 4.0) : Vector(3, -1.0);
    std::vector<Vector> vector_values(2, is_source ? Vector(2, 5.0) : Vector(2, -1.0));
    Matrix matrix_value = is_source ? Matrix(2, 2, 6.0) : Matrix(2, 2, -1.0);
    std::string string_value = is_source ? "kratos" : "      ";

    std::vector<DataCommunicatorRequest> requests;
    requests.push_back(mpi_world_communicator.IBroadcast(int_value, source_rank));
    requests.push_back(mpi_world_communicator.IBroadcast(double_values, source_rank));
    requests.push_back(mpi_world_communicator.IBroadcast(vector_value, source_rank));
    requests.push_back(mpi_world_communicator.IBroadcast(vector_values, source_rank));
    requests.push_back(mpi_world_communicator.IBroadcast(matrix_value, source_rank));
    requests.push_back(mpi_world_communicator.IBroadcast(string_value, source_rank));
    DataCommunicatorRequest::WaitAll(requests);

    KRATOS_EXPECT_EQ(int_value, 7);
    KRATOS_EXPECT_EQ(double_values[0], 1.5);
    KRATOS_EXPECT_EQ(double_values[1], 2.5);
    KRATOS_EXPECT_VECTOR_EQ(vector_value, Vector(3, 4.0));
    KRATOS_EXPECT_VECTOR_EQ(vector_values[0], Vector(2, 5.0));
    KRATOS_EXPECT_VECTOR_EQ(vector_values[1], Vector(2, 5.0));
    KRATOS_EXPECT_MATRIX_EQ(matrix_value, Matrix(2, 2, 6.0));
    KRATOS_EXPECT_EQ(string_value, "kratos");
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorIAllReduceInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    int sum = 0, min = 0, max = 0;
    const std::vector<int> local_values{world_rank, -world_rank};
    std::vector<int> sums(2), mins(2), maxs(2);

    std::vector<DataCommunicatorRequest> requests;
    requests.push_back(mpi_world_communicator.ISumAll(world_rank, sum));
    requests.push_back(mpi_world_communicator.IMinAll(world_rank, min));
    requests.push_back(mpi_world_communicator.IMaxAll(world_rank, max));
    requests.push_back(mpi_world_communicator.ISumAll(local_values, sums));
    requests.push_back(mpi_world_communicator.IMinAll(local_values, mins));
    requests.push_back(mpi_world_communicator.IMaxAll(local_values, maxs));
    DataCommunicatorRequest::WaitAll(requests);

    const int expected_sum = world_size * (world_size - 1) / 2;
    KRATOS_EXPECT_EQ(sum, expected_sum);
    KRATOS_EXPECT_EQ(min, 0);
    KRATOS_EXPECT_EQ(max, world_size - 1);
    KRATOS_EXPECT_EQ(sums, (std::vector<int>{expected_sum, -expected_sum}));
    KRATOS_EXPECT_EQ(mins, (std::vector<int>{0, 1 - world_size}));
    KRATOS_EXPECT_EQ(maxs, (std::vector<int>{world_size - 1, 0}));
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorIAllReduceDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const double world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    double sum = 0.0, min = 0.0, max = 0.0;
    auto sum_request = mpi_world_communicator.ISumAll(0.5 * world_rank, sum);
    auto min_request = mpi_world_communicator.IMinAll(0.5 * world_rank, min);
    auto max_request = mpi_world_communicator.IMaxAll(0.5 * world_rank, max);
    sum_request.Wait();
    min_request.Wait();
    max_request.Wait();

    KRATOS_EXPECT_NEAR(sum, 0.25 * world_size * (world_size - 1), 1e-12);
    KRATOS_EXPECT_NEAR(min, 0.0, 1e-12);
    KRATOS_EXPECT_NEAR(max, 0.5 * (world_size - 1), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorIAllReduceArray1d, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const double world_rank = mpi_world_communicator.Rank();
    const double world_size = mpi_world_communicator.Size();

    const array_1d<double, 3> local_value{world_rank, -world_rank, 1.0};
    array_1d<double, 3> sum, min, max;

    std::vector<DataCommunicatorRequest> requests;
    requests.push_back(mpi_world_communicator.ISumAll(local_value, sum));
    requests.push_back(mpi_world_communicator.IMinAll(local_value, min));
    requests.push_back(mpi_world_communicator.IMaxAll(local_value, max));
    DataCommunicatorRequest::WaitAll(requests);

    const double expected_sum = 0.5 * world_size * (world_size - 1);
    KRATOS_EXPECT_VECTOR_NEAR(sum, (array_1d<double, 3>{expected_sum, -expected_sum, world_size}), 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(min, (array_1d<double, 3>{0.0, 1.0 - world_size, 1.0}), 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(max, (array_1d<double, 3>{world_size - 1.0, 0.0, 1.0}), 1e-12);
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorIAllReduceVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const double world_rank = mpi_world_communicator.Rank();
    const double world_size = mpi_world_communicator.Size();

    const std::vector<Vector> local_values{Vector(2, world_rank), Vector(2, -world_rank)};
    // Output buffers must have the right shape in advance
    std::vector<Vector> sums(2, Vector(2)), mins(2, Vector(2)), maxs(2, Vector(2));

    std::vector<DataCommunicatorRequest> requests;
    requests.push_back(mpi_world_communicator.ISumAll(local_values, sums));
    requests.push_back(mpi_world_communicator.IMinAll(local_values, mins));
    requests.push_back(mpi_world_communicator.IMaxAll(local_values, maxs));
    DataCommunicatorRequest::WaitAll(requests);

    const double expected_sum = 0.5 * world_size * (world_size - 1);
    KRATOS_EXPECT_VECTOR_NEAR(sums[0], Vector(2, expected_sum), 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(sums[1], Vector(2, -expected_sum), 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(mins[0], Vector(2, 0.0), 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(mins[1], Vector(2, 1.0 - world_size), 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(maxs[0], Vector(2, world_size - 1.0), 1e-12);
    KRATOS_EXPECT_VECTOR_NEAR(maxs[1], Vector(2, 0.0), 1e-12);
}

} // namespace Kratos::Testing
