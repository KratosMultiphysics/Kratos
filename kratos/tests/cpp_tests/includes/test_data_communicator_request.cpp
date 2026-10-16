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

// Project includes
#include "testing/testing.h"
#include "includes/data_communicator.h"
#include "includes/data_communicator_request.h"

namespace Kratos::Testing
{

namespace
{

/// Dummy pending operation, completing after a given number of calls to Test().
class DummyRequestImplementation final : public DataCommunicatorRequest::Implementation
{
public:
    DummyRequestImplementation(int& rNumberOfWaits, const int NumberOfTestsToComplete)
        : mrNumberOfWaits(rNumberOfWaits),
          mNumberOfTestsToComplete(NumberOfTestsToComplete)
    {
    }

    void Wait() override
    {
        ++mrNumberOfWaits;
    }

    bool Test() override
    {
        return --mNumberOfTestsToComplete <= 0;
    }

private:
    int& mrNumberOfWaits;
    int mNumberOfTestsToComplete;
};

}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorRequestDefaultIsCompleted, KratosCoreFastSuite)
{
    DataCommunicatorRequest request;
    KRATOS_EXPECT_TRUE(request.IsCompleted());
    KRATOS_EXPECT_TRUE(request.Test());
    request.Wait();
    KRATOS_EXPECT_TRUE(request.IsCompleted());
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorRequestWaitAndTest, KratosCoreFastSuite)
{
    int number_of_waits = 0;

    DataCommunicatorRequest wait_request(Kratos::make_unique<DummyRequestImplementation>(number_of_waits, 1));
    KRATOS_EXPECT_FALSE(wait_request.IsCompleted());
    wait_request.Wait();
    KRATOS_EXPECT_TRUE(wait_request.IsCompleted());
    wait_request.Wait(); // no-op on a completed request
    KRATOS_EXPECT_EQ(number_of_waits, 1);

    DataCommunicatorRequest test_request(Kratos::make_unique<DummyRequestImplementation>(number_of_waits, 3));
    KRATOS_EXPECT_FALSE(test_request.Test());
    KRATOS_EXPECT_FALSE(test_request.Test());
    KRATOS_EXPECT_TRUE(test_request.Test());
    KRATOS_EXPECT_TRUE(test_request.IsCompleted());
    KRATOS_EXPECT_EQ(number_of_waits, 1);

    // Moving transfers the pending operation
    DataCommunicatorRequest moved_from(Kratos::make_unique<DummyRequestImplementation>(number_of_waits, 1));
    DataCommunicatorRequest moved_to(std::move(moved_from));
    KRATOS_EXPECT_FALSE(moved_to.IsCompleted());

    std::vector<DataCommunicatorRequest> requests;
    requests.push_back(std::move(moved_to));
    requests.emplace_back(Kratos::make_unique<DummyRequestImplementation>(number_of_waits, 1));
    requests.emplace_back();
    DataCommunicatorRequest::WaitAll(requests);
    for (const auto& r_request : requests) {
        KRATOS_EXPECT_TRUE(r_request.IsCompleted());
    }
    KRATOS_EXPECT_EQ(number_of_waits, 3);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSerialPointToPointAsync, KratosCoreFastSuite)
{
    DataCommunicator serial_communicator;

    const int send_value = 3;
    const std::vector<double> send_values{1.0, 2.0};
    const std::string send_string("kratos");
    KRATOS_EXPECT_TRUE(serial_communicator.ISend(send_value, 0).IsCompleted());
    KRATOS_EXPECT_TRUE(serial_communicator.ISend(send_values, 0, 1).IsCompleted());
    KRATOS_EXPECT_TRUE(serial_communicator.ISend(send_string, 0, 2).IsCompleted());

    KRATOS_EXPECT_EXCEPTION_IS_THROWN(serial_communicator.ISend(send_value, 1),
        "Communication between different ranks is not possible with a serial DataCommunicator.");

    int recv_value = 0;
    std::vector<double> recv_values(2);
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(serial_communicator.IRecv(recv_value, 0),
        "Calling serial DataCommunicator::IRecv, which has no meaningful return.");
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(serial_communicator.IRecv(recv_values, 0),
        "Calling serial DataCommunicator::IRecv, which has no meaningful return.");
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSerialCollectiveAsync, KratosCoreFastSuite)
{
    DataCommunicator serial_communicator;

    KRATOS_EXPECT_TRUE(serial_communicator.IBarrier().IsCompleted());

    double broadcast_value = 2.0;
    std::vector<int> broadcast_values{1, 2, 3};
    KRATOS_EXPECT_TRUE(serial_communicator.IBroadcast(broadcast_value, 0).IsCompleted());
    KRATOS_EXPECT_TRUE(serial_communicator.IBroadcast(broadcast_values, 0).IsCompleted());
    KRATOS_EXPECT_EQ(broadcast_value, 2.0);
    KRATOS_EXPECT_EQ(broadcast_values, (std::vector<int>{1, 2, 3}));

    const double local_value = 1.5;
    double sum = 0.0, min = 0.0, max = 0.0;
    KRATOS_EXPECT_TRUE(serial_communicator.ISumAll(local_value, sum).IsCompleted());
    KRATOS_EXPECT_TRUE(serial_communicator.IMinAll(local_value, min).IsCompleted());
    KRATOS_EXPECT_TRUE(serial_communicator.IMaxAll(local_value, max).IsCompleted());
    KRATOS_EXPECT_EQ(sum, local_value);
    KRATOS_EXPECT_EQ(min, local_value);
    KRATOS_EXPECT_EQ(max, local_value);

    const std::vector<int> local_values{-1, 4};
    std::vector<int> sums(2), mins(2), maxs(2);
    KRATOS_EXPECT_TRUE(serial_communicator.ISumAll(local_values, sums).IsCompleted());
    KRATOS_EXPECT_TRUE(serial_communicator.IMinAll(local_values, mins).IsCompleted());
    KRATOS_EXPECT_TRUE(serial_communicator.IMaxAll(local_values, maxs).IsCompleted());
    KRATOS_EXPECT_EQ(sums, local_values);
    KRATOS_EXPECT_EQ(mins, local_values);
    KRATOS_EXPECT_EQ(maxs, local_values);
}

} // namespace Kratos::Testing
