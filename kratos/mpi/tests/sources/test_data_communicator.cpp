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

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorRankAndSize, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    KRATOS_CHECK_EQUAL(serial_communicator.Rank(), 0);
    KRATOS_CHECK_EQUAL(serial_communicator.Size(), 1);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorInquiryChecks, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    KRATOS_CHECK_EQUAL(serial_communicator.IsDefinedOnThisRank(), true);
    KRATOS_CHECK_EQUAL(serial_communicator.IsNullOnThisRank(), false);
    KRATOS_CHECK_EQUAL(serial_communicator.IsDistributed(), false);
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
            KRATOS_CHECK_EQUAL(output[i], local[i]);
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

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Sum(local, wrong_size_global, root),
        "Input error in call to DataCommunicator::Sum"
    );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSumDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
            KRATOS_CHECK_EQUAL(output[i], local[i]);
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

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Sum(local, wrong_size_global, root),
        "Input error in call to DataCommunicator::Sum"
    );
    #endif
}

// Min ////////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMinInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
            KRATOS_CHECK_EQUAL(output[i], local[i]);
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

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Min(local, wrong_size_global, root),
        "Input error in call to DataCommunicator::Min"
    );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMinDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
            KRATOS_CHECK_EQUAL(output[i], local[i]);
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

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Min(local, wrong_size_global, root),
        "Input error in call to DataCommunicator::Min"
    );
    #endif
}

// Max ////////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommuniactorMaxInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
            KRATOS_CHECK_EQUAL(output[i], local[i]);
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

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Max(local, wrong_size_global, root),
        "Input error in call to DataCommunicator::Max"
    );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorMaxDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
            KRATOS_CHECK_EQUAL(output[i], local[i]);
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

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Max(local, wrong_size_global, root),
        "Input error in call to DataCommunicator::Max"
    );
    #endif
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
        KRATOS_CHECK_EQUAL(output[i], local[i]);
    }

    // return buffer version
    std::vector<int> returned_result = serial_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.SumAll(local, wrong_size_global),
        "Input error in call to DataCommunicator::SumAll"
    );
    #endif
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
        KRATOS_CHECK_EQUAL(output[i], local[i]);
    }

    // return buffer version
    std::vector<double> returned_result = serial_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.SumAll(local, wrong_size_global),
        "Input error in call to DataCommunicator::SumAll"
    );
    #endif
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
        KRATOS_CHECK_EQUAL(output[i], local[i]);
    }

    // return buffer version
    std::vector<int> returned_result = serial_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.MinAll(local, wrong_size_global),
        "Input error in call to DataCommunicator::MinAll"
    );
    #endif
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
        KRATOS_CHECK_EQUAL(output[i], local[i]);
    }

    // return buffer version
    std::vector<double> returned_result = serial_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.MinAll(local, wrong_size_global),
        "Input error in call to DataCommunicator::MinAll"
    );
    #endif
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
        KRATOS_CHECK_EQUAL(output[i], local[i]);
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
        KRATOS_CHECK_EQUAL(output[i], local[i]);
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
        KRATOS_CHECK_EQUAL(output[i], local[i]);
    }

    // return buffer version
    std::vector<int> returned_result = serial_communicator.ScanSum(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.ScanSum(local, wrong_size_global),
        "Input error in call to DataCommunicator::ScanSum"
    );
    #endif
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
        KRATOS_CHECK_EQUAL(output[i], local[i]);
    }

    // return buffer version
    std::vector<double> returned_result = serial_communicator.ScanSum(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], local[i]);
    }

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_global{-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.ScanSum(local, wrong_size_global),
        "Input error in call to DataCommunicator::ScanSum"
    );
    #endif
}

// SendRecv ///////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSendRecvInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();

    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();

    // Serial communication is only well defined for the trivial case.
    int send_rank = 0;
    int recv_rank = 0;

    std::vector<int> send_buffer = {world_rank, world_rank};
    std::vector<int> recv_buffer = {-1, -1};

    // two-buffer version
    serial_communicator.SendRecv(send_buffer, send_rank, 0, recv_buffer, recv_rank, 0);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }

    // return version
    std::vector<int> return_buffer = serial_communicator.SendRecv(send_buffer, send_rank, recv_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer[i]);
    }

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_recv = {-1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.SendRecv(send_buffer, send_rank, 0, wrong_size_recv, recv_rank, 0),
        "Input error in call to DataCommunicator::SendRecv"
    );
    #endif

    // remote calls are not supported
    if (world_size > 2) {
        send_rank = world_rank + 1 == world_size ? 0 : world_rank + 1;
        recv_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.SendRecv(send_buffer, send_rank, 0, recv_buffer, recv_rank, 0),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.SendRecv(send_buffer, send_rank, recv_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSendRecvDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();

    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();

    // Serial communication is only well defined for the trivial case.
    int send_rank = 0;
    int recv_rank = 0;

    std::vector<double> send_buffer = {2.0*world_rank, 2.0*world_rank};
    std::vector<double> recv_buffer = {-1.0, -1.0};

    // two-buffer version
    serial_communicator.SendRecv(send_buffer, send_rank, 0, recv_buffer, recv_rank, 0);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }

    // return version
    std::vector<double> return_buffer = serial_communicator.SendRecv(send_buffer, send_rank, recv_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer[i]);
    }

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_recv = {-1.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.SendRecv(send_buffer, send_rank, 0, wrong_size_recv, recv_rank, 0),
        "Input error in call to DataCommunicator::SendRecv"
    );
    #endif

    // remote calls are not supported
    if (world_size > 2) {
        send_rank = world_rank + 1 == world_size ? 0 : world_rank + 1;
        recv_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.SendRecv(send_buffer, send_rank, 0, recv_buffer, recv_rank, 0),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.SendRecv(send_buffer, send_rank, recv_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorSendRecvString, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();

    const int world_size = r_world.Size();
    const int world_rank = r_world.Rank();

    // Serial communication is only well defined for the trivial case.
    int send_rank = 0;
    int recv_rank = 0;

    std::string send_buffer("Hello world!");
    // The output is only needed to be of the same size, I initialize it to check it later.
    std::string recv_buffer("************");

    // two-buffer version
    serial_communicator.SendRecv(send_buffer, send_rank, 0, recv_buffer, recv_rank, 0);
    KRATOS_CHECK_C_STRING_EQUAL(recv_buffer.c_str(), "Hello world!");

    // return version
    std::string return_buffer = serial_communicator.SendRecv(send_buffer, send_rank, recv_rank);
    KRATOS_CHECK_C_STRING_EQUAL(recv_buffer.c_str(), "Hello world!");

    #ifdef KRATOS_DEBUG
    std::string wrong_size_recv("*");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.SendRecv(send_buffer, send_rank, 0, wrong_size_recv, recv_rank, 0),
        "Input error in call to DataCommunicator::SendRecv"
    );
    #endif

    // remote calls are not supported
    if (world_size > 2) {
        send_rank = world_rank + 1 == world_size ? 0 : world_rank + 1;
        recv_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.SendRecv(send_buffer, send_rank, 0, recv_buffer, recv_rank, 0),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.SendRecv(send_buffer, send_rank, recv_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }
}

// Broadcast //////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorBroadcastInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
    const int world_rank = r_world.Rank();
    const int send_rank = 0;

    int send = 1 + world_rank;
    serial_communicator.Broadcast(send,send_rank);
    KRATOS_CHECK_EQUAL(send, 1 + world_rank);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorBroadcastDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
    const int world_rank = r_world.Rank();
    const int send_rank = 0;

    double send = 1.0 + world_rank;
    serial_communicator.Broadcast(send,send_rank);
    KRATOS_CHECK_EQUAL(send, 1.0 + world_rank);
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorBroadcastIntVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
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

    // the serial version of scatter only works for the trivial case (from 0 to 0)
    std::vector<int> send_buffer = {1, 1};
    std::vector<int> recv_buffer = {-1, -1};
    int send_rank = 0;

    // two-buffer version
    serial_communicator.Scatter(send_buffer, recv_buffer, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }

    // return version
    std::vector<int> return_buffer = serial_communicator.Scatter(send_buffer, send_rank);
    KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer[i]);
    }

    // remote calls are not supported
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
    const int world_size = r_world.Size();

    if (world_size > 1) {
        send_rank = world_size - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.Scatter(send_buffer, recv_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.Scatter(send_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_recv = {-1, -1, -1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Scatter(send_buffer, wrong_size_recv, send_rank),
        "Input error in call to DataCommunicator::Scatter"
    );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScatterDoubleVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // the serial version of scatter only works for the trivial case (from 0 to 0)
    std::vector<double> send_buffer = {2.0, 2.0};
    std::vector<double> recv_buffer = {-1.0, -1.0};
    int send_rank = 0;

    // two-buffer version
    serial_communicator.Scatter(send_buffer, recv_buffer, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }

    // return version
    std::vector<double> return_buffer = serial_communicator.Scatter(send_buffer, send_rank);
    KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer[i]);
    }

    // remote calls are not supported
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
    const int world_size = r_world.Size();

    if (world_size > 1) {
        send_rank = world_size - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.Scatter(send_buffer, recv_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.Scatter(send_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_recv = {-1.0, -1.0, -1.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Scatter(send_buffer, wrong_size_recv, send_rank),
        "Input error in call to DataCommunicator::Scatter"
    );
    #endif
}

// Scatterv ///////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScattervInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // the serial version of scatterv only works for the trivial case (from 0 to 0)
    std::vector<int> send_buffer_single = {1, 1};

    std::vector<std::vector<int>> send_buffer_multiple = {send_buffer_single};
    std::vector<int> send_offsets = {0};
    std::vector<int> send_counts = {2};

    std::vector<int> recv_buffer = {-1, -1};
    int send_rank = 0;

    // two-buffer version
    serial_communicator.Scatterv(send_buffer_single, send_counts, send_offsets, recv_buffer, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer_single[i]);
    }

    // return version
    std::vector<int> return_buffer = serial_communicator.Scatterv(send_buffer_multiple, send_rank);
    KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer_single.size());
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer_single[i]);
    }

    // remote calls are not supported
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
    const int world_size = r_world.Size();

    if (world_size > 1) {
        send_rank = world_size - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.Scatterv(send_buffer_single, send_counts, send_offsets, recv_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.Scatterv(send_buffer_multiple, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_recv = {-1, -1, -1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Scatterv(send_buffer_single, send_counts, send_offsets, wrong_size_recv, send_rank),
        "Input error in call to DataCommunicator::Scatterv"
    );
    std::vector<int> wrong_counts = {2, 3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Scatterv(send_buffer_single, wrong_counts, send_offsets, recv_buffer, send_rank),
        "Input error in call to DataCommunicator::Scatterv"
    );
    std::vector<int> wrong_offsets = {0, 1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Scatterv(send_buffer_single, send_counts, wrong_offsets, recv_buffer, send_rank),
        "Input error in call to DataCommunicator::Scatterv"
    );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorScattervDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // the serial version of scatterv only works for the trivial case (from 0 to 0)
    std::vector<double> send_buffer_single = {2.0, 2.0};

    std::vector<std::vector<double>> send_buffer_multiple = {send_buffer_single};
    std::vector<int> send_offsets = {0};
    std::vector<int> send_counts = {2};

    std::vector<double> recv_buffer = {-1.0, -1.0};
    int send_rank = 0;

    // two-buffer version
    serial_communicator.Scatterv(send_buffer_single, send_counts, send_offsets, recv_buffer, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer_single[i]);
    }

    // return version
    std::vector<double> return_buffer = serial_communicator.Scatterv(send_buffer_multiple, send_rank);
    KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer_single.size());
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer_single[i]);
    }

    // remote calls are not supported
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
    const int world_size = r_world.Size();

    if (world_size > 1) {
        send_rank = world_size - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.Scatterv(send_buffer_single, send_counts, send_offsets, recv_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.Scatterv(send_buffer_multiple, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_recv = {-1.0, -1.0, -1.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Scatterv(send_buffer_single, send_counts, send_offsets, wrong_size_recv, send_rank),
        "Input error in call to DataCommunicator::Scatterv"
    );
    std::vector<int> wrong_counts = {2, 3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Scatterv(send_buffer_single, wrong_counts, send_offsets, recv_buffer, send_rank),
        "Input error in call to DataCommunicator::Scatterv"
    );
    std::vector<int> wrong_offsets = {0, 1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Scatterv(send_buffer_single, send_counts, wrong_offsets, recv_buffer, send_rank),
        "Input error in call to DataCommunicator::Scatterv"
    );
    #endif
}

// Gather /////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorGatherInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // the serial version of scatter only works for the trivial case (from 0 to 0)
    std::vector<int> send_buffer = {1, 1};
    std::vector<int> recv_buffer = {-1, -1};
    int send_rank = 0;

    // two-buffer version
    serial_communicator.Gather(send_buffer, recv_buffer, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }

    // return version
    std::vector<int> return_buffer = serial_communicator.Gather(send_buffer, send_rank);
    KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer[i]);
    }

    // remote calls are not supported
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
    const int world_size = r_world.Size();

    if (world_size > 1) {
        send_rank = world_size - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.Gather(send_buffer, recv_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.Gather(send_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_recv = {-1, -1, -1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Gather(send_buffer, wrong_size_recv, send_rank),
        "Input error in call to DataCommunicator::Gather"
    );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorGatherDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // the serial version of scatter only works for the trivial case (from 0 to 0)
    std::vector<double> send_buffer = {2.0, 2.0};
    std::vector<double> recv_buffer = {-1.0, -1.0};
    int send_rank = 0;

    // two-buffer version
    serial_communicator.Gather(send_buffer, recv_buffer, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }

    // return version
    std::vector<double> return_buffer = serial_communicator.Gather(send_buffer, send_rank);
    KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer[i]);
    }

    // remote calls are not supported
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
    const int world_size = r_world.Size();

    if (world_size > 1) {
        send_rank = world_size - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.Gather(send_buffer, recv_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.Gather(send_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_recv = {-1.0, -1.0, -1.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Gather(send_buffer, wrong_size_recv, send_rank),
        "Input error in call to DataCommunicator::Gather"
    );
    #endif
}

// Gatherv ////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorGathervInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // the serial version of scatterv only works for the trivial case (from 0 to 0)
    std::vector<int> send_buffer = {1, 1};

    std::vector<int> recv_offsets = {0};
    std::vector<int> recv_counts = {2};

    std::vector<int> recv_buffer = {-1, -1};
    int send_rank = 0;

    // two-buffer version
    serial_communicator.Gatherv(send_buffer, recv_buffer, recv_counts, recv_offsets, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }

    // return version
    std::vector<std::vector<int>> return_buffer = serial_communicator.Gatherv(send_buffer, send_rank);
    KRATOS_CHECK_EQUAL(return_buffer.size(), 1);
    KRATOS_CHECK_EQUAL(return_buffer[0].size(), send_buffer.size());
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[0][i], send_buffer[i]);
    }

    // remote calls are not supported
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
    const int world_size = r_world.Size();

    if (world_size > 1) {
        send_rank = world_size - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.Gatherv(send_buffer, recv_buffer, recv_counts, recv_offsets, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.Gatherv(send_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_recv = {-1, -1, -1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Gatherv(send_buffer, wrong_size_recv, recv_counts, recv_offsets, send_rank),
        "Input error in call to DataCommunicator::Gatherv"
    );
    std::vector<int> wrong_counts = {2, 3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Gatherv(send_buffer, recv_buffer, wrong_counts, recv_offsets, send_rank),
        "Input error in call to DataCommunicator::Gatherv"
    );
    std::vector<int> wrong_offsets = {0, 1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Gatherv(send_buffer, recv_buffer, recv_counts, wrong_offsets, send_rank),
        "Input error in call to DataCommunicator::Gatherv"
    );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorGathervDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // the serial version of scatterv only works for the trivial case (from 0 to 0)
    std::vector<double> send_buffer = {2.0, 2.0};

    std::vector<int> recv_offsets = {0};
    std::vector<int> recv_counts = {2};

    std::vector<double> recv_buffer = {-1.0, -1.0};
    int send_rank = 0;

    // two-buffer version
    serial_communicator.Gatherv(send_buffer, recv_buffer, recv_counts, recv_offsets, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }

    // return version
    std::vector<std::vector<double>> return_buffer = serial_communicator.Gatherv(send_buffer, send_rank);
    KRATOS_CHECK_EQUAL(return_buffer.size(), 1);
    KRATOS_CHECK_EQUAL(return_buffer[0].size(), send_buffer.size());
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[0][i], send_buffer[i]);
    }

    // remote calls are not supported
    const DataCommunicator& r_world = DataCommunicator::GetDefault();
    const int world_size = r_world.Size();

    if (world_size > 1) {
        send_rank = world_size - 1;

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            serial_communicator.Gatherv(send_buffer, recv_buffer, recv_counts, recv_offsets, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );

        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            return_buffer = serial_communicator.Gatherv(send_buffer, send_rank),
            "Communication between different ranks is not possible with a serial DataCommunicator."
        );
    }

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_recv = {-1.0, -1.0, -1.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Gatherv(send_buffer, wrong_size_recv, recv_counts, recv_offsets, send_rank),
        "Input error in call to DataCommunicator::Gatherv"
    );
    std::vector<int> wrong_counts = {2, 3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Gatherv(send_buffer, recv_buffer, wrong_counts, recv_offsets, send_rank),
        "Input error in call to DataCommunicator::Gatherv"
    );
    std::vector<int> wrong_offsets = {0, 1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.Gatherv(send_buffer, recv_buffer, recv_counts, wrong_offsets, send_rank),
        "Input error in call to DataCommunicator::Gatherv"
    );
    #endif
}

// AllGather //////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorAllGatherInt, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // the serial version of scatter only works for the trivial case (from 0 to 0)
    const std::vector<int> send_buffer = {1, 1};
    std::vector<int> recv_buffer = {-1, -1};

    // two-buffer version
    serial_communicator.AllGather(send_buffer, recv_buffer);
    for (unsigned int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }

    // return buffer version
    std::vector<int> return_buffer = serial_communicator.AllGather(send_buffer);
    KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
    for (unsigned int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer[i]);
    }

    #ifdef KRATOS_DEBUG
    std::vector<int> wrong_size_recv = {-1, -1, -1};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.AllGather(send_buffer, wrong_size_recv),
        "Input error in call to DataCommunicator::AllGather"
    );
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(DataCommunicatorAllGatherDouble, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    // the serial version of scatter only works for the trivial case (from 0 to 0)
    const std::vector<double> send_buffer = {2.0, 2.0};
    std::vector<double> recv_buffer = {-1.0, -1.0};

    // two-buffer version
    serial_communicator.AllGather(send_buffer, recv_buffer);
    for (unsigned int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], send_buffer[i]);
    }

    // return buffer version
    std::vector<double> return_buffer = serial_communicator.AllGather(send_buffer);
    KRATOS_CHECK_EQUAL(return_buffer.size(), send_buffer.size());
    for (unsigned int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], send_buffer[i]);
    }

    #ifdef KRATOS_DEBUG
    std::vector<double> wrong_size_recv = {-1.0, -1.0, -1.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        serial_communicator.AllGather(send_buffer, wrong_size_recv),
        "Input error in call to DataCommunicator::AllGather"
    );
    #endif
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
