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
#include "mpi/includes/mpi_data_communicator.h"

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorRankAndSize, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;

    KRATOS_CHECK_EQUAL(serial_communicator.Rank(), 0);
    KRATOS_CHECK_EQUAL(serial_communicator.Size(), 1);

    MPIDataCommunicator mpi_self_communicator(MPI_COMM_SELF);

    KRATOS_CHECK_EQUAL(mpi_self_communicator.Rank(), 0);
    KRATOS_CHECK_EQUAL(mpi_self_communicator.Size(), 1);

    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);

    KRATOS_CHECK_EQUAL(mpi_world_communicator.Rank(), world_rank);
    KRATOS_CHECK_EQUAL(mpi_world_communicator.Size(), world_size);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPICommRetrieval, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_self_communicator(MPI_COMM_SELF);
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    KRATOS_CHECK_EQUAL(MPIDataCommunicator::GetMPICommunicator(serial_communicator), MPI_COMM_SELF);
    KRATOS_CHECK_EQUAL(MPIDataCommunicator::GetMPICommunicator(mpi_self_communicator), MPI_COMM_SELF);
    KRATOS_CHECK_EQUAL(MPIDataCommunicator::GetMPICommunicator(mpi_world_communicator), MPI_COMM_WORLD);
    KRATOS_CHECK_NOT_EQUAL(MPIDataCommunicator::GetMPICommunicator(mpi_world_communicator), MPI_COMM_SELF);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorFromKratosComponents, KratosMPICoreFastSuite)
{
    // This should work always
    KRATOS_CHECK_EQUAL(KratosComponents<DataCommunicator>::Has("Serial"), true);
    const DataCommunicator& r_serial = KratosComponents<DataCommunicator>::Get("Serial");
    KRATOS_CHECK_EQUAL(r_serial.IsDistributed(), false);
    // This assumes running Kratos with mpi (this should be the case, since this test should be auto-disabled in serial runs)
    KRATOS_CHECK_EQUAL(KratosComponents<DataCommunicator>::Has("World"), true);
    const DataCommunicator& r_world = KratosComponents<DataCommunicator>::Get("World");
    KRATOS_CHECK_EQUAL(r_world.IsDistributed(), true);
}

// Sum ////////////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorSumIntegralTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const T world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    T local = 1;
    T result = mpi_world_communicator.Sum(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, world_size);
    }

    #ifdef KRATOS_DEBUG
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Sum(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Sum(local, -1),"is not a valid rank.");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumIntegralTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumIntegralTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumIntegralTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    double local = 2.0;
    double result = mpi_world_communicator.Sum(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, 2.0*world_size);
    }

    #ifdef KRATOS_DEBUG
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Sum(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Sum(local, -1),"is not a valid rank.");
    #endif
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumArray1d, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    constexpr int root = 0;
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    array_1d<double,3> local;
    local[0] = -1.0;
    local[1] =  0.0;
    local[2] =  1.0;

    array_1d<double,3> result = mpi_world_communicator.Sum(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result[0], -1.0*world_size);
        KRATOS_CHECK_EQUAL(result[1],  0.0);
        KRATOS_CHECK_EQUAL(result[2],  1.0*world_size);
    }

    #ifdef KRATOS_DEBUG
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Sum(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Sum(local, -1),"is not a valid rank.");
    #endif
}

namespace {
template<typename T> void MPIDataCommunicatorSumIntegralTypeVectorTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const T world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    std::vector<T> local{1, 1};
    std::vector<T> output{999, 999};

    // two-buffer version
    mpi_world_communicator.Sum(local, output, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(output[i], world_size);
        }
    }

    // return buffer version
    std::vector<T> returned_result = mpi_world_communicator.Sum(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(returned_result[i], world_size);
        }
    }

    #ifdef KRATOS_DEBUG
    // One of the inputs has a different size
    if (world_size > 1)
    {
        if (world_rank == 0) {
            local.resize(3);
            local = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Sum(local, output, root),"Input error in call to MPI_Reduce");
    }
    // Input size != output size
    std::vector<T> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Sum(local_vector_wrong_size, output, root),"Error:");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumIntegralTypeVectorTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumIntegralTypeVectorTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumLongUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumIntegralTypeVectorTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumDoubleVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    std::vector<double> local{2.0, 2.0};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    mpi_world_communicator.Sum(local, output, root);
    if (world_rank == root)
    {
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(output[i], 2.0*world_size);
        }
    }

    // return buffer version
    std::vector<double> returned_result = mpi_world_communicator.Sum(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(returned_result[i], 2.0*world_size);
        }
    }

    #ifdef KRATOS_DEBUG
    // One of the inputs has a different size
    if (world_size > 1)
    {
        if (world_rank == 0) {
            local.resize(3);
            local = {1.0,2.0,3.0};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Sum(local, output, root),"Input error in call to MPI_Reduce");
    }
    // Input size != output size
    std::vector<double> local_vector_wrong_size{1.0,2.0,3.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Sum(local_vector_wrong_size, output, root),"Error:");
    #endif
}

// Min ////////////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorMinIntegralTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    constexpr int root = 0;

    T local = world_rank;
    T result = mpi_world_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, 0);
    }

    #ifdef KRATOS_DEBUG
    const int world_size = mpi_world_communicator.Size();
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Min(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Min(local, -1),"is not a valid rank.");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinIntegralTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinIntegralTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinIntegralTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    constexpr int root = 0;

    double local = 2.0*world_rank;
    double result = mpi_world_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, 0.0);
    }

    #ifdef KRATOS_DEBUG
    const int world_size = mpi_world_communicator.Size();
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Min(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Min(local, -1),"is not a valid rank.");
    #endif
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinArray1d, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    constexpr int root = 0;
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    array_1d<double,3> local;
    local[0] = -1.0*world_rank;
    local[1] =  0.0;
    local[2] =  1.0*world_rank;

    array_1d<double,3> result = mpi_world_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result[0], -1.0*(world_size-1));
        KRATOS_CHECK_EQUAL(result[1],  0.0);
        KRATOS_CHECK_EQUAL(result[2],  0.0);
    }

    #ifdef KRATOS_DEBUG
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Min(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Min(local, -1),"is not a valid rank.");
    #endif
}

namespace {
template<typename T> void MPIDataCommunicatorMinIntegralVectorTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    constexpr int root = 0;

    std::vector<T> local{(T)world_rank, 0};
    std::vector<T> output{999, 999};

    // two-buffer version
    mpi_world_communicator.Min(local, output, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(output[0], 0);
        KRATOS_CHECK_EQUAL(output[1], 0);
    }

    // return buffer version
    std::vector<T> returned_result = mpi_world_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        KRATOS_CHECK_EQUAL(returned_result[0], 0);
        KRATOS_CHECK_EQUAL(returned_result[1], 0);
    }

    #ifdef KRATOS_DEBUG
    const int world_size = mpi_world_communicator.Size();
    // One of the inputs has a different size
    if (world_size > 1)
    {
        if (world_rank == 0) {
            local.resize(3);
            local = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Min(local, output, root),"Input error in call to MPI_Reduce");
    }
    // Input size != output size
    std::vector<T> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Min(local_vector_wrong_size, output, root),"Error:");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinIntegralVectorTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinIntegralVectorTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinLongUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinIntegralVectorTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinDoubleVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    std::vector<double> local{2.0*world_rank, -2.0*world_rank};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    mpi_world_communicator.Min(local, output, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(output[0], 0);
        KRATOS_CHECK_EQUAL(output[1], -2.0*(world_size-1));
    }

    // return buffer version
    std::vector<double> returned_result = mpi_world_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        KRATOS_CHECK_EQUAL(returned_result[0], 0.0);
        KRATOS_CHECK_EQUAL(returned_result[1], -2.0*(world_size-1));
    }

    #ifdef KRATOS_DEBUG
    // One of the inputs has a different size
    if (world_size > 1)
    {
        if (world_rank == 0) {
            local.resize(3);
            local = {1.0,2.0,3.0};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Min(local, output, root),"Input error in call to MPI_Reduce");
    }
    // Input size != output size
    std::vector<double> local_vector_wrong_size{1.0,2.0,3.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Min(local_vector_wrong_size, output, root),"Error:");
    #endif
}


// Max ////////////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorMaxIntegralTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    T local = world_rank;
    T result = mpi_world_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, (T)(world_size-1));
    }

    #ifdef KRATOS_DEBUG
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Max(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Max(local, -1),"is not a valid rank.");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxIntegralTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxIntegralTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxIntegralTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    double local = 2.0*world_rank;
    double result = mpi_world_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, 2.0*(world_size-1));
    }

    #ifdef KRATOS_DEBUG
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Max(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Max(local, -1),"is not a valid rank.");
    #endif
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxArray1d, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    constexpr int root = 0;
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    array_1d<double,3> local;
    local[0] = -1.0*world_rank;
    local[1] =  0.0;
    local[2] =  1.0*world_rank;

    array_1d<double,3> result = mpi_world_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result[0], 0.0);
        KRATOS_CHECK_EQUAL(result[1], 0.0);
        KRATOS_CHECK_EQUAL(result[2], 1.0*(world_size-1));
    }

    #ifdef KRATOS_DEBUG
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Max(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Max(local, -1),"is not a valid rank.");
    #endif
}

namespace {
template<typename T> void MPIDataCommunicatorMaxIntegralTypeVectorTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    std::vector<T> local{(T)world_rank, 0};
    std::vector<T> output{999, 999};

    // two-buffer version
    mpi_world_communicator.Max(local, output, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(output[0], (T)(world_size-1));
        KRATOS_CHECK_EQUAL(output[1], 0);
    }

    // return buffer version
    std::vector<T> returned_result = mpi_world_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        KRATOS_CHECK_EQUAL(returned_result[0], (T)(world_size-1));
        KRATOS_CHECK_EQUAL(returned_result[1], 0);
    }

    #ifdef KRATOS_DEBUG
    // One of the inputs has a different size
    if (world_size > 1)
    {
        if (world_rank == 0) {
            local.resize(3);
            local = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Max(local, output, root),"Input error in call to MPI_Reduce");
    }
    // Input size != output size
    std::vector<T> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Max(local_vector_wrong_size, output, root),"Error:");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxIntegralTypeVectorTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxIntegralTypeVectorTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxLongUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxIntegralTypeVectorTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxDoubleVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    std::vector<double> local{2.0*world_rank, -2.0*world_rank};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    mpi_world_communicator.Max(local, output, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(output[0], 2.0*(world_size-1));
        KRATOS_CHECK_EQUAL(output[1], 0.0);
    }

    // return buffer version
    std::vector<double> returned_result = mpi_world_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        KRATOS_CHECK_EQUAL(returned_result[0], 2.0*(world_size-1));
        KRATOS_CHECK_EQUAL(returned_result[1], 0.0);
    }

    #ifdef KRATOS_DEBUG
    // One of the inputs has a different size
    if (world_size > 1)
    {
        if (world_rank == 0) {
            local.resize(3);
            local = {1.0,2.0,3.0};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Max(local, output, root),"Input error in call to MPI_Reduce");
    }
    // Input size != output size
    std::vector<double> local_vector_wrong_size{1.0,2.0,3.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Max(local_vector_wrong_size, output, root),"Error:");
    #endif
}

// SumAll /////////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorSumAllIntegralTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();

    T local = 1;
    T result = mpi_world_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(result, (T)world_size);
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumAllIntegralTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumAllIntegralTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumAllIntegralTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();

    double local = 2.0;
    double result = mpi_world_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(result, 2.0*world_size);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllArray1d, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();
    array_1d<double,3> local;
    local[0] = -1.0;
    local[1] =  0.0;
    local[2] =  1.0;

    array_1d<double,3> result = mpi_world_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(result[0], -1.0*world_size);
    KRATOS_CHECK_EQUAL(result[1],  0.0);
    KRATOS_CHECK_EQUAL(result[2],  1.0*world_size);
}

namespace {
template<typename T> void MPIDataCommunicatorSumAllIntegralTypeVectorTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();

    std::vector<T> local{1, 1};
    std::vector<T> output{0, 0};

    const T expected = world_size;

    // two-buffer version
    mpi_world_communicator.SumAll(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], expected);
    }

    // return buffer version
    std::vector<T> returned_result = mpi_world_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], expected);
    }

    #ifdef KRATOS_DEBUG
    if (world_size > 1)
    {
        // One of the inputs has a different size
        if (mpi_world_communicator.Rank() == 0) {
            local.resize(3);
            local = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.SumAll(local, output),"Input error in call to MPI_Allreduce");
    }
    // Input size != output size
    std::vector<T> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.SumAll(local_vector_wrong_size, output),"Input error in call to MPI_Allreduce");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumAllIntegralTypeVectorTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumAllIntegralTypeVectorTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllUnsignedLongIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSumAllIntegralTypeVectorTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllDoubleVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();

    std::vector<double> local{2.0, 2.0};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    mpi_world_communicator.SumAll(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], 2.0*world_size);
    }

    // return buffer version
    std::vector<double> returned_result = mpi_world_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], 2.0*world_size);
    }

    #ifdef KRATOS_DEBUG
    if (world_size > 1)
    {
        // One of the inputs has a different size
        if (mpi_world_communicator.Rank() == 0) {
            local.resize(3);
            local = {1.0,2.0,3.0};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.SumAll(local, output),"Input error in call to MPI_Allreduce");
    }
    // Input size != output size
    std::vector<double> local_vector_wrong_size{1.0,2.0,3.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.SumAll(local_vector_wrong_size, output),"Input error in call to MPI_Allreduce");
    #endif
}

// MinAll /////////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorMinAllIntegralTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();

    T local = world_rank;
    T result = mpi_world_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(result, 0);
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinAllIntegralTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinAllIntegralTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinAllIntegralTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();

    double local = 2.0*world_rank;
    double result = mpi_world_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(result, 0.0);
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllArray1d, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    array_1d<double,3> local;
    local[0] = -1.0*world_rank;
    local[1] =  0.0;
    local[2] =  1.0*world_rank;

    array_1d<double,3> result = mpi_world_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(result[0], -1.0*(world_size-1));
    KRATOS_CHECK_EQUAL(result[1],  0.0);
    KRATOS_CHECK_EQUAL(result[2],  0.0);
}

namespace {
template<typename T> void MPIDataCommunicatorMinAllIntegralTypeVectorTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();

    std::vector<T> local{(T)world_rank, 0};
    std::vector<T> output{999, 999};

    // two-buffer version
    mpi_world_communicator.MinAll(local, output);
    KRATOS_CHECK_EQUAL(output[0], 0);
    KRATOS_CHECK_EQUAL(output[1], 0);

    // return buffer version
    std::vector<T> returned_result = mpi_world_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    KRATOS_CHECK_EQUAL(returned_result[0], 0);
    KRATOS_CHECK_EQUAL(returned_result[1], 0);

    #ifdef KRATOS_DEBUG
    const int world_size = mpi_world_communicator.Size();
    if (world_size > 1)
    {
        // One of the inputs has a different size
        if (mpi_world_communicator.Rank() == 0) {
            local.resize(3);
            local = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.MinAll(local, output),"Input error in call to MPI_Allreduce");
    }
    // Input size != output size
    std::vector<T> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.MinAll(local_vector_wrong_size, output),"Input error in call to MPI_Allreduce");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinAllIntegralTypeVectorTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinAllIntegralTypeVectorTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllLongUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMinAllIntegralTypeVectorTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllDoubleVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    std::vector<double> local{2.0*world_rank, -2.0*world_rank};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    mpi_world_communicator.MinAll(local, output);
    KRATOS_CHECK_EQUAL(output[0], 0);
    KRATOS_CHECK_EQUAL(output[1], -2.0*(world_size-1));

    // return buffer version
    std::vector<double> returned_result = mpi_world_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    KRATOS_CHECK_EQUAL(returned_result[0], 0.0);
    KRATOS_CHECK_EQUAL(returned_result[1], -2.0*(world_size-1));

    #ifdef KRATOS_DEBUG
    if (world_size > 1)
    {
        // One of the inputs has a different size
        if (mpi_world_communicator.Rank() == 0) {
            local.resize(3);
            local = {1.0,2.0,3.0};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.MinAll(local, output),"Input error in call to MPI_Allreduce");
    }
    // Input size != output size
    std::vector<double> local_vector_wrong_size{1.0,2.0,3.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.MinAll(local_vector_wrong_size, output),"Input error in call to MPI_Allreduce");
    #endif
}

// MaxAll /////////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorMaxAllIntegralTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    T local = world_rank;
    T result = mpi_world_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(result, (T)(world_size-1));
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxAllIntegralTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxAllIntegralTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxAllIntegralTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    double local = 2.0*world_rank;
    double result = mpi_world_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(result, 2.0*(world_size-1));
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllArray1d, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    array_1d<double,3> local;
    local[0] = -1.0*world_rank;
    local[1] =  0.0;
    local[2] =  1.0*world_rank;

    array_1d<double,3> result = mpi_world_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(result[0], 0.0);
    KRATOS_CHECK_EQUAL(result[1], 0.0);
    KRATOS_CHECK_EQUAL(result[2], 1.0*(world_size-1));
}

namespace {
template<typename T> void MPIDataCommunicatorMaxAllIntegralTypeVectorTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    std::vector<T> local{(T)world_rank, 0};
    std::vector<T> output{999, 999};

    const T expected = world_size - 1;

    // two-buffer version
    mpi_world_communicator.MaxAll(local, output);
    KRATOS_CHECK_EQUAL(output[0], expected);
    KRATOS_CHECK_EQUAL(output[1], 0);

    // return buffer version
    std::vector<T> returned_result = mpi_world_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    KRATOS_CHECK_EQUAL(returned_result[0], expected);
    KRATOS_CHECK_EQUAL(returned_result[1], 0);

    #ifdef KRATOS_DEBUG
    if (world_size > 1)
    {
        // One of the inputs has a different size
        if (mpi_world_communicator.Rank() == 0) {
            local.resize(3);
            local = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.MaxAll(local, output),"Input error in call to MPI_Allreduce");
    }
    // Input size != output size
    std::vector<T> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.MaxAll(local_vector_wrong_size, output),"Input error in call to MPI_Allreduce");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxAllIntegralTypeVectorTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxAllIntegralTypeVectorTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllLongUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorMaxAllIntegralTypeVectorTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllDoubleVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    std::vector<double> local{2.0*world_rank, -2.0*world_rank};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    mpi_world_communicator.MaxAll(local, output);
    KRATOS_CHECK_EQUAL(output[0], 2.0*(world_size-1));
    KRATOS_CHECK_EQUAL(output[1], 0.0);

    // return buffer version
    std::vector<double> returned_result = mpi_world_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    KRATOS_CHECK_EQUAL(returned_result[0], 2.0*(world_size-1));
    KRATOS_CHECK_EQUAL(returned_result[1], 0.0);

    #ifdef KRATOS_DEBUG
    if (world_size > 1)
    {
        // One of the inputs has a different size
        if (mpi_world_communicator.Rank() == 0) {
            local.resize(3);
            local = {1.0,2.0,3.0};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.MaxAll(local, output),"Input error in call to MPI_Allreduce");
    }
    // Input size != output size
    std::vector<double> local_vector_wrong_size{1.0,2.0,3.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.MaxAll(local_vector_wrong_size, output),"Input error in call to MPI_Allreduce");
    #endif
}

// ScanSum ////////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorScanSumIntegralTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    int rank = mpi_world_communicator.Rank();

    T local_total = 1;
    T partial_sum = mpi_world_communicator.ScanSum(local_total);
    KRATOS_CHECK_EQUAL(partial_sum, (T)(rank + 1));
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScanSumInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScanSumIntegralTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScanSumUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScanSumIntegralTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScanSumLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScanSumIntegralTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScanSumDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    int rank = mpi_world_communicator.Rank();

    double local_total = 2.0;
    double partial_sum = mpi_world_communicator.ScanSum(local_total);
    KRATOS_CHECK_EQUAL(partial_sum, 2.0*(rank + 1));
}

namespace {
template<typename T> void MPIDataCommunicatorScanSumIntegralVectorTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    int rank = mpi_world_communicator.Rank();

    std::vector<T> local_total{1,1};
    std::vector<T> output{0, 0};

    const T expected = rank + 1;

    // two-buffer version
    mpi_world_communicator.ScanSum(local_total, output);
    for (unsigned int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], expected);
    }

    // return buffer version
    std::vector<T> partial_sum = mpi_world_communicator.ScanSum(local_total);
    KRATOS_CHECK_EQUAL(partial_sum.size(), 2);
    KRATOS_CHECK_EQUAL(partial_sum[0], expected);
    KRATOS_CHECK_EQUAL(partial_sum[1], expected);

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the inputs has a different size
        if (rank == 0) {
            local_total.resize(3);
            local_total = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.ScanSum(local_total, partial_sum),"Input error in call to MPI_Scan");
    }
    // Input size != output size
    std::vector<T> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.ScanSum(local_vector_wrong_size, partial_sum),"Input error in call to MPI_Scan");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScanSumVectorInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScanSumIntegralVectorTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScanSumVectorUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScanSumIntegralVectorTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScanSumVectorLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScanSumIntegralVectorTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScanSumVectorDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    int rank = mpi_world_communicator.Rank();

    std::vector<double> local_total{2.0,2.0};
    std::vector<double> output{-1.0, -1.0};

    // two-buffer version
    mpi_world_communicator.ScanSum(local_total, output);
    for (unsigned int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], 2.0*(rank + 1));
    }

    // return buffer version
    std::vector<double> partial_sum = mpi_world_communicator.ScanSum(local_total);
    KRATOS_CHECK_EQUAL(partial_sum.size(), 2);
    KRATOS_CHECK_EQUAL(partial_sum[0], 2.0*(rank + 1));
    KRATOS_CHECK_EQUAL(partial_sum[1], 2.0*(rank + 1));

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the inputs has a different size
        if (rank == 0) {
            local_total.resize(3);
            local_total = {1.0,2.0,3.0};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.ScanSum(local_total, partial_sum),"Input error in call to MPI_Scan");
    }
    // Input size != output size
    std::vector<double> local_vector_wrong_size{1.0,2.0,3.0};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.ScanSum(local_vector_wrong_size, partial_sum),"Input error in call to MPI_Scan");
    #endif
}

// SendRecv ///////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorSendRecvIntegralTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_rank + 1 == world_size ? 0 : world_rank + 1;
    const int recv_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;

    T send_value(world_rank);
    T recv_value(999);
    std::vector<T> send_buffer{send_value, send_value};
    std::vector<T> recv_buffer{999, 999};

    if (world_size > 1)
    {
        const T expected_recv = world_rank > 0 ? world_rank - 1 : world_size - 1;

        // value two-buffer version
        mpi_world_communicator.SendRecv(send_value, send_rank, 0, recv_value, recv_rank, 0);
        KRATOS_CHECK_EQUAL(recv_value, expected_recv);

        // value return version
        T return_value = mpi_world_communicator.SendRecv(send_value, send_rank, 0, recv_rank, 0);
        KRATOS_CHECK_EQUAL(return_value, expected_recv);

        // vector two-buffer version
        mpi_world_communicator.SendRecv(send_buffer, send_rank, 0, recv_buffer, recv_rank, 0);

        // vector return version
        std::vector<T> return_buffer = mpi_world_communicator.SendRecv(send_buffer, send_rank, recv_rank);

        KRATOS_CHECK_EQUAL(return_buffer.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(recv_buffer[i], expected_recv);
            KRATOS_CHECK_EQUAL(return_buffer[i], expected_recv);
        }
    }
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSendRecvInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSendRecvIntegralTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSendRecvUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSendRecvIntegralTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSendRecvLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorSendRecvIntegralTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSendRecvDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_rank + 1 == world_size ? 0 : world_rank + 1;
    const int recv_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;

    double send_value(2.0*world_rank);
    double recv_value(-1.0);
    std::vector<double> send_buffer{2.0*world_rank, 2.0*world_rank};
    std::vector<double> recv_buffer{-1.0, -1.0};

    if (world_size > 1)
    {
        const double expected_recv = world_rank > 0 ? 2.0*(world_rank - 1) : 2.0*(world_size - 1);

        // value two-buffer version
        mpi_world_communicator.SendRecv(send_value, send_rank, 0, recv_value, recv_rank, 0);
        KRATOS_CHECK_EQUAL(recv_value, expected_recv);

        // value return version
        double return_value = mpi_world_communicator.SendRecv(send_value, send_rank, 0, recv_rank, 0);
        KRATOS_CHECK_EQUAL(return_value, expected_recv);

        // two-buffer version
        mpi_world_communicator.SendRecv(send_buffer, send_rank, 0, recv_buffer, recv_rank, 0);

        // return version
        std::vector<double> return_buffer = mpi_world_communicator.SendRecv(send_buffer, send_rank, recv_rank);

        KRATOS_CHECK_EQUAL(return_buffer.size(), 2);
        for (int i = 0; i < 2; i++)
        {
            KRATOS_CHECK_EQUAL(recv_buffer[i], expected_recv);
            KRATOS_CHECK_EQUAL(return_buffer[i], expected_recv);
        }
    }
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorSendRecvString, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_rank + 1 == world_size ? 0 : world_rank + 1;
    const int recv_rank = world_rank == 0 ? world_size - 1 : world_rank - 1;

    std::string send_buffer("Hello world!");
    std::string recv_buffer;
    recv_buffer.resize(send_buffer.size()); // here we assume both send and recv buffers have the same size

    if (world_size > 1)
    {
        // two-buffer version
        mpi_world_communicator.SendRecv(send_buffer, send_rank, 0, recv_buffer, recv_rank, 0);

        // return version
        std::string return_buffer = mpi_world_communicator.SendRecv(send_buffer, send_rank, recv_rank);

        KRATOS_CHECK_EQUAL(return_buffer.size(), 12);
        KRATOS_CHECK_EQUAL(recv_buffer, send_buffer);
        KRATOS_CHECK_EQUAL(return_buffer, send_buffer);
    }
}

// Broadcast //////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorBroadcastIntegralTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_size-1;

    T send = world_rank == send_rank ? 1 : 0;

    mpi_world_communicator.Broadcast(send,send_rank);
    KRATOS_CHECK_EQUAL(send, 1);
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorBroadcastInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorBroadcastIntegralTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorBroadcastUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorBroadcastIntegralTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorBroadcastLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorBroadcastIntegralTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorBroadcastDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_size-1;

    double send = world_rank == send_rank ? 2.0 : 0.0;

    mpi_world_communicator.Broadcast(send,send_rank);
    KRATOS_CHECK_EQUAL(send, 2.0);
}

namespace {
template<typename T> void MPIDataCommunicatorBroadcastIntegralTypeVectorTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_size-1;

    std::vector<T> send = world_rank == send_rank ? std::vector<T>{1, 1} : std::vector<T>{0, 0};

    mpi_world_communicator.Broadcast(send,send_rank);
    for (unsigned int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(send[i], 1);
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        if (world_rank == 0)
        {
            send.resize(3);
            send = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Broadcast(send, send_rank),"Input error in call to MPI_Bcast");
    }
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorBroadcastIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorBroadcastIntegralTypeVectorTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorBroadcastUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorBroadcastIntegralTypeVectorTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorBroadcastLongUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorBroadcastIntegralTypeVectorTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorBroadcastDoubleVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_size-1;

    std::vector<double> send = world_rank == send_rank ? std::vector<double>{2.0, 2.0} : std::vector<double>{0.0, 0.0};

    mpi_world_communicator.Broadcast(send,send_rank);
    for (unsigned int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(send[i], 2.0);
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        if (world_rank == 0)
        {
            send.resize(3);
            send = {1.0,2.0,3.0};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Broadcast(send, send_rank),"Input error in call to MPI_Bcast");
    }
    #endif
}

// Scatter ////////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorScatterIntegralTypeVectorTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = 0;

    std::vector<T> send_buffer(0);
    std::vector<T> recv_buffer{0, 0};

    if (world_rank == send_rank)
    {
        send_buffer.resize(2*world_size);
        for (int i = 0; i < 2*world_size; i++)
        {
            send_buffer[i] = 1;
        }
    }

    // two-buffer version
    mpi_world_communicator.Scatter(send_buffer, recv_buffer, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], 1);
    }

    // return version
    std::vector<T> return_buffer = mpi_world_communicator.Scatter(send_buffer, send_rank);
    KRATOS_CHECK_EQUAL(return_buffer.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], 1);
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        std::vector<T> wrong_recv = {999, 999};
        if (world_rank == 0)
        {
            recv_buffer.resize(3);
            recv_buffer = {999, 999, 999};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Scatter(send_buffer, recv_buffer, send_rank),"Error");
    }
    // send rank has wrong size
    std::vector<T> wrong_send(0);
    if (world_rank == send_rank)
    {
        wrong_send = {1};
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Scatter(send_buffer, recv_buffer, send_rank),"Error");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScatterIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScatterIntegralTypeVectorTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScatterUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScatterIntegralTypeVectorTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScatterLongUnsignedIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScatterIntegralTypeVectorTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScatterDoubleVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = 0;

    std::vector<double> send_buffer(0);
    std::vector<double> recv_buffer{0.0, 0.0};

    if (world_rank == send_rank)
    {
        send_buffer.resize(2*world_size);
        for (int i = 0; i < 2*world_size; i++)
        {
            send_buffer[i] = 2.0;
        }
    }

    // two-buffer version
    mpi_world_communicator.Scatter(send_buffer, recv_buffer, send_rank);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], 2.0);
    }

    // return version
    std::vector<double> return_buffer = mpi_world_communicator.Scatter(send_buffer, send_rank);
    KRATOS_CHECK_EQUAL(return_buffer.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(return_buffer[i], 2.0);
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        std::vector<double> wrong_recv = {-1.0, -1.0};
        if (world_rank == 0)
        {
            recv_buffer.resize(3);
            recv_buffer = {-1.0,-1.0,-1.0};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Scatter(send_buffer, recv_buffer, send_rank),"Error");
    }
    // send rank has wrong size
    std::vector<double> wrong_send(0);
    if (world_rank == send_rank)
    {
        wrong_send = {1.0};
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Scatter(send_buffer, recv_buffer, send_rank),"Error");
    #endif
}

// Scatterv ///////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorScattervIntegralTypeVectorTest()
{
    /* send message for ints is {0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, ...} (max 6 values per rank)
     * read only first <rank> values of the message per rank (up to 5 values per rank)
     * message containing doubles is double of message containing ints
     */
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int send_rank = world_size-1;

    std::vector<T> send_buffer(0);
    std::vector<int> send_sizes(0);
    std::vector<int> send_offsets(0);

    auto make_message_size = [](int rank) { return rank < 5 ? rank : 5; };
    auto make_message_distance = [](int rank, int padding) {
        return rank < 5 ? ((rank-1)*rank)/2 + rank*padding : rank*(5+padding) - 15;
    };

    const int recv_size = make_message_size(world_rank);
    std::vector<T> recv_buffer(recv_size, 999);

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
    mpi_world_communicator.Scatterv(send_buffer, send_sizes, send_offsets, recv_buffer, send_rank);

    for (int i = 0; i < recv_size; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], (T)world_rank);
    }

    // return buffer version

    // pre-process: prepare input vector-of-vectors
    std::vector<std::vector<T>> scatterv_message;
    if (world_rank == send_rank)
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
    std::vector<T> result = mpi_world_communicator.Scatterv(scatterv_message, send_rank);
    for (int i = 0; i < recv_size; i++)
    {
        KRATOS_CHECK_EQUAL(result[i], (T)world_rank);
    }

    #ifdef KRATOS_DEBUG
    // send sizes do not match
    std::vector<int> wrong_send_sizes = send_sizes;
    if (world_rank == send_rank) {
        wrong_send_sizes[0] += 1;
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Scatterv(send_buffer, wrong_send_sizes, send_offsets, recv_buffer, send_rank),
        "Error");

    // sent message is too large
    if (world_size > 1) // This test should be skipped when running with one rank, since recv_buffer.size() is 0 in that case
    {
        std::vector<T> wrong_recv_message;
        if (world_rank == send_rank)
        {
            wrong_recv_message.resize(recv_buffer.size()-1);
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            mpi_world_communicator.Scatterv(send_buffer, send_sizes, send_offsets, wrong_recv_message, send_rank),
            "Error");
    }

    // sent offsets overflow
    std::vector<int> wrong_send_offsets = send_offsets;
    if (world_rank == send_rank)
    {
        wrong_send_offsets[world_size - 1] += 5;
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Scatterv(send_buffer, send_sizes, wrong_send_offsets, recv_buffer, send_rank),
        "Error");

    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScattervInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScattervIntegralTypeVectorTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScattervUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScattervIntegralTypeVectorTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScattervLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorScattervIntegralTypeVectorTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorScattervDouble, KratosMPICoreFastSuite)
{
    /* send message for ints is {0, 1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4, ...} (max 6 values per rank)
     * read only first <rank> values of the message per rank (up to 5 values per rank)
     * message containing doubles is double of message containing ints
     */
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
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
    mpi_world_communicator.Scatterv(send_buffer, send_sizes, send_offsets, recv_buffer, send_rank);

    for (int i = 0; i < recv_size; i++)
    {
        KRATOS_CHECK_EQUAL(recv_buffer[i], 2.0*world_rank);
    }

    // return buffer version

    // pre-process: prepare input vector-of-vectors
    std::vector<std::vector<double>> scatterv_message;
    if (world_rank == send_rank)
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
    std::vector<double> result = mpi_world_communicator.Scatterv(scatterv_message, send_rank);
    for (int i = 0; i < recv_size; i++)
    {
        KRATOS_CHECK_EQUAL(result[i], 2.0*world_rank);
    }

    #ifdef KRATOS_DEBUG
    // send sizes do not match
    std::vector<int> wrong_send_sizes = send_sizes;
    if (world_rank == send_rank) {
        wrong_send_sizes[0] += 1;
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Scatterv(send_buffer, wrong_send_sizes, send_offsets, recv_buffer, send_rank),
        "Error");

    // sent message is too large
    if (world_size > 1) // This test should be skipped when running with one rank, since recv_buffer.size() is 0 in that case
    {
        std::vector<double> wrong_recv_message;
        if (world_rank == send_rank)
        {
            wrong_recv_message.resize(recv_buffer.size()-1);
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            mpi_world_communicator.Scatterv(send_buffer, send_sizes, send_offsets, wrong_recv_message, send_rank),
            "Error");
    }

    // sent offsets overflow
    std::vector<int> wrong_send_offsets = send_offsets;
    if (world_rank == send_rank)
    {
        wrong_send_offsets[world_size - 1] += 5;
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Scatterv(send_buffer, send_sizes, wrong_send_offsets, recv_buffer, send_rank),
        "Error");

    #endif
}

// Gather /////////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorGatherIntegralTypeVectorTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int recv_rank = 0;

    std::vector<T> send_buffer{(T)world_rank, (T)world_rank};
    std::vector<T> recv_buffer;

    if (world_rank == recv_rank)
    {
        recv_buffer = std::vector<T>(2*world_size, 999);
    }

    // two-buffer version
    mpi_world_communicator.Gather(send_buffer, recv_buffer, recv_rank);

    if (world_rank == recv_rank)
    {
        for (int rank = 0; rank < world_size; rank++)
        {
            for (int j = 2*rank; j < 2*rank+2; j++)
            {
                KRATOS_CHECK_EQUAL(recv_buffer[j], (T)rank);
            }
        }
    }

    // return buffer version
    std::vector<T> return_buffer = mpi_world_communicator.Gather(send_buffer, recv_rank);
    if (world_rank == recv_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), static_cast<unsigned int>(2*world_size));
        for (int rank = 0; rank < world_size; rank++)
        {
            for (int j = 2*rank; j < 2*rank+2; j++)
            {
                KRATOS_CHECK_EQUAL(return_buffer[j], (T)rank);
            }
        }
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        std::vector<T> wrong_buffer{(T)world_rank, (T)world_rank};
        if (world_rank == 0)
        {
            wrong_buffer.push_back((T)world_rank);
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Gather(wrong_buffer, recv_buffer, recv_rank),"Error");
    }
    // recv rank has wrong size
    recv_buffer.push_back(0);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Gather(send_buffer, recv_buffer, recv_rank),"Error");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorGatherInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorGatherIntegralTypeVectorTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorGatherUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorGatherIntegralTypeVectorTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorGatherLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorGatherIntegralTypeVectorTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorGatherDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int recv_rank = 0;

    std::vector<double> send_buffer{2.0*world_rank, 2.0*world_rank};
    std::vector<double> recv_buffer;

    if (world_rank == recv_rank)
    {
        recv_buffer = std::vector<double>(2*world_size, -1);
    }

    // two-buffer version
    mpi_world_communicator.Gather(send_buffer, recv_buffer, recv_rank);

    if (world_rank == recv_rank)
    {
        for (int rank = 0; rank < world_size; rank++)
        {
            for (int j = 2*rank; j < 2*rank+2; j++)
            {
                KRATOS_CHECK_EQUAL(recv_buffer[j], 2.0*rank);
            }
        }
    }

    // return buffer version
    std::vector<double> return_buffer = mpi_world_communicator.Gather(send_buffer, recv_rank);
    if (world_rank == recv_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), static_cast<unsigned int>(2*world_size));
        for (int rank = 0; rank < world_size; rank++)
        {
            for (int j = 2*rank; j < 2*rank+2; j++)
            {
                KRATOS_CHECK_EQUAL(return_buffer[j], 2.0*rank);
            }
        }
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        std::vector<double> wrong_buffer{2.0*world_rank, 2.0*world_rank};
        if (world_rank == 0)
        {
            wrong_buffer.push_back(2.0*world_rank);
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Gather(wrong_buffer, recv_buffer, recv_rank),"Error");
    }
    // recv rank has wrong size
    recv_buffer.push_back(0.0);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Gather(send_buffer, recv_buffer, recv_rank),"Error");
    #endif
}

// Gatherv ////////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorGathervIntegralTypeVectorTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
    const int recv_rank = world_size-1;

    auto make_message_size = [](int rank) { return rank < 5 ? rank : 5; };
    auto make_message_distance = [](int rank, int padding) {
        return rank < 5 ? ((rank-1)*rank)/2 + rank*padding : rank*(5+padding) - 15;
    };

    const int send_size = make_message_size(world_rank);
    std::vector<T> send_buffer(send_size, (T)world_rank);

    // two-buffer version
    const int message_padding = 1;
    const int recv_size = make_message_distance(world_size, message_padding);
    std::vector<T> recv_buffer(0);
    std::vector<int> recv_sizes(0);
    std::vector<int> recv_offsets(0);

    if (world_rank == recv_rank)
    {
        recv_buffer.resize(recv_size, 999);
        recv_sizes.resize(world_size);
        recv_offsets.resize(world_size);
        for (int rank = 0; rank < world_size; rank++)
        {
            recv_sizes[rank] = make_message_size(rank);
            recv_offsets[rank] = make_message_distance(rank, message_padding);
        }
    }

    mpi_world_communicator.Gatherv(send_buffer, recv_buffer, recv_sizes, recv_offsets, recv_rank);

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
                KRATOS_CHECK_EQUAL(recv_buffer[i], (T)rank);
            }
            // ...followed by the expected padding.
            for (int i = recv_offset + recv_size; i < recv_offset + recv_size + message_padding; i++)
            {
                KRATOS_CHECK_EQUAL(recv_buffer[i], 999);
            }

        }
    }

    // return buffer version
    std::vector<std::vector<T>> return_buffer = mpi_world_communicator.Gatherv(send_buffer, recv_rank);

    if (world_rank == recv_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), static_cast<unsigned int>(world_size));
        for (int rank = 0; rank < world_size; rank++)
        {
            unsigned int expected_size = make_message_size(rank);
            KRATOS_CHECK_EQUAL(return_buffer[rank].size(), expected_size);
            for (unsigned int i = 0; i < expected_size; i++)
            {
                KRATOS_CHECK_EQUAL(return_buffer[rank][i], (T)rank);
            }
            // no padding in return version
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
            send_buffer, recv_buffer, wrong_recv_sizes, recv_offsets, recv_rank),
            "Error");

    // recv message is too small
    std::vector<T> wrong_recv_message;
    if (world_rank == recv_size)
    {
        wrong_recv_message.resize(recv_buffer.size()-1);
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Gatherv(send_buffer, wrong_recv_message, recv_sizes, recv_offsets, recv_rank),
        "Error");
    // sent offsets overflow
    std::vector<int> wrong_recv_offsets = recv_offsets;
    if (world_rank == recv_rank)
    {
        wrong_recv_offsets[world_size - 1] += 5;
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Gatherv(send_buffer, recv_buffer, recv_sizes, wrong_recv_offsets, recv_rank),
        "Error");

    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorGathervInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorGathervIntegralTypeVectorTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorGathervUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorGathervIntegralTypeVectorTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorGathervLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorGathervIntegralTypeVectorTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorGathervDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();
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

    mpi_world_communicator.Gatherv(send_buffer, recv_buffer, recv_sizes, recv_offsets, recv_rank);

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
                KRATOS_CHECK_EQUAL(recv_buffer[i], 2.0*rank);
            }
            // ...followed by the expected padding.
            for (int i = recv_offset + recv_size; i < recv_offset + recv_size + message_padding; i++)
            {
                KRATOS_CHECK_EQUAL(recv_buffer[i], -1.0);
            }

        }
    }

    // return buffer version
    std::vector<std::vector<double>> return_buffer = mpi_world_communicator.Gatherv(send_buffer, recv_rank);

    if (world_rank == recv_rank)
    {
        KRATOS_CHECK_EQUAL(return_buffer.size(), static_cast<unsigned int>(world_size));
        for (int rank = 0; rank < world_size; rank++)
        {
            unsigned int expected_size = make_message_size(rank);
            KRATOS_CHECK_EQUAL(return_buffer[rank].size(), expected_size);
            for (unsigned int i = 0; i < expected_size; i++)
            {
                KRATOS_CHECK_EQUAL(return_buffer[rank][i], 2.0*rank);
            }
            // no padding in return version
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
            send_buffer, recv_buffer, wrong_recv_sizes, recv_offsets, recv_rank),
            "Error");

    // recv message is too small
    std::vector<double> wrong_recv_message;
    if (world_rank == recv_size)
    {
        wrong_recv_message.resize(recv_buffer.size()-1);
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Gatherv(send_buffer, wrong_recv_message, recv_sizes, recv_offsets, recv_rank),
        "Error");

    // sent offsets overflow
    std::vector<int> wrong_recv_offsets = recv_offsets;
    if (world_rank == recv_rank)
    {
        wrong_recv_offsets[world_size - 1] += 5;
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Gatherv(send_buffer, recv_buffer, recv_sizes, wrong_recv_offsets, recv_rank),
        "Error");

    #endif
}

// AllGather //////////////////////////////////////////////////////////////////

namespace {
template<typename T> void MPIDataCommunicatorAllGatherIntegralTypeTest()
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();

    std::vector<T> send_buffer{(T)world_rank, (T)world_rank};

    // two-buffer version
    std::vector<T> recv_buffer(2*world_size, 999);
    mpi_world_communicator.AllGather(send_buffer, recv_buffer);

    for (int rank = 0; rank < world_size; rank++)
    {
        for (int j = 2*rank; j < 2*rank+2; j++)
        {
            KRATOS_CHECK_EQUAL(recv_buffer[j], (T)rank);
        }
    }

    // return buffer version
    std::vector<T> return_buffer = mpi_world_communicator.AllGather(send_buffer);

    KRATOS_CHECK_EQUAL(return_buffer.size(), static_cast<unsigned int>(2*world_size));
    for (int rank = 0; rank < world_size; rank++)
    {
        for (int j = 2*rank; j < 2*rank+2; j++)
        {
            KRATOS_CHECK_EQUAL(return_buffer[j], (T)rank);
        }
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        std::vector<T> wrong_send{999, 999};
        if (world_rank == 0)
        {
            wrong_send.push_back(999);
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            mpi_world_communicator.AllGather(wrong_send, recv_buffer),
            "Input error in call to MPI_Allgather");
    }
    // recv rank has wrong size
    recv_buffer.push_back(999);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.AllGather(send_buffer, recv_buffer),
        "Input error in call to MPI_Allgather");
    #endif
}
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorAllGatherInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorAllGatherIntegralTypeTest<int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorAllGatherUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorAllGatherIntegralTypeTest<unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorAllGatherLongUnsignedInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicatorAllGatherIntegralTypeTest<long unsigned int>();
}

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorAllGatherDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    const int world_size = mpi_world_communicator.Size();
    const int world_rank = mpi_world_communicator.Rank();

    std::vector<double> send_buffer{2.0*world_rank, 2.0*world_rank};

    // two-buffer version
    std::vector<double> recv_buffer(2*world_size, -1.0);
    mpi_world_communicator.AllGather(send_buffer, recv_buffer);

    for (int rank = 0; rank < world_size; rank++)
    {
        for (int j = 2*rank; j < 2*rank+2; j++)
        {
            KRATOS_CHECK_EQUAL(recv_buffer[j], 2.0*rank);
        }
    }

    // return buffer version
    std::vector<double> return_buffer = mpi_world_communicator.AllGather(send_buffer);

    KRATOS_CHECK_EQUAL(return_buffer.size(), static_cast<unsigned int>(2*world_size));
    for (int rank = 0; rank < world_size; rank++)
    {
        for (int j = 2*rank; j < 2*rank+2; j++)
        {
            KRATOS_CHECK_EQUAL(return_buffer[j], 2.0*rank);
        }
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has a different size
        std::vector<double> wrong_send{-1.0, -1.0};
        if (world_rank == 0)
        {
            wrong_send.push_back(-1.0);
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(
            mpi_world_communicator.AllGather(wrong_send, recv_buffer),
            "Input error in call to MPI_Allgather");
    }
    // recv rank has wrong size
    recv_buffer.push_back(-1);
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.AllGather(send_buffer, recv_buffer),
        "Input error in call to MPI_Allgather");
    #endif
}

// Error broadcasting methods /////////////////////////////////////////////////

KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPIDataCommunicatorErrorBroadcasting, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int rank = mpi_world_communicator.Rank();
    const int size = mpi_world_communicator.Size();

    auto broadcast_if_true_test = [&mpi_world_communicator](){
        KRATOS_ERROR_IF(mpi_world_communicator.BroadcastErrorIfTrue(true,0) )
        << "Something went wrong in rank 0." << std::endl;
    };
    auto broadcast_if_false_test = [&mpi_world_communicator](){
        KRATOS_ERROR_IF_NOT(mpi_world_communicator.BroadcastErrorIfFalse(false,0) )
        << "Something went wrong in rank 0." << std::endl;
    };

    std::stringstream broadcast_message;
    if (rank == 0)
    {
        broadcast_message << "Something went wrong in rank 0.";
    }
    else
    {
        broadcast_message << "Stopping because of error in rank 0.";
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        broadcast_if_true_test(),
        broadcast_message.str()
    );

    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        broadcast_if_false_test(),
        broadcast_message.str()
    );

    auto true_on_any_rank_test = [&mpi_world_communicator, &rank](int fail_rank) {
        KRATOS_ERROR_IF(mpi_world_communicator.ErrorIfTrueOnAnyRank(rank == fail_rank))
        << "Something went wrong on rank " << rank << "." << std::endl;
    };

    auto false_on_any_rank_test = [&mpi_world_communicator, &rank](int fail_rank) {
        KRATOS_ERROR_IF_NOT(mpi_world_communicator.ErrorIfFalseOnAnyRank(rank < fail_rank))
        << "Something went wrong on rank " << rank << "." << std::endl;
    };

    std::stringstream error_on_any_rank_message;
    if (rank == size - 1)
    {
        error_on_any_rank_message << "Something went wrong on rank " << rank << ".";
    }
    else
    {
        error_on_any_rank_message << "Stopping because an error was detected on a different rank.";
    }
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        true_on_any_rank_test(size-1),
        error_on_any_rank_message.str()
    );
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        false_on_any_rank_test(size-1),
        error_on_any_rank_message.str()
    );
}

}
}
