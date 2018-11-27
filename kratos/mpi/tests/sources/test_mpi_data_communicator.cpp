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

#include "testing/testing.h"

namespace Kratos {

namespace Testing {

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorRankAndSize, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPICommRetrieval, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_self_communicator(MPI_COMM_SELF);
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    KRATOS_CHECK_EQUAL(MPIEnvironment::GetMPICommunicator(serial_communicator), MPI_COMM_SELF);
    KRATOS_CHECK_EQUAL(MPIEnvironment::GetMPICommunicator(mpi_self_communicator), MPI_COMM_SELF);
    KRATOS_CHECK_EQUAL(MPIEnvironment::GetMPICommunicator(mpi_world_communicator), MPI_COMM_WORLD);
    KRATOS_CHECK_NOT_EQUAL(MPIEnvironment::GetMPICommunicator(mpi_world_communicator), MPI_COMM_SELF);
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorFromKratosComponents, KratosMPICoreFastSuite)
{
    // This should work always
    KRATOS_CHECK_EQUAL(KratosComponents<DataCommunicator>::Has("Serial"), true);
    const DataCommunicator& r_serial = KratosComponents<DataCommunicator>::Get("Serial");
    KRATOS_CHECK_EQUAL(r_serial.IsDistributed(), false);
    // This assumes running Kratos with mpi (this should be the case, since this test's suite is part of the MPI core)
    KRATOS_CHECK_EQUAL(KratosComponents<DataCommunicator>::Has("World"), true);
    const DataCommunicator& r_world = KratosComponents<DataCommunicator>::Get("World");
    KRATOS_CHECK_EQUAL(r_world.IsDistributed(), true);
}

// Sum ////////////////////////////////////////////////////////////////////////

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorSumInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    int local = 1;
    int result = mpi_world_communicator.Sum(local, root);
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorSumDouble, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorSumArray1d, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    std::vector<int> local{1, 1};
    std::vector<int> output{-1, -1};

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
    std::vector<int> returned_result = mpi_world_communicator.Sum(local, root);
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
    std::vector<int> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Sum(local_vector_wrong_size, output, root),"Error:");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumDoubleVector, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMinInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    int local = world_rank;
    int result = mpi_world_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, 0);
    }

    #ifdef KRATOS_DEBUG
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Min(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Min(local, -1),"is not a valid rank.");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMinDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    double local = 2.0*world_rank;
    double result = mpi_world_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, 0.0);
    }

    #ifdef KRATOS_DEBUG
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Min(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Min(local, -1),"is not a valid rank.");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMinArray1d, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    std::vector<int> local{world_rank, -world_rank};
    std::vector<int> output{-1, -1};

    // two-buffer version
    mpi_world_communicator.Min(local, output, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(output[0], 0);
        KRATOS_CHECK_EQUAL(output[1], -(world_size-1));
    }

    // return buffer version
    std::vector<int> returned_result = mpi_world_communicator.Min(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        KRATOS_CHECK_EQUAL(returned_result[0], 0);
        KRATOS_CHECK_EQUAL(returned_result[1], -(world_size-1));
    }

    #ifdef KRATOS_DEBUG
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
    std::vector<int> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Min(local_vector_wrong_size, output, root),"Error:");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinDoubleVector, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMaxInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    int local = world_rank;
    int result = mpi_world_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(result, world_size-1);
    }

    #ifdef KRATOS_DEBUG
    // passing invalid rank as argument
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Max(local, world_size),"is not a valid rank.");
    KRATOS_CHECK_EXCEPTION_IS_THROWN(
        mpi_world_communicator.Max(local, -1),"is not a valid rank.");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMaxDouble, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMaxArray1d, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();
    constexpr int root = 0;

    std::vector<int> local{world_rank, -world_rank};
    std::vector<int> output{-1, -1};

    // two-buffer version
    mpi_world_communicator.Max(local, output, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(output[0], world_size-1);
        KRATOS_CHECK_EQUAL(output[1], 0);
    }

    // return buffer version
    std::vector<int> returned_result = mpi_world_communicator.Max(local, root);
    if (world_rank == root)
    {
        KRATOS_CHECK_EQUAL(returned_result.size(), 2);
        KRATOS_CHECK_EQUAL(returned_result[0], world_size-1);
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
    std::vector<int> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.Max(local_vector_wrong_size, output, root),"Error:");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxDoubleVector, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorSumAllInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();

    int local = 1;
    int result = mpi_world_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(result, world_size);
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorSumAllDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();

    double local = 2.0;
    double result = mpi_world_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(result, 2.0*world_size);
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorSumAllArray1d, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_size = mpi_world_communicator.Size();

    std::vector<int> local{1, 1};
    std::vector<int> output{-1, -1};

    // two-buffer version
    mpi_world_communicator.SumAll(local, output);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(output[i], world_size);
    }

    // return buffer version
    std::vector<int> returned_result = mpi_world_communicator.SumAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(returned_result[i], world_size);
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
    std::vector<int> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.SumAll(local_vector_wrong_size, output),"Input error in call to MPI_Allreduce");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorSumAllDoubleVector, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMinAllInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();

    int local = world_rank;
    int result = mpi_world_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(result, 0);
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMinAllDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();

    double local = 2.0*world_rank;
    double result = mpi_world_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(result, 0.0);
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMinAllArray1d, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    std::vector<int> local{world_rank, -world_rank};
    std::vector<int> output{-1, -1};

    // two-buffer version
    mpi_world_communicator.MinAll(local, output);
    KRATOS_CHECK_EQUAL(output[0], 0);
    KRATOS_CHECK_EQUAL(output[1], -(world_size-1));

    // return buffer version
    std::vector<int> returned_result = mpi_world_communicator.MinAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    KRATOS_CHECK_EQUAL(returned_result[0], 0);
    KRATOS_CHECK_EQUAL(returned_result[1], -(world_size-1));

    #ifdef KRATOS_DEBUG
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
    std::vector<int> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.MinAll(local_vector_wrong_size, output),"Input error in call to MPI_Allreduce");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorMinAllDoubleVector, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMaxAllInt, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    int local = world_rank;
    int result = mpi_world_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(result, world_size-1);
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMaxAllDouble, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    double local = 2.0*world_rank;
    double result = mpi_world_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(result, 2.0*(world_size-1));
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommuniactorMaxAllArray1d, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllIntVector, KratosMPICoreFastSuite)
{
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);
    const int world_rank = mpi_world_communicator.Rank();
    const int world_size = mpi_world_communicator.Size();

    std::vector<int> local{world_rank, -world_rank};
    std::vector<int> output{-1, -1};

    // two-buffer version
    mpi_world_communicator.MaxAll(local, output);
    KRATOS_CHECK_EQUAL(output[0], world_size-1);
    KRATOS_CHECK_EQUAL(output[1], 0);

    // return buffer version
    std::vector<int> returned_result = mpi_world_communicator.MaxAll(local);
    KRATOS_CHECK_EQUAL(returned_result.size(), 2);
    KRATOS_CHECK_EQUAL(returned_result[0], world_size-1);
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
    std::vector<int> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.MaxAll(local_vector_wrong_size, output),"Input error in call to MPI_Allreduce");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorMaxAllDoubleVector, KratosMPICoreFastSuite)
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


KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorScanSum, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorScanSumVector, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

    std::vector<int> local_total_int{1,1};
    std::vector<int> partial_sum_int{-1,-1};
    std::vector<double> local_total_double{2.0, 2.0};
    std::vector<double> partial_sum_double{-1.0,-1.0};

    // local version: do nothing
    serial_communicator.ScanSum(local_total_int, partial_sum_int);
    serial_communicator.ScanSum(local_total_double, partial_sum_double);
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(partial_sum_int[i], -1);
        KRATOS_CHECK_EQUAL(partial_sum_double[i], -1.0);
    }

    // MPI version
    mpi_world_communicator.ScanSum(local_total_int, partial_sum_int);
    mpi_world_communicator.ScanSum(local_total_double, partial_sum_double);
    int mpi_world_rank = mpi_world_communicator.Rank();
    for (int i = 0; i < 2; i++)
    {
        KRATOS_CHECK_EQUAL(partial_sum_int[i], mpi_world_rank + 1);
        KRATOS_CHECK_EQUAL(partial_sum_double[i], 2.0*(mpi_world_rank + 1));
    }

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the inputs has a different size
        if (mpi_world_rank == 0) {
            local_total_int.resize(3);
            local_total_int = {1,2,3};
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.ScanSum(local_total_int, partial_sum_int),"Input error in call to MPI_Scan");
    }
    // Input size != output size
    std::vector<int> local_vector_wrong_size{1,2,3};
    KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.ScanSum(local_vector_wrong_size, partial_sum_int),"Input error in call to MPI_Scan");
    #endif
}

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorSendRecv, KratosMPICoreFastSuite)
{
    DataCommunicator serial_communicator;
    MPIDataCommunicator mpi_world_communicator(MPI_COMM_WORLD);

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

    #ifdef KRATOS_DEBUG
    if (mpi_world_communicator.Size() > 1)
    {
        // One of the ranks has the wrong source/destination
        int wrong_send_rank = send_rank;
        if (mpi_world_communicator.Rank() == 0) {
            wrong_send_rank = 2;
        }
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.SendRecv(send_buffer_int, wrong_send_rank, recv_buffer_int, recv_rank),"Error:");

        // Input size != output size
        std::vector<int> local_vector_wrong_size{1,2,3};
        KRATOS_CHECK_EXCEPTION_IS_THROWN(mpi_world_communicator.SendRecv(local_vector_wrong_size, send_rank, recv_buffer_int, recv_rank),"Input error in call to MPI_Sendrecv");
    }
    #endif
}


KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorBroadcast, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorBroadcastVector, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorScatter, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorScatterv, KratosMPICoreFastSuite)
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


KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorGather, KratosMPICoreFastSuite)
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


KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorGatherv, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorAllGather, KratosMPICoreFastSuite)
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

KRATOS_TEST_CASE_IN_SUITE(MPIDataCommunicatorErrorBroadcasting, KratosMPICoreFastSuite)
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