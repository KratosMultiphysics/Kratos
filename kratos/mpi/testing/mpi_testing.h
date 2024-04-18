//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//
//

#pragma once

// System includes

// External includes
#include <mpi.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// Project includes
#include "testing/testing.h"
#include "includes/parallel_environment.h"

// Parallel Extension
#include "mpi/includes/mpi_expect.h"
#include "mpi/includes/mpi_communicator.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"

namespace Kratos::Testing 
{

/*
 * This Fixture creates a new kernel instance for kratos, so the test is able to interact with the database.
 * Its called this way to that all tests belong to a existing kernel fixture
*/
class KratosMpiTestEnv : public ::testing::Environment 
{
    protected:
        ~KratosMpiTestEnv() override {}

        void SetUp() override {
            // Define the World DataCommunicator as a wrapper for MPI_COMM_WORLD and make it the default.
            ParallelEnvironment::RegisterDataCommunicator("World", MPIDataCommunicator::Create(MPI_COMM_WORLD), ParallelEnvironment::MakeDefault);

            // Register the MPICommunicator to be used as factory for the communicator.
            ParallelEnvironment::RegisterCommunicatorFactory<const std::string>([](ModelPart& rModelPart, const std::string& rDataCommunicatorName) -> Communicator::UniquePointer {
                KRATOS_ERROR_IF_NOT(ParallelEnvironment::HasDataCommunicator(rDataCommunicatorName)) << "Asking for an unregistered \'" << rDataCommunicatorName <<  "\' data communicator." << std::endl;
                const auto& r_data_communicator = ParallelEnvironment::GetDataCommunicator(rDataCommunicatorName);
                KRATOS_ERROR_IF_NOT(r_data_communicator.IsDistributed()) << "Trying to create an MPI communicator with the non-distributed \'" << rDataCommunicatorName << "\'. data communicator" << std::endl;
                return Kratos::make_unique<MPICommunicator>(&(rModelPart.GetNodalSolutionStepVariablesList()), r_data_communicator);
            });

            // Register the MPICommunicator to be used as factory for the communicator.
            ParallelEnvironment::RegisterCommunicatorFactory<const DataCommunicator>([](ModelPart& rModelPart, const DataCommunicator& rDataCommunicator) -> Communicator::UniquePointer {
                KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Trying to create an MPI communicator with a non-distributed data communicator." << std::endl;
                return Kratos::make_unique<MPICommunicator>(&(rModelPart.GetNodalSolutionStepVariablesList()), rDataCommunicator);
            });

            // Register the ParallelFillCommunicator to be used as factory for the parallel communicators fill.
            ParallelEnvironment::RegisterFillCommunicatorFactory<const std::string>([](ModelPart& rModelPart, const std::string& rDataCommunicatorName) -> FillCommunicator::Pointer {
                KRATOS_ERROR_IF_NOT(ParallelEnvironment::HasDataCommunicator(rDataCommunicatorName)) << "Asking for an unregistered \'" << rDataCommunicatorName <<  "\' data communicator." << std::endl;
                const auto& r_data_communicator = ParallelEnvironment::GetDataCommunicator(rDataCommunicatorName);
                KRATOS_ERROR_IF_NOT(r_data_communicator.IsDistributed()) << "Trying to create an MPI communicator with the non-distributed \'" << rDataCommunicatorName << "\'. data communicator" << std::endl;
                return Kratos::make_shared<ParallelFillCommunicator>(rModelPart, r_data_communicator);
            });

            // Register the ParallelFillCommunicator to be used as factory for the parallel communicators fill.
            ParallelEnvironment::RegisterFillCommunicatorFactory<const DataCommunicator>([](ModelPart& rModelPart, const DataCommunicator& rDataCommunicator) -> FillCommunicator::Pointer {
                KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Trying to create an MPI communicator with a non-distributed data communicator." << std::endl;
                return Kratos::make_shared<ParallelFillCommunicator>(rModelPart, rDataCommunicator);
            });
        }

        void TearDown() override {
            int has_failure = ::testing::Test::HasFailure();

            // Synchronize the failre status
            MPI_Allreduce(MPI_IN_PLACE, &has_failure, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

            // If any of the processes has issued a failure, fail the whole tests
            if (has_failure) {
                KRATOS_FAIL();
            }
        }
};

/*
 * Suite for the mpi testing environment (mKernel(true))
*/
class KratosMPICoreFastSuite : public ::testing::Test 
{
    protected:
        KratosMPICoreFastSuite(): mKernel(true) {}
        ~KratosMPICoreFastSuite() {}

    	Kratos::Kernel mKernel;
};

/*
 * Initializes the parallel testing environment. This is usefull for other tests depending on a parallel environment.
*/
class MPIGTestMain {
    public:
        static int InitializeMPIKernel(int argc, char* argv[]) {
            // Initialize MPI
            // int err = MPI_Init(&argc, &argv);
            MPI_Init(&argc, &argv);

            // Initialize the tests
            ::testing::InitGoogleTest(&argc, argv);

            // Get the size and rank 
            // TODO: we should ask this from the KratosMpiTestEnv, not from MPI directly, but for now it will have to suffice.
            int rank = 0;
            int size = 0;

            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);

            // Remove the default listener
            testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
            auto default_printer = listeners.Release(listeners.default_result_printer());

            // Create a configurable listener
            Kratos::Testing::ConfigurableEventListener *listener = new Kratos::Testing::ConfigurableEventListener(default_printer);

            // Set the listener configuration (by default false to all if rank != 0)
            if (rank != 0) {
                listener->showStart = false;
                listener->showIterations = false;
                listener->showEnvironment = false;
                listener->showTestCases = false;
                listener->showTestNames = false;
                listener->showSuccesses = false;
                listener->showInlineFailures = false;
                listener->showResult = false;
                listener->showEnd = false;
            }

            // Add the MPI environment to the test 
            ::testing::AddGlobalTestEnvironment(new Kratos::Testing::KratosMpiTestEnv);

            // Add our listener
            listeners.Append(listener);

            // Finalize MPI
            MPI_Finalize();

            // Run the tests
            return RUN_ALL_TESTS();
        }
};

} // namespace Kratos::Testing
