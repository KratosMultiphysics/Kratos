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
 * ConfigurableEventListener provides a configurable event listener for the test output. 
 * In Kratos this is used to remove the output from the tests
 * not executed by the main rank(0)
 * Inspiration from: From: https://gist.github.com/elliotchance/8215283
*/
class ConfigurableEventListener : public ::testing::TestEventListener
{
    protected:
        ::testing::TestEventListener * eventListener;
        
    public:

        /**
        * Swtiches to enable or disable the output of the test. Separated in different sections:
        * - showStart: Show the start of the test program.
        * - showIterations: Show the start of an iteration of the test program.
        * - showTestCases: Show the names of each test case.
        * - showTestNames: Show the names of each test.
        * - showSuccesses: Show each success.
        * - showInlineFailures: Show each failure as it occurs. You will also see it at the bottom after the full suite is run.
        * - showEnvironment: Show the setup of the global environment.
        * - showResult: Show the results of the test program.
        * - showEnd: Show the end of the test program.
        */

        /// Show the start of the test program. 
        bool showStart;

        /// Show the start of an iteration of the test program. 
        bool showIterations;

        /// Show the start of the test program. 
        bool showTestCases;
        
        /// Show the start of the test program. 
        bool showTestNames;
        
        /// Show the start of the test program. 
        bool showSuccesses;
        
        /// Show the start of the test program. 
        bool showInlineFailures;
        
        /// Show the start of the test program. 
        bool showEnvironment;

        /// Show the start of the test program. 
        bool showResult;

        /// Show the start of the test program. 
        bool showEnd;
        
        explicit ConfigurableEventListener(::testing::TestEventListener* theEventListener) : eventListener(theEventListener)
        {
            showStart = true;
            showIterations = true;
            showTestCases = true;
            showTestNames = true;
            showSuccesses = true;
            showInlineFailures = true;
            showEnvironment = true;
            showResult = true;
            showEnd = true;
        }
        
        virtual ~ConfigurableEventListener()
        {
            delete eventListener;
        }
        
        virtual void OnTestProgramStart(const ::testing::UnitTest& unit_test) override
        {
            if(showStart) {
                eventListener->OnTestProgramStart(unit_test);
            }
        }
        
        virtual void OnTestIterationStart(const ::testing::UnitTest& unit_test, int iteration) override
        {
            if(showIterations) {
                eventListener->OnTestIterationStart(unit_test, iteration);
            }
        }
        
        virtual void OnEnvironmentsSetUpStart(const ::testing::UnitTest& unit_test) override
        {
            if(showEnvironment) {
                eventListener->OnEnvironmentsSetUpStart(unit_test);
            }
        }
        
        virtual void OnEnvironmentsSetUpEnd(const ::testing::UnitTest& unit_test) override
        {
            if(showEnvironment) {
                eventListener->OnEnvironmentsSetUpEnd(unit_test);
            }
        }
        
        virtual void OnTestCaseStart(const ::testing::TestCase& test_case) override
        {
            if(showTestCases) {
                eventListener->OnTestCaseStart(test_case);
            }
        }
        
        virtual void OnTestStart(const ::testing::TestInfo& test_info) override
        {
            if(showTestNames) {
                eventListener->OnTestStart(test_info);
            }
        }
        
        virtual void OnTestPartResult(const ::testing::TestPartResult& result) override
        {
            if(showResult) {
                eventListener->OnTestPartResult(result);
            } 
        }
        
        virtual void OnTestEnd(const ::testing::TestInfo& test_info) override
        {
            if((showInlineFailures && test_info.result()->Failed()) || (showSuccesses && !test_info.result()->Failed())) {
                eventListener->OnTestEnd(test_info);
            }
        }
        
        virtual void OnTestCaseEnd(const ::testing::TestCase& test_case) override
        {
            if(showTestCases) {
                eventListener->OnTestCaseEnd(test_case);
            }
        }
        
        virtual void OnEnvironmentsTearDownStart(const ::testing::UnitTest& unit_test) override
        {
            if(showEnvironment) {
                eventListener->OnEnvironmentsTearDownStart(unit_test);
            }
        }
        
        virtual void OnEnvironmentsTearDownEnd(const ::testing::UnitTest& unit_test) override
        {
            if(showEnvironment) {
                eventListener->OnEnvironmentsTearDownEnd(unit_test);
            }
        }
        
        virtual void OnTestIterationEnd(const ::testing::UnitTest& unit_test, int iteration) override
        {
            if(showIterations) {
                eventListener->OnTestIterationEnd(unit_test, iteration);
            }
        }
        
        virtual void OnTestProgramEnd(const ::testing::UnitTest& unit_test) override
        {
            if(showEnd) {
                eventListener->OnTestProgramEnd(unit_test);
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
class GTestMain {
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

            std::cout << "Initializing GTEST MPI environment with " << size << " ranks (" << rank << ")" << std::endl;
            
            // Add the MPI environment to the test 
            ::testing::AddGlobalTestEnvironment(new Kratos::Testing::KratosMpiTestEnv);

            // Add our listener
            listeners.Append(listener);

            // Run the tests
            return RUN_ALL_TESTS();

            // Finalize MPI
            // err = MPI_Finalize();
            MPI_Finalize();
        }
};

}
