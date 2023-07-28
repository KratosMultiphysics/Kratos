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
            int argc = 0;
            char** argv;

            int err = MPI_Init(&argc, &argv);

            // Define the World DataCommunicator as a wrapper for MPI_COMM_WORLD and make it the default.
            Kratos::ParallelEnvironment::RegisterDataCommunicator("World", MPIDataCommunicator::Create(MPI_COMM_WORLD), ParallelEnvironment::MakeDefault);

            // Register the MPICommunicator to be used as factory for the communicator.
            Kratos::ParallelEnvironment::RegisterCommunicatorFactory<const std::string>([](ModelPart& rModelPart, const std::string& rDataCommunicatorName)->Communicator::UniquePointer {
                KRATOS_ERROR_IF_NOT(ParallelEnvironment::HasDataCommunicator(rDataCommunicatorName)) << "Asking for an unregistered \'" << rDataCommunicatorName <<  "\' data communicator." << std::endl;
                const auto& r_data_communicator = ParallelEnvironment::GetDataCommunicator(rDataCommunicatorName);
                KRATOS_ERROR_IF_NOT(r_data_communicator.IsDistributed()) << "Trying to create an MPI communicator with the non-distributed \'" << rDataCommunicatorName << "\'. data communicator" << std::endl;
                return Kratos::make_unique<MPICommunicator>(&(rModelPart.GetNodalSolutionStepVariablesList()), r_data_communicator);
            });
            // Register the MPICommunicator to be used as factory for the communicator.
            Kratos::ParallelEnvironment::RegisterCommunicatorFactory<const DataCommunicator>([](ModelPart& rModelPart, const DataCommunicator& rDataCommunicator)->Communicator::UniquePointer {
                KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Trying to create an MPI communicator with a non-distributed data communicator." << std::endl;
                return Kratos::make_unique<MPICommunicator>(&(rModelPart.GetNodalSolutionStepVariablesList()), rDataCommunicator);
            });

            // Register the ParallelFillCommunicator to be used as factory for the parallel communicators fill.
            Kratos::ParallelEnvironment::RegisterFillCommunicatorFactory<const std::string>([](ModelPart& rModelPart, const std::string& rDataCommunicatorName)->FillCommunicator::Pointer {
                KRATOS_ERROR_IF_NOT(ParallelEnvironment::HasDataCommunicator(rDataCommunicatorName)) << "Asking for an unregistered \'" << rDataCommunicatorName <<  "\' data communicator." << std::endl;
                const auto& r_data_communicator = ParallelEnvironment::GetDataCommunicator(rDataCommunicatorName);
                KRATOS_ERROR_IF_NOT(r_data_communicator.IsDistributed()) << "Trying to create an MPI communicator with the non-distributed \'" << rDataCommunicatorName << "\'. data communicator" << std::endl;
                return Kratos::make_shared<ParallelFillCommunicator>(rModelPart, r_data_communicator);
            });
            // Register the ParallelFillCommunicator to be used as factory for the parallel communicators fill.
            Kratos::ParallelEnvironment::RegisterFillCommunicatorFactory<const DataCommunicator>([](ModelPart& rModelPart, const DataCommunicator& rDataCommunicator)->FillCommunicator::Pointer {
                KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Trying to create an MPI communicator with a non-distributed data communicator." << std::endl;
                return Kratos::make_shared<ParallelFillCommunicator>(rModelPart, rDataCommunicator);
            });

            KRATOS_EXPECT_FALSE(err);
        }

        void TearDown() override {
            int err = MPI_Finalize();

            KRATOS_EXPECT_FALSE(err);
        }
};

class KratosMPICoreFastSuite : public ::testing::Test 
{
    protected:
        KratosMPICoreFastSuite(): mKernel(true) {} // 
        ~KratosMPICoreFastSuite() {}

        Kratos::Kernel mKernel;
};

}
