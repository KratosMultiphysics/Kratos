//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/parallel_environment.h"

// Parallel Extension
#include "mpi/includes/mpi_expect.h"
#include "mpi/includes/mpi_communicator.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "mpi/tests/test_utilities/mpi_test_environment.h"

namespace Kratos::Testing 
{

void KratosMpiTestEnv::SetUp() {
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

void KratosMpiTestEnv::TearDown() {
    // Create a serial communicator to avoid the destructor of the MPI communicator to be called after the MPI_Finalize function has been called.
    ParallelEnvironment::RegisterDataCommunicator("ExitComm", DataCommunicator::Create(), ParallelEnvironment::MakeDefault);

    // Unregister the World DataCommunicator
    // TODO: It should be possible to obtain a list of all data_communicators and unregister them all (see test: - PointerCommunicatorPartialPartitions)
    ParallelEnvironment::UnregisterDataCommunicator("World");
}

} // namespace Kratos::Testing