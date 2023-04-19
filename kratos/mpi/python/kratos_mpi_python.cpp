//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//                   Pooyan Dadvand
//

// External includes
#include <pybind11/pybind11.h>
#include "mpi.h"

// Module includes
#include "mpi/includes/mpi_communicator.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "add_mpi_communicator_to_python.h"
#include "add_mpi_data_communicator_to_python.h"
#include "add_mpi_utilities_to_python.h"
#include "add_mpi_debug_utilities_to_python.h"
#include "add_distributed_sparse_matrices_to_python.h"
#include "includes/parallel_environment.h"
#include "add_mpi_search_strategies_to_python.h"

namespace Kratos::Python {

void InitializeMPIParallelRun()
{
    // Define the World DataCommunicator as a wrapper for MPI_COMM_WORLD and make it the default.
    ParallelEnvironment::RegisterDataCommunicator("World", MPIDataCommunicator::Create(MPI_COMM_WORLD), ParallelEnvironment::MakeDefault);

    // Register the MPICommunicator to be used as factory for the communicator.
    ParallelEnvironment::RegisterCommunicatorFactory<const std::string>([](ModelPart& rModelPart, const std::string& rDataCommunicatorName)->Communicator::UniquePointer{
        KRATOS_ERROR_IF_NOT(ParallelEnvironment::HasDataCommunicator(rDataCommunicatorName)) << "Asking for an unregistered \'" << rDataCommunicatorName <<  "\' data communicator." << std::endl;
        const auto& r_data_communicator = ParallelEnvironment::GetDataCommunicator(rDataCommunicatorName);
        KRATOS_ERROR_IF_NOT(r_data_communicator.IsDistributed()) << "Trying to create an MPI communicator with the non-distributed \'" << rDataCommunicatorName << "\'. data communicator" << std::endl;
        return Kratos::make_unique<MPICommunicator>(&(rModelPart.GetNodalSolutionStepVariablesList()), r_data_communicator);
    });
    // Register the MPICommunicator to be used as factory for the communicator.
    ParallelEnvironment::RegisterCommunicatorFactory<const DataCommunicator>([](ModelPart& rModelPart, const DataCommunicator& rDataCommunicator)->Communicator::UniquePointer{
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Trying to create an MPI communicator with a non-distributed data communicator." << std::endl;
        return Kratos::make_unique<MPICommunicator>(&(rModelPart.GetNodalSolutionStepVariablesList()), rDataCommunicator);
    });

    // Register the ParallelFillCommunicator to be used as factory for the parallel communicators fill.
    ParallelEnvironment::RegisterFillCommunicatorFactory<const std::string>([](ModelPart& rModelPart, const std::string& rDataCommunicatorName)->FillCommunicator::Pointer{
        KRATOS_ERROR_IF_NOT(ParallelEnvironment::HasDataCommunicator(rDataCommunicatorName)) << "Asking for an unregistered \'" << rDataCommunicatorName <<  "\' data communicator." << std::endl;
        const auto& r_data_communicator = ParallelEnvironment::GetDataCommunicator(rDataCommunicatorName);
        KRATOS_ERROR_IF_NOT(r_data_communicator.IsDistributed()) << "Trying to create an MPI communicator with the non-distributed \'" << rDataCommunicatorName << "\'. data communicator" << std::endl;
        return Kratos::make_shared<ParallelFillCommunicator>(rModelPart, r_data_communicator);
    });
    // Register the ParallelFillCommunicator to be used as factory for the parallel communicators fill.
    ParallelEnvironment::RegisterFillCommunicatorFactory<const DataCommunicator>([](ModelPart& rModelPart, const DataCommunicator& rDataCommunicator)->FillCommunicator::Pointer{
        KRATOS_ERROR_IF_NOT(rDataCommunicator.IsDistributed()) << "Trying to create an MPI communicator with a non-distributed data communicator." << std::endl;
        return Kratos::make_shared<ParallelFillCommunicator>(rModelPart, rDataCommunicator);
    });
}

PYBIND11_MODULE(KratosMPI, m)
{
    namespace py = pybind11;

    m.def("InitializeMPIParallelRun",&InitializeMPIParallelRun,"Initialitze MPI and set up Kratos for a parallel run.");

    AddMPICommunicatorToPython(m);
    AddMPIDataCommunicatorToPython(m);
    AddMPIUtilitiesToPython(m);
    AddMPIDebugUtilitiesToPython(m);
    AddDistributedSparseMatricesToPython(m);
    AddMPISearchStrategiesToPython(m);
}

}
