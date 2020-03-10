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
#include "mpi/includes/mpi_data_communicator.h"
#include "add_mpi_communicator_to_python.h"
#include "add_mpi_data_communicator_to_python.h"
#include "add_mpi_utilities_to_python.h"

#include "includes/parallel_environment.h"

namespace Kratos {
namespace Python {

void InitializeMPIParallelRun()
{
    // Define the World DataCommunicator as a wrapper for MPI_COMM_WORLD and make it the default.
    ParallelEnvironment::RegisterDataCommunicator("World", MPIDataCommunicator::Create(MPI_COMM_WORLD), ParallelEnvironment::MakeDefault);
}

PYBIND11_MODULE(KratosMPI, m)
{
    namespace py = pybind11;

    m.def("InitializeMPIParallelRun",&InitializeMPIParallelRun,"Initialitze MPI and set up Kratos for a parallel run.");

    AddMPICommunicatorToPython(m);
    AddMPIDataCommunicatorToPython(m);
    AddMPIUtilitiesToPython(m);
}

}
}