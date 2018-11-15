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

// System includes
#include "includes/define_python.h"
#include "includes/kratos_version.h"

// External includes
#include <pybind11/pybind11.h>

// Module includes
#include "mpi/mpi_environment.h"
#include "add_mpi_communicator_to_python.h"
#include "add_mpi_data_communicator_to_python.h"
#include "add_mpi_utilities_to_python.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosMPI, m)
{
    namespace py = pybind11;

    // Initialize MPI when loading this module
    MPIEnvironment::Initialize();
    // Configure Kratos::ParallelEnvironment for MPI runs on module load
    MPIEnvironment::InitializeKratosParallelEnvironment();

    // Define a callback to finalize MPI on module cleanup
    auto cleanup_callback = []() {
        MPIEnvironment::Finalize();
    };

    m.add_object("_cleanup", py::capsule(cleanup_callback));

    AddMPICommunicatorToPython(m);
    AddMPIDataCommunicatorToPython(m);
    AddMPIUtilitiesToPython(m);
}

}
}