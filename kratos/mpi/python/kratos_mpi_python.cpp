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

// Module includes
#include "mpi/mpi_environment.h"
#include "add_mpi_communicator_to_python.h"
#include "add_mpi_data_communicator_to_python.h"
#include "add_mpi_utilities_to_python.h"

#include "includes/parallel_environment.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosMPI, m)
{
    namespace py = pybind11;

    MPIEnvironment& mpi_environment = MPIEnvironment::Instance();

    // Initialize MPI when loading this module
    mpi_environment.Initialize();

    // Configure Kratos::ParallelEnvironment for MPI runs on module load
    mpi_environment.InitializeKratosParallelEnvironment();

    AddMPICommunicatorToPython(m);
    AddMPIDataCommunicatorToPython(m);
    AddMPIUtilitiesToPython(m);
}

}
}