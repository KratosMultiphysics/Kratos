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

void InitializeMPIParallelRun()
{
    MPIEnvironment& mpi_environment = MPIEnvironment::Instance();

    // Initialize MPI
    mpi_environment.Initialize();

    // SetUp the default communicator to MPI comm world
    mpi_environment.InitializeKratosParallelEnvironment();
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