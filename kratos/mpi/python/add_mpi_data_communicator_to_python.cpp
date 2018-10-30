//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/data_communicator.h"
#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos {

namespace Python {

void AddMPIDataCommunicatorToPython(pybind11::module &m)
{
    namespace py = pybind11;

    py::class_<MPIDataCommunicator, MPIDataCommunicator::Pointer, DataCommunicator>(m,"MPIDataCommunicator");
}

} // namespace Python.

} // Namespace Kratos
