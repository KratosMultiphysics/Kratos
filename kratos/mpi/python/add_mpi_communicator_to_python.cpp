//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
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
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "add_mpi_communicator_to_python.h"
#include "mpi/includes/mpi_communicator.h"

namespace Kratos {
namespace Python {

void AddMPICommunicatorToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<MPICommunicator, Communicator>(m,"MPICommunicator")
    .def("__str__", PrintObject<MPICommunicator>);
    ;
}

} // namespace Python
} // namespace Kratos

