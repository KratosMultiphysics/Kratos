//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include <pybind11/pybind11.h>
#include "includes/model_part.h"

// Project includes
#include "includes/define_python.h"
#include "includes/mpi_communicator.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

void  AddTrilinosCommunicatorToPython(pybind11::module& m)
{
    py::class_<MPICommunicator,Communicator>(m,"MPICommunicator")
    .def("__str__", PrintObject<MPICommunicator>)
    ;
}
}  // namespace Python.

} // namespace Kratos

