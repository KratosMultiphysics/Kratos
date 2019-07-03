//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Philipp Bucher
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "co_simulation_application.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosCoSimulationApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosCoSimulationApplication,
        KratosCoSimulationApplication::Pointer,
        KratosApplication>(m, "KratosCoSimulationApplication")
        .def(py::init<>())
        ;
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
