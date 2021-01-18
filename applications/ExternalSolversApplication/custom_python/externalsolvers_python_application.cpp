//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "externalsolvers_application.h"
#include "custom_python/add_linear_solvers_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosExternalSolversApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosExternalSolversApplication,
        KratosExternalSolversApplication::Pointer,
        KratosApplication >(m,"KratosExternalSolversApplication")
        .def(py::init<>())
        ;

    AddLinearSolversToPython(m);
}

}  // namespace Python.
}  // namespace Kratos.

#endif // KRATOS_PYTHON defined
