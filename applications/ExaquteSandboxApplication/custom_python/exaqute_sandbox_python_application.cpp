//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "exaqute_sandbox_application.h"
#include "exaqute_sandbox_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosExaquteSandboxApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosExaquteSandboxApplication,
        KratosExaquteSandboxApplication::Pointer,
        KratosApplication>(m, "KratosExaquteSandboxApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomProcessesToPython(m);

    //registering variables in python
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AVERAGED_DIVERGENCE )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DIVERGENCE )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITY_H1_SEMINORM )

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
