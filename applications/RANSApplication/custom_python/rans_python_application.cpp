//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya (https://github.com/sunethwarna)
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "rans_application.h"
#include "rans_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosRANSApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosRANSApplication,
        KratosRANSApplication::Pointer,
        KratosApplication>(m, "KratosRANSApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    // incompressible potential flow specific variables
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VELOCITY_POTENTIAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_POTENTIAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_IS_INLET )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_IS_OUTLET )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RANS_IS_STRUCTURE )
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
