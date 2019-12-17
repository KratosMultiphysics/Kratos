//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:   Raul Bravo
//
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "rom_application.h"
#include "rom_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosRomApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosRomApplication,
        KratosRomApplication::Pointer,
        KratosApplication>(m, "KratosRomApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, AUX_ID )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ROM_BASIS )

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
