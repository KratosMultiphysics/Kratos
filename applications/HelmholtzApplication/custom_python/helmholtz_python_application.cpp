//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "helmholtz_application.h"
#include "helmholtz_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosHelmholtzApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosHelmholtzApplication,
        KratosHelmholtzApplication::Pointer,
        KratosApplication>(m, "KratosHelmholtzApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);

    //registering variables in python
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_DIRECTION )
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, HELMHOLTZ_POISSON_RATIO )
      KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, HELMHOLTZ_VARS )
      KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, HELMHOLTZ_SOURCE )

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
