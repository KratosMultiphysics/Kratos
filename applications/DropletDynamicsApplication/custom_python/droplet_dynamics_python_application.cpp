//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "droplet_dynamics_application.h"
#include "droplet_dynamics_application_variables.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosDropletDynamicsApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosDropletDynamicsApplication,
        KratosDropletDynamicsApplication::Pointer,
        KratosApplication>(m, "KratosDropletDynamicsApplication")
        .def(py::init<>())
        ;

    AddCustomUtilitiesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,EXT_INT_FORCE)

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
