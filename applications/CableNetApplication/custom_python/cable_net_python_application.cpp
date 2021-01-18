//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "cable_net_application.h"
#include "cable_net_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosCableNetApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosCableNetApplication,
        KratosCableNetApplication::Pointer,
        KratosApplication>(m, "KratosCableNetApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    //registering variables in python

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SPRING_DEFORMATION_EMPIRICAL_POLYNOMIAL )


}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
