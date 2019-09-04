//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author1 Fullname
//                   Author2 Fullname
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "mor_application.h"
#include "mor_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_solvers_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosMORApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosMORApplication,
        KratosMORApplication::Pointer,
        KratosApplication>(m, "KratosMORApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomSolversToPython(m);

    //registering variables in python
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, FREQUENCY )

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
