//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Tosi
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define.h"
#include "multilevel_monte_carlo_application.h"
#include "multilevel_monte_carlo_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_statistics_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosMultilevelMonteCarloApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosMultilevelMonteCarloApplication,
        KratosMultilevelMonteCarloApplication::Pointer,
        KratosApplication>(m, "KratosMultilevelMonteCarloApplication")
        .def(py::init<>())
        ;

    // AddCustomStrategiesToPython(m);
    // AddCustomUtilitiesToPython(m);
    AddCustomStatisticsToPython(m);

    //registering variables in python

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_SUM_1 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_SUM_2 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_SUM_3 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_SUM_4 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_SUM_5 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_SUM_6 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_SUM_7 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_SUM_8 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_SUM_9 )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, POWER_SUM_10 )
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
