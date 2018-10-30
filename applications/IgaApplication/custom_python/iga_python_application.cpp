/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
*/

#if defined(KRATOS_PYTHON)

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "iga_application.h"
#include "iga_application_variables.h"
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosIgaApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosIgaApplication, KratosIgaApplication::Pointer, 
        KratosApplication>(m, "KratosIgaApplication")
        .def(py::init<>())
    ;

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NURBS_CONTROL_POINT_WEIGHT)
    
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, COORDINATES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, TANGENTS)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CROSS_AREA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESTRESS_CAUCHY)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHAPE_FUNCTION_VALUES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHAPE_FUNCTION_LOCAL_DERIVATIVES)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAYLEIGH_ALPHA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAYLEIGH_BETA)

    AddCustomUtilitiesToPython(m);
}

} // namespace Python
} // namespace Kratos

#endif // defined(KRATOS_PYTHON)
