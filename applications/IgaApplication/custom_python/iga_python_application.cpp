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
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosIgaApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosIgaApplication, KratosIgaApplication::Pointer,
        KratosApplication>(m, "KratosIgaApplication")
        .def(py::init<>())
    ;

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CROSS_AREA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESTRESS_CAUCHY)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAYLEIGH_ALPHA)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, RAYLEIGH_BETA)

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, POINT_LOAD)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, LINE_LOAD)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, SURFACE_LOAD)

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, DEAD_LOAD)

    KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(m, PK2_STRESS)
    KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(m, CAUCHY_STRESS)
    KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(m, CAUCHY_STRESS_TOP)
    KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(m, CAUCHY_STRESS_BOTTOM)
    KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(m, MEMBRANE_FORCE)
    KRATOS_REGISTER_IN_PYTHON_SYMMETRIC_2D_TENSOR_VARIABLE_WITH_COMPONENTS(m, INTERNAL_MOMENT)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHEAR_FORCE_1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SHEAR_FORCE_2)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PENALTY_FACTOR)

    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);
}

} // namespace Python
} // namespace Kratos

#endif // defined(KRATOS_PYTHON)
