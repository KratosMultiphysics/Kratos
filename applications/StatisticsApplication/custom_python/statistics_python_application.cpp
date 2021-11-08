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
#include "custom_python/add_custom_methods_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "includes/define_python.h"
#include "statistics_application.h"
#include "statistics_application_variables.h"

namespace Kratos
{
namespace Python
{
PYBIND11_MODULE(KratosStatisticsApplication, m)
{
    namespace py = pybind11;

    py::class_<KratosStatisticsApplication, KratosStatisticsApplication::Pointer, KratosApplication>(
        m, "KratosStatisticsApplication")
        .def(py::init<>());

    AddCustomUtilitiesToPython(m);
    AddCustomMethodsToPython(m);

    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VECTOR_3D_SUM)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VECTOR_3D_MEAN)
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VECTOR_3D_VARIANCE)

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VECTOR_3D_NORM)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_NORM)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_SUM)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_MEAN)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCALAR_VARIANCE)
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
