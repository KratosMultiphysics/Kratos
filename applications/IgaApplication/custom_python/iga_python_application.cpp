/*
//  KRATOS .___  ________    _____
//         |   |/  _____/   /  _  \
//         |   /   \  ___  /  /_\  \
//         |   \    \_\  \/    |    \
//         |___|\______  /\____|__  /
//                     \/         \/  Application
//
//  License: BSD License
//           Kratos default license: kratos/license.txt
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

    AddCustomUtilitiesToPython(m);

    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, NURBS_CONTROLPOINT_WEIGHT)
}

} // namespace Python
} // namespace Kratos

#endif // defined(KRATOS_PYTHON)
