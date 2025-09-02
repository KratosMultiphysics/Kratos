
// System includes

#if defined(KRATOS_PYTHON)

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"

// Application includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/uvec_utilities.h"
#include "railway_application.h"

namespace Kratos::Python
{

namespace py = pybind11;

PYBIND11_MODULE(KratosRailwayApplication, m)
{
    py::class_<KratosRailwayApplication, KratosRailwayApplication::Pointer, KratosApplication>(
        m, "KratosRailwayApplication")
        .def(py::init<>());

    AddCustomUtilitiesToPython(m);
    
}

} // namespace Kratos::Python.

#endif // KRATOS_PYTHON defined
