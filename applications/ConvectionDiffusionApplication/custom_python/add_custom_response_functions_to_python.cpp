// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Jordi Cotela
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_response_functions_to_python.h"

// Response Functions
#include "custom_response_functions/local_temperature_average_response_function.h"


namespace Kratos {
namespace Python {

void  AddCustomResponseFunctionsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Response Functions
    py::class_<LocalTemperatureAverageResponseFunction, LocalTemperatureAverageResponseFunction::Pointer, AdjointResponseFunction>
        (m, "LocalTemperatureAverageResponseFunction")
        .def(py::init<Parameters, ModelPart&>());

}

}  // namespace Python.
} // Namespace Kratos