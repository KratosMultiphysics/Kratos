//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_responses/measurement_residual_response_function.h"

// Include base h
#include "custom_python/add_custom_responses_to_python.h"

namespace Kratos::Python {

void AddCustomResponsesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto responses_module = m.def_submodule("Responses");

    py::class_<MeasurementResidualResponseFunction, MeasurementResidualResponseFunction::Pointer, AdjointResponseFunction>(responses_module, "MeasurementResidualResponseFunction")
        .def(py::init<const double>(), py::arg("p_coefficient"))
        .def("AddSensor", &MeasurementResidualResponseFunction::AddSensor, py::arg("sensor"))
        .def("Clear", &MeasurementResidualResponseFunction::Clear)
        .def("GetSensorsList", &MeasurementResidualResponseFunction::GetSensorsList)
        .def("__str__", PrintObject<MeasurementResidualResponseFunction>)
        ;

}

} // namespace Kratos::Python
