//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//                   Ihar Antonau
//                   Fabian Meister
//

// System includes

// External includes
#include <pybind11/stl.h>

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_utilities/response/sensor_localization_response_utils.h"

// Include base h
#include "custom_python/add_custom_response_utilities_to_python.h"

namespace Kratos::Python {

void AddCustomResponseUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<SensorLocalizationResponseUtils, SensorLocalizationResponseUtils::Pointer>(m, "SensorLocalizationResponseUtils")
        .def(py::init<ModelPart&, const double>(), py::arg("sensor_model_part"), py::arg("p_coefficient"))
        .def("CalculateValue", &SensorLocalizationResponseUtils::CalculateValue, py::arg("list_of_masks"))
        .def("CalculateGradient", &SensorLocalizationResponseUtils::CalculateGradient, py::arg("list_of_masks"))
        ;
}

} // namespace Kratos::Python
