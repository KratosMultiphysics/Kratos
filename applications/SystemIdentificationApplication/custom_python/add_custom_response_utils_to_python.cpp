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
#include "custom_utilities/response/sensor_coverage_response_utils.h"

// Include base h
#include "custom_python/add_custom_response_utils_to_python.h"

namespace Kratos::Python {

void AddCustomResponseUtilsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto responses_module = m.def_submodule("Responses");

    responses_module.def_submodule("SensorCoverageResponseUtils")
        .def("CalculateValue", &SensorCoverageResponseUtils::CalculateValue, py::arg("sensor_mask_status"))
        .def("CalculateGradient", &SensorCoverageResponseUtils::CalculateGradient, py::arg("sensor_mask_status"))
        ;
}

} // namespace Kratos::Python
