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
#include "custom_utilities/response/sensor_distance_summation_response_utils.h"
#include "custom_utilities/response/sensor_isolation_response_utils.h"

// Include base h
#include "custom_python/add_custom_response_utils_to_python.h"

namespace Kratos::Python {

void AddCustomResponseUtilsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    m.def_submodule("SensorDistanceSummationResponseUtils")
        .def("CalculateValue", &SensorDistanceSummationResponseUtils::CalculateValue, py::arg("model_part"), py::arg("distance_matrix"))
        .def("CalculateGradient", &SensorDistanceSummationResponseUtils::CalculateGradient, py::arg("model_part"), py::arg("distance_matrix"))
        ;

    m.def_submodule("SensorIsolationResponseUtils")
        .def("CalculateValue", &SensorIsolationResponseUtils::CalculateValue, py::arg("model_part"), py::arg("radius"), py::arg("distance_matrix"))
        .def("CalculateGradient", &SensorIsolationResponseUtils::CalculateGradient, py::arg("model_part"), py::arg("radius"), py::arg("distance_matrix"))
        ;
}

} // namespace Kratos::Python
