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
#include "custom_sensors/measurement_residual_response_function.h"
#include "custom_utilities/response/sensor_coverage_response_utils.h"
#include "custom_utilities/response/sensor_localization_response_utils.h"
#include "custom_utilities/response/sensor_isolation_response_utils.h"
#include "custom_utilities/response/sensor_inverse_distance_summation_response_utils.h"

// Include base h
#include "custom_python/add_custom_response_utils_to_python.h"

namespace Kratos::Python {

void AddCustomResponseUtilsToPython(pybind11::module& m)
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

    responses_module.def_submodule("SensorCoverageResponseUtils")
        .def("CalculateValue", &SensorCoverageResponseUtils::CalculateValue, py::arg("sensor_mask_status"))
        .def("CalculateGradient", &SensorCoverageResponseUtils::CalculateGradient, py::arg("sensor_mask_status"))
        ;

    py::class_<SensorLocalizationResponseUtils, SensorLocalizationResponseUtils::Pointer>(responses_module, "SensorLocalizationResponseUtils")
        .def(py::init<SensorMaskStatusKDTree::Pointer, const double, const double, const double>(), py::arg("sensor_mask_kd_tree"), py::arg("minimum_cluster_size_ratio"), py::arg("p_coefficient"), py::arg("allowed_dissimilarity"))
        .def("CalculateValue", &SensorLocalizationResponseUtils::CalculateValue)
        .def("CalculateGradient", &SensorLocalizationResponseUtils::CalculateGradient)
        .def("GetClusterSizes", &SensorLocalizationResponseUtils::GetClusterSizes)
        ;

    responses_module.def_submodule("SensorIsolationResponseUtils")
        .def("CalculateValue", &SensorIsolationResponseUtils::CalculateValue, py::arg("model_part"), py::arg("radius"), py::arg("distance_matrix"))
        .def("CalculateGradient", &SensorIsolationResponseUtils::CalculateGradient, py::arg("model_part"), py::arg("radius"), py::arg("distance_matrix"))
        ;

    responses_module.def_submodule("SensorInverseDistanceSummationResponseUtils")
        .def("CalculateValue", &SensorInverseDistanceSummationResponseUtils::CalculateValue, py::arg("model_part"), py::arg("p"), py::arg("distance_matrix"))
        .def("CalculateGradient", &SensorInverseDistanceSummationResponseUtils::CalculateGradient, py::arg("model_part"), py::arg("p"), py::arg("distance_matrix"))
        ;
}

} // namespace Kratos::Python
