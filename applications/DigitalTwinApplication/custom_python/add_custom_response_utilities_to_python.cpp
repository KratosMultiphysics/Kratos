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
#include "custom_utilities/response/sensor_isolation_response_utils.h"
#include "custom_utilities/response/sensor_distance_p_norm_utils.h"
#include "custom_utilities/response/sensor_cosine_distance_utils.h"

// Include base h
#include "custom_python/add_custom_response_utilities_to_python.h"

namespace Kratos::Python {

void AddCustomResponseUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<SensorLocalizationResponseUtils, SensorLocalizationResponseUtils::Pointer>(m, "SensorLocalizationResponseUtils")
        .def(py::init<ModelPart&, const std::vector<ContainerExpression<ModelPart::ElementsContainerType>::Pointer>&, const double>(), py::arg("sensor_model_part"), py::arg("sensor_masks_list"), py::arg("p_coefficient"))
        .def("CalculateValue", &SensorLocalizationResponseUtils::CalculateValue)
        .def("CalculateGradient", &SensorLocalizationResponseUtils::CalculateGradient)
        ;

    py::class_<SensorIsolationResponseUtils, SensorIsolationResponseUtils::Pointer>(m, "SensorIsolationResponseUtils")
        .def(py::init<ModelPart&, const IndexType, const double>(), py::arg("sensor_model_part"), py::arg("max_number_of_neighbours"), py::arg("radius"))
        .def("Initialize", &SensorIsolationResponseUtils::Initialize)
        .def("CalculateValue", &SensorIsolationResponseUtils::CalculateValue)
        .def("CalculateGradient", &SensorIsolationResponseUtils::CalculateGradient)
        ;

    py::class_<SensorDistancePNormResponseUtils, SensorDistancePNormResponseUtils::Pointer>(m, "SensorDistancePNormResponseUtils")
        .def(py::init<ModelPart&, const double>(), py::arg("sensor_model_part"), py::arg("p_coefficient"))
        .def("Initialize", &SensorDistancePNormResponseUtils::Initialize)
        .def("CalculateValue", &SensorDistancePNormResponseUtils::CalculateValue)
        .def("CalculateGradient", &SensorDistancePNormResponseUtils::CalculateGradient)
        ;

    py::class_<SensorCosineDistanceResponseUtils, SensorCosineDistanceResponseUtils::Pointer>(m, "SensorCosineDistanceResponseUtils")
        .def(py::init<ModelPart&, const double>(), py::arg("sensor_model_part"), py::arg("p_coefficient"))
        .def("Initialize", &SensorCosineDistanceResponseUtils::Initialize, py::arg("masks_list"))
        .def("CalculateValue", &SensorCosineDistanceResponseUtils::CalculateValue)
        .def("CalculateGradient", &SensorCosineDistanceResponseUtils::CalculateGradient)
        ;
}

} // namespace Kratos::Python
