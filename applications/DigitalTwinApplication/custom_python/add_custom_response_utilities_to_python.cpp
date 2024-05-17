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
#include "custom_utilities/response/sensor_distance_boltzmann_operator_utils.h"
#include "custom_utilities/response/sensor_cosine_distance_utils.h"
#include "custom_utilities/response/sensor_coverage_response_utils.h"
#include "custom_utilities/response/sensor_sensitivity_boltzmann_operator_response.h"

// Include base h
#include "custom_python/add_custom_response_utilities_to_python.h"

namespace Kratos::Python {

void AddCustomResponseUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<SensorLocalizationResponseUtils, SensorLocalizationResponseUtils::Pointer>(m, "SensorLocalizationResponseUtils")
        .def(py::init<SensorMaskStatusKDTree<ModelPart::ElementsContainerType>::Pointer, const double, const double>(), py::arg("sensor_mask_kd_tree"), py::arg("beta"), py::arg("allowed_dissimilarity"))
        .def("CalculateValue", &SensorLocalizationResponseUtils::CalculateValue)
        .def("CalculateGradient", &SensorLocalizationResponseUtils::CalculateGradient)
        .def("GetClusterSizes", &SensorLocalizationResponseUtils::GetClusterSizes)
        ;

    py::class_<SensorIsolationResponseUtils, SensorIsolationResponseUtils::Pointer>(m, "SensorIsolationResponseUtils")
        .def(py::init<ModelPart&, const IndexType, const double, const double>(), py::arg("sensor_model_part"), py::arg("max_number_of_neighbours"), py::arg("radius"), py::arg("beta"))
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

    py::class_<SensorDistanceBoltzmannOperatorResponseUtils, SensorDistanceBoltzmannOperatorResponseUtils::Pointer>(m, "SensorDistanceBoltzmannOperatorResponseUtils")
        .def(py::init<ModelPart&, const double, const double>(), py::arg("sensor_model_part"), py::arg("p"), py::arg("beta"))
        .def("Initialize", &SensorDistanceBoltzmannOperatorResponseUtils::Initialize)
        .def("CalculateValue", &SensorDistanceBoltzmannOperatorResponseUtils::CalculateValue)
        .def("CalculateGradient", &SensorDistanceBoltzmannOperatorResponseUtils::CalculateGradient)
        ;

    py::class_<SensorCosineDistanceResponseUtils, SensorCosineDistanceResponseUtils::Pointer>(m, "SensorCosineDistanceResponseUtils")
        .def(py::init<ModelPart&, const double>(), py::arg("sensor_model_part"), py::arg("p_coefficient"))
        .def("Initialize", &SensorCosineDistanceResponseUtils::Initialize, py::arg("masks_list"))
        .def("CalculateValue", &SensorCosineDistanceResponseUtils::CalculateValue)
        .def("CalculateGradient", &SensorCosineDistanceResponseUtils::CalculateGradient)
        ;

    m.def_submodule("SensorCoverageResponseUtils")
        .def("CalculateValue", &SensorCoverageResponseUtils::CalculateValue<ModelPart::ConditionsContainerType>, py::arg("sensor_mask_status"))
        .def("CalculateGradient", &SensorCoverageResponseUtils::CalculateGradient<ModelPart::ConditionsContainerType>, py::arg("sensor_mask_status"))
        .def("CalculateValue", &SensorCoverageResponseUtils::CalculateValue<ModelPart::ElementsContainerType>, py::arg("sensor_mask_status"))
        .def("CalculateGradient", &SensorCoverageResponseUtils::CalculateGradient<ModelPart::ElementsContainerType>, py::arg("sensor_mask_status"))
        ;

    py::class_<SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::NodesContainerType>, SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::NodesContainerType>::Pointer>(m, "SensorNodalSensitivityBoltzmannOperatorResponseUtils")
        .def(py::init<
            ModelPart&,
            const std::vector<typename ContainerExpression<ModelPart::NodesContainerType>::Pointer>&,
            const double>(),
            py::arg("sensor_model_part"),
            py::arg("sensor_sensitivity_distributions"),
            py::arg("beta"))
        .def("CalculateValue", &SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::NodesContainerType>::CalculateValue)
        .def("CalculateGradient", &SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::NodesContainerType>::CalculateGradient)
        ;

    py::class_<SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::ConditionsContainerType>, SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::ConditionsContainerType>::Pointer>(m, "SensorConditionSensitivityBoltzmannOperatorResponseUtils")
        .def(py::init<
            ModelPart&,
            const std::vector<typename ContainerExpression<ModelPart::ConditionsContainerType>::Pointer>&,
            const double>(),
            py::arg("sensor_model_part"),
            py::arg("sensor_sensitivity_distributions"),
            py::arg("beta"))
        .def("CalculateValue", &SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::ConditionsContainerType>::CalculateValue)
        .def("CalculateGradient", &SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::ConditionsContainerType>::CalculateGradient)
        ;

    py::class_<SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::ElementsContainerType>, SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::ElementsContainerType>::Pointer>(m, "SensorElementSensitivityBoltzmannOperatorResponseUtils")
        .def(py::init<
            ModelPart&,
            const std::vector<typename ContainerExpression<ModelPart::ElementsContainerType>::Pointer>&,
            const double>(),
            py::arg("sensor_model_part"),
            py::arg("sensor_sensitivity_distributions"),
            py::arg("beta"))
        .def("CalculateValue", &SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::ElementsContainerType>::CalculateValue)
        .def("CalculateGradient", &SensorSensitivityBoltzmannOperatorResponseUtils<ModelPart::ElementsContainerType>::CalculateGradient)
        ;
}

} // namespace Kratos::Python
