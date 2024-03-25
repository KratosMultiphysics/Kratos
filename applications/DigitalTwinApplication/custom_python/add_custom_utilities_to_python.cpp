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
#include "custom_utilities/control_utils.h"
#include "custom_utilities/smooth_clamper.h"
#include "custom_utilities/mask_utils.h"
#include "custom_utilities/sensor_utils.h"

// Include base h
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos::Python {

template<class TContainerType>
void AddSmoothClamper(
    pybind11::module& m,
    const std::string& rName)
{
    namespace py = pybind11;

    using smooth_clamper_type = SmoothClamper<TContainerType>;
    py::class_<smooth_clamper_type, typename smooth_clamper_type::Pointer>(m, (rName + "SmoothClamper").c_str())
        .def(py::init<const double, const double>(), py::arg("min"), py::arg("max"))
        .def("Clamp", &smooth_clamper_type::Clamp, py::arg("x_expression"))
        .def("ClampDerivative", &smooth_clamper_type::ClampDerivative, py::arg("x_expression"))
        .def("InverseClamp", &smooth_clamper_type::InverseClamp, py::arg("y_expression"))
        ;
}

template<class TContainerType>
void AddMaskUtilsToPython(
    pybind11::module& m)
{
    namespace py = pybind11;

    std::string lower_prefix;
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        lower_prefix = "nodal";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        lower_prefix = "condition";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        lower_prefix = "element";
    }

    m.def("GetMaskSize", &MaskUtils::GetMaskSize<TContainerType>, py::arg((lower_prefix + "_mask_expression").c_str()), py::arg("required_minimum_redundancy") = 1);
    m.def("GetMask", py::overload_cast<const ContainerExpression<TContainerType>&>(&MaskUtils::GetMask<TContainerType>), py::arg((lower_prefix + "_scalar_expression").c_str()));
    m.def("GetMask", py::overload_cast<const ContainerExpression<TContainerType>&, const double>(&MaskUtils::GetMask<TContainerType>), py::arg((lower_prefix + "_scalar_expression").c_str()), py::arg("threshold"));
    m.def("Union", &MaskUtils::Union<TContainerType>, py::arg((lower_prefix + "_mask_1_expression").c_str()), py::arg((lower_prefix + "_mask_2_expression").c_str()), py::arg("required_minimum_redundancy") = 1);
    m.def("Intersect", &MaskUtils::Intersect<TContainerType>, py::arg((lower_prefix + "_mask_1_expression").c_str()), py::arg((lower_prefix + "_mask_2_expression").c_str()), py::arg("required_minimum_redundancy") = 1);
    m.def("Subtract", &MaskUtils::Subtract<TContainerType>, py::arg((lower_prefix + "_mask_1_expression").c_str()), py::arg((lower_prefix + "_mask_2_expression").c_str()), py::arg("required_minimum_redundancy") = 1);
    m.def("Scale", &MaskUtils::Scale<TContainerType>, py::arg((lower_prefix + "_scalar_expression").c_str()), py::arg((lower_prefix + "_mask_expression").c_str()), py::arg("required_minimum_redundancy") = 1);
    m.def("ClusterMasks", &MaskUtils::ClusterMasks<TContainerType>, py::arg(("list_of_" + lower_prefix + "_mask_expressions").c_str()), py::arg("required_minimum_redundancy") = 1);
    m.def("GetMasksDividingReferenceMask", &MaskUtils::GetMasksDividingReferenceMask<TContainerType>, py::arg(("reference_" + lower_prefix + "_mask_expression").c_str()), py::arg(("list_of_" + lower_prefix + "_mask_expressions").c_str()), py::arg("required_minimum_redundancy") = 1);
}

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto control_utils = m.def_submodule("ControlUtils");
    control_utils.def("AssignEquivalentProperties", &ControlUtils::AssignEquivalentProperties<ModelPart::ConditionsContainerType>, py::arg("source_conditions"), py::arg("destination_conditions"));
    control_utils.def("AssignEquivalentProperties", &ControlUtils::AssignEquivalentProperties<ModelPart::ElementsContainerType>, py::arg("source_elements"), py::arg("destination_elements"));

    AddSmoothClamper<ModelPart::NodesContainerType>(m, "Node");
    AddSmoothClamper<ModelPart::ConditionsContainerType>(m, "Condition");
    AddSmoothClamper<ModelPart::ElementsContainerType>(m, "Element");

    auto mask_utils = m.def_submodule("MaskUtils");
    AddMaskUtilsToPython<ModelPart::NodesContainerType>(mask_utils);
    AddMaskUtilsToPython<ModelPart::ConditionsContainerType>(mask_utils);
    AddMaskUtilsToPython<ModelPart::ElementsContainerType>(mask_utils);

    m.def_submodule("SensorUtils")
        .def("AddSensors", &SensorUtils::AddSensors, py::arg("model_part"), py::arg("list_of_sensors"))
        .def("SetSensor", &SensorUtils::SetSensor, py::arg("node"), py::arg("sensor"))
        .def("GetSensor", &SensorUtils::GetSensor, py::arg("node"))
        .def("AssignSensorIds", &SensorUtils::AssignSensorIds, py::arg("list_of_sensors"))
        .def("GetMostDistanced", &SensorUtils::GetMostDistanced, py::arg("origin_sensors"), py::arg("test_sensors"))
        ;
}

} // namespace Kratos::Python
