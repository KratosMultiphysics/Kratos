//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         SystemIdentificationApplication/license.txt
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
#include "custom_utilities/mask_utils.h"
#include "custom_utilities/smooth_clamper.h"
#include "custom_utilities/sensor_utils.h"
#include "custom_utilities/distance_matrix.h"

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
        .def("ProjectForward", py::overload_cast<const ContainerExpression<TContainerType>&>(&smooth_clamper_type::ProjectForward, py::const_), py::arg("x_expression"))
        .def("CalculateForwardProjectionGradient", py::overload_cast<const ContainerExpression<TContainerType>&>(&smooth_clamper_type::CalculateForwardProjectionGradient, py::const_), py::arg("x_expression"))
        .def("ProjectBackward", py::overload_cast<const ContainerExpression<TContainerType>&>(&smooth_clamper_type::ProjectBackward, py::const_), py::arg("y_expression"))
        ;
}

template<class TContainerType>
void AddMaskUtilsToPython(pybind11::module& m)
{
    namespace py = pybind11;

    std::string lower_prefix, upper_prefix;
    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        lower_prefix = "nodal";
        upper_prefix = "Nodal";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        lower_prefix = "condition";
        upper_prefix = "Condition";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        lower_prefix = "element";
        upper_prefix = "Element";
    }

    m.def("GetMaskSize", &MaskUtils::GetMaskSize<TContainerType>, py::arg((lower_prefix + "_mask_expression").c_str()), py::arg("required_minimum_redundancy") = 1);
    m.def("GetMask", py::overload_cast<const ContainerExpression<TContainerType>&>(&MaskUtils::GetMask<TContainerType>), py::arg((lower_prefix + "_scalar_expression").c_str()));
    m.def("GetMask", py::overload_cast<const ContainerExpression<TContainerType>&, const double>(&MaskUtils::GetMask<TContainerType>), py::arg((lower_prefix + "_scalar_expression").c_str()), py::arg("threshold"));
    m.def("GetMaskThreshold", &MaskUtils::GetMaskThreshold<TContainerType>, py::arg((lower_prefix + "_scalar_expression").c_str()));
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
    control_utils.def("ClipContainerExpression", &ControlUtils::ClipContainerExpression<ModelPart::NodesContainerType>, py::arg("nodal_expression"), py::arg("min"), py::arg("max"));
    control_utils.def("ClipContainerExpression", &ControlUtils::ClipContainerExpression<ModelPart::ConditionsContainerType>, py::arg("condition_expression"), py::arg("min"), py::arg("max"));
    control_utils.def("ClipContainerExpression", &ControlUtils::ClipContainerExpression<ModelPart::ElementsContainerType>, py::arg("element_expression"), py::arg("min"), py::arg("max"));

    AddSmoothClamper<ModelPart::NodesContainerType>(m, "Node");
    AddSmoothClamper<ModelPart::ConditionsContainerType>(m, "Condition");
    AddSmoothClamper<ModelPart::ElementsContainerType>(m, "Element");

    auto mask_utils = m.def_submodule("MaskUtils");
    AddMaskUtilsToPython<ModelPart::NodesContainerType>(mask_utils);
    AddMaskUtilsToPython<ModelPart::ConditionsContainerType>(mask_utils);
    AddMaskUtilsToPython<ModelPart::ElementsContainerType>(mask_utils);

    auto sensor_utils = m.def_submodule("SensorUtils");
    sensor_utils.def("IsPointInGeometry", &SensorUtils::IsPointInGeometry, py::arg("point"), py::arg("geometry"));
    sensor_utils.def("CreateSensorView", &SensorUtils::CreateSensorView, py::arg("sensor"), py::arg("expression_name"));

    py::class_<DistanceMatrix, DistanceMatrix::Pointer>(m, "DistanceMatrix")
        .def(py::init<>())
        .def("Update", &DistanceMatrix::Update, py::arg("values_container_expression"))
        .def("GetDistance", py::overload_cast<const IndexType, const IndexType>(&DistanceMatrix::GetDistance, py::const_), py::arg("index_i"), py::arg("index_j"))
        .def("GetEntriesSize", &DistanceMatrix::GetEntriesSize)
        .def("GetNumberOfItems", &DistanceMatrix::GetNumberOfItems)
        ;
}

} // namespace Kratos::Python
