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
#include "custom_utilities/mask_utils.h"
#include "custom_utilities/control_utils.h"
#include "custom_utilities/smooth_clamper.h"
#include "custom_utilities/distance_matrix.h"
#include "custom_utilities/sensor_mask_status.h"
#include "custom_utilities/sensor_mask_status_kd_tree.h"

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

    py::class_<DistanceMatrix, DistanceMatrix::Pointer>(m, "DistanceMatrix")
        .def(py::init<>())
        .def("Update", &DistanceMatrix::Update, py::arg("values_container_expression"))
        .def("GetDistance", py::overload_cast<const IndexType, const IndexType>(&DistanceMatrix::GetDistance, py::const_), py::arg("index_i"), py::arg("index_j"))
        .def("GetEntriesSize", &DistanceMatrix::GetEntriesSize)
        .def("GetNumberOfItems", &DistanceMatrix::GetNumberOfItems)
        ;

    auto mask_utils = m.def_submodule("MaskUtils");
    AddMaskUtilsToPython<ModelPart::NodesContainerType>(mask_utils);
    AddMaskUtilsToPython<ModelPart::ConditionsContainerType>(mask_utils);
    AddMaskUtilsToPython<ModelPart::ElementsContainerType>(mask_utils);

    py::class_<SensorMaskStatus, SensorMaskStatus::Pointer>(m, "SensorMaskStatus")
        .def(py::init<const ModelPart&, const SensorMaskStatus::MasksListType&>(), py::arg("sensor_model_part"), py::arg("sensor_masks_list"))
        .def("GetMasks", &SensorMaskStatus::GetMasks)
        .def("GetMaskStatuses", &SensorMaskStatus::GetMaskStatuses)
        .def("GetSensorModelPart", &SensorMaskStatus::GetSensorModelPart)
        .def("Update", &SensorMaskStatus::Update)
        ;

    py::class_<SensorMaskStatusKDTree, SensorMaskStatusKDTree::Pointer>(m, "SensorMaskStatusKDTree")
        .def(py::init<SensorMaskStatus::Pointer, const IndexType>(), py::arg("sensor_mask_status"), py::arg("leaf_max_size"))
        .def("RadiusSearch", [](const SensorMaskStatusKDTree& rSelf, const Vector& rQueryPoint, const double Radius)
             {
                std::vector<nanoflann::ResultItem<unsigned int, double>> result;
                rSelf.RadiusSearch(rQueryPoint, Radius, result);

                std::vector<std::pair<unsigned int, double>> std_result;
                std_result.resize(result.size());
                std::transform(result.begin(), result.end(), std_result.begin(), [](const auto& rInput)
                               { return std::make_pair(rInput.first, rInput.second); });
                return std_result;

             }, py::arg("query_point"), py::arg("radius"))
        .def("GetSensorMaskStatus", &SensorMaskStatusKDTree::GetSensorMaskStatus)
        .def("Update", &SensorMaskStatusKDTree::Update);
}

} // namespace Kratos::Python
