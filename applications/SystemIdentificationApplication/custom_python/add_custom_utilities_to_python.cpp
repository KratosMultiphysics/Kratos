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
#include "custom_utilities/sensor_mask_status.h"
#include "custom_utilities/sensor_mask_status_kd_tree.h"
#include "custom_utilities/point_data_interpolation_utils.h"

// Include base h
#include "custom_python/add_custom_utilities_to_python.h"

namespace Kratos::Python {

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto control_utils = m.def_submodule("ControlUtils");
    control_utils.def("AssignEquivalentProperties", &ControlUtils::AssignEquivalentProperties<ModelPart::ConditionsContainerType>, py::arg("source_conditions"), py::arg("destination_conditions"));
    control_utils.def("AssignEquivalentProperties", &ControlUtils::AssignEquivalentProperties<ModelPart::ElementsContainerType>, py::arg("source_elements"), py::arg("destination_elements"));
    control_utils.def("GetIntegrationPoints", [](const ModelPart::ConditionsContainerType& rContainer){ std::vector<Point> result; ControlUtils::GetIntegrationPoints(result, rContainer); return result; }, py::arg("container"));
    control_utils.def("GetIntegrationPoints", [](const ModelPart::ElementsContainerType& rContainer){ std::vector<Point> result; ControlUtils::GetIntegrationPoints(result, rContainer); return result; }, py::arg("container"));
    control_utils.def("GetIntegrationPointAreas", [](const ModelPart::ConditionsContainerType& rContainer){ std::vector<double> result; ControlUtils::GetIntegrationPointAreas(result, rContainer); return result; }, py::arg("container"));
    control_utils.def("GetIntegrationPointAreas", [](const ModelPart::ElementsContainerType& rContainer){ std::vector<double> result; ControlUtils::GetIntegrationPointAreas(result, rContainer); return result; }, py::arg("container"));

    py::class_<SmoothClamper, typename SmoothClamper::Pointer>(m, "SmoothClamper")
        .def(py::init<const double, const double>(), py::arg("min"), py::arg("max"))
        .def("ProjectForward", py::overload_cast<const TensorAdaptor<double>&>(&SmoothClamper::ProjectForward, py::const_), py::arg("x_tensor_adaptor"))
        .def("CalculateForwardProjectionGradient", py::overload_cast<const TensorAdaptor<double>&>(&SmoothClamper::CalculateForwardProjectionGradient, py::const_), py::arg("x_tensor_adaptor"))
        .def("ProjectBackward", py::overload_cast<const TensorAdaptor<double>&>(&SmoothClamper::ProjectBackward, py::const_), py::arg("y_tensor_adaptor"))
        ;

    auto mask_utils = m.def_submodule("MaskUtils");
    mask_utils.def("GetMaskSize", &MaskUtils::GetMaskSize, py::arg("mask_tensor_adaptor"), py::arg("required_minimum_redundancy") = 1);
    mask_utils.def("GetMask", py::overload_cast<const TensorAdaptor<double>&>(&MaskUtils::GetMask), py::arg("scalar_tensor_adaptor"));
    mask_utils.def("GetMask", py::overload_cast<const TensorAdaptor<double>&, const double>(&MaskUtils::GetMask), py::arg("scalar_tensor_adaptor"), py::arg("threshold"));
    mask_utils.def("GetMaskThreshold", &MaskUtils::GetMaskThreshold, py::arg("scalar_tensor_adaptor"));
    mask_utils.def("Union", &MaskUtils::Union, py::arg("mask_1_tensor_adaptor"), py::arg("mask_2_tensor_adaptor"), py::arg("required_minimum_redundancy") = 1);
    mask_utils.def("Intersect", &MaskUtils::Intersect, py::arg("mask_1_tensor_adaptor"), py::arg("mask_2_tensor_adaptor"), py::arg("required_minimum_redundancy") = 1);
    mask_utils.def("Subtract", &MaskUtils::Subtract, py::arg("mask_1_tensor_adaptor"), py::arg("mask_2_tensor_adaptor"), py::arg("required_minimum_redundancy") = 1);
    mask_utils.def("Scale", &MaskUtils::Scale, py::arg("scalar_tensor_adaptor"), py::arg("mask_tensor_adaptor"), py::arg("required_minimum_redundancy") = 1);
    mask_utils.def("ClusterMasks", &MaskUtils::ClusterMasks, py::arg("list_of_mask_tensor_adaptors"), py::arg("required_minimum_redundancy") = 1);
    mask_utils.def("FillModelPartUsingClusterMask", &MaskUtils::FillModelPartUsingClusterMask, py::arg("model_part"), py::arg("cluster_mask"), py::arg("required_minimum_redundancy") = 1);
    mask_utils.def("GetMasksDividingReferenceMask", &MaskUtils::GetMasksDividingReferenceMask, py::arg("reference_mask_tensor_adaptor"), py::arg("list_of_mask_tensor_adaptors"), py::arg("required_minimum_redundancy") = 1);
    mask_utils.def("SmoothenMasks", &MaskUtils::SmoothenMasks, py::arg("list_of_mask_tensor_adaptors"), py::arg("radius"), py::arg("maximum_number_of_leafs") = 10);

    auto sensor_utils = m.def_submodule("SensorUtils");
    sensor_utils.def("IsPointInGeometry", &SensorUtils::IsPointInGeometry, py::arg("point"), py::arg("geometry"));

    py::class_<DistanceMatrix, DistanceMatrix::Pointer>(m, "DistanceMatrix")
        .def(py::init<>())
        .def("Update", &DistanceMatrix::Update, py::arg("values_container_tensor_adaptor"))
        .def("GetDistance", py::overload_cast<const IndexType, const IndexType>(&DistanceMatrix::GetDistance, py::const_), py::arg("index_i"), py::arg("index_j"))
        .def("GetEntriesSize", &DistanceMatrix::GetEntriesSize)
        .def("GetNumberOfItems", &DistanceMatrix::GetNumberOfItems)
        .def("__str__", PrintObject<DistanceMatrix>)
        ;

    py::class_<SensorMaskStatus, SensorMaskStatus::Pointer>(m, "SensorMaskStatus")
        .def(py::init<ModelPart&, const std::vector<TensorAdaptor<double>::Pointer>&, const IndexType>(), py::arg("sensor_model_part"), py::arg("sensor_masks_list"), py::arg("echo_level"))
        .def("GetMasks", &SensorMaskStatus::GetMasks)
        .def("GetMaskStatuses", &SensorMaskStatus::GetMaskStatuses)
        .def("GetSensorModelPart", &SensorMaskStatus::GetSensorModelPart)
        .def("GetMaskContainer", &SensorMaskStatus::pGetMaskContainer)
        .def("Update", &SensorMaskStatus::Update)
        ;

    py::class_<SensorMaskStatusKDTree, SensorMaskStatusKDTree::Pointer>(m, "SensorMaskStatusKDTree")
        .def(py::init<SensorMaskStatus::Pointer, const IndexType, const IndexType>(), py::arg("sensor_mask_status"), py::arg("leaf_max_size"), py::arg("echo_level"))
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

    using element_point_data_interpolation_utils = PointDataInterpolationUtils<ModelPart::ElementType>;
    py::class_<element_point_data_interpolation_utils, element_point_data_interpolation_utils::Pointer>(m, "ElementPointDataInterpolationUtils")
        .def(py::init<ModelPart&>(), py::arg("model_part"))
        .def("CalculateInterpolatedNodalValues", [](const element_point_data_interpolation_utils& rSelf, const Variable<double>& rVariable) { std::vector<double> output; rSelf.CalculateInterpolatedNodalValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedNodalValues", [](const element_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 3>>& rVariable) { std::vector<array_1d<double, 3>> output; rSelf.CalculateInterpolatedNodalValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedNodalValues", [](const element_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 4>>& rVariable) { std::vector<array_1d<double, 4>> output; rSelf.CalculateInterpolatedNodalValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedNodalValues", [](const element_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 6>>& rVariable) { std::vector<array_1d<double, 6>> output; rSelf.CalculateInterpolatedNodalValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedNodalValues", [](const element_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 9>>& rVariable) { std::vector<array_1d<double, 9>> output; rSelf.CalculateInterpolatedNodalValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedEntityValues", [](const element_point_data_interpolation_utils& rSelf, const Variable<double>& rVariable) { std::vector<double> output; rSelf.CalculateInterpolatedEntityValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedEntityValues", [](const element_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 3>>& rVariable) { std::vector<array_1d<double, 3>> output; rSelf.CalculateInterpolatedEntityValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedEntityValues", [](const element_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 4>>& rVariable) { std::vector<array_1d<double, 4>> output; rSelf.CalculateInterpolatedEntityValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedEntityValues", [](const element_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 6>>& rVariable) { std::vector<array_1d<double, 6>> output; rSelf.CalculateInterpolatedEntityValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedEntityValues", [](const element_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 9>>& rVariable) { std::vector<array_1d<double, 9>> output; rSelf.CalculateInterpolatedEntityValues(output, rVariable); return output;}, py::arg("variable"))
        .def("UpdatePoints", &element_point_data_interpolation_utils::UpdatePoints, py::arg("list_of_points"))
        ;

    using condition_point_data_interpolation_utils = PointDataInterpolationUtils<ModelPart::ConditionType>;
    py::class_<condition_point_data_interpolation_utils, condition_point_data_interpolation_utils::Pointer>(m, "ConditionPointDataInterpolationUtils")
        .def(py::init<ModelPart&>(), py::arg("model_part"))
        .def("CalculateInterpolatedNodalValues", [](const condition_point_data_interpolation_utils& rSelf, const Variable<double>& rVariable) { std::vector<double> output; rSelf.CalculateInterpolatedNodalValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedNodalValues", [](const condition_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 3>>& rVariable) { std::vector<array_1d<double, 3>> output; rSelf.CalculateInterpolatedNodalValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedNodalValues", [](const condition_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 4>>& rVariable) { std::vector<array_1d<double, 4>> output; rSelf.CalculateInterpolatedNodalValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedNodalValues", [](const condition_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 6>>& rVariable) { std::vector<array_1d<double, 6>> output; rSelf.CalculateInterpolatedNodalValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedNodalValues", [](const condition_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 9>>& rVariable) { std::vector<array_1d<double, 9>> output; rSelf.CalculateInterpolatedNodalValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedEntityValues", [](const condition_point_data_interpolation_utils& rSelf, const Variable<double>& rVariable) { std::vector<double> output; rSelf.CalculateInterpolatedEntityValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedEntityValues", [](const condition_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 3>>& rVariable) { std::vector<array_1d<double, 3>> output; rSelf.CalculateInterpolatedEntityValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedEntityValues", [](const condition_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 4>>& rVariable) { std::vector<array_1d<double, 4>> output; rSelf.CalculateInterpolatedEntityValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedEntityValues", [](const condition_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 6>>& rVariable) { std::vector<array_1d<double, 6>> output; rSelf.CalculateInterpolatedEntityValues(output, rVariable); return output;}, py::arg("variable"))
        .def("CalculateInterpolatedEntityValues", [](const condition_point_data_interpolation_utils& rSelf, const Variable<array_1d<double, 9>>& rVariable) { std::vector<array_1d<double, 9>> output; rSelf.CalculateInterpolatedEntityValues(output, rVariable); return output;}, py::arg("variable"))
        .def("UpdatePoints", &condition_point_data_interpolation_utils::UpdatePoints, py::arg("list_of_points"))
        ;
}

} // namespace Kratos::Python
