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

void AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto control_utils = m.def_submodule("ControlUtils");
    control_utils.def("AssignEquivalentProperties", &ControlUtils::AssignEquivalentProperties<ModelPart::ConditionsContainerType>, py::arg("source_conditions"), py::arg("destination_conditions"));
    control_utils.def("AssignEquivalentProperties", &ControlUtils::AssignEquivalentProperties<ModelPart::ElementsContainerType>, py::arg("source_elements"), py::arg("destination_elements"));

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
    mask_utils.def("GetMasksDividingReferenceMask", &MaskUtils::GetMasksDividingReferenceMask, py::arg("reference_mask_tensor_adaptor"), py::arg("list_of_mask_tensor_adaptors"), py::arg("required_minimum_redundancy") = 1);

    auto sensor_utils = m.def_submodule("SensorUtils");
    sensor_utils.def("IsPointInGeometry", &SensorUtils::IsPointInGeometry, py::arg("point"), py::arg("geometry"));

    py::class_<DistanceMatrix, DistanceMatrix::Pointer>(m, "DistanceMatrix")
        .def(py::init<>())
        .def("Update", &DistanceMatrix::Update, py::arg("values_container_expression"))
        .def("GetDistance", py::overload_cast<const IndexType, const IndexType>(&DistanceMatrix::GetDistance, py::const_), py::arg("index_i"), py::arg("index_j"))
        .def("GetEntriesSize", &DistanceMatrix::GetEntriesSize)
        .def("GetNumberOfItems", &DistanceMatrix::GetNumberOfItems)
        ;
}

} // namespace Kratos::Python
