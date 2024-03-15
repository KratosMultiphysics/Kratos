//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// External includes

// Project includes

// Application includes
#include "custom_utilities/filtering/explicit_damping_utils.h"
#include "custom_utilities/filtering/explicit_filter_utils.h"

// Include base h
#include "add_custom_filters_to_python.h"

namespace Kratos {
namespace Python {

namespace Detail
{

template <class TContainerType>
void AddExplicitFilterUtils(
    pybind11::module& m,
    const std::string& rName)
{
    namespace py = pybind11;

    py::class_<ExplicitFilterUtils<TContainerType>, typename ExplicitFilterUtils<TContainerType>::Pointer>(m, rName.c_str())
        .def(py::init<const ModelPart&, const std::string&, const std::size_t, const std::size_t>(), py::arg("model_part"), py::arg("kernel_function_type"), py::arg("max_number_of_neighbours"), py::arg("echo_level"))
        .def("SetFilterRadius", &ExplicitFilterUtils<TContainerType>::SetFilterRadius, py::arg("filter_radius"))
        .def("SetDampingCoefficients", &ExplicitFilterUtils<TContainerType>::SetDampingCoefficients, py::arg("damping_coefficients"))
        .def("ForwardFilterField", &ExplicitFilterUtils<TContainerType>::ForwardFilterField, py::arg("mesh_independent_control_space_field"))
        .def("BackwardFilterField", &ExplicitFilterUtils<TContainerType>::BackwardFilterField, py::arg("physical_space_mesh_independent_gradient"))
        .def("BackwardFilterIntegratedField", &ExplicitFilterUtils<TContainerType>::BackwardFilterIntegratedField, py::arg("physical_space_mesh_dependent_gradient"))
        .def("GetFilterRadius", &ExplicitFilterUtils<TContainerType>::GetFilterRadius)
        .def("GetDampingCoefficients", &ExplicitFilterUtils<TContainerType>::GetDampingCoefficients)
        .def("GetIntegrationWeights", &ExplicitFilterUtils<TContainerType>::GetIntegrationWeights, py::arg("integration_weight_field"))
        .def("Update", &ExplicitFilterUtils<TContainerType>::Update)
        .def("__str__", &ExplicitFilterUtils<TContainerType>::Info)
        ;
}

} // namespace Detail

void AddCustomFiltersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto explicit_damping_utils = m.def_submodule("ExplicitDampingUtils");
    explicit_damping_utils.def("ComputeDampingCoefficientsBasedOnNearestEntity", &ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntity<ModelPart::NodesContainerType>,
        py::arg("damping_radius_nodal_expression"),
        py::arg("lists_of_model_parts_for_components"),
        py::arg("shape_of_the_damping_variable"),
        py::arg("damping_function_type"),
        py::arg("bucket_size"));
    explicit_damping_utils.def("ComputeDampingCoefficientsBasedOnNearestEntity", &ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntity<ModelPart::ConditionsContainerType>,
        py::arg("damping_radius_condition_expression"),
        py::arg("lists_of_model_parts_for_components"),
        py::arg("shape_of_the_damping_variable"),
        py::arg("damping_function_type"),
        py::arg("bucket_size"));
    explicit_damping_utils.def("ComputeDampingCoefficientsBasedOnNearestEntity", &ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntity<ModelPart::ElementsContainerType>,
        py::arg("damping_radius_element_expression"),
        py::arg("lists_of_model_parts_for_components"),
        py::arg("shape_of_the_damping_variable"),
        py::arg("damping_function_type"),
        py::arg("bucket_size"));

    Detail::AddExplicitFilterUtils<ModelPart::NodesContainerType>(m, "NodalExplicitFilterUtils");
    Detail::AddExplicitFilterUtils<ModelPart::ConditionsContainerType>(m, "ConditionExplicitFilterUtils");
    Detail::AddExplicitFilterUtils<ModelPart::ElementsContainerType>(m, "ElementExplicitFilterUtils");
}

} // namespace Python.
} // namespace Kratos.
