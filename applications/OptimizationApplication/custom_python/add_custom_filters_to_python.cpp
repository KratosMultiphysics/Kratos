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
        .def(py::init<const ModelPart&, const std::string&, const std::size_t>(), py::arg("model_part"), py::arg("kernel_function_type"), py::arg("max_number_of_neighbours"))
        .def("SetFilterRadius", &ExplicitFilterUtils<TContainerType>::SetFilterRadius, py::arg("filter_radius"))
        .def("SetDampingCoefficients", &ExplicitFilterUtils<TContainerType>::SetDampingCoefficients, py::arg("damping_coefficients"))
        .def("FilterField", &ExplicitFilterUtils<TContainerType>::FilterField, py::arg("unfiltered_field"))
        .def("FilterIntegratedField", &ExplicitFilterUtils<TContainerType>::FilterIntegratedField, py::arg("unfiltered_field"))
        .def("UnfilterField", &ExplicitFilterUtils<TContainerType>::UnfilterField, py::arg("filtered_field"))
        .def("UnfilterIntegratedField", &ExplicitFilterUtils<TContainerType>::UnfilterIntegratedField, py::arg("filtered_field"))
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

    auto filter_utils = m.def_submodule("ExplicitDampingUtils");
    filter_utils.def("GetComponentWiseDampedModelParts", &ExplicitDampingUtils::GetComponentWiseDampedModelParts,
        py::arg("model"),
        py::arg("parameters"),
        py::arg("number_of_components"),
        py::return_value_policy::reference);
    filter_utils.def("ComputeDampingCoefficientsBasedOnNearestEntity", &ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntity<ModelPart::NodesContainerType>,
        py::arg("damping_radius_nodal_expression"),
        py::arg("lists_of_model_parts_for_components"),
        py::arg("shape_of_the_damping_variable"),
        py::arg("damping_function_type"),
        py::arg("bucket_size"));
    filter_utils.def("ComputeDampingCoefficientsBasedOnNearestEntity", &ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntity<ModelPart::ConditionsContainerType>,
        py::arg("damping_radius_condition_expression"),
        py::arg("lists_of_model_parts_for_components"),
        py::arg("shape_of_the_damping_variable"),
        py::arg("damping_function_type"),
        py::arg("bucket_size"));
    filter_utils.def("ComputeDampingCoefficientsBasedOnNearestEntity", &ExplicitDampingUtils::ComputeDampingCoefficientsBasedOnNearestEntity<ModelPart::ElementsContainerType>,
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
