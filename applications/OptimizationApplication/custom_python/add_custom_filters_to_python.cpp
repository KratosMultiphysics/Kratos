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
#include "custom_utilities/filtering/explicit_filter_utils.h"
#include "custom_utilities/filtering/explicit_damping.h"
#include "custom_utilities/filtering/neareset_entity_explicit_damping.h"
#include "custom_utilities/filtering/integrated_neareset_entity_explicit_damping.h"

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
        .def(py::init<ModelPart&, const std::string&, const std::size_t, const std::size_t>(), py::arg("model_part"), py::arg("kernel_function_type"), py::arg("max_number_of_leaf_items"), py::arg("echo_level"))
        .def("SetRadius", &ExplicitFilterUtils<TContainerType>::SetRadius, py::arg("filter_radius"))
        .def("SetDamping", &ExplicitFilterUtils<TContainerType>::SetDamping, py::arg("damping"))
        .def("ForwardFilterField", &ExplicitFilterUtils<TContainerType>::ForwardFilterField, py::arg("mesh_independent_control_space_field"))
        .def("BackwardFilterField", &ExplicitFilterUtils<TContainerType>::BackwardFilterField, py::arg("physical_space_mesh_independent_gradient"))
        .def("BackwardFilterIntegratedField", &ExplicitFilterUtils<TContainerType>::BackwardFilterIntegratedField, py::arg("physical_space_mesh_dependent_gradient"))
        .def("GetRadius", &ExplicitFilterUtils<TContainerType>::GetRadius)
        .def("GetIntegrationWeights", &ExplicitFilterUtils<TContainerType>::GetIntegrationWeights, py::arg("integration_weight_field"))
        .def("CalculateMatrix", &ExplicitFilterUtils<TContainerType>::CalculateMatrix, py::arg("output_filtering_matrix"))
        .def("Update", &ExplicitFilterUtils<TContainerType>::Update)
        .def("__str__", &ExplicitFilterUtils<TContainerType>::Info)
        ;
}

template<class TContainerType>
void AddExplicitDamping(
    pybind11::module& m,
    const std::string& rSuffix)
{
    namespace py = pybind11;

    using explicit_damping_type = ExplicitDamping<TContainerType>;
    py::class_<explicit_damping_type, typename explicit_damping_type::Pointer>(m, (rSuffix + "ExplicitDamping").c_str())
        .def("SetRadius", &explicit_damping_type::SetRadius, py::arg("radius_expression"))
        .def("GetRadius", &explicit_damping_type::GetRadius)
        .def("GetStride", &explicit_damping_type::GetStride)
        .def("GetDampedModelParts", &explicit_damping_type::GetDampedModelParts)
        .def("CalculateMatrix", &explicit_damping_type::CalculateMatrix, py::arg("output_matrix"), py::arg("component_index"))
        .def("Update", &explicit_damping_type::Update)
        ;

    using nearest_entity_explicit_damping_type = NearestEntityExplicitDamping<TContainerType>;
    py::class_<nearest_entity_explicit_damping_type, typename nearest_entity_explicit_damping_type::Pointer, explicit_damping_type>(m, ("Nearest" + rSuffix + "ExplicitDamping").c_str())
        .def(py::init<Model&, Parameters, const IndexType>(), py::arg("model"), py::arg("parameters"), py::arg("stride"))
        ;

    using integrated_nearest_entity_explicit_damping_type = IntegratedNearestEntityExplicitDamping<TContainerType>;
    py::class_<integrated_nearest_entity_explicit_damping_type, typename integrated_nearest_entity_explicit_damping_type::Pointer, explicit_damping_type>(m, ("IntegratedNearest" + rSuffix + "ExplicitDamping").c_str())
        .def(py::init<Model&, Parameters, const IndexType>(), py::arg("model"), py::arg("parameters"), py::arg("stride"))
        ;
}


} // namespace Detail

void AddCustomFiltersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // Add damping interface
    Detail::AddExplicitDamping<ModelPart::NodesContainerType>(m, "Node");
    Detail::AddExplicitDamping<ModelPart::ConditionsContainerType>(m, "Condition");
    Detail::AddExplicitDamping<ModelPart::ElementsContainerType>(m, "Element");

    Detail::AddExplicitFilterUtils<ModelPart::NodesContainerType>(m, "NodeExplicitFilterUtils");
    Detail::AddExplicitFilterUtils<ModelPart::ConditionsContainerType>(m, "ConditionExplicitFilterUtils");
    Detail::AddExplicitFilterUtils<ModelPart::ElementsContainerType>(m, "ElementExplicitFilterUtils");

}

} // namespace Python.
} // namespace Kratos.
