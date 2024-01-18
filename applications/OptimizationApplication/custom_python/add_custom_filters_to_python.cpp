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
#include "custom_utilities/filtering/filter.h"
#include "custom_utilities/filtering/explicit_filter.h"
#include "custom_utilities/filtering/entity_node_entity_filter.h"

// Include base h
#include "add_custom_filters_to_python.h"

namespace Kratos {
namespace Python {

namespace Detail
{

template <class TContainerType>
void AddExplicitFilter(
    pybind11::module& m,
    const std::string& rName)
{
    namespace py = pybind11;

    py::class_<ExplicitFilter<TContainerType>, typename ExplicitFilter<TContainerType>::Pointer, Filter<TContainerType>>(m, rName.c_str())
        .def(py::init<const ModelPart&, const std::string&, const std::size_t>(), py::arg("model_part"), py::arg("kernel_function_type"), py::arg("max_number_of_neighbours"))
        .def(py::init<const ModelPart&, const ModelPart&, const std::string&, const std::string&, const std::size_t>(), py::arg("model_part"), py::arg("fixed_model_part"), py::arg("kernel_function_type"), py::arg("damping_function_type"), py::arg("max_number_of_neighbours"))
        .def("SetFilterRadius", &ExplicitFilter<TContainerType>::SetFilterRadius, py::arg("filter_radius"))
        .def("GetIntegrationWeights", &ExplicitFilter<TContainerType>::GetIntegrationWeights, py::arg("integration_weight_field"))
        ;
}

} // namespace Detail

void AddCustomFiltersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using nodal_filter_type = Filter<ModelPart::NodesContainerType>;
    py::class_<nodal_filter_type, nodal_filter_type::Pointer>(m, "NodalFilter")
        .def("FilterField", &nodal_filter_type::FilterField, py::arg("unfiltered_field"))
        .def("FilterIntegratedField", &nodal_filter_type::FilterIntegratedField, py::arg("filtered_field"))
        .def("Update", &nodal_filter_type::Update)
        .def("__str__", &nodal_filter_type::Info)
        ;

    using condition_filter_type = Filter<ModelPart::ConditionsContainerType>;
    py::class_<condition_filter_type, condition_filter_type::Pointer>(m, "ConditionFilter")
        .def("FilterField", &condition_filter_type::FilterField, py::arg("unfiltered_field"))
        .def("FilterIntegratedField", &condition_filter_type::FilterIntegratedField, py::arg("filtered_field"))
        .def("Update", &condition_filter_type::Update)
        .def("__str__", &condition_filter_type::Info)
        ;

    using element_filter_type = Filter<ModelPart::ElementsContainerType>;
    py::class_<element_filter_type, element_filter_type::Pointer>(m, "ElementFilter")
        .def("FilterField", &element_filter_type::FilterField, py::arg("unfiltered_field"))
        .def("FilterIntegratedField", &element_filter_type::FilterIntegratedField, py::arg("filtered_field"))
        .def("Update", &element_filter_type::Update)
        .def("__str__", &element_filter_type::Info)
        ;

    py::class_<EntityNodeEntityFilter<ModelPart::ConditionsContainerType>, EntityNodeEntityFilter<ModelPart::ConditionsContainerType>::Pointer, Filter<ModelPart::ConditionsContainerType>>(m, "ConditionNodeConditionFilter")
        .def(py::init<ModelPart&>(), py::arg("model_part"))
        ;

    py::class_<EntityNodeEntityFilter<ModelPart::ElementsContainerType>, EntityNodeEntityFilter<ModelPart::ElementsContainerType>::Pointer, Filter<ModelPart::ElementsContainerType>>(m, "ElementNodeElementFilter")
        .def(py::init<ModelPart&>(), py::arg("model_part"))
        ;

    Detail::AddExplicitFilter<ModelPart::NodesContainerType>(m, "NodalExplicitFilter");
    Detail::AddExplicitFilter<ModelPart::ConditionsContainerType>(m, "ConditionExplicitFilter");
    Detail::AddExplicitFilter<ModelPart::ElementsContainerType>(m, "ElementExplicitFilter");

}

} // namespace Python.
} // namespace Kratos.
