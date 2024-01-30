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
#include "custom_utilities/filtering/explicit_filter.h"

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

    py::class_<ExplicitFilter<TContainerType>, typename ExplicitFilter<TContainerType>::Pointer>(m, rName.c_str())
        .def(py::init<const ModelPart&, const std::string&, const std::size_t>(), py::arg("model_part"), py::arg("kernel_function_type"), py::arg("max_number_of_neighbours"))
        .def(py::init<const ModelPart&, const ModelPart&, const std::string&, const std::string&, const std::size_t>(), py::arg("model_part"), py::arg("fixed_model_part"), py::arg("kernel_function_type"), py::arg("damping_function_type"), py::arg("max_number_of_neighbours"))
        .def("SetFilterRadius", &ExplicitFilter<TContainerType>::SetFilterRadius, py::arg("filter_radius"))
        .def("FilterField", &ExplicitFilter<TContainerType>::FilterField, py::arg("unfiltered_field"))
        .def("FilterIntegratedField", &ExplicitFilter<TContainerType>::FilterIntegratedField, py::arg("filtered_field"))
        .def("GetIntegrationWeights", &ExplicitFilter<TContainerType>::GetIntegrationWeights, py::arg("integration_weight_field"))
        .def("Update", &ExplicitFilter<TContainerType>::Update)
        .def("__str__", &ExplicitFilter<TContainerType>::Info)
        ;
}

} // namespace Detail

void AddCustomFiltersToPython(pybind11::module& m)
{
    Detail::AddExplicitFilter<ModelPart::NodesContainerType>(m, "NodalExplicitFilter");
    Detail::AddExplicitFilter<ModelPart::ConditionsContainerType>(m, "ConditionExplicitFilter");
    Detail::AddExplicitFilter<ModelPart::ElementsContainerType>(m, "ElementExplicitFilter");

}

} // namespace Python.
} // namespace Kratos.
