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
#include <pybind11/stl.h>

// Project includes

// Application includes
#include "custom_utilities/filtering/filter.h"
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

    py::class_<ExplicitFilter<TContainerType>, typename ExplicitFilter<TContainerType>::Pointer, Filter<TContainerType>>(m, rName.c_str())
        .def(py::init<Model&, Parameters>(), py::arg("model"), py::arg("settings"), py::doc("help"))
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
        .def("Initialize", &nodal_filter_type::Initialize)
        .def("Check", &nodal_filter_type::Check)
        .def("Finalize", &nodal_filter_type::Finalize)
        .def("Update", &nodal_filter_type::Update)
        .def("__str__", &nodal_filter_type::Info)
        ;

    using condition_filter_type = Filter<ModelPart::ConditionsContainerType>;
    py::class_<condition_filter_type, condition_filter_type::Pointer>(m, "ConditionFilter")
        .def("FilterField", &condition_filter_type::FilterField, py::arg("unfiltered_field"))
        .def("FilterIntegratedField", &condition_filter_type::FilterIntegratedField, py::arg("filtered_field"))
        .def("Initialize", &condition_filter_type::Initialize)
        .def("Check", &condition_filter_type::Check)
        .def("Finalize", &condition_filter_type::Finalize)
        .def("Update", &condition_filter_type::Update)
        .def("__str__", &condition_filter_type::Info)
        ;

    using element_filter_type = Filter<ModelPart::ElementsContainerType>;
    py::class_<element_filter_type, element_filter_type::Pointer>(m, "ElementFilter")
        .def("FilterField", &element_filter_type::FilterField, py::arg("unfiltered_field"))
        .def("FilterIntegratedField", &element_filter_type::FilterIntegratedField, py::arg("filtered_field"))
        .def("Initialize", &element_filter_type::Initialize)
        .def("Check", &element_filter_type::Check)
        .def("Finalize", &element_filter_type::Finalize)
        .def("Update", &element_filter_type::Update)
        .def("__str__", &element_filter_type::Info)
        ;


    Detail::AddExplicitFilter<ModelPart::NodesContainerType>(m, "NodalExplicitFilter");
    Detail::AddExplicitFilter<ModelPart::ConditionsContainerType>(m, "ConditionExplicitFilter");
    Detail::AddExplicitFilter<ModelPart::ElementsContainerType>(m, "ElementExplicitFilter");

}

} // namespace Python.
} // namespace Kratos.
