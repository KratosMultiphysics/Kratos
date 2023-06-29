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
#include "custom_utilities/filtering/explicit_vertex_morphing_filter.h"

// Include base h
#include "add_custom_filters_to_python.h"

namespace Kratos {
namespace Python {

namespace Detail
{

template <class TContainerType>
void AddExplicitVertexMorphingFilter(
    pybind11::module& m,
    const std::string& rName)
{
    namespace py = pybind11;

    py::class_<ExplicitVertexMorphingFilter<TContainerType>, typename ExplicitVertexMorphingFilter<TContainerType>::Pointer>(m, rName.c_str())
        .def(py::init<const ModelPart&, const std::string&, const std::size_t>(), py::arg("model_part"), py::arg("kernel_function_type"), py::arg("max_number_of_neighbours"))
        .def("SetFilterRadius", &ExplicitVertexMorphingFilter<TContainerType>::SetFilterRadius, py::arg("filter_radius"))
        .def("FilterField", &ExplicitVertexMorphingFilter<TContainerType>::FilterField, py::arg("unfiltered_field"))
        .def("FilterIntegratedField", &ExplicitVertexMorphingFilter<TContainerType>::FilterIntegratedField, py::arg("filtered_field"))
        .def("Update", &ExplicitVertexMorphingFilter<TContainerType>::Update)
        .def("__str__", &ExplicitVertexMorphingFilter<TContainerType>::Info)
        ;
}

} // namespace Detail

void AddCustomFiltersToPython(pybind11::module& m)
{
    Detail::AddExplicitVertexMorphingFilter<ModelPart::NodesContainerType>(m, "NodalExplicitVertexMorphingFilter");
    Detail::AddExplicitVertexMorphingFilter<ModelPart::ConditionsContainerType>(m, "ConditionExplicitVertexMorphingFilter");
    Detail::AddExplicitVertexMorphingFilter<ModelPart::ElementsContainerType>(m, "ElementExplicitVertexMorphingFilter");

}

} // namespace Python.
} // namespace Kratos.
