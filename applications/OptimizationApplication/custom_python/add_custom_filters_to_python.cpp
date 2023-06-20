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

void AddCustomFiltersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<ExplicitVertexMorphingFilter<ModelPart::NodesContainerType>, ExplicitVertexMorphingFilter<ModelPart::NodesContainerType>::Pointer>(m, "NodalExplicitVertexMorphingFilter")
        .def(py::init<const ModelPart&, const std::string&>())
        .def("SetFilterRadius", &ExplicitVertexMorphingFilter<ModelPart::NodesContainerType>::SetFilterRadius)
        .def("FilterField", &ExplicitVertexMorphingFilter<ModelPart::NodesContainerType>::FilterField)
        .def("UnFilterField", &ExplicitVertexMorphingFilter<ModelPart::NodesContainerType>::UnFilterField)
        .def("__str__", &ExplicitVertexMorphingFilter<ModelPart::NodesContainerType>::Info)
        ;

    py::class_<ExplicitVertexMorphingFilter<ModelPart::ConditionsContainerType>, ExplicitVertexMorphingFilter<ModelPart::ConditionsContainerType>::Pointer>(m, "ConditionExplicitVertexMorphingFilter")
        .def(py::init<const ModelPart&, const std::string&>())
        .def("SetFilterRadius", &ExplicitVertexMorphingFilter<ModelPart::ConditionsContainerType>::SetFilterRadius)
        .def("FilterField", &ExplicitVertexMorphingFilter<ModelPart::ConditionsContainerType>::FilterField)
        .def("UnFilterField", &ExplicitVertexMorphingFilter<ModelPart::ConditionsContainerType>::UnFilterField)
        .def("__str__", &ExplicitVertexMorphingFilter<ModelPart::ConditionsContainerType>::Info)
        ;

    py::class_<ExplicitVertexMorphingFilter<ModelPart::ElementsContainerType>, ExplicitVertexMorphingFilter<ModelPart::ElementsContainerType>::Pointer>(m, "ElementExplicitVertexMorphingFilter")
        .def(py::init<const ModelPart&, const std::string&>())
        .def("SetFilterRadius", &ExplicitVertexMorphingFilter<ModelPart::ElementsContainerType>::SetFilterRadius)
        .def("FilterField", &ExplicitVertexMorphingFilter<ModelPart::ElementsContainerType>::FilterField)
        .def("UnFilterField", &ExplicitVertexMorphingFilter<ModelPart::ElementsContainerType>::UnFilterField)
        .def("__str__", &ExplicitVertexMorphingFilter<ModelPart::ElementsContainerType>::Info)
        ;
}

} // namespace Python.
} // namespace Kratos.
