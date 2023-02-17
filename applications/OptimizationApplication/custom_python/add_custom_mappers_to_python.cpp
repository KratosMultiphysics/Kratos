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
#include "custom_utilities/mappers/container_variable_data_mapper.h"
#include "custom_utilities/mappers/vertex_morphing_container_variable_data_mapper.h"

// Include base h
#include "add_custom_mappers_to_python.h"

namespace Kratos {
namespace Python {

void  AddCustomMappersToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using nodal_container_variable_data_mapper = ContainerVariableDataMapper<ModelPart::NodesContainerType>;
    py::class_<nodal_container_variable_data_mapper, typename nodal_container_variable_data_mapper::Pointer>(m, "NodalContainerVariableDataMapper")
        .def("Update", &nodal_container_variable_data_mapper::Update)
        .def("Map", &nodal_container_variable_data_mapper::Map, py::arg("origin_container_data"), py::arg("destination_container_data"))
        .def("InverseMap", &nodal_container_variable_data_mapper::InverseMap, py::arg("origin_container_data"), py::arg("destination_container_data"))
        .def("__str__", &nodal_container_variable_data_mapper::Info)
        ;

    using condition_container_variable_data_mapper = ContainerVariableDataMapper<ModelPart::ConditionsContainerType>;
    py::class_<condition_container_variable_data_mapper, typename condition_container_variable_data_mapper::Pointer>(m, "ConditionContainerVariableDataMapper")
        .def("Update", &condition_container_variable_data_mapper::Update)
        .def("Map", &condition_container_variable_data_mapper::Map, py::arg("origin_container_data"), py::arg("destination_container_data"))
        .def("InverseMap", &condition_container_variable_data_mapper::InverseMap, py::arg("origin_container_data"), py::arg("destination_container_data"))
        .def("__str__", &condition_container_variable_data_mapper::Info)
        ;

    using element_container_variable_data_mapper = ContainerVariableDataMapper<ModelPart::ElementsContainerType>;
    py::class_<element_container_variable_data_mapper, typename element_container_variable_data_mapper::Pointer>(m, "ElementContainerVariableDataMapper")
        .def("Update", &element_container_variable_data_mapper::Update)
        .def("Map", &element_container_variable_data_mapper::Map, py::arg("origin_container_data"), py::arg("destination_container_data"))
        .def("InverseMap", &element_container_variable_data_mapper::InverseMap, py::arg("origin_container_data"), py::arg("destination_container_data"))
        .def("__str__", &element_container_variable_data_mapper::Info)
        ;

    using vertex_morphing_nodal_container_variable_data_mapper = VertexMorphingContainerVariableDataMapper<ModelPart::NodesContainerType>;
    py::class_<vertex_morphing_nodal_container_variable_data_mapper, typename vertex_morphing_nodal_container_variable_data_mapper::Pointer, nodal_container_variable_data_mapper>(m, "VertexMorphingNodalContainerVariableDataMapper")
        .def(py::init<ModelPart&, ModelPart&, Parameters>(), py::arg("origin_model_part"), py::arg("destination_model_part"), py::arg("parameters"))
        ;

    using vertex_morphing_condition_container_variable_data_mapper = VertexMorphingContainerVariableDataMapper<ModelPart::ConditionsContainerType>;
    py::class_<vertex_morphing_condition_container_variable_data_mapper, typename vertex_morphing_condition_container_variable_data_mapper::Pointer, condition_container_variable_data_mapper>(m, "VertexMorphingConditionContainerVariableDataMapper")
        .def(py::init<ModelPart&, ModelPart&, Parameters>(), py::arg("origin_model_part"), py::arg("destination_model_part"), py::arg("parameters"))
        ;

    using vertex_morphing_element_container_variable_data_mapper = VertexMorphingContainerVariableDataMapper<ModelPart::ElementsContainerType>;
    py::class_<vertex_morphing_element_container_variable_data_mapper, typename vertex_morphing_element_container_variable_data_mapper::Pointer, element_container_variable_data_mapper>(m, "VertexMorphingElementContainerVariableDataMapper")
        .def(py::init<ModelPart&, ModelPart&, Parameters>(), py::arg("origin_model_part"), py::arg("destination_model_part"), py::arg("parameters"))
        ;

}

}  // namespace Python.
} // Namespace Kratos

