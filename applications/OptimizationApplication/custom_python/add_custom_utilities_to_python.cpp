//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

// System includes

// External includes
#include <pybind11/stl.h>
#include <pybind11/operators.h>

// Project includes
#include "containers/model.h"
#include "includes/model_part.h"
#include "includes/data_communicator.h"
#include "spaces/ublas_space.h"

// Application includes
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/container_variable_data_holder/container_variable_data_holder_base.h"
#include "custom_utilities/container_variable_data_holder/container_variable_data_holder.h"
#include "custom_utilities/container_variable_data_holder_utils.h"

#include "custom_utilities/mappers/container_variable_data_mapper.h"
#include "custom_utilities/mappers/vertex_morphing_container_variable_data_mapper.h"

// Include base h
#include "add_custom_utilities_to_python.h"

#define PYBIND11_DETAILED_ERROR_MESSAGES

namespace Kratos {
namespace Python {

template<class TContainerType>
void AddContainerVariableDataHolderBaseTypeToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_variable_data_holder_base = ContainerVariableDataHolderBase<TContainerType>;
    py::class_<container_variable_data_holder_base, typename container_variable_data_holder_base::Pointer>(m, rName.c_str())
        .def("CopyDataFrom", &container_variable_data_holder_base::CopyDataFrom, py::arg("origin_container_data"))
        .def("GetModelPart", py::overload_cast<>(&container_variable_data_holder_base::GetModelPart), py::return_value_policy::reference)
        .def("GetContainer", py::overload_cast<>(&container_variable_data_holder_base::GetContainer), py::return_value_policy::reference)
        .def("__str__", &container_variable_data_holder_base::Info)
        ;
}

template<class TContainerType, class TContainerIO>
void AddContainerVariableDataHolderTypeToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_type = ContainerVariableDataHolder<TContainerType, TContainerIO>;
    py::class_<container_type, typename container_type::Pointer, ContainerVariableDataHolderBase<TContainerType>>(m, rName.c_str())
        .def(py::init<ModelPart&>(), py::arg("model_part"), py::doc("Creates a new container data object with model_part."))
        .def(py::init<const container_type&>(), py::arg("other_container_data_to_copy_from"), py::doc("Creates a new same type container data object by copying data from other_container_data_to_copy_from."))
        .def(py::init<const typename container_type::BaseType&>(), py::arg("other_container_data_to_copy_from"), py::doc("Creates a new destination type container data object by copying data from compatible other_container_data_to_copy_from."))
        .def("AssignDataToContainerVariable", &container_type::template AssignDataToContainerVariable<double>, py::arg("scalar_variable"))
        .def("AssignDataToContainerVariable", &container_type::template AssignDataToContainerVariable<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("ReadDataFromContainerVariable", &container_type::template ReadDataFromContainerVariable<double>, py::arg("scalar_variable"))
        .def("ReadDataFromContainerVariable", &container_type::template ReadDataFromContainerVariable<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("SetDataForContainerVariable", &container_type::template SetDataForContainerVariable<double>, py::arg("scalar_variable"), py::arg("scalar_value"))
        .def("SetDataForContainerVariable", &container_type::template SetDataForContainerVariable<array_1d<double, 3>>, py::arg("Array3_variable"), py::arg("Array3_value"))
        .def("SetDataForContainerVariableToZero", &container_type::template SetDataForContainerVariableToZero<double>, py::arg("scalar_variable"))
        .def("SetDataForContainerVariableToZero", &container_type::template SetDataForContainerVariableToZero<array_1d<double, 3>>, py::arg("Array3_variable"))
        .def("Clone", &container_type::Clone)
        .def(py::self +  py::self)
        .def(py::self += py::self)
        .def(py::self +  float())
        .def(py::self += float())
        .def(py::self -  py::self)
        .def(py::self -= py::self)
        .def(py::self -  float())
        .def(py::self -= float())
        .def(py::self *  float())
        .def(py::self *= float())
        .def(py::self /  float())
        .def(py::self /= float())
        .def("__pow__", &container_type::operator^)
        .def("__ipow__", &container_type::operator^=)
        .def("__neg__", [](container_type& rSelf) { return rSelf.operator*(-1.0); })
        ;
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

    using SparseMatrixType = SparseSpaceType::MatrixType;

    AddContainerVariableDataHolderBaseTypeToPython<ModelPart::NodesContainerType>(m, "NodalContainerVariableDataHolderBase");
    AddContainerVariableDataHolderBaseTypeToPython<ModelPart::ConditionsContainerType>(m, "ConditionContainerVariableDataHolderBase");
    AddContainerVariableDataHolderBaseTypeToPython<ModelPart::ElementsContainerType>(m, "ElementContainerVariableDataHolderBase");

    AddContainerVariableDataHolderTypeToPython<ModelPart::NodesContainerType, HistoricalContainerDataIO>(m, "HistoricalContainerVariableDataHolder");
    AddContainerVariableDataHolderTypeToPython<ModelPart::NodesContainerType, NonHistoricalContainerDataIO>(m, "NodalContainerVariableDataHolder");
    AddContainerVariableDataHolderTypeToPython<ModelPart::ConditionsContainerType, NonHistoricalContainerDataIO>(m, "ConditionContainerVariableDataHolder");
    AddContainerVariableDataHolderTypeToPython<ModelPart::ElementsContainerType, NonHistoricalContainerDataIO>(m, "ElementContainerVariableDataHolder");
    AddContainerVariableDataHolderTypeToPython<ModelPart::ConditionsContainerType, PropertiesContainerDataIO>(m, "ConditionPropertiesContainerVariableDataHolder");
    AddContainerVariableDataHolderTypeToPython<ModelPart::ElementsContainerType, PropertiesContainerDataIO>(m, "ElementPropertiesContainerVariableDataHolder");

    py::class_<ContainerVariableDataHolderUtils>(m, "ContainerVariableDataHolderUtils")
        .def_static("NormInf", &ContainerVariableDataHolderUtils::NormInf<ModelPart::NodesContainerType>, py::arg("container_data"))
        .def_static("NormInf", &ContainerVariableDataHolderUtils::NormInf<ModelPart::ConditionsContainerType>, py::arg("container_data"))
        .def_static("NormInf", &ContainerVariableDataHolderUtils::NormInf<ModelPart::ElementsContainerType>, py::arg("container_data"))
        .def_static("EntityMaxNormL2", &ContainerVariableDataHolderUtils::EntityMaxNormL2<ModelPart::NodesContainerType>, py::arg("container_data"))
        .def_static("EntityMaxNormL2", &ContainerVariableDataHolderUtils::EntityMaxNormL2<ModelPart::ConditionsContainerType>, py::arg("container_data"))
        .def_static("EntityMaxNormL2", &ContainerVariableDataHolderUtils::EntityMaxNormL2<ModelPart::ElementsContainerType>, py::arg("container_data"))
        .def_static("InnerProduct", &ContainerVariableDataHolderUtils::InnerProduct<ModelPart::NodesContainerType>, py::arg("container_data_1"), py::arg("container_data_2"))
        .def_static("InnerProduct", &ContainerVariableDataHolderUtils::InnerProduct<ModelPart::ConditionsContainerType>, py::arg("container_data_1"), py::arg("container_data_2"))
        .def_static("InnerProduct", &ContainerVariableDataHolderUtils::InnerProduct<ModelPart::ElementsContainerType>, py::arg("container_data_1"), py::arg("container_data_2"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&, const Matrix&, const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::NodesContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>&, const Matrix&, const ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::ConditionsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>&, const Matrix&, const ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::ElementsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&, const SparseMatrixType&, const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::NodesContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>&, const SparseMatrixType&, const ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::ConditionsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>&, const SparseMatrixType&, const ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::ElementsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("Transpose", py::overload_cast<SparseMatrixType&,const SparseMatrixType&>(&ContainerVariableDataHolderUtils::Transpose), py::arg("output_matrix"), py::arg("input_matrix"))
        .def_static("Transpose", py::overload_cast<Matrix&,const Matrix&>(&ContainerVariableDataHolderUtils::Transpose), py::arg("output_matrix"), py::arg("input_matrix"))
        ;

    py::class_<OptimizationUtils >(m, "OptimizationUtils")
        .def_static("IsSameModelPart", &OptimizationUtils::IsSameModelPart)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,double>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,double>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def_static("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def_static("AreAllEntitiesOfSameGeometryType", [](ModelPart::ConditionsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def_static("AreAllEntitiesOfSameGeometryType", [](ModelPart::ElementsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def_static("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ConditionsContainerType>)
        .def_static("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ElementsContainerType>)
        ;

    // add mappers
    using nodal_container_variable_data_mapper = ContainerVariableDataMapper<ModelPart::NodesContainerType>;
    py::class_<nodal_container_variable_data_mapper, typename nodal_container_variable_data_mapper::Pointer>(m, "NodalContainerVariableDataMapper")
        .def(py::init<>())
        .def("Update", &nodal_container_variable_data_mapper::Update)
        .def("Map", &nodal_container_variable_data_mapper::Map, py::arg("origin_container_data"), py::arg("destination_container_data"))
        .def("InverseMap", &nodal_container_variable_data_mapper::InverseMap, py::arg("origin_container_data"), py::arg("destination_container_data"))
        .def("__str__", &nodal_container_variable_data_mapper::Info)
        ;

    using condition_container_variable_data_mapper = ContainerVariableDataMapper<ModelPart::ConditionsContainerType>;
    py::class_<condition_container_variable_data_mapper, typename condition_container_variable_data_mapper::Pointer>(m, "ConditionContainerVariableDataMapper")
        .def(py::init<>())
        .def("Update", &condition_container_variable_data_mapper::Update)
        .def("Map", &condition_container_variable_data_mapper::Map, py::arg("origin_container_data"), py::arg("destination_container_data"))
        .def("InverseMap", &condition_container_variable_data_mapper::InverseMap, py::arg("origin_container_data"), py::arg("destination_container_data"))
        .def("__str__", &condition_container_variable_data_mapper::Info)
        ;

    using element_container_variable_data_mapper = ContainerVariableDataMapper<ModelPart::ElementsContainerType>;
    py::class_<element_container_variable_data_mapper, typename element_container_variable_data_mapper::Pointer>(m, "ElementContainerVariableDataMapper")
        .def(py::init<>())
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

