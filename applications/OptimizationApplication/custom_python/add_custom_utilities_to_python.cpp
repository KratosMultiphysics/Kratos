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
#include <pybind11/operators.h>

// Project includes

// Application includes
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/container_variable_data_holder/container_variable_data_holder.h"
#include "custom_utilities/container_variable_data_holder/collective_variable_data_holder.h"
#include "custom_utilities/container_variable_data_holder_utils.h"
#include "custom_utilities/controls/simp_utils.h"

// Include base h
#include "add_custom_response_utilities_to_python.h"

namespace Kratos {
namespace Python {

template<class TContainerType>
void AddContainerVariableDataHolderBaseTypeToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_variable_data_holder_base = ContainerVariableDataHolderBase<TContainerType>;
    py::class_<container_variable_data_holder_base, typename container_variable_data_holder_base::Pointer>(m, rName.c_str())
        .def("CopyDataFrom", &container_variable_data_holder_base::CopyDataFrom, py::arg("origin_container_data"))
        .def("GetDataDimension", &container_variable_data_holder_base::GetDataDimension)
        .def("GetModelPart", py::overload_cast<>(&container_variable_data_holder_base::GetModelPart), py::return_value_policy::reference)
        .def("GetContainer", py::overload_cast<>(&container_variable_data_holder_base::GetContainer), py::return_value_policy::reference)
        .def("PrintData", &container_variable_data_holder_base::PrintData)
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
        .def("CloneWithDataInitializedToZero", &container_type::CloneWithDataInitializedToZero)
        .def(py::self +  py::self)
        .def(py::self += py::self)
        .def(py::self +  float())
        .def(py::self += float())
        .def(py::self -  py::self)
        .def(py::self -= py::self)
        .def(py::self -  float())
        .def(py::self -= float())
        .def(py::self *  py::self)
        .def(py::self *= py::self)
        .def(py::self *  float())
        .def(py::self *= float())
        .def(py::self /  py::self)
        .def(py::self /= py::self)
        .def(py::self /  float())
        .def(py::self /= float())
        .def("__pow__", py::overload_cast<const container_type&>(&container_type::operator^, py::const_))
        .def("__ipow__", py::overload_cast<const container_type&>(&container_type::operator^=))
        .def("__pow__", py::overload_cast<const double>(&container_type::operator^, py::const_))
        .def("__ipow__", py::overload_cast<const double>(&container_type::operator^=))
        .def("__neg__", [](container_type& rSelf) { return rSelf.operator*(-1.0); })
        ;
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

    using SparseMatrixType = SparseSpaceType::MatrixType;

    py::class_<OptimizationUtils >(m, "OptimizationUtils")
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
        .def_static("GetVariableDimension", &OptimizationUtils::GetVariableDimension<double>)
        .def_static("GetVariableDimension", &OptimizationUtils::GetVariableDimension<array_1d<double, 3>>)
        .def_static("CopySolutionStepVariablesList", &OptimizationUtils::CopySolutionStepVariablesList)
        ;

    AddContainerVariableDataHolderBaseTypeToPython<ModelPart::NodesContainerType>(m, "NodalContainerVariableDataHolderBase");
    AddContainerVariableDataHolderBaseTypeToPython<ModelPart::ConditionsContainerType>(m, "ConditionContainerVariableDataHolderBase");
    AddContainerVariableDataHolderBaseTypeToPython<ModelPart::ElementsContainerType>(m, "ElementContainerVariableDataHolderBase");

    AddContainerVariableDataHolderTypeToPython<ModelPart::NodesContainerType, HistoricalContainerDataIO>(m, "HistoricalContainerVariableDataHolder");
    AddContainerVariableDataHolderTypeToPython<ModelPart::NodesContainerType, NonHistoricalContainerDataIO>(m, "NodalContainerVariableDataHolder");
    AddContainerVariableDataHolderTypeToPython<ModelPart::ConditionsContainerType, NonHistoricalContainerDataIO>(m, "ConditionContainerVariableDataHolder");
    AddContainerVariableDataHolderTypeToPython<ModelPart::ElementsContainerType, NonHistoricalContainerDataIO>(m, "ElementContainerVariableDataHolder");
    AddContainerVariableDataHolderTypeToPython<ModelPart::ConditionsContainerType, PropertiesContainerDataIO>(m, "ConditionPropertiesContainerVariableDataHolder");
    AddContainerVariableDataHolderTypeToPython<ModelPart::ElementsContainerType, PropertiesContainerDataIO>(m, "ElementPropertiesContainerVariableDataHolder");

    py::class_<CollectiveVariableDataHolder, CollectiveVariableDataHolder::Pointer>(m, "CollectiveVariableDataHolder")
        .def(py::init<>())
        .def(py::init<const CollectiveVariableDataHolder&>())
        .def(py::init<const std::vector<CollectiveVariableDataHolder::ContainerVariableDataHolderPointerVariantType>&>())
        .def("AddVariableDataHolder", &CollectiveVariableDataHolder::AddVariableDataHolder)
        .def("AddCollectiveVariableDataHolder", &CollectiveVariableDataHolder::AddCollectiveVariableDataHolder)
        .def("ClearVariableDataHolders", &CollectiveVariableDataHolder::ClearVariableDataHolders)
        .def("GetVariableDataHolders", py::overload_cast<>(&CollectiveVariableDataHolder::GetVariableDataHolders))
        .def("IsCompatibleWith", &CollectiveVariableDataHolder::IsCompatibleWith)
        .def("Clone", &CollectiveVariableDataHolder::Clone)
        .def("CloneWithDataInitializedToZero", &CollectiveVariableDataHolder::CloneWithDataInitializedToZero)
        .def(py::self +  py::self)
        .def(py::self += py::self)
        .def(py::self +  float())
        .def(py::self += float())
        .def(py::self -  py::self)
        .def(py::self -= py::self)
        .def(py::self -  float())
        .def(py::self -= float())
        .def(py::self *  py::self)
        .def(py::self *= py::self)
        .def(py::self *  float())
        .def(py::self *= float())
        .def(py::self /  py::self)
        .def(py::self /= py::self)
        .def(py::self /  float())
        .def(py::self /= float())
        .def("__pow__", py::overload_cast<const CollectiveVariableDataHolder&>(&CollectiveVariableDataHolder::operator^, py::const_))
        .def("__ipow__", py::overload_cast<const CollectiveVariableDataHolder&>(&CollectiveVariableDataHolder::operator^=))
        .def("__pow__", py::overload_cast<const double>(&CollectiveVariableDataHolder::operator^, py::const_))
        .def("__ipow__", py::overload_cast<const double>(&CollectiveVariableDataHolder::operator^=))
        .def("__neg__", [](CollectiveVariableDataHolder& rSelf) { return rSelf.operator*(-1.0); })
        ;

    py::class_<ContainerVariableDataHolderUtils>(m, "ContainerVariableDataHolderUtils")
        .def_static("NormInf", &ContainerVariableDataHolderUtils::NormInf<ModelPart::NodesContainerType>, py::arg("container_data"))
        .def_static("NormInf", &ContainerVariableDataHolderUtils::NormInf<ModelPart::ConditionsContainerType>, py::arg("container_data"))
        .def_static("NormInf", &ContainerVariableDataHolderUtils::NormInf<ModelPart::ElementsContainerType>, py::arg("container_data"))
        .def_static("NormInf", [](const CollectiveVariableDataHolder& rContainer) { return ContainerVariableDataHolderUtils::NormInf(rContainer); }, py::arg("container_data"))
        .def_static("NormL2", &ContainerVariableDataHolderUtils::NormL2<ModelPart::NodesContainerType>, py::arg("container_data"))
        .def_static("NormL2", &ContainerVariableDataHolderUtils::NormL2<ModelPart::ConditionsContainerType>, py::arg("container_data"))
        .def_static("NormL2", &ContainerVariableDataHolderUtils::NormL2<ModelPart::ElementsContainerType>, py::arg("container_data"))
        .def_static("NormL2", [](const CollectiveVariableDataHolder& rContainer) { return ContainerVariableDataHolderUtils::NormL2(rContainer); }, py::arg("container_data"))
        .def_static("EntityMaxNormL2", &ContainerVariableDataHolderUtils::EntityMaxNormL2<ModelPart::NodesContainerType>, py::arg("container_data"))
        .def_static("EntityMaxNormL2", &ContainerVariableDataHolderUtils::EntityMaxNormL2<ModelPart::ConditionsContainerType>, py::arg("container_data"))
        .def_static("EntityMaxNormL2", &ContainerVariableDataHolderUtils::EntityMaxNormL2<ModelPart::ElementsContainerType>, py::arg("container_data"))
        .def_static("InnerProduct", &ContainerVariableDataHolderUtils::InnerProduct<ModelPart::NodesContainerType>, py::arg("container_data_1"), py::arg("container_data_2"))
        .def_static("InnerProduct", &ContainerVariableDataHolderUtils::InnerProduct<ModelPart::ConditionsContainerType>, py::arg("container_data_1"), py::arg("container_data_2"))
        .def_static("InnerProduct", &ContainerVariableDataHolderUtils::InnerProduct<ModelPart::ElementsContainerType>, py::arg("container_data_1"), py::arg("container_data_2"))
        .def_static("InnerProduct", [](const CollectiveVariableDataHolder& rV1, const CollectiveVariableDataHolder& rV2) { return ContainerVariableDataHolderUtils::InnerProduct(rV1, rV2); }, py::arg("container_data_1"), py::arg("container_data_2"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&, const Matrix&, const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::NodesContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>&, const Matrix&, const ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::ConditionsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>&, const Matrix&, const ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::ElementsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&, const SparseMatrixType&, const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::NodesContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>&, const SparseMatrixType&, const ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::ConditionsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("ProductWithEntityMatrix", py::overload_cast<ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>&, const SparseMatrixType&, const ContainerVariableDataHolderBase<ModelPart::ElementsContainerType>&>(&ContainerVariableDataHolderUtils::ProductWithEntityMatrix<ModelPart::ElementsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def_static("Transpose", py::overload_cast<SparseMatrixType&,const SparseMatrixType&>(&ContainerVariableDataHolderUtils::Transpose), py::arg("output_matrix"), py::arg("input_matrix"))
        .def_static("Transpose", py::overload_cast<Matrix&,const Matrix&>(&ContainerVariableDataHolderUtils::Transpose), py::arg("output_matrix"), py::arg("input_matrix"))
        .def_static("ComputeNumberOfNeighbourEntities", &ContainerVariableDataHolderUtils::ComputeNumberOfNeighbourEntities<ModelPart::ConditionsContainerType>)
        .def_static("ComputeNumberOfNeighbourEntities", &ContainerVariableDataHolderUtils::ComputeNumberOfNeighbourEntities<ModelPart::ElementsContainerType>)
        .def_static("MapContainerVariableDataHolderToNodalVariableDataHolder", &ContainerVariableDataHolderUtils::MapContainerVariableDataHolderToNodalVariableDataHolder<ModelPart::ConditionsContainerType>)
        .def_static("MapContainerVariableDataHolderToNodalVariableDataHolder", &ContainerVariableDataHolderUtils::MapContainerVariableDataHolderToNodalVariableDataHolder<ModelPart::ElementsContainerType>)
        .def_static("MapNodalVariableDataHolderToContainerVariableDataHolder", &ContainerVariableDataHolderUtils::MapNodalVariableDataHolderToContainerVariableDataHolder<ModelPart::ConditionsContainerType>)
        .def_static("MapNodalVariableDataHolderToContainerVariableDataHolder", &ContainerVariableDataHolderUtils::MapNodalVariableDataHolderToContainerVariableDataHolder<ModelPart::ElementsContainerType>)
        .def_static("ComputeVariableDataHolderProductWithEntityMatrix", &ContainerVariableDataHolderUtils::ComputeVariableDataHolderProductWithEntityMatrix<ModelPart::ConditionsContainerType>)
        .def_static("ComputeVariableDataHolderProductWithEntityMatrix", &ContainerVariableDataHolderUtils::ComputeVariableDataHolderProductWithEntityMatrix<ModelPart::ElementsContainerType>)
        ;

    py::class_<SimpUtils>(m, "SimpUtils")
        .def_static("ProjectContainerVariableDataHolderForward", &SimpUtils::ProjectContainerVariableDataHolderForward<ModelPart::NodesContainerType>)
        .def_static("ProjectContainerVariableDataHolderBackward", &SimpUtils::ProjectContainerVariableDataHolderBackward<ModelPart::NodesContainerType>)
        .def_static("ProjectContainerVariableDataHolderDerivative", &SimpUtils::ProjectContainerVariableDataHolderDerivative<ModelPart::NodesContainerType>)
        .def_static("ProjectContainerVariableDataHolderForward", &SimpUtils::ProjectContainerVariableDataHolderForward<ModelPart::ConditionsContainerType>)
        .def_static("ProjectContainerVariableDataHolderBackward", &SimpUtils::ProjectContainerVariableDataHolderBackward<ModelPart::ConditionsContainerType>)
        .def_static("ProjectContainerVariableDataHolderDerivative", &SimpUtils::ProjectContainerVariableDataHolderDerivative<ModelPart::ConditionsContainerType>)
        .def_static("ProjectContainerVariableDataHolderForward", &SimpUtils::ProjectContainerVariableDataHolderForward<ModelPart::ElementsContainerType>)
        .def_static("ProjectContainerVariableDataHolderBackward", &SimpUtils::ProjectContainerVariableDataHolderBackward<ModelPart::ElementsContainerType>)
        .def_static("ProjectContainerVariableDataHolderDerivative", &SimpUtils::ProjectContainerVariableDataHolderDerivative<ModelPart::ElementsContainerType>)
        ;
}

}  // namespace Python.
} // Namespace Kratos

