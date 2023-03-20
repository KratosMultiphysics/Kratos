//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl,
//                   Suneth Warnakulasuriya
//

// System includes

// External includes
#include <pybind11/stl.h>
#include <pybind11/operators.h>

// Project includes

// Application includes
#include "custom_utilities/geometrical/symmetry_utility.h"
#include "custom_utilities/geometrical/model_part_utils.h"
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/container_variable_data/container_data_io.h"
#include "custom_utilities/container_variable_data/container_variable_data.h"
#include "custom_utilities/container_variable_data/specialized_container_variable_data.h"
#include "custom_utilities/container_variable_data_utils.h"

// Include base h
#include "add_custom_response_utilities_to_python.h"

namespace Kratos {
namespace Python {

template<class TContainerType>
void AddContainerVariableDataToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_variable_data_holder_base = ContainerVariableData<TContainerType>;
    py::class_<container_variable_data_holder_base, typename container_variable_data_holder_base::Pointer>(m, rName.c_str())
        .def("CopyDataFrom", &container_variable_data_holder_base::CopyDataFrom, py::arg("origin_container_data"))
        .def("GetDataDimension", &container_variable_data_holder_base::GetDataDimension)
        .def("GetModelPart", py::overload_cast<>(&container_variable_data_holder_base::GetModelPart), py::return_value_policy::reference)
        .def("GetContainer", py::overload_cast<>(&container_variable_data_holder_base::GetContainer), py::return_value_policy::reference)
        .def("PrintData", &container_variable_data_holder_base::PrintData)
        .def("__str__", &container_variable_data_holder_base::Info)
        ;
}

template<class TContainerType, class TContainerDataIOTag>
void AddSpecializedContainerVariableDataToPython(pybind11::module& m, const std::string& rName)
{
    namespace py = pybind11;

    using container_type = SpecializedContainerVariableData<TContainerType, ContainerDataIO<TContainerDataIOTag>>;
    py::class_<container_type, typename container_type::Pointer, ContainerVariableData<TContainerType>>(m, rName.c_str())
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
        .def("__pow__", [](container_type& rSelf, const container_type& rInput) { container_type result(rSelf.GetModelPart()); ContainerVariableDataUtils::Pow(result, rSelf, rInput); return result; })
        .def("__ipow__", [](container_type& rSelf, const container_type& rInput) { ContainerVariableDataUtils::Pow(rSelf, rSelf, rInput); return rSelf; })
        .def("__pow__", [](container_type& rSelf, const double Value) { container_type result(rSelf.GetModelPart()); ContainerVariableDataUtils::Pow(result, rSelf, Value); return result; })
        .def("__ipow__", [](container_type& rSelf, const double Value) { ContainerVariableDataUtils::Pow(rSelf, rSelf, Value); return rSelf; })
        .def("__neg__", [](container_type& rSelf) { return rSelf.operator*(-1.0); })
        ;
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

    using SparseMatrixType = SparseSpaceType::MatrixType;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

    using SparseMatrixType = SparseSpaceType::MatrixType;

    using namespace pybind11::literals;

    py::class_<SymmetryUtility >(m, "SymmetryUtility")
        .def(py::init<std::string, ModelPart&, Parameters>())
        .def("Initialize", &SymmetryUtility::Initialize)
        .def("Update", &SymmetryUtility::Update)
        .def("ApplyOnVectorField", &SymmetryUtility::ApplyOnVectorField)
        .def("ApplyOnScalarField", &SymmetryUtility::ApplyOnScalarField)
        ;

    m.def_submodule("ModelPartUtils")
        .def("GetModelPartsWithCommonReferenceEntities", [](
            const std::vector<ModelPart*>& rEvaluatedModelPartsList,
            const std::vector<ModelPart*>& rReferenceModelPartsList,
            const bool AreNodesConsidered,
            const bool AreConditionsConsidered,
            const bool AreElementsConsidered,
            const bool AreParentsConsidered,
            const IndexType EchoLevel){
                const auto& r_model_parts = ModelPartUtils::GetModelPartsWithCommonReferenceEntities(
                    rEvaluatedModelPartsList,
                    rReferenceModelPartsList,
                    AreNodesConsidered,
                    AreConditionsConsidered,
                    AreElementsConsidered,
                    AreParentsConsidered,
                    EchoLevel);

                auto pylist = py::list();
                for (auto p_model_part : r_model_parts) {
                    auto pyobj = py::cast(*p_model_part, py::return_value_policy::reference);
                    pylist.append(pyobj);
                }
                return pylist;
            },
            "examined_model_parts_list"_a,
            "reference_model_parts"_a,
            "are_nodes_considered"_a,
            "are_conditions_considered"_a,
            "are_elements_considered"_a,
            "are_parents_considered"_a,
            "echo_level"_a = 0)
        .def("RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList", &ModelPartUtils::RemoveModelPartsWithCommonReferenceEntitiesBetweenReferenceListAndExaminedList,
            "model_parts_list"_a)
        ;

    m.def_submodule("OptimizationUtils")
        .def("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,double>)
        .def("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def("IsVariableExistsInAllContainerProperties", &OptimizationUtils::IsVariableExistsInAllContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, double>)
        .def("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,double>)
        .def("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ConditionsContainerType, array_1d<double, 3>>)
        .def("IsVariableExistsInAtLeastOneContainerProperties", &OptimizationUtils::IsVariableExistsInAtLeastOneContainerProperties<ModelPart::ElementsContainerType,array_1d<double, 3>>)
        .def("AreAllEntitiesOfSameGeometryType", [](ModelPart::ConditionsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def("AreAllEntitiesOfSameGeometryType", [](ModelPart::ElementsContainerType& rContainer, const DataCommunicator& rDataCommunicator) { return OptimizationUtils::GetContainerEntityGeometryType(rContainer, rDataCommunicator) != GeometryData::KratosGeometryType::Kratos_generic_type; } )
        .def("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ConditionsContainerType>)
        .def("CreateEntitySpecificPropertiesForContainer", &OptimizationUtils::CreateEntitySpecificPropertiesForContainer<ModelPart::ElementsContainerType>)
        .def("GetVariableDimension", &OptimizationUtils::GetVariableDimension<double>)
        .def("GetVariableDimension", &OptimizationUtils::GetVariableDimension<array_1d<double, 3>>)
        .def("CopySolutionStepVariablesList", &OptimizationUtils::CopySolutionStepVariablesList)
        ;

    AddContainerVariableDataToPython<ModelPart::NodesContainerType>(m, "NodalVariableData");
    AddContainerVariableDataToPython<ModelPart::ConditionsContainerType>(m, "ConditionVariableData");
    AddContainerVariableDataToPython<ModelPart::ElementsContainerType>(m, "ElementVariableData");

    AddSpecializedContainerVariableDataToPython<ModelPart::NodesContainerType, ContainerDataIOTags::Historical>(m, "HistoricalVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::NodesContainerType, ContainerDataIOTags::NonHistorical>(m, "NodalNonHistoricalVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::ConditionsContainerType, ContainerDataIOTags::NonHistorical>(m, "ConditionNonHistoricalVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::ElementsContainerType, ContainerDataIOTags::NonHistorical>(m, "ElementNonHistoricalVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::ConditionsContainerType, ContainerDataIOTags::Properties>(m, "ConditionPropertiesVariableData");
    AddSpecializedContainerVariableDataToPython<ModelPart::ElementsContainerType, ContainerDataIOTags::Properties>(m, "ElementPropertiesVariableData");

    m.def_submodule("ContainerVariableDataUtils")
        .def("NormInf", &ContainerVariableDataUtils::NormInf<ModelPart::NodesContainerType>, py::arg("container_data"))
        .def("NormInf", &ContainerVariableDataUtils::NormInf<ModelPart::ConditionsContainerType>, py::arg("container_data"))
        .def("NormInf", &ContainerVariableDataUtils::NormInf<ModelPart::ElementsContainerType>, py::arg("container_data"))
        .def("NormL2", &ContainerVariableDataUtils::NormL2<ModelPart::NodesContainerType>, py::arg("container_data"))
        .def("NormL2", &ContainerVariableDataUtils::NormL2<ModelPart::ConditionsContainerType>, py::arg("container_data"))
        .def("NormL2", &ContainerVariableDataUtils::NormL2<ModelPart::ElementsContainerType>, py::arg("container_data"))
        .def("EntityMaxNormL2", &ContainerVariableDataUtils::EntityMaxNormL2<ModelPart::NodesContainerType>, py::arg("container_data"))
        .def("EntityMaxNormL2", &ContainerVariableDataUtils::EntityMaxNormL2<ModelPart::ConditionsContainerType>, py::arg("container_data"))
        .def("EntityMaxNormL2", &ContainerVariableDataUtils::EntityMaxNormL2<ModelPart::ElementsContainerType>, py::arg("container_data"))
        .def("InnerProduct", &ContainerVariableDataUtils::InnerProduct<ModelPart::NodesContainerType>, py::arg("container_data_1"), py::arg("container_data_2"))
        .def("InnerProduct", &ContainerVariableDataUtils::InnerProduct<ModelPart::ConditionsContainerType>, py::arg("container_data_1"), py::arg("container_data_2"))
        .def("InnerProduct", &ContainerVariableDataUtils::InnerProduct<ModelPart::ElementsContainerType>, py::arg("container_data_1"), py::arg("container_data_2"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerVariableData<ModelPart::NodesContainerType>&, const Matrix&, const ContainerVariableData<ModelPart::NodesContainerType>&>(&ContainerVariableDataUtils::ProductWithEntityMatrix<ModelPart::NodesContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerVariableData<ModelPart::ConditionsContainerType>&, const Matrix&, const ContainerVariableData<ModelPart::ConditionsContainerType>&>(&ContainerVariableDataUtils::ProductWithEntityMatrix<ModelPart::ConditionsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerVariableData<ModelPart::ElementsContainerType>&, const Matrix&, const ContainerVariableData<ModelPart::ElementsContainerType>&>(&ContainerVariableDataUtils::ProductWithEntityMatrix<ModelPart::ElementsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerVariableData<ModelPart::NodesContainerType>&, const SparseMatrixType&, const ContainerVariableData<ModelPart::NodesContainerType>&>(&ContainerVariableDataUtils::ProductWithEntityMatrix<ModelPart::NodesContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerVariableData<ModelPart::ConditionsContainerType>&, const SparseMatrixType&, const ContainerVariableData<ModelPart::ConditionsContainerType>&>(&ContainerVariableDataUtils::ProductWithEntityMatrix<ModelPart::ConditionsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerVariableData<ModelPart::ElementsContainerType>&, const SparseMatrixType&, const ContainerVariableData<ModelPart::ElementsContainerType>&>(&ContainerVariableDataUtils::ProductWithEntityMatrix<ModelPart::ElementsContainerType>), py::arg("output_container_data"), py::arg("matrix_with_entity_size"), py::arg("input_container_data_for_multiplication"))
        .def("Transpose", py::overload_cast<SparseMatrixType&,const SparseMatrixType&>(&ContainerVariableDataUtils::Transpose), py::arg("output_matrix"), py::arg("input_matrix"))
        .def("Transpose", py::overload_cast<Matrix&,const Matrix&>(&ContainerVariableDataUtils::Transpose), py::arg("output_matrix"), py::arg("input_matrix"))
        .def("ComputeNumberOfNeighbourConditions", &ContainerVariableDataUtils::ComputeNumberOfNeighbourEntities<ModelPart::ConditionsContainerType>, py::arg("output_nodal_container_data"))
        .def("ComputeNumberOfNeighbourElements", &ContainerVariableDataUtils::ComputeNumberOfNeighbourEntities<ModelPart::ElementsContainerType>, py::arg("output_nodal_container_data"))
        .def("MapContainerVariableDataToNodalVariableData", &ContainerVariableDataUtils::MapContainerVariableDataToNodalVariableData<ModelPart::ConditionsContainerType>, py::arg("output_nodal_container_data"), py::arg("input_container_data_to_map"), py::arg("nodal_container_data_with_neighbour_entities_for_input_container_data"))
        .def("MapContainerVariableDataToNodalVariableData", &ContainerVariableDataUtils::MapContainerVariableDataToNodalVariableData<ModelPart::ElementsContainerType>, py::arg("output_nodal_container_data"), py::arg("input_container_data_to_map"), py::arg("nodal_container_data_with_neighbour_entities_for_input_container_data"))
        .def("MapNodalVariableDataToContainerVariableData", &ContainerVariableDataUtils::MapNodalVariableDataToContainerVariableData<ModelPart::ConditionsContainerType>, py::arg("output_container_data"), py::arg("input_nodal_container_data_to_map"))
        .def("MapNodalVariableDataToContainerVariableData", &ContainerVariableDataUtils::MapNodalVariableDataToContainerVariableData<ModelPart::ElementsContainerType>, py::arg("output_container_data"), py::arg("input_nodal_container_data_to_map"))
        .def("ComputeVariableDataProductWithEntityMatrix", &ContainerVariableDataUtils::ComputeVariableDataProductWithEntityMatrix<ModelPart::ConditionsContainerType>, py::arg("output_nodal_container_data"), py::arg("input_nodal_values_container_data"), py::arg("matrix_variable"), py::arg("entities"))
        .def("ComputeVariableDataProductWithEntityMatrix", &ContainerVariableDataUtils::ComputeVariableDataProductWithEntityMatrix<ModelPart::ElementsContainerType>, py::arg("output_nodal_container_data"), py::arg("input_nodal_values_container_data"), py::arg("matrix_variable"), py::arg("entities"))
        ;
}

}  // namespace Python.
} // Namespace Kratos

