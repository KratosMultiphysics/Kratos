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

// Project includes
#include "python/add_container_expression_to_python_utils.h"
#include "containers/container_expression/specialized_container_expression.h"

// Application includes
#include "custom_utilities/geometrical/symmetry_utility.h"
#include "custom_utilities/geometrical/model_part_utils.h"
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/container_properties_data_io.h"
#include "custom_utilities/collective_expressions.h"
#include "custom_utilities/container_expression_utils.h"

// Include base h
#include "add_custom_response_utilities_to_python.h"

namespace Kratos {
namespace Python {

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

    auto sub_module = m.def_submodule("ContainerExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::ConditionsContainerType, ContainerDataIOTags::Properties>(sub_module, "ConditionPropertiesExpression");
    AddSpecializedContainerExpressionToPython<ModelPart::ElementsContainerType, ContainerDataIOTags::Properties>(sub_module, "ElementPropertiesExpression");

    py::class_<CollectiveExpressions, CollectiveExpressions::Pointer>(sub_module, "CollectiveExpressions")
        .def(py::init<>())
        .def(py::init<const CollectiveExpressions&>())
        .def(py::init<const std::vector<CollectiveExpressions::CollectiveExpressionType>&>())
        .def("Add", py::overload_cast<const CollectiveExpressions::CollectiveExpressionType&>(&CollectiveExpressions::Add))
        .def("Add", py::overload_cast<const CollectiveExpressions&>(&CollectiveExpressions::Add))
        .def("Clear", &CollectiveExpressions::Clear)
        .def("GetContainerExpressions", py::overload_cast<>(&CollectiveExpressions::GetContainerExpressions))
        .def("IsCompatibleWith", &CollectiveExpressions::IsCompatibleWith)
        .def("Clone", &CollectiveExpressions::Clone)
        .def("SetToZero", &CollectiveExpressions::SetToZero)
        .def("__add__", [](const CollectiveExpressions& rSelf, const CollectiveExpressions& rOther) { return rSelf + rOther; })
        .def("__iadd__", [](CollectiveExpressions& rSelf, const CollectiveExpressions& rOther) { rSelf = rSelf + rOther; return rSelf; })
        .def("__add__", [](const CollectiveExpressions& rSelf, const double Value) { return rSelf + Value; })
        .def("__iadd__", [](CollectiveExpressions& rSelf, const double Value) { rSelf = rSelf + Value; return rSelf; })
        .def("__sub__", [](const CollectiveExpressions& rSelf, const CollectiveExpressions& rOther) { return rSelf - rOther; })
        .def("__isub__", [](CollectiveExpressions& rSelf, const CollectiveExpressions& rOther) { rSelf = rSelf - rOther; return rSelf; })
        .def("__sub__", [](const CollectiveExpressions& rSelf, const double Value) { return rSelf - Value; })
        .def("__isub__", [](CollectiveExpressions& rSelf, const double Value) { rSelf = rSelf - Value; return rSelf; })
        .def("__mul__", [](const CollectiveExpressions& rSelf, const CollectiveExpressions& rOther) { return rSelf * rOther; })
        .def("__imul__", [](CollectiveExpressions& rSelf, const CollectiveExpressions& rOther) { rSelf = rSelf * rOther; return rSelf; })
        .def("__mul__", [](const CollectiveExpressions& rSelf, const double Value) { return rSelf * Value; })
        .def("__imul__", [](CollectiveExpressions& rSelf, const double Value) { rSelf = rSelf * Value; return rSelf; })
        .def("__truediv__", [](const CollectiveExpressions& rSelf, const CollectiveExpressions& rOther) { return rSelf / rOther; })
        .def("__itruediv__", [](CollectiveExpressions& rSelf, const CollectiveExpressions& rOther) { rSelf = rSelf / rOther; return rSelf; })
        .def("__truediv__", [](const CollectiveExpressions& rSelf, const double Value) { return rSelf / Value; })
        .def("__itruediv__", [](CollectiveExpressions& rSelf, const double Value) { rSelf = rSelf / Value; return rSelf; })
        .def("__pow__", [](CollectiveExpressions& rSelf, const CollectiveExpressions& rInput) { CollectiveExpressions result(rSelf); result = rSelf.Pow(rInput); return result; })
        .def("__ipow__", [](CollectiveExpressions& rSelf, const CollectiveExpressions& rInput) { rSelf = rSelf.Pow(rInput); return rSelf; })
        .def("__pow__", [](CollectiveExpressions& rSelf, const double Value) { CollectiveExpressions result(rSelf); result = rSelf.Pow(Value); return result; })
        .def("__ipow__", [](CollectiveExpressions& rSelf, const double Value) { rSelf = rSelf.Pow(Value); return rSelf; })
        .def("__neg__", [](CollectiveExpressions& rSelf) { return rSelf.operator*(-1.0); })
        ;

    m.def_submodule("ContainerExpressionUtils")
        .def("NormInf", &ContainerExpressionUtils::NormInf<ModelPart::NodesContainerType>, py::arg("container_expression"))
        .def("NormInf", &ContainerExpressionUtils::NormInf<ModelPart::ConditionsContainerType>, py::arg("container_expression"))
        .def("NormInf", &ContainerExpressionUtils::NormInf<ModelPart::ElementsContainerType>, py::arg("container_expression"))
        .def("NormInf", [](const CollectiveExpressions& rData){ return ContainerExpressionUtils::NormInf(rData); }, py::arg("collective_expressions"))
        .def("NormL2", &ContainerExpressionUtils::NormL2<ModelPart::NodesContainerType>, py::arg("container_expression"))
        .def("NormL2", &ContainerExpressionUtils::NormL2<ModelPart::ConditionsContainerType>, py::arg("container_expression"))
        .def("NormL2", &ContainerExpressionUtils::NormL2<ModelPart::ElementsContainerType>, py::arg("container_expression"))
        .def("NormL2", [](const CollectiveExpressions& rData){ return ContainerExpressionUtils::NormL2(rData); }, py::arg("collective_expressions"))
        .def("EntityMaxNormL2", &ContainerExpressionUtils::EntityMaxNormL2<ModelPart::NodesContainerType>, py::arg("container_expression"))
        .def("EntityMaxNormL2", &ContainerExpressionUtils::EntityMaxNormL2<ModelPart::ConditionsContainerType>, py::arg("container_expression"))
        .def("EntityMaxNormL2", &ContainerExpressionUtils::EntityMaxNormL2<ModelPart::ElementsContainerType>, py::arg("container_expression"))
        .def("InnerProduct", &ContainerExpressionUtils::InnerProduct<ModelPart::NodesContainerType>, py::arg("container_expression_1"), py::arg("container_expression_2"))
        .def("InnerProduct", &ContainerExpressionUtils::InnerProduct<ModelPart::ConditionsContainerType>, py::arg("container_expression_1"), py::arg("container_expression_2"))
        .def("InnerProduct", &ContainerExpressionUtils::InnerProduct<ModelPart::ElementsContainerType>, py::arg("container_expression_1"), py::arg("container_expression_2"))
        .def("InnerProduct", [](const CollectiveExpressions& rData1, const CollectiveExpressions& rData2){ return ContainerExpressionUtils::InnerProduct(rData1, rData2); }, py::arg("collective_expressions_1"), py::arg("collective_expressions_2"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerExpression<ModelPart::NodesContainerType>&, const Matrix&, const ContainerExpression<ModelPart::NodesContainerType>&>(&ContainerExpressionUtils::ProductWithEntityMatrix<ModelPart::NodesContainerType>), py::arg("output_container_expression"), py::arg("matrix_with_entity_size"), py::arg("input_container_expression_for_multiplication"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerExpression<ModelPart::ConditionsContainerType>&, const Matrix&, const ContainerExpression<ModelPart::ConditionsContainerType>&>(&ContainerExpressionUtils::ProductWithEntityMatrix<ModelPart::ConditionsContainerType>), py::arg("output_container_expression"), py::arg("matrix_with_entity_size"), py::arg("input_container_expression_for_multiplication"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerExpression<ModelPart::ElementsContainerType>&, const Matrix&, const ContainerExpression<ModelPart::ElementsContainerType>&>(&ContainerExpressionUtils::ProductWithEntityMatrix<ModelPart::ElementsContainerType>), py::arg("output_container_expression"), py::arg("matrix_with_entity_size"), py::arg("input_container_expression_for_multiplication"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerExpression<ModelPart::NodesContainerType>&, const SparseMatrixType&, const ContainerExpression<ModelPart::NodesContainerType>&>(&ContainerExpressionUtils::ProductWithEntityMatrix<ModelPart::NodesContainerType>), py::arg("output_container_expression"), py::arg("matrix_with_entity_size"), py::arg("input_container_expression_for_multiplication"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerExpression<ModelPart::ConditionsContainerType>&, const SparseMatrixType&, const ContainerExpression<ModelPart::ConditionsContainerType>&>(&ContainerExpressionUtils::ProductWithEntityMatrix<ModelPart::ConditionsContainerType>), py::arg("output_container_expression"), py::arg("matrix_with_entity_size"), py::arg("input_container_expression_for_multiplication"))
        .def("ProductWithEntityMatrix", py::overload_cast<ContainerExpression<ModelPart::ElementsContainerType>&, const SparseMatrixType&, const ContainerExpression<ModelPart::ElementsContainerType>&>(&ContainerExpressionUtils::ProductWithEntityMatrix<ModelPart::ElementsContainerType>), py::arg("output_container_expression"), py::arg("matrix_with_entity_size"), py::arg("input_container_expression_for_multiplication"))
        .def("Transpose", py::overload_cast<SparseMatrixType&,const SparseMatrixType&>(&ContainerExpressionUtils::Transpose), py::arg("output_matrix"), py::arg("input_matrix"))
        .def("Transpose", py::overload_cast<Matrix&,const Matrix&>(&ContainerExpressionUtils::Transpose), py::arg("output_matrix"), py::arg("input_matrix"))
        .def("ComputeNumberOfNeighbourConditions", &ContainerExpressionUtils::ComputeNumberOfNeighbourEntities<ModelPart::ConditionsContainerType>, py::arg("output_nodal_container_expression"))
        .def("ComputeNumberOfNeighbourElements", &ContainerExpressionUtils::ComputeNumberOfNeighbourEntities<ModelPart::ElementsContainerType>, py::arg("output_nodal_container_expression"))
        .def("MapContainerVariableToNodalVariable", &ContainerExpressionUtils::MapContainerVariableToNodalVariable<ModelPart::ConditionsContainerType>, py::arg("output_nodal_container_expression"), py::arg("input_container_expression_to_map"), py::arg("neighbour_container_for_nodes"))
        .def("MapContainerVariableToNodalVariable", &ContainerExpressionUtils::MapContainerVariableToNodalVariable<ModelPart::ElementsContainerType>, py::arg("output_nodal_container_expression"), py::arg("input_container_expression_to_map"), py::arg("neighbour_container_for_nodes"))
        .def("MapNodalVariableToContainerVariable", &ContainerExpressionUtils::MapNodalVariableToContainerVariable<ModelPart::ConditionsContainerType>, py::arg("output_container_expression"), py::arg("input_nodal_container_expression_to_map"))
        .def("MapNodalVariableToContainerVariable", &ContainerExpressionUtils::MapNodalVariableToContainerVariable<ModelPart::ElementsContainerType>, py::arg("output_container_expression"), py::arg("input_nodal_container_expression_to_map"))
        .def("ComputeNodalVariableProductWithEntityMatrix", &ContainerExpressionUtils::ComputeNodalVariableProductWithEntityMatrix<ModelPart::ConditionsContainerType>, py::arg("output_nodal_container_expression"), py::arg("input_nodal_values_container_expression"), py::arg("matrix_variable"), py::arg("entities"))
        .def("ComputeNodalVariableProductWithEntityMatrix", &ContainerExpressionUtils::ComputeNodalVariableProductWithEntityMatrix<ModelPart::ElementsContainerType>, py::arg("output_nodal_container_expression"), py::arg("input_nodal_values_container_expression"), py::arg("matrix_variable"), py::arg("entities"))
        ;
}

}  // namespace Python.
} // Namespace Kratos

