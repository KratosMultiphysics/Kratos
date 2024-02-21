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
#include <pybind11/numpy.h>

// Project includes
#include "python/numpy_utils.h"

// Application includes
#include "custom_utilities/collective_expression.h"
#include "custom_utilities/collective_expression_io.h"
#include "custom_utilities/container_expression_utils.h"
#include "custom_utilities/geometrical/model_part_utils.h"
#include "custom_utilities/geometrical/symmetry_utility.h"
#include "custom_utilities/implicit_filter_utils.h"
#include "custom_utilities/optimization_utils.h"
#include "custom_utilities/properties_variable_expression_io.h"

// Include base h
#include "add_custom_response_utilities_to_python.h"

namespace Kratos {
namespace Python {

namespace Detail
{

template<class TCArrayData>
void CollectiveExpressionFromPythonArray(
    CollectiveExpression& rCollectiveExpression,
    const pybind11::array_t<double>& rData,
    const std::vector<std::vector<int>>& rListOfShapes,
    TCArrayData rCArrayData)
{
    KRATOS_TRY

    auto r_container_expressions = rCollectiveExpression.GetContainerExpressions();

    KRATOS_ERROR_IF(rData.ndim() == 0) << "Passed data is not compatible.\n";

    KRATOS_ERROR_IF(r_container_expressions.size() != rListOfShapes.size())
        << "Number of container expressions in the Collective expression and list of shapes size mismatch. "
        << "[ Number of container expressions = "
        << r_container_expressions.size()
        << ", list of shapes size = " << rListOfShapes.size() << " ].\n";

    int required_number_of_values = 0;
    for (IndexType i = 0; i < r_container_expressions.size(); ++i) {
        required_number_of_values += std::visit([&rListOfShapes, i](auto& pContainerExpression){
            return pContainerExpression->GetContainer().size() * std::accumulate(rListOfShapes[i].begin(), rListOfShapes[i].end(), 1, [](const int V1, const int V2) { return V1 * V2; });
        }, r_container_expressions[i]);
    }

    KRATOS_ERROR_IF_NOT(required_number_of_values == rData.size())
        << "Required number of values for the specified ContainerExpressions "
           "and shapes mismatch with the values present in the rData array. [ "
           "required number of values = "
        << required_number_of_values << ", rData.size() = " << rData.size()
        << ", shapes = " << rListOfShapes
        << ", collective expression = " << rCollectiveExpression
        << " ].\n";

    // create c style double pointer from std::vector list for shapes.
    int const** p_list_of_shapes = new int const*[rListOfShapes.size()];
    std::transform(rListOfShapes.begin(), rListOfShapes.end(), p_list_of_shapes,
                   [](const auto& rShape) { return rShape.data(); });

    // create c style double pointer for sizes of std::vector list.
    int* p_list_of_shape_dimensions =  new int[rListOfShapes.size()];
    std::transform(rListOfShapes.begin(), rListOfShapes.end(), p_list_of_shape_dimensions,
                   [](const auto& rShape) { return rShape.size(); });

    // create c style double pointer for number of entities in each container expression.
    int* p_list_of_number_of_entities_in_container =  new int[rListOfShapes.size()];
    std::transform(r_container_expressions.begin(), r_container_expressions.end(), p_list_of_number_of_entities_in_container,
                   [](const auto& pContainerExpressionVariant) { return std::visit([](const auto pContainerExpression) { return pContainerExpression->GetContainer().size(); }, pContainerExpressionVariant); });

    if constexpr(std::is_same_v<TCArrayData, double*>) {
        CollectiveExpressionIO::Move(rCollectiveExpression, rCArrayData, p_list_of_number_of_entities_in_container,
                p_list_of_shapes, p_list_of_shape_dimensions, rListOfShapes.size());
    } else {
        CollectiveExpressionIO::Read(rCollectiveExpression, rCArrayData, p_list_of_number_of_entities_in_container,
                p_list_of_shapes, p_list_of_shape_dimensions, rListOfShapes.size());
    }

    // delete allocated memories
    delete[] p_list_of_shapes;
    delete[] p_list_of_shape_dimensions;
    delete[] p_list_of_number_of_entities_in_container;

    KRATOS_CATCH("");
}

} // namespace Detail

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
        .def("LogModelPartStatus", &ModelPartUtils::LogModelPartStatus, py::arg("model_part"), py::arg("status_to_log"))
        .def("GetModelPartStatusLog", &ModelPartUtils::GetModelPartStatusLog, py::arg("model_part"))
        .def("CheckModelPartStatus", &ModelPartUtils::CheckModelPartStatus, py::arg("model_part"), py::arg("status_to_check"))
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

    // Add collective expression to python
    pybind11::class_<CollectiveExpression, CollectiveExpression::Pointer>(m, "CollectiveExpression")
        .def(pybind11::init<>())
        .def(pybind11::init<const std::vector<CollectiveExpression::CollectiveExpressionType>&>())
        .def("Add", pybind11::overload_cast<const CollectiveExpression::CollectiveExpressionType&>(&CollectiveExpression::Add))
        .def("Add", pybind11::overload_cast<const CollectiveExpression&>(&CollectiveExpression::Add))
        .def("Clear", &CollectiveExpression::Clear)
        .def("Evaluate", [](const CollectiveExpression& rSelf){
            const auto& r_container_expressions = rSelf.GetContainerExpressions();
            if (r_container_expressions.size() > 0) {
                bool is_same_item_shape = true;
                IndexType number_of_entities = 0;
                auto current_shape =  std::visit([](auto& pContainerExpression) {return pContainerExpression->GetItemShape();}, r_container_expressions[0]);

                for (const auto& r_container_expression : r_container_expressions) {
                    is_same_item_shape = is_same_item_shape && std::visit([&current_shape, &number_of_entities](auto& pContainerExpression){
                        number_of_entities += pContainerExpression->GetContainer().size();
                        return pContainerExpression->GetItemShape() == current_shape;
                    }, r_container_expression);
                }

                // if all the container expressions does not have the same shape, make the output numpy array scalar.
                if (!is_same_item_shape) {
                    current_shape.clear();
                    number_of_entities = rSelf.GetCollectiveFlattenedDataSize();
                }

                auto array = AllocateNumpyArray<double>(number_of_entities, current_shape);
                CollectiveExpressionIO::Write(rSelf, array.mutable_data(), array.size());
                return array;
            } else {
                return AllocateNumpyArray<double>(0, {});
            }
        })
        .def("GetCollectiveFlattenedDataSize", &CollectiveExpression::GetCollectiveFlattenedDataSize)
        .def("GetContainerExpressions", pybind11::overload_cast<>(&CollectiveExpression::GetContainerExpressions))
        .def("IsCompatibleWith", &CollectiveExpression::IsCompatibleWith)
        .def("Clone", &CollectiveExpression::Clone)
        .def("__add__", [](const CollectiveExpression& rSelf, const CollectiveExpression& rOther) { return rSelf + rOther; })
        .def("__iadd__", [](CollectiveExpression& rSelf, const CollectiveExpression& rOther) { rSelf = rSelf + rOther; return rSelf; })
        .def("__add__", [](const CollectiveExpression& rSelf, const double Value) { return rSelf + Value; })
        .def("__iadd__", [](CollectiveExpression& rSelf, const double Value) { rSelf = rSelf + Value; return rSelf; })
        .def("__sub__", [](const CollectiveExpression& rSelf, const CollectiveExpression& rOther) { return rSelf - rOther; })
        .def("__isub__", [](CollectiveExpression& rSelf, const CollectiveExpression& rOther) { rSelf = rSelf - rOther; return rSelf; })
        .def("__sub__", [](const CollectiveExpression& rSelf, const double Value) { return rSelf - Value; })
        .def("__isub__", [](CollectiveExpression& rSelf, const double Value) { rSelf = rSelf - Value; return rSelf; })
        .def("__mul__", [](const CollectiveExpression& rSelf, const CollectiveExpression& rOther) { return rSelf * rOther; })
        .def("__imul__", [](CollectiveExpression& rSelf, const CollectiveExpression& rOther) { rSelf = rSelf * rOther; return rSelf; })
        .def("__mul__", [](const CollectiveExpression& rSelf, const double Value) { return rSelf * Value; })
        .def("__imul__", [](CollectiveExpression& rSelf, const double Value) { rSelf = rSelf * Value; return rSelf; })
        .def("__truediv__", [](const CollectiveExpression& rSelf, const CollectiveExpression& rOther) { return rSelf / rOther; })
        .def("__itruediv__", [](CollectiveExpression& rSelf, const CollectiveExpression& rOther) { rSelf = rSelf / rOther; return rSelf; })
        .def("__truediv__", [](const CollectiveExpression& rSelf, const double Value) { return rSelf / Value; })
        .def("__itruediv__", [](CollectiveExpression& rSelf, const double Value) { rSelf = rSelf / Value; return rSelf; })
        .def("__pow__", [](CollectiveExpression& rSelf, const CollectiveExpression& rInput) { CollectiveExpression result; result = ContainerExpressionUtils::Pow(rSelf, rInput); return result; })
        .def("__ipow__", [](CollectiveExpression& rSelf, const CollectiveExpression& rInput) { rSelf = ContainerExpressionUtils::Pow(rSelf, rInput); return rSelf; })
        .def("__pow__", [](CollectiveExpression& rSelf, const double Value) { CollectiveExpression result; result = ContainerExpressionUtils::Pow(rSelf, Value); return result; })
        .def("__ipow__", [](CollectiveExpression& rSelf, const double Value) { rSelf = ContainerExpressionUtils::Pow(rSelf, Value); return rSelf; })
        .def("__neg__", [](CollectiveExpression& rSelf) { return rSelf * -1.0; })
        .def("__str__", &CollectiveExpression::Info)
        ;



    m.def_submodule("ExpressionUtils")
        .def("Collapse", &ContainerExpressionUtils::Collapse, py::arg("collective_expressions"))
        .def("Abs", &ContainerExpressionUtils::Abs, py::arg("collective_expressions"))
        .def("EntityMin", &ContainerExpressionUtils::EntityMin, py::arg("collective_expressions"))
        .def("EntityMax", &ContainerExpressionUtils::EntityMax, py::arg("collective_expressions"))
        .def("EntitySum", &ContainerExpressionUtils::EntitySum, py::arg("collective_expressions"))
        .def("Sum", &ContainerExpressionUtils::Sum, py::arg("collective_expressions"))
        .def("NormInf", &ContainerExpressionUtils::NormInf, py::arg("collective_expressions"))
        .def("NormL2", &ContainerExpressionUtils::NormL2, py::arg("collective_expressions"))
        .def("Pow", py::overload_cast<const CollectiveExpression&, const double>(&ContainerExpressionUtils::Pow), py::arg("collective_expression"), py::arg("power_coeff"))
        .def("Pow", py::overload_cast<const CollectiveExpression&, const CollectiveExpression&>(&ContainerExpressionUtils::Pow), py::arg("collective_expression"), py::arg("power_coeff_collective_expression"))
        .def("Scale", py::overload_cast<const CollectiveExpression&, const double>(&ContainerExpressionUtils::Scale), py::arg("collective_expression"), py::arg("scaling_coeff"))
        .def("Scale", py::overload_cast<const CollectiveExpression&, const CollectiveExpression&>(&ContainerExpressionUtils::Scale), py::arg("collective_expression"), py::arg("scaling_coeff_collective_expression"))
        .def("InnerProduct", &ContainerExpressionUtils::InnerProduct, py::arg("collective_expressions_1"), py::arg("collective_expressions_2"))
        .def("EntityMaxNormL2", &ContainerExpressionUtils::EntityMaxNormL2<ModelPart::NodesContainerType>, py::arg("container_expression"))
        .def("EntityMaxNormL2", &ContainerExpressionUtils::EntityMaxNormL2<ModelPart::ConditionsContainerType>, py::arg("container_expression"))
        .def("EntityMaxNormL2", &ContainerExpressionUtils::EntityMaxNormL2<ModelPart::ElementsContainerType>, py::arg("container_expression"))
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

    auto collective_expression_io = m.def_submodule("CollectiveExpressionIO");
    py::class_<CollectiveExpressionIO::HistoricalVariable, CollectiveExpressionIO::HistoricalVariable::Pointer>(collective_expression_io, "HistoricalVariable")
        .def(py::init<const CollectiveExpressionIO::VariableType&>(), py::arg("variable"));
    py::class_<CollectiveExpressionIO::NonHistoricalVariable, CollectiveExpressionIO::NonHistoricalVariable::Pointer>(collective_expression_io, "NonHistoricalVariable")
        .def(py::init<const CollectiveExpressionIO::VariableType&>(), py::arg("variable"));
    py::class_<CollectiveExpressionIO::PropertiesVariable, CollectiveExpressionIO::PropertiesVariable::Pointer>(collective_expression_io, "PropertiesVariable")
        .def(py::init<const CollectiveExpressionIO::VariableType&>(), py::arg("variable"));

    collective_expression_io.def("Read", [](CollectiveExpression& rCExpression, const CollectiveExpressionIO::ContainerVariableType& rContainerVariable){ CollectiveExpressionIO::Read(rCExpression, rContainerVariable); }, py::arg("collective_expression"), py::arg("variable_container"));
    collective_expression_io.def("Read", [](CollectiveExpression& rCExpression, const std::vector<CollectiveExpressionIO::ContainerVariableType>& rContainerVariable){ CollectiveExpressionIO::Read(rCExpression, rContainerVariable); }, py::arg("collective_expression"), py::arg("list_of_variable_containers"));
    collective_expression_io.def("Read", [](CollectiveExpression& rCExpression, const py::array_t<double>& rData) {
        KRATOS_ERROR_IF(rData.ndim() == 0) << "Passed data is not compatible.\n";

        // get the shape
        std::vector<int> current_shape;
        for (int i = 1; i < rData.ndim(); ++i) {
            current_shape.push_back(rData.shape(i));
        }

        std::vector<std::vector<int>> list_of_shapes(rCExpression.GetContainerExpressions().size(), current_shape);

        Detail::CollectiveExpressionFromPythonArray(rCExpression, rData, list_of_shapes, rData.data());
    }, py::arg("collective_expression"), py::arg("numpy_data_array").noconvert());
    collective_expression_io.def("Read", [](CollectiveExpression& rCExpression, const py::array_t<double>& rData, const std::vector<std::vector<int>>& rListOfShapes) { Detail::CollectiveExpressionFromPythonArray(rCExpression, rData, rListOfShapes, rData.data()); }, py::arg("collective_expression"), py::arg("numpy_data_array").noconvert(), py::arg("list_of_shapes"));
    collective_expression_io.def("Move", [](CollectiveExpression& rCExpression, py::array_t<double>& rData) {
        KRATOS_ERROR_IF(rData.ndim() == 0) << "Passed data is not compatible.\n";

        // get the shape
        std::vector<int> current_shape;
        for (int i = 1; i < rData.ndim(); ++i) {
            current_shape.push_back(rData.shape(i));
        }

        std::vector<std::vector<int>> list_of_shapes(rCExpression.GetContainerExpressions().size(), current_shape);

        Detail::CollectiveExpressionFromPythonArray(rCExpression, rData, list_of_shapes, rData.mutable_data());
    }, py::arg("collective_expression"), py::arg("numpy_data_array").noconvert());
    collective_expression_io.def("Move", [](CollectiveExpression& rCExpression, py::array_t<double>& rData, const std::vector<std::vector<int>>& rListOfShapes) { Detail::CollectiveExpressionFromPythonArray(rCExpression, rData, rListOfShapes, rData.mutable_data()); }, py::arg("collective_expression"), py::arg("numpy_data_array").noconvert(), py::arg("list_of_shapes"));
    collective_expression_io.def("Write", [](CollectiveExpression& rCExpression, const CollectiveExpressionIO::ContainerVariableType& rContainerVariable){ CollectiveExpressionIO::Write(rCExpression, rContainerVariable); }, py::arg("collective_expression"), py::arg("variable_container"));
    collective_expression_io.def("Write", [](CollectiveExpression& rCExpression, const std::vector<CollectiveExpressionIO::ContainerVariableType>& rContainerVariable){ CollectiveExpressionIO::Write(rCExpression, rContainerVariable); }, py::arg("collective_expression"), py::arg("list_of_variable_containers"));

    m.def_submodule("ImplicitFilterUtils")
        .def("CalculateNodeNeighbourCount", &ImplicitFilterUtils::CalculateNodeNeighbourCount, py::arg("input_model_part"))
        .def("SetBulkRadiusForShapeFiltering", &ImplicitFilterUtils::SetBulkRadiusForShapeFiltering, py::arg("input_model_part"))
        .def("AssignProperties", &ImplicitFilterUtils::AssignProperties, py::arg("model_part"), py::arg("properties_parameters"))
        ;


    auto properties_variable_expression_io = m.def_submodule("PropertiesVariableExpressionIO");
    properties_variable_expression_io.def("Read", &PropertiesVariableExpressionIO::Read<ModelPart::ConditionsContainerType>, py::arg("condition_container_expression"), py::arg("variable"));
    properties_variable_expression_io.def("Read", &PropertiesVariableExpressionIO::Read<ModelPart::ElementsContainerType>, py::arg("element_container_expression"), py::arg("variable"));
    properties_variable_expression_io.def("Check", &PropertiesVariableExpressionIO::Check<ModelPart::ConditionsContainerType>, py::arg("condition_container_expression"), py::arg("variable"));
    properties_variable_expression_io.def("Check", &PropertiesVariableExpressionIO::Check<ModelPart::ElementsContainerType>, py::arg("element_container_expression"), py::arg("variable"));
    properties_variable_expression_io.def("Write", &PropertiesVariableExpressionIO::Write<ModelPart::ConditionsContainerType>, py::arg("condition_container_expression"), py::arg("variable"));
    properties_variable_expression_io.def("Write", &PropertiesVariableExpressionIO::Write<ModelPart::ElementsContainerType>, py::arg("element_container_expression"), py::arg("variable"));

    py::class_<PropertiesVariableExpressionIO::Input, PropertiesVariableExpressionIO::Input::Pointer, ExpressionInput>(properties_variable_expression_io, "Input")
        .def(py::init<const ModelPart&,
                      const PropertiesVariableExpressionIO::VariableType&,
                      Globals::DataLocation>(),
             py::arg("model_part"),
             py::arg("variable"),
             py::arg("data_location"))
        ;

    py::class_<PropertiesVariableExpressionIO::Output, PropertiesVariableExpressionIO::Output::Pointer, ExpressionOutput>(properties_variable_expression_io, "Output")
        .def(py::init<ModelPart&,
                      const PropertiesVariableExpressionIO::VariableType&,
                      Globals::DataLocation>(),
             py::arg("model_part"),
             py::arg("variable"),
             py::arg("data_location"))
        ;
}

}  // namespace Python.
} // Namespace Kratos

