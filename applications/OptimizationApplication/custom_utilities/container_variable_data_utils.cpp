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
#include <cmath>
#include <variant>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"
#include "utilities/atomic_utilities.h"

// Application includes
#include "custom_utilities/container_variable_data/container_data_io.h"
#include "custom_utilities/container_variable_data/specialized_container_variable_data.h"
#include "custom_utilities/container_variable_data/variable_expression_data_io.h"
#include "custom_utilities/container_variable_data/expressions/expression.h"
#include "custom_utilities/container_variable_data/expressions/literal/literal_expression.h"
#include "custom_utilities/container_variable_data/expressions/literal/literal_flat_expression.h"
#include "custom_utilities/container_variable_data/expressions/binary/binary_expression.h"
#include "custom_utilities/container_variable_data/expressions/binary/binary_weighted_multiplication_expression.h"
#include "optimization_application_variables.h"

// Include base h
#include "container_variable_data_utils.h"

namespace Kratos
{

namespace ContainerVariableDataHolderUtilsHelper
{
using VariableVariantType = std::variant<const Variable<double>*, const Variable<array_1d<double, 3>>*>;
using VariablePairVariantType = std::variant<
                                    std::pair<const Variable<double>*, const Variable<double>*>,
                                    std::pair<const Variable<array_1d<double, 3>>*, const Variable<array_1d<double, 3>>*>>;

VariableVariantType GetTemporaryVariable(const std::vector<IndexType>& rShape)
{
    if (rShape.size() == 0) {
        return &TEMPORARY_SCALAR_VARIABLE_1;
    } else if (rShape == std::vector<IndexType>{3}) {
        return &TEMPORARY_ARRAY3_VARIABLE_1;
    } else {
        KRATOS_ERROR << "Unsupported data shape = "
                     << rShape << ". Only scalar and array3 data shapes are supported for temporary variable retrieval.\n";
    }
    return &TEMPORARY_SCALAR_VARIABLE_1;
}

VariablePairVariantType GetTemporaryVariable1And2(const std::vector<IndexType>& rShape)
{
    if (rShape.size() == 0) {
        return std::make_pair(&TEMPORARY_SCALAR_VARIABLE_1, &TEMPORARY_SCALAR_VARIABLE_2);
    } else if (rShape == std::vector<IndexType>{3}) {
        return std::make_pair(&TEMPORARY_ARRAY3_VARIABLE_2, &TEMPORARY_ARRAY3_VARIABLE_2);
    } else {
        KRATOS_ERROR << "Unsupported data shape = "
                     << rShape << ". Only scalar and array3 data shapes are supported for temporary variable retrieval.\n";
    }
    return std::make_pair(&TEMPORARY_SCALAR_VARIABLE_1, &TEMPORARY_SCALAR_VARIABLE_2);
}

template<class TDataType>
typename VariableExpressionDataIO<TDataType>::Pointer GetVariableExpressionDataIO(
    const Variable<TDataType>& rVaraible,
    const std::vector<IndexType>& rShape)
{
    return VariableExpressionDataIO<TDataType>::Create(rShape);
}

void ComputeMatrixExpressionProduct(
    LiteralFlatExpression& rOutputExpression,
    const Matrix& rMatrix,
    const LiteralFlatExpression& rInputExpression,
    const IndexType NumberOfEntities,
    const IndexType ExpressionLocalSize)
{
    KRATOS_TRY

    const IndexType row_local_size = rMatrix.size1() / NumberOfEntities;
    const IndexType col_local_size = rMatrix.size2() / NumberOfEntities;

    for (IndexType i = 0; i < NumberOfEntities; ++i) {
        const IndexType matrix_row_begin = i * row_local_size;
        const IndexType output_data_begin = i * ExpressionLocalSize;

        IndexType output_component_index;

        // fill the output expression with the values from the matrix vector multiplication.
        for (output_component_index = 0; output_component_index < row_local_size; ++output_component_index) {
            for (IndexType j = 0; j < NumberOfEntities; ++j) {
                const IndexType matrix_col_begin = j * col_local_size;
                const IndexType input_data_begin = j * ExpressionLocalSize;

                for (IndexType input_component_index = 0; input_component_index < col_local_size; ++input_component_index) {
                    rOutputExpression.SetData(
                        output_data_begin, output_component_index,
                        rMatrix(matrix_row_begin + output_component_index,
                                matrix_col_begin + input_component_index) *
                            rInputExpression.Evaluate(input_data_begin, input_component_index));
                }
            }
        }

        // set the rest of the components in 2D to zero because, the expression
        // local size is 3 even for 2D in array3 case.
        for (; output_component_index < ExpressionLocalSize; ++output_component_index) {
            rOutputExpression.SetData(output_data_begin, output_component_index, 0.0);
        }
    }

    KRATOS_CATCH("");
}
} // namespace ContainerVariableDataHolderUtilsHelper

template<class TContainerType>
void ContainerVariableDataUtils::Pow(
    ContainerVariableData<TContainerType>& rOutputContainer,
    const ContainerVariableData<TContainerType>& rInputContainer1,
    const ContainerVariableData<TContainerType>& rInputContainer2)
{
    KRATOS_ERROR_IF(&rOutputContainer.GetModelPart() != &rInputContainer1.GetModelPart())
        << "Output container " << rOutputContainer.GetModelPart().FullName()
        << " mismatch with input container 1 "
        << rInputContainer1.GetModelPart().FullName() << ".\n";

    KRATOS_ERROR_IF(&rOutputContainer.GetModelPart() != &rInputContainer2.GetModelPart())
        << "Output container " << rOutputContainer.GetModelPart().FullName()
        << " mismatch with input container 2 "
        << rInputContainer2.GetModelPart().FullName() << ".\n";

    rOutputContainer.SetExpression(BinaryExpression<BinaryOperations::Power>::Create(rInputContainer1.pGetExpression(), rInputContainer2.pGetExpression()));
}

template<class TContainerType>
void ContainerVariableDataUtils::Pow(
    ContainerVariableData<TContainerType>& rOutputContainer,
    const ContainerVariableData<TContainerType>& rInputContainer1,
    const double Value)
{
    KRATOS_ERROR_IF(&rOutputContainer.GetModelPart() != &rInputContainer1.GetModelPart())
        << "Output container " << rOutputContainer.GetModelPart().FullName()
        << " mismatch with input container 1 "
        << rInputContainer1.GetModelPart().FullName() << ".\n";

    rOutputContainer.SetExpression(BinaryExpression<BinaryOperations::Power>::Create(rInputContainer1.pGetExpression(), LiteralExpression<double>::Create(Value)));
}

template<class TContainerType>
void ContainerVariableDataUtils::WeightedProduct(
    ContainerVariableData<TContainerType>& rOutputContainer,
    const ContainerVariableData<TContainerType>& rInputContainer,
    const ContainerVariableData<TContainerType>& rWeightsContainer)
{
    KRATOS_ERROR_IF(&rOutputContainer.GetModelPart() != &rInputContainer.GetModelPart())
        << "Output container " << rOutputContainer.GetModelPart().FullName()
        << " mismatch with input container "
        << rInputContainer.GetModelPart().FullName() << ".\n";

    KRATOS_ERROR_IF(&rOutputContainer.GetModelPart() != &rWeightsContainer.GetModelPart())
        << "Output container " << rOutputContainer.GetModelPart().FullName()
        << " mismatch with weights container "
        << rWeightsContainer.GetModelPart().FullName() << ".\n";

    rOutputContainer.SetExpression(BinaryWeightedMultiplicationExpression::Create(rInputContainer.pGetExpression(), rWeightsContainer.pGetExpression()));
}

template<class TContainerType>
double ContainerVariableDataUtils::EntityMaxNormL2(const ContainerVariableData<TContainerType>& rContainer)
{
    if (rContainer.GetLocalSize() == 0) {
        return 0.0;
    }

    const auto& r_expression = rContainer.GetExpression();
    const IndexType local_size = rContainer.GetLocalSize();
    const IndexType number_of_entities = rContainer.GetContainer().size();

    return std::sqrt(rContainer.GetModelPart().GetCommunicator().GetDataCommunicator().MaxAll(IndexPartition<IndexType>(number_of_entities).for_each<MaxReduction<double>>([&r_expression, local_size](const IndexType EntityIndex) {
        const IndexType local_data_begin_index = EntityIndex * local_size;
        double value = 0.0;
        for (IndexType i = 0; i < local_size; ++i) {
            value += std::pow(r_expression.Evaluate(local_data_begin_index, i), 2);
        }
        return value;
    })));
}

template<class TContainerType>
double ContainerVariableDataUtils::NormInf(const ContainerVariableData<TContainerType>& rContainer)
{
    const auto& r_expression = rContainer.GetExpression();
    const IndexType local_size = rContainer.GetLocalSize();
    const IndexType number_of_entities = rContainer.GetContainer().size();

    return rContainer.GetModelPart().GetCommunicator().GetDataCommunicator().MaxAll(IndexPartition<IndexType>(number_of_entities).for_each<MaxReduction<double>>([&r_expression, local_size](const IndexType EntityIndex) {
        const IndexType local_data_begin_index = EntityIndex * local_size;
        double value = 0.0;
        for (IndexType i = 0; i < local_size; ++i) {
            value = std::max(value, r_expression.Evaluate(local_data_begin_index, i));
        }
        return value;
    }));
}

template<class TContainerType>
double ContainerVariableDataUtils::NormL2(
    const ContainerVariableData<TContainerType>& rContainer)
{
    const auto& r_expression = rContainer.GetExpression();
    const IndexType local_size = rContainer.GetLocalSize();
    const IndexType number_of_entities = rContainer.GetContainer().size();

    const double local_l2_norm_square = IndexPartition<IndexType>(number_of_entities).for_each<SumReduction<double>>([&r_expression, local_size](const IndexType EntityIndex) {
        const IndexType local_data_begin_index = EntityIndex * local_size;
        double value = 0.0;
        for (IndexType i = 0; i < local_size; ++i) {
            value += std::pow(r_expression.Evaluate(local_data_begin_index, i), 2);
        }
        return value;
    });

    return std::sqrt(rContainer.GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(local_l2_norm_square));
}

template<class TContainerType>
double ContainerVariableDataUtils::InnerProduct(
    const ContainerVariableData<TContainerType>& rContainer1,
    const ContainerVariableData<TContainerType>& rContainer2)
{
    const auto& r_expression_1 = rContainer1.GetExpression();
    const IndexType local_size_1 = rContainer1.GetLocalSize();
    const IndexType number_of_entities_1 = rContainer1.GetContainer().size();

    const auto& r_expression_2 = rContainer2.GetExpression();
    const IndexType local_size_2 = rContainer2.GetLocalSize();
    const IndexType number_of_entities_2 = rContainer2.GetContainer().size();

    KRATOS_ERROR_IF(local_size_1 != local_size_2)
        << "Data dimension mismatch in InnerProduct calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rContainer1 << "\n"
        << "   Container 2: " << rContainer2 << "\n";

    KRATOS_ERROR_IF(number_of_entities_1 != number_of_entities_2)
        << "Number of entities mismatch in InnerProduct calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rContainer1 << "\n"
        << "   Container 2: " << rContainer2 << "\n";

    KRATOS_ERROR_IF(&rContainer1.GetModelPart() != &rContainer2.GetModelPart())
        << "Model part mismatch in InnerProduct calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rContainer1 << "\n"
        << "   Container 2: " << rContainer2 << "\n";

    return rContainer1.GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(IndexPartition<IndexType>(number_of_entities_1).for_each<SumReduction<double>>([&r_expression_1, &r_expression_2, local_size_1](const IndexType EntityIndex) {
        const IndexType local_data_begin_index = EntityIndex * local_size_1;
        double value = 0.0;
        for (IndexType i = 0; i < local_size_1; ++i) {
            value += r_expression_1.Evaluate(local_data_begin_index, i) * r_expression_2.Evaluate(local_data_begin_index, i);
        }
        return value;
    }));
}

template<class TContainerType>
void ContainerVariableDataUtils::ProductWithEntityMatrix(
    ContainerVariableData<TContainerType>& rOutput,
    const SparseMatrixType& rMatrix,
    const ContainerVariableData<TContainerType>& rInput)
{
    KRATOS_ERROR_IF(rInput.GetModelPart().IsDistributed() ||
                    rOutput.GetModelPart().IsDistributed())
        << "ProductWithEntityMatrix does not support MPI yet.\n";

    const IndexType number_of_output_entities = rOutput.GetContainer().size();
    const IndexType number_of_input_entities = rInput.GetContainer().size();

    KRATOS_ERROR_IF_NOT(number_of_input_entities == rMatrix.size2())
        << "Input container size and matrix size2 mismatch. [ Container size = "
        << number_of_input_entities << ", Matrix.size1() = " << rMatrix.size2()
        << " ]. Followings are the given containers:"
        << "\n\tInput data container : " << rInput
        << "\n\tOutput data container: " << rOutput << "\n\t";

    KRATOS_ERROR_IF_NOT(number_of_output_entities == rMatrix.size1())
        << "Output container size and matrix size1 mismatch. [ Container size = "
        << number_of_output_entities << ", Matrix.size1() = " << rMatrix.size1()
        << " ]. Followings are the given containers:"
        << "\n\tInput data container : " << rInput
        << "\n\tOutput data container: " << rOutput << "\n\t";

    auto p_flat_data_expression = LiteralFlatExpression::Create(number_of_output_entities, rInput.GetShape());
    rOutput.SetExpression(p_flat_data_expression);

    const IndexType local_size = rInput.GetLocalSize();
    const auto& r_input_expression = rInput.GetExpression();
    auto& r_output_expression = *p_flat_data_expression;

    const double* a_values = rMatrix.value_data().begin();
    const IndexType* a_row_indices = rMatrix.index1_data().begin();
    const IndexType* a_col_indices = rMatrix.index2_data().begin();

    IndexPartition<IndexType>(rMatrix.size1()).for_each([a_values, &r_input_expression, &r_output_expression, local_size, a_row_indices, a_col_indices](const IndexType i) {
        const IndexType col_begin = a_row_indices[i];
        const IndexType col_end = a_row_indices[i + 1];
        const IndexType local_data_begin_index = i * local_size;

        for (IndexType d = 0; d < local_size; ++d) {
            double result = 0.0;
            for (IndexType j = col_begin; j < col_end; ++j) {
                result += a_values[j] * r_input_expression.Evaluate(a_col_indices[j] * local_size, d);
            }
            r_output_expression.SetData(local_data_begin_index, d, result);
        }
    });
}

template<class TContainerType>
void ContainerVariableDataUtils::ProductWithEntityMatrix(
    ContainerVariableData<TContainerType>& rOutput,
    const Matrix& rMatrix,
    const ContainerVariableData<TContainerType>& rInput)
{
    KRATOS_ERROR_IF(rInput.GetModelPart().IsDistributed() ||
                    rOutput.GetModelPart().IsDistributed())
        << "ProductWithEntityMatrix does not support MPI yet.\n";

    const IndexType number_of_output_entities = rOutput.GetContainer().size();
    const IndexType number_of_input_entities = rInput.GetContainer().size();

    KRATOS_ERROR_IF_NOT(number_of_input_entities == rMatrix.size2())
        << "Input container size and matrix size2 mismatch. [ Container size = "
        << number_of_input_entities << ", Matrix.size1() = " << rMatrix.size2()
        << " ]. Followings are the given containers:"
        << "\n\tInput data container : " << rInput
        << "\n\tOutput data container: " << rOutput << "\n\t";

    KRATOS_ERROR_IF_NOT(number_of_output_entities == rMatrix.size1())
        << "Output container size and matrix size1 mismatch. [ Container size = "
        << number_of_output_entities << ", Matrix.size1() = " << rMatrix.size1()
        << " ]. Followings are the given containers:"
        << "\n\tInput data container : " << rInput
        << "\n\tOutput data container: " << rOutput << "\n\t";

    auto p_flat_data_expression = LiteralFlatExpression::Create(number_of_output_entities, rInput.GetShape());
    rOutput.SetExpression(p_flat_data_expression);

    const auto& r_input_expression = rInput.GetExpression();
    const IndexType local_size = rInput.GetLocalSize();
    auto& r_output_expression = *p_flat_data_expression;

    IndexPartition<IndexType>(rMatrix.size1()).for_each([&rMatrix, &r_input_expression, &r_output_expression, local_size](const IndexType i) {
        const IndexType local_data_begin_index = i * local_size;
        for (IndexType d = 0; d < local_size; ++d) {
            double result = 0.0;
            for (IndexType j = 0; j < rMatrix.size2(); ++j) {
                result += rMatrix(i, j) * r_input_expression.Evaluate(j * local_size, d);;
            }
            r_output_expression.SetData(local_data_begin_index, d, result);
        }
    });
}

void ContainerVariableDataUtils::Transpose(
    Matrix& rOutput,
    const Matrix& rInput)
{
    KRATOS_ERROR_IF(&rOutput == &rInput)
        << "Please provide a matrix different to input matrix for output.";

    if (rOutput.size1() != rInput.size2() || rOutput.size2() != rInput.size1()) {
        rOutput.resize(rInput.size2(), rInput.size1(), false);
    }

    IndexPartition<IndexType>(rInput.size1()).for_each([&rInput, &rOutput](const IndexType i){
        for (IndexType j = 0; j < rInput.size2(); ++j) {
            rOutput(j, i) = rInput(i, j);
        }
    });
}

void ContainerVariableDataUtils::Transpose(
    SparseMatrixType& rOutput,
    const SparseMatrixType& rInput)
{
    if (rOutput.size1() != rInput.size2() || rOutput.size2() != rInput.size1()) {
        rOutput.resize(rInput.size2(), rInput.size1(), false);
    }

    const double* a_values = rInput.value_data().begin();
    const IndexType* a_row_indices = rInput.index1_data().begin();
    const IndexType* a_col_indices = rInput.index2_data().begin();

    for (IndexType i = 0; i < rInput.size1(); ++i) {
        const IndexType col_begin = a_row_indices[i];
        const IndexType col_end = a_row_indices[i + 1];

        for (IndexType j = col_begin; j < col_end; ++j) {
            rOutput.insert_element(a_col_indices[j], i, a_values[j]);
        }
    }
}

template<class TContainerType>
void ContainerVariableDataUtils::ComputeNumberOfNeighbourEntities(
    ContainerVariableData<ModelPart::NodesContainerType>& rOutput)
{
    // clear the neighbour count storing variable
    VariableUtils().SetNonHistoricalVariableToZero(TEMPORARY_SCALAR_VARIABLE_1, rOutput.GetModelPart().Nodes());

    // create a dummy copy input data container to access its nodes and modify data
    SpecializedContainerVariableData<TContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>> dummy_input_container(rOutput.GetModelPart());

    block_for_each(dummy_input_container.GetContainer(), [](auto& rEntity) {
        auto& r_geometry = rEntity.GetGeometry();
        for (auto& r_node : r_geometry) {
            AtomicAdd(r_node.GetValue(TEMPORARY_SCALAR_VARIABLE_1), 1.0);
        }
    });

    rOutput.GetModelPart().GetCommunicator().AssembleNonHistoricalData(TEMPORARY_SCALAR_VARIABLE_1);

    // now read in the nodal data
    SpecializedContainerVariableData<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>> dummy_read(rOutput.GetModelPart());
    dummy_read.ReadData(TEMPORARY_SCALAR_VARIABLE_1);

    // now fill the rOutput
    rOutput.CopyDataFrom(dummy_read);
}

template<class TContainerType>
void ContainerVariableDataUtils::MapContainerVariableToNodalVariable(
    ContainerVariableData<ModelPart::NodesContainerType>& rOutput,
    const ContainerVariableData<TContainerType>& rInput,
    const ContainerVariableData<ModelPart::NodesContainerType>& rNeighbourEntities)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(&rOutput.GetModelPart() != &rInput.GetModelPart())
        << "Output container and input container model parts mismatch. "
           "Followings are the container details:"
        << "\n\tOutput container: " << rOutput
        << "\n\tInput container : " << rInput << "\n";

    KRATOS_ERROR_IF(&rOutput.GetModelPart() != &rNeighbourEntities.GetModelPart())
        << "Output container and neighbour entities container model parts mismatch. "
           "Followings are the container details:"
        << "\n\tOutput container            : " << rOutput
        << "\n\tNeighbour entities container: " << rNeighbourEntities << "\n";

    KRATOS_ERROR_IF(rNeighbourEntities.GetLocalSize() != 1)
        << "Neighbour entities container can only have data with dimensionality = 1"
        << "\n\tNeighbour entities container:" << rNeighbourEntities << "\n";

    // reset temporary variables
    std::visit([&rOutput](auto&& p_variable) {
        VariableUtils().SetNonHistoricalVariableToZero(*p_variable, rOutput.GetModelPart().Nodes());
    }, ContainerVariableDataHolderUtilsHelper::GetTemporaryVariable(rInput.GetShape()));

    // copy number of neighbours
    SpecializedContainerVariableData<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>> dummy_weights(rNeighbourEntities);

    // assign dummy weights to nodes
    dummy_weights.AssignData(TEMPORARY_SCALAR_VARIABLE_2);

    // create a dummy copy input data container to access its nodes and modify data
    SpecializedContainerVariableData<TContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>> dummy_input_container(rOutput.GetModelPart());

    auto& r_container = dummy_input_container.GetContainer();
    auto& r_communicator = rOutput.GetModelPart().GetCommunicator();
    const IndexType number_of_entities = r_container.size();
    const auto& r_expression = rInput.GetExpression();

    // now distribute the entity values to nodes
    std::visit([&r_communicator, &r_container, &r_expression, number_of_entities](auto&& p_variable) {
        // now create the variable_expression_data_io
        auto p_variable_expression_data_io = ContainerVariableDataHolderUtilsHelper::GetVariableExpressionDataIO(*p_variable, r_expression.GetShape());

        IndexPartition<IndexType>(number_of_entities).for_each(p_variable->Zero(), [&p_variable_expression_data_io, &p_variable, &r_container, &r_expression](const IndexType EntityIndex, auto& rValue) {
            p_variable_expression_data_io->Assign(rValue, r_expression, EntityIndex);

            auto p_entity = (r_container.begin() + EntityIndex);
            auto& r_geometry = p_entity->GetGeometry();
            for (auto& r_node : r_geometry) {
                r_node.SetLock();
                r_node.GetValue(*p_variable) += rValue / r_node.GetValue(TEMPORARY_SCALAR_VARIABLE_2);
                r_node.UnSetLock();
            }
        });

        // now assemble data at nodes
        r_communicator.AssembleNonHistoricalData(*p_variable);

    }, ContainerVariableDataHolderUtilsHelper::GetTemporaryVariable(rInput.GetShape()));

    // now read in the nodal data
    SpecializedContainerVariableData<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>> dummy_read(rOutput.GetModelPart());
    std::visit([&dummy_read](auto&& p_variable) {
        dummy_read.ReadData(*p_variable);
    }, ContainerVariableDataHolderUtilsHelper::GetTemporaryVariable(rInput.GetShape()));

    // now copy back the data to the output container
    rOutput.CopyDataFrom(dummy_read);

    KRATOS_CATCH("");
}

template<class TContainerType>
void ContainerVariableDataUtils::MapNodalVariableToContainerVariable(
    ContainerVariableData<TContainerType>& rOutput,
    const ContainerVariableData<ModelPart::NodesContainerType>& rInput)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(&rOutput.GetModelPart() != &rInput.GetModelPart())
        << "Output container and input container model parts mismatch. "
           "Followings are the container details:"
        << "\n\tOutput container: " << rOutput
        << "\n\tInput container : " << rInput << "\n";

    std::visit([&rOutput, &rInput](auto&& p_variable) {
        // now create the variable_expression_data_io
        auto p_variable_expression_data_io = ContainerVariableDataHolderUtilsHelper::GetVariableExpressionDataIO(*p_variable, rInput.GetShape());

        // create a dummy with the model part and nodes
        SpecializedContainerVariableData<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>> dummy(rOutput.GetModelPart());
        dummy.CopyDataFrom(rInput);

        const auto& r_output_container = rOutput.GetContainer();
        const IndexType number_of_entities = r_output_container.size();

        // create output expression
        auto p_flat_data_expression = LiteralFlatExpression::Create(number_of_entities, rInput.GetShape());
        auto& r_flat_data_expression = *p_flat_data_expression;
        rOutput.SetExpression(p_flat_data_expression);

        // first assign input data to nodes
        dummy.AssignData(*p_variable);

        // compute the entity valeus.
        IndexPartition<IndexType>(number_of_entities).for_each([&p_variable_expression_data_io, &r_output_container, &r_flat_data_expression, &p_variable](const IndexType EntityIndex) {
            const auto p_entity = (r_output_container.begin() + EntityIndex);
            const auto& r_geometry = p_entity->GetGeometry();
            const IndexType number_of_nodes = r_geometry.size();

            auto value = r_geometry[0].GetValue(*p_variable);

            for (IndexType i = 1; i < number_of_nodes; ++i) {
                value += r_geometry[i].GetValue(*p_variable);
            }
            p_variable_expression_data_io->Read(r_flat_data_expression, EntityIndex, value / number_of_nodes);
        });
    }, ContainerVariableDataHolderUtilsHelper::GetTemporaryVariable(rInput.GetShape()));

    KRATOS_CATCH("");
}

template<class TContainerType>
void ContainerVariableDataUtils::ComputeNodalVariableProductWithEntityMatrix(
    ContainerVariableData<ModelPart::NodesContainerType>& rOutput,
    const ContainerVariableData<ModelPart::NodesContainerType>& rNodalValues,
    const Variable<Matrix>& rMatrixVariable,
    TContainerType& rEntities)
{
    KRATOS_TRY

    using tls_type = std::tuple<Vector, Vector, Matrix>;

    KRATOS_ERROR_IF(&rOutput.GetModelPart() != &rNodalValues.GetModelPart())
        << "Output container and input container model parts mismatch. "
           "Followings are the container details:"
        << "\n\tOutput container: " << rOutput
        << "\n\tNodal container : " << rNodalValues << "\n";

    SpecializedContainerVariableData<TContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>> dummy_container(rOutput.GetModelPart());

    KRATOS_ERROR_IF(dummy_container.GetContainer().size() != rEntities.size())
        << "Provided entities container size mismatch with output container variable data holder size. "
        << "[ Provided entities size = " << rEntities.size()
        << ", output data container entities size = " << rOutput.GetContainer().size() << " ].\n";

    const IndexType local_size = rNodalValues.GetLocalSize();

    std::visit([&rOutput, &rNodalValues, &rMatrixVariable, &rEntities, &dummy_container, local_size](auto&& p_variable_pair) {
        const auto& r_input_variable = *std::get<0>(p_variable_pair);
        const auto& r_output_variable = *std::get<1>(p_variable_pair);

        // now create the variable_expression_data_io
        auto p_variable_expression_data_io = ContainerVariableDataHolderUtilsHelper::GetVariableExpressionDataIO(r_input_variable, rNodalValues.GetExpression().GetShape());

        // assign nodal values
        SpecializedContainerVariableData<ModelPart::NodesContainerType, ContainerDataIO<ContainerDataIOTags::NonHistorical>> nodal_write(rNodalValues);
        nodal_write.AssignData(r_input_variable);

        // clear the output variable
        VariableUtils().SetNonHistoricalVariableToZero(r_output_variable, rOutput.GetModelPart().Nodes());

        const auto& process_info = rOutput.GetModelPart().GetProcessInfo();

        block_for_each(rEntities, tls_type(), [&p_variable_expression_data_io, &process_info, &rMatrixVariable, &r_input_variable, &r_output_variable, &rNodalValues, local_size](auto& rEntity, tls_type& rTLS) {
            Matrix& r_entity_matrix = std::get<2>(rTLS);

            auto& r_geometry = rEntity.GetGeometry();

            const IndexType number_of_nodes = r_geometry.size();

            LiteralFlatExpression input_expression(number_of_nodes, rNodalValues.GetExpression().GetShape());
            LiteralFlatExpression output_expression(number_of_nodes, rNodalValues.GetExpression().GetShape());

            for (IndexType i = 0; i < number_of_nodes; ++i) {
                p_variable_expression_data_io->Read(input_expression, i, r_geometry[i].GetValue(r_input_variable));
            }

            rEntity.Calculate(rMatrixVariable, r_entity_matrix, process_info);

            ContainerVariableDataHolderUtilsHelper::ComputeMatrixExpressionProduct(
                output_expression, r_entity_matrix, input_expression,
                number_of_nodes, local_size);

            // assign the vector to entity nodal values
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                auto& r_node = r_geometry[i];
                r_node.SetLock();
                p_variable_expression_data_io->Assign(r_node.GetValue(r_output_variable), output_expression, i);
                r_node.UnSetLock();
            }
        });

        // assemble data for MPI
        dummy_container.GetModelPart().GetCommunicator().AssembleNonHistoricalData(r_output_variable);

        // now read the data
        nodal_write.ReadData(r_output_variable);
        rOutput.CopyDataFrom(nodal_write);

    }, ContainerVariableDataHolderUtilsHelper::GetTemporaryVariable1And2(rNodalValues.GetShape()));

    KRATOS_CATCH("");
}

// template instantiations
#define KRATOS_INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ContainerType)                                                                                                                                                          \
    template void ContainerVariableDataUtils::Pow(ContainerVariableData<ContainerType>&, const ContainerVariableData<ContainerType>&, const ContainerVariableData<ContainerType>&);                                                  \
    template void ContainerVariableDataUtils::Pow(ContainerVariableData<ContainerType>&, const ContainerVariableData<ContainerType>&, const double);                                                                                 \
    template void ContainerVariableDataUtils::WeightedProduct(ContainerVariableData<ContainerType>&, const ContainerVariableData<ContainerType>&, const ContainerVariableData<ContainerType>&);                                      \
    template double ContainerVariableDataUtils::EntityMaxNormL2(const ContainerVariableData<ContainerType>&);                                                                                                                        \
    template double ContainerVariableDataUtils::NormInf(const ContainerVariableData<ContainerType>&);                                                                                                                                \
    template double ContainerVariableDataUtils::NormL2(const ContainerVariableData<ContainerType>&);                                                                                                                                 \
    template double ContainerVariableDataUtils::InnerProduct(const ContainerVariableData<ContainerType>&, const ContainerVariableData<ContainerType>&);                                                                              \
    template void ContainerVariableDataUtils::ProductWithEntityMatrix(ContainerVariableData<ContainerType>&, const typename UblasSpace<double, CompressedMatrix, Vector>::MatrixType&, const ContainerVariableData<ContainerType>&); \
    template void ContainerVariableDataUtils::ProductWithEntityMatrix(ContainerVariableData<ContainerType>&, const Matrix&, const ContainerVariableData<ContainerType>&);

#define KRATOS_INSTANTIATE_NON_NODAL_UTILITY_METHOD_FOR_CONTAINER_TYPE(ContainerType)                                                                                                                                                                                \
    template void ContainerVariableDataUtils::ComputeNumberOfNeighbourEntities<ContainerType>(ContainerVariableData<ModelPart::NodesContainerType>&);                                                                                                                \
    template void ContainerVariableDataUtils::MapContainerVariableToNodalVariable(ContainerVariableData<ModelPart::NodesContainerType>&, const ContainerVariableData<ContainerType>&, const ContainerVariableData<ModelPart::NodesContainerType>&);          \
    template void ContainerVariableDataUtils::MapNodalVariableToContainerVariable(ContainerVariableData<ContainerType>&, const ContainerVariableData<ModelPart::NodesContainerType>&);                                                                       \
    template void ContainerVariableDataUtils::ComputeNodalVariableProductWithEntityMatrix(ContainerVariableData<ModelPart::NodesContainerType>&, const ContainerVariableData<ModelPart::NodesContainerType>&, const Variable<Matrix>& rMatrixVariable, ContainerType&);


KRATOS_INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::NodesContainerType)
KRATOS_INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ConditionsContainerType)
KRATOS_INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ElementsContainerType)

KRATOS_INSTANTIATE_NON_NODAL_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ConditionsContainerType)
KRATOS_INSTANTIATE_NON_NODAL_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ElementsContainerType)

#undef KRATOS_INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE
#undef KRATOS_INSTANTIATE_NON_NODAL_UTILITY_METHOD_FOR_CONTAINER_TYPE

}