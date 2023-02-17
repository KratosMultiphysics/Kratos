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
#include "custom_utilities/container_variable_data_holder/container_variable_data_holder.h"
#include "optimization_application_variables.h"

// Include base h
#include "container_variable_data_holder_utils.h"

namespace Kratos
{

namespace ContainerVariableDataHolderUtilsHelper
{
using VariableVariantType = std::variant<const Variable<double>*, const Variable<array_1d<double, 3>>*>;
using VariablePairVariantType = std::variant<
                                    std::pair<const Variable<double>*, const Variable<double>*>,
                                    std::pair<const Variable<array_1d<double, 3>>*, const Variable<array_1d<double, 3>>*>>;

VariableVariantType GetTemporaryVariable(const IndexType DataDimension)
{
    switch (DataDimension) {
        case 1: return &TEMPORARY_SCALAR_VARIABLE_1;
        case 2:
        case 3: return &TEMPORARY_ARRAY3_VARIABLE_1;
        default: KRATOS_ERROR << "Unsupported data dimension = " << DataDimension << ". Only 2 and 3 data dimension is supported for temporary variable retrieval.\n";
    }

    return &TEMPORARY_SCALAR_VARIABLE_1;
}

VariablePairVariantType GetTemporaryVariable1And2(const IndexType DataDimension)
{
    switch (DataDimension) {
        case 1: return std::make_pair(&TEMPORARY_SCALAR_VARIABLE_1, &TEMPORARY_SCALAR_VARIABLE_2);
        case 2:
        case 3: return std::make_pair(&TEMPORARY_ARRAY3_VARIABLE_2, &TEMPORARY_ARRAY3_VARIABLE_2);
        default: KRATOS_ERROR << "Unsupported data dimension = " << DataDimension << ". Only 2 and 3 data dimension is supported for temporary variable retrieval.\n";
    }

    return std::make_pair(&TEMPORARY_SCALAR_VARIABLE_1, &TEMPORARY_SCALAR_VARIABLE_2);
}

double GetEntityData(
    const Vector& rData,
    const IndexType EntityStartIndex,
    const IndexType DataDimension,
    const Variable<double>& rVariable)
{
    return rData[EntityStartIndex];
}

array_1d<double, 3> GetEntityData(
    const Vector& rData,
    const IndexType EntityStartIndex,
    const IndexType DataDimension,
    const Variable<array_1d<double, 3>>& rVariable)
{
    array_1d<double, 3> result = ZeroVector(3);
    for (IndexType i = 0; i < DataDimension; ++i) {
        result[i] = rData[EntityStartIndex * DataDimension + i];
    }
    return result;
}

void AddToEntityData(
    Vector& rData,
    const IndexType EntityStartIndex,
    const IndexType DataDimension,
    const double& rValue)
{
    switch (DataDimension) {
        case 1:
            rData[EntityStartIndex] += rValue;
            break;
        default: KRATOS_ERROR << "Unsupported data dimension. Only 1 data dimension is supported for temporary variable entity data setting.\n";
    }
}

void AddToEntityData(
    Vector& rData,
    const IndexType EntityStartIndex,
    const IndexType DataDimension,
    const array_1d<double, 3>& rValue)
{
    switch (DataDimension) {
        case 2:
        case 3:
            for (IndexType i = 0; i < DataDimension; ++i) {
                rData[EntityStartIndex * DataDimension + i] += rValue[i];
            }
            break;
        default: KRATOS_ERROR << "Unsupported data dimension. Only 2 and 3 data dimension is supported for temporary variable entity data setting.\n";
    }
}
}

template<class TContainerType>
double ContainerVariableDataHolderUtils::EntityMaxNormL2(const ContainerVariableDataHolderBase<TContainerType>& rContainer)
{
    if (rContainer.GetDataDimension() == 0) {
        return 0.0;
    }

    const auto& r_data = rContainer.GetData();
    const IndexType data_dimension = rContainer.GetDataDimension();

    const IndexType number_of_entities = r_data.size() / data_dimension;

    return std::sqrt(rContainer.GetModelPart().GetCommunicator().GetDataCommunicator().MaxAll(IndexPartition<IndexType>(number_of_entities).for_each<MaxReduction<double>>([&](const IndexType Index) {
        double value = 0.0;
        const IndexType local_start = Index * data_dimension;
        for (IndexType i = 0; i < data_dimension; ++i) {
            value += std::pow(r_data[local_start + i], 2);
        }
        return value;
    })));
}

template<class TContainerType>
double ContainerVariableDataHolderUtils::NormInf(const ContainerVariableDataHolderBase<TContainerType>& rContainer)
{
    const auto& r_data = rContainer.GetData();
    return rContainer.GetModelPart().GetCommunicator().GetDataCommunicator().MaxAll(IndexPartition<IndexType>(r_data.size()).for_each<MaxReduction<double>>([&](const IndexType Index) {
        return r_data[Index];
    }));
}

double ContainerVariableDataHolderUtils::NormInf(
    const CollectiveVariableDataHolder& rContainer)
{
    double max_norm = 0.0;
    for (const auto& p_variable_data_container : rContainer.GetVariableDataHolders()) {
        std::visit([&](auto&& v) { max_norm = std::max(max_norm, NormInf(*v));}, p_variable_data_container);
    }

    return max_norm;
}

template<class TContainerType>
double ContainerVariableDataHolderUtils::NormL2(
    const ContainerVariableDataHolderBase<TContainerType>& rContainer)
{
    const double local_l2_norm_square = IndexPartition<IndexType>(rContainer.GetData().size()).for_each<SumReduction<double>>([&](const IndexType Index) {
        return std::pow(rContainer.GetData()[Index], 2);
    });

    return std::sqrt(rContainer.GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(local_l2_norm_square));
}

double ContainerVariableDataHolderUtils::NormL2(
    const CollectiveVariableDataHolder& rContainer)
{
    double l2_norm_square = 0.0;
    for (const auto& p_variable_data_container : rContainer.GetVariableDataHolders()) {
        std::visit([&](auto&& v) { l2_norm_square += std::pow(NormL2(*v), 2);}, p_variable_data_container);
    }

    return std::sqrt(l2_norm_square);
}

template<class TContainerType>
double ContainerVariableDataHolderUtils::InnerProduct(
    const ContainerVariableDataHolderBase<TContainerType>& rContainer1,
    const ContainerVariableDataHolderBase<TContainerType>& rContainer2)
{
    const auto& r_data_1 = rContainer1.GetData();
    const auto& r_data_2 = rContainer2.GetData();

    KRATOS_ERROR_IF(r_data_1.size() != r_data_2.size())
        << "Data size mismatch in InnerProduct calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rContainer1 << "\n"
        << "   Container 2: " << rContainer2 << "\n";

    KRATOS_ERROR_IF(rContainer1.GetModelPart() != rContainer2.GetModelPart())
        << "Model part mismatch in InnerProduct calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rContainer1 << "\n"
        << "   Container 2: " << rContainer2 << "\n";

    return rContainer1.GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(IndexPartition<IndexType>(r_data_1.size()).for_each<SumReduction<double>>([&](const IndexType Index) {
        return r_data_1[Index] * r_data_2[Index];
    }));
}

double ContainerVariableDataHolderUtils::InnerProduct(
    const CollectiveVariableDataHolder& rContainer1,
    const CollectiveVariableDataHolder& rContainer2)
{
    KRATOS_ERROR_IF_NOT(rContainer1.IsCompatibleWith(rContainer2))
        << "Unsupported collective variable data holders provided for \"+\" operation."
        << "\nLeft operand : " << rContainer1 << "\nRight operand: " << rContainer2 << std::endl;

    double inner_product_value = 0.0;
    for (IndexType i = 0; i < rContainer1.GetVariableDataHolders().size(); ++i) {
        std::visit([&](auto&& v) {
            using v_type = std::decay_t<decltype(v)>;
            inner_product_value += InnerProduct(*v, *std::get<v_type>(rContainer2.GetVariableDataHolders()[i]));
        }, rContainer1.GetVariableDataHolders()[i]);
    }
    return inner_product_value;
}

template<class TContainerType>
void ContainerVariableDataHolderUtils::ProductWithEntityMatrix(
    ContainerVariableDataHolderBase<TContainerType>& rOutput,
    const SparseMatrixType& rMatrix,
    const ContainerVariableDataHolderBase<TContainerType>& rInput)
{
    KRATOS_ERROR_IF(rInput.GetModelPart().IsDistributed() ||
                    rOutput.GetModelPart().IsDistributed())
        << "ProductWithEntityMatrix does not support MPI yet.\n";

    KRATOS_ERROR_IF_NOT(rInput.GetContainer().size() == rMatrix.size2())
        << "Input container size and matrix size2 mismatch. [ Container size = "
        << rOutput.GetContainer().size() << ", Matrix.size1() = " << rMatrix.size2()
        << "Followings are the given containers:"
        << "\n\tInput data container : " << rInput
        << "\n\tOutput data container: " << rOutput << "\n\t";

    KRATOS_ERROR_IF_NOT(rOutput.GetContainer().size() == rMatrix.size1())
        << "Output container size and matrix size1 mismatch. [ Container size = "
        << rOutput.GetContainer().size() << ", Matrix.size1() = " << rMatrix.size1()
        << "Followings are the given containers:"
        << "\n\tInput data container : " << rInput
        << "\n\tOutput data container: " << rOutput << "\n\t";

    const IndexType data_dimension = rInput.GetDataDimension();

    rOutput.SetDataToZero(data_dimension);

    auto& r_output_data = rOutput.GetData();
    const auto& r_input_data = rInput.GetData();

    const double* a_values = rMatrix.value_data().begin();
    const IndexType* a_row_indices = rMatrix.index1_data().begin();
    const IndexType* a_col_indices = rMatrix.index2_data().begin();

    IndexPartition<IndexType>(rMatrix.size1()).for_each([&](const IndexType i) {
        const IndexType col_begin = a_row_indices[i];
        const IndexType col_end = a_row_indices[i + 1];

        for (IndexType d = 0; d < data_dimension; ++d) {
            auto& r_value = r_output_data[i * data_dimension + d];
            for (IndexType j = col_begin; j < col_end; ++j) {
                r_value += a_values[j] * r_input_data[a_col_indices[j] * data_dimension + d];
            }
        }
    });
}

template<class TContainerType>
void ContainerVariableDataHolderUtils::ProductWithEntityMatrix(
    ContainerVariableDataHolderBase<TContainerType>& rOutput,
    const Matrix& rMatrix,
    const ContainerVariableDataHolderBase<TContainerType>& rInput)
{
    KRATOS_ERROR_IF(rInput.GetModelPart().IsDistributed() ||
                    rOutput.GetModelPart().IsDistributed())
        << "ProductWithEntityMatrix does not support MPI yet.\n";

    KRATOS_ERROR_IF_NOT(rInput.GetContainer().size() == rMatrix.size2())
        << "Input container size and matrix size2 mismatch. [ Container size = "
        << rOutput.GetContainer().size() << ", Matrix.size1() = " << rMatrix.size2()
        << "Followings are the given containers:"
        << "\n\tInput data container : " << rInput
        << "\n\tOutput data container: " << rOutput << "\n\t";

    KRATOS_ERROR_IF_NOT(rOutput.GetContainer().size() == rMatrix.size1())
        << "Output container size and matrix size1 mismatch. [ Container size = "
        << rOutput.GetContainer().size() << ", Matrix.size1() = " << rMatrix.size1()
        << "Followings are the given containers:"
        << "\n\tInput data container : " << rInput
        << "\n\tOutput data container: " << rOutput << "\n\t";

    const IndexType data_dimension = rInput.GetDataDimension();

    rOutput.SetDataToZero(data_dimension);

    auto& r_output_data = rOutput.GetData();
    const auto& r_input_data = rInput.GetData();

    IndexPartition<IndexType>(rMatrix.size1()).for_each([&](const IndexType i) {
        for (IndexType d = 0; d < data_dimension; ++d) {
            auto& r_value = r_output_data[i * data_dimension + d];
            for (IndexType j = 0; j < rMatrix.size2(); ++j) {
                r_value += rMatrix(i, j) * r_input_data[j * data_dimension + d];
            }
        }
    });
}

void ContainerVariableDataHolderUtils::Transpose(
    Matrix& rOutput,
    const Matrix& rInput)
{
    if (rOutput.size1() != rInput.size2() || rOutput.size2() != rInput.size1()) {
        rOutput.resize(rInput.size2(), rInput.size1(), false);
    }

    IndexPartition<IndexType>(rInput.size1()).for_each([&](const IndexType i){
        for (IndexType j = 0; j < rInput.size2(); ++j) {
            rOutput(j, i) = rInput(i, j);
        }
    });
}

void ContainerVariableDataHolderUtils::Transpose(
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
void ContainerVariableDataHolderUtils::ComputeNumberOfNeighbourEntities(
    ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rOutput,
    const ContainerVariableDataHolderBase<TContainerType>& rInput)
{
    KRATOS_ERROR_IF(&rOutput.GetModelPart() != &rInput.GetModelPart())
        << "Output container and input container model parts mismatch. "
           "Followings are the container details:"
        << "\n\tOutput container: " << rOutput
        << "\n\tInput container : " << rInput << "\n";

    // clear the neighbour count storing variable
    VariableUtils().SetNonHistoricalVariableToZero(TEMPORARY_SCALAR_VARIABLE_1, rOutput.GetModelPart().Nodes());

    // create a dummy copy input data container to access its nodes and modify data
    ContainerVariableDataHolder<TContainerType, NonHistoricalContainerDataIO> dummy_input_container(rOutput.GetModelPart());

    block_for_each(dummy_input_container.GetContainer(), [&](auto& rEntity) {
        auto& r_geometry = rEntity.GetGeometry();
        for (auto& r_node : r_geometry) {
            AtomicAdd(r_node.GetValue(TEMPORARY_SCALAR_VARIABLE_1), 1.0);
        }
    });

    rOutput.GetModelPart().GetCommunicator().AssembleNonHistoricalData(TEMPORARY_SCALAR_VARIABLE_1);

    // now read in the nodal data
    ContainerVariableDataHolder<ModelPart::NodesContainerType, NonHistoricalContainerDataIO> dummy_read(rOutput.GetModelPart());
    dummy_read.ReadDataFromContainerVariable(TEMPORARY_SCALAR_VARIABLE_1);

    // now fill the rOutput
    rOutput.CopyDataFrom(dummy_read);
}

template<class TContainerType>
void ContainerVariableDataHolderUtils::MapContainerVariableDataHolderToNodalVariableDataHolder(
    ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rOutput,
    const ContainerVariableDataHolderBase<TContainerType>& rInput,
    const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rNeighbourEntities)
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

    KRATOS_ERROR_IF(rNeighbourEntities.GetDataDimension() != 1)
        << "Neighbour entities container can only have data with dimensionality = 1"
        << "\n\tNeighbour entities container:" << rNeighbourEntities << "\n";

    // reset temporary variables
    std::visit([&](auto&& p_variable) {
        VariableUtils().SetNonHistoricalVariableToZero(*p_variable, rOutput.GetModelPart().Nodes());
    }, ContainerVariableDataHolderUtilsHelper::GetTemporaryVariable(rInput.GetDataDimension()));

    // copy number of neighbours
    ContainerVariableDataHolder<ModelPart::NodesContainerType, NonHistoricalContainerDataIO> dummy_weights(rNeighbourEntities);

    // assign dummy weights to nodes
    dummy_weights.AssignDataToContainerVariable(TEMPORARY_SCALAR_VARIABLE_2);

    // create a dummy copy input data container to access its nodes and modify data
    ContainerVariableDataHolder<TContainerType, NonHistoricalContainerDataIO> dummy_input_container(rOutput.GetModelPart());

    auto& r_input_data_container = dummy_input_container.GetContainer();
    const auto& r_input_data = rInput.GetData();

    // now distribute the entity values to nodes
    std::visit([&](auto&& p_variable) {
        IndexPartition<IndexType>(r_input_data_container.size()).for_each([&](const IndexType EntityIndex) {
            auto r_entity = (r_input_data_container.begin() + EntityIndex);
            auto& r_geometry = r_entity->GetGeometry();
            const auto& entity_value = ContainerVariableDataHolderUtilsHelper::GetEntityData(r_input_data, EntityIndex, rInput.GetDataDimension(), *p_variable);

            for (auto& r_node : r_geometry) {
                r_node.SetLock();
                r_node.GetValue(*p_variable) += entity_value / r_node.GetValue(TEMPORARY_SCALAR_VARIABLE_2);
                r_node.UnSetLock();
            }
        });
    }, ContainerVariableDataHolderUtilsHelper::GetTemporaryVariable(rInput.GetDataDimension()));

    // now read in the nodal data
    ContainerVariableDataHolder<ModelPart::NodesContainerType, NonHistoricalContainerDataIO> dummy_read(rOutput.GetModelPart());
    std::visit([&](auto&& p_variable) {
        dummy_read.ReadDataFromContainerVariable(*p_variable);
    }, ContainerVariableDataHolderUtilsHelper::GetTemporaryVariable(rInput.GetDataDimension()));

    // now copy back the data to the output container
    rOutput.CopyDataFrom(dummy_read);

    KRATOS_CATCH("");
}

template<class TContainerType>
void ContainerVariableDataHolderUtils::MapNodalVariableDataHolderToContainerVariableDataHolder(
    ContainerVariableDataHolderBase<TContainerType>& rOutput,
    const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rInput)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(&rOutput.GetModelPart() != &rInput.GetModelPart())
        << "Output container and input container model parts mismatch. "
           "Followings are the container details:"
        << "\n\tOutput container: " << rOutput
        << "\n\tInput container : " << rInput << "\n";

    // create a dummy with the model part and nodes
    ContainerVariableDataHolder<ModelPart::NodesContainerType, NonHistoricalContainerDataIO> dummy(rOutput.GetModelPart());
    dummy.CopyDataFrom(rInput);

    std::visit([&](auto&& p_variable) {
        // first assign input data to nodes
        dummy.AssignDataToContainerVariable(*p_variable);

        rOutput.SetDataToZero(rInput.GetDataDimension());

        // compute the entity valeus.
        IndexPartition<IndexType>(rOutput.GetContainer().size()).for_each([&](const IndexType EntityIndex) {
            const auto r_entity = (rOutput.GetContainer().begin() + EntityIndex);
            const auto& r_geometry = r_entity->GetGeometry();
            const IndexType number_of_nodes = r_geometry.size();

            for (const auto& r_node : r_geometry) {
                ContainerVariableDataHolderUtilsHelper::AddToEntityData(rOutput.GetData(), EntityIndex, rInput.GetDataDimension(), r_node.GetValue(*p_variable) / number_of_nodes);
            }
        });
    }, ContainerVariableDataHolderUtilsHelper::GetTemporaryVariable(rInput.GetDataDimension()));

    KRATOS_CATCH("");
}

template<class TContainerType>
void ContainerVariableDataHolderUtils::ComputeVariableDataHolderProductWithEntityMatrix(
    ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rOutput,
    const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>& rNodalValues,
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

    ContainerVariableDataHolder<TContainerType, NonHistoricalContainerDataIO> dummy_container(rOutput.GetModelPart());

    KRATOS_ERROR_IF(dummy_container.GetContainer().size() != rEntities.size())
        << "Provided entities container size mismatch with output container variable data holder size. "
        << "[ Provided entities size = " << rEntities.size()
        << ", output data container entities size = " << rOutput.GetContainer().size() << " ].\n";

    const IndexType data_dimension = rNodalValues.GetDataDimension();

    std::visit([&](auto&& p_variable_pair) {
        const auto& r_input_variable = *std::get<0>(p_variable_pair);
        const auto& r_output_variable = *std::get<1>(p_variable_pair);

        // assign nodal values
        ContainerVariableDataHolder<ModelPart::NodesContainerType, NonHistoricalContainerDataIO> nodal_write(rNodalValues);
        nodal_write.AssignDataToContainerVariable(r_input_variable);

        // clear the output variable
        VariableUtils().SetNonHistoricalVariableToZero(r_output_variable, rOutput.GetModelPart().Nodes());

        block_for_each(rEntities, tls_type(), [&](auto& rEntity, tls_type& rTLS) {
            Vector& r_input_values = std::get<0>(rTLS);
            Vector& r_output_values = std::get<1>(rTLS);
            Matrix& r_entity_matrix = std::get<2>(rTLS);
            auto& r_geometry = rEntity.GetGeometry();

            const IndexType number_of_nodes = r_geometry.size();
            const IndexType local_size = data_dimension * number_of_nodes;

            // initialize data vector, matrix sizes
            if (r_input_values.size() != local_size) {
                r_input_values.resize(local_size, false);
                r_output_values.resize(local_size, false);
            }

            noalias(r_input_values) = ZeroVector(local_size);
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                ContainerVariableDataHolderUtilsHelper::AddToEntityData(r_input_values, i, data_dimension, r_geometry[i].GetValue(r_input_variable));
            }

            rEntity.Calculate(rMatrixVariable, r_entity_matrix, rOutput.GetModelPart().GetProcessInfo());
            KRATOS_ERROR_IF(r_entity_matrix.size1() != local_size ||
                            r_entity_matrix.size2() != local_size)
                << "Entity matrix with entity id " << rEntity.Id() << " for "
                << rMatrixVariable.Name() << " is not having required size of ("
                << local_size << ", " << local_size
                << "). [ Obtained matrix: " << r_entity_matrix << " ].\n";

            noalias(r_output_values) = prod(r_entity_matrix, r_input_values);

            // assign the vector to entity nodal values
            for (IndexType i = 0; i < number_of_nodes; ++i) {
                auto& r_node = r_geometry[i];
                r_node.SetLock();
                r_node.GetValue(r_output_variable) += ContainerVariableDataHolderUtilsHelper::GetEntityData(r_output_values, i, data_dimension, r_output_variable);
                r_node.UnSetLock();
            }
        });

        // assemble data for MPI
        dummy_container.GetModelPart().GetCommunicator().AssembleNonHistoricalData(r_output_variable);

        // now read the data
        nodal_write.ReadDataFromContainerVariable(r_output_variable);
        rOutput.CopyDataFrom(nodal_write);

    }, ContainerVariableDataHolderUtilsHelper::GetTemporaryVariable1And2(data_dimension));

    KRATOS_CATCH("");
}

// template instantiations
#define INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ContainerType)                                                                                                                                                                                           \
    template double ContainerVariableDataHolderUtils::EntityMaxNormL2(const ContainerVariableDataHolderBase<ContainerType>&);                                                                                                                                  \
    template double ContainerVariableDataHolderUtils::NormInf(const ContainerVariableDataHolderBase<ContainerType>&);                                                                                                                                          \
    template double ContainerVariableDataHolderUtils::NormL2(const ContainerVariableDataHolderBase<ContainerType>&);                                                                                                                                           \
    template double ContainerVariableDataHolderUtils::InnerProduct(const ContainerVariableDataHolderBase<ContainerType>&, const ContainerVariableDataHolderBase<ContainerType>&);                                                                              \
    template void ContainerVariableDataHolderUtils::ProductWithEntityMatrix(ContainerVariableDataHolderBase<ContainerType>&, const typename UblasSpace<double, CompressedMatrix, Vector>::MatrixType&, const ContainerVariableDataHolderBase<ContainerType>&); \
    template void ContainerVariableDataHolderUtils::ProductWithEntityMatrix(ContainerVariableDataHolderBase<ContainerType>&, const Matrix&, const ContainerVariableDataHolderBase<ContainerType>&);

// template void ContainerVariableDataHolderUtils::ComputeNumberOfNeighbourEntities(ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&, const ContainerVariableDataHolderBase<ModelPart::ConditionsContainerType>&);

#define INSTANTIATE_NON_NODAL_UTILITY_METHOD_FOR_CONTAINER_TYPE(ContainerType)                                                                                                                                                                                                                              \
    template void ContainerVariableDataHolderUtils::ComputeNumberOfNeighbourEntities(ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&, const ContainerVariableDataHolderBase<ContainerType>&);                                                                                                 \
    template void ContainerVariableDataHolderUtils::MapContainerVariableDataHolderToNodalVariableDataHolder(ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&, const ContainerVariableDataHolderBase<ContainerType>&, const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&);   \
    template void ContainerVariableDataHolderUtils::MapNodalVariableDataHolderToContainerVariableDataHolder(ContainerVariableDataHolderBase<ContainerType>&, const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&);                                                                          \
    template void ContainerVariableDataHolderUtils::ComputeVariableDataHolderProductWithEntityMatrix(ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&, const ContainerVariableDataHolderBase<ModelPart::NodesContainerType>&, const Variable<Matrix>& rMatrixVariable, ContainerType&);


INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::NodesContainerType)
INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ConditionsContainerType)
INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ElementsContainerType)

INSTANTIATE_NON_NODAL_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ConditionsContainerType)
INSTANTIATE_NON_NODAL_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ElementsContainerType)

#undef INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE
#undef INSTANTIATE_NON_NODAL_UTILITY_METHOD_FOR_CONTAINER_TYPE

}