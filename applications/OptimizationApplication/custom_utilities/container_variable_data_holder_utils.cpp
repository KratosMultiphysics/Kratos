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
#include <cmath>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "container_variable_data_holder_utils.h"

namespace Kratos
{

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

    KRATOS_ERROR_IF_NOT(OptimizationUtils::IsSameModelPart(rContainer1.GetModelPart(), rContainer2.GetModelPart()))
        << "Model part mismatch in InnerProduct calculation. "
        << "Followings are the given containers: \n"
        << "   Container 1: " << rContainer1 << "\n"
        << "   Container 2: " << rContainer2 << "\n";

    return rContainer1.GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(IndexPartition<IndexType>(r_data_1.size()).for_each<SumReduction<double>>([&](const IndexType Index) {
        return r_data_1[Index] * r_data_2[Index];
    }));
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
        << "\n\tOutput data contaienr: " << rOutput << "\n\t";

    KRATOS_ERROR_IF_NOT(rOutput.GetContainer().size() == rMatrix.size1())
        << "Output container size and matrix size1 mismatch. [ Container size = "
        << rOutput.GetContainer().size() << ", Matrix.size1() = " << rMatrix.size1()
        << "Followings are the given containers:"
        << "\n\tInput data container : " << rInput
        << "\n\tOutput data contaienr: " << rOutput << "\n\t";

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
        << "\n\tOutput data contaienr: " << rOutput << "\n\t";

    KRATOS_ERROR_IF_NOT(rOutput.GetContainer().size() == rMatrix.size1())
        << "Output container size and matrix size1 mismatch. [ Container size = "
        << rOutput.GetContainer().size() << ", Matrix.size1() = " << rMatrix.size1()
        << "Followings are the given containers:"
        << "\n\tInput data container : " << rInput
        << "\n\tOutput data contaienr: " << rOutput << "\n\t";

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

// template instantiations
#define INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ContainerType)                                                                                                                                                                           \
    template double ContainerVariableDataHolderUtils::EntityMaxNormL2(const ContainerVariableDataHolderBase<ContainerType>&);                                                                                                                  \
    template double ContainerVariableDataHolderUtils::NormInf(const ContainerVariableDataHolderBase<ContainerType>&);                                                                                                                          \
    template double ContainerVariableDataHolderUtils::InnerProduct(const ContainerVariableDataHolderBase<ContainerType>&, const ContainerVariableDataHolderBase<ContainerType>&);                                                              \
    template void ContainerVariableDataHolderUtils::ProductWithEntityMatrix(ContainerVariableDataHolderBase<ContainerType>&, const typename UblasSpace<double, CompressedMatrix, Vector>::MatrixType&, const ContainerVariableDataHolderBase<ContainerType>&); \
    template void ContainerVariableDataHolderUtils::ProductWithEntityMatrix(ContainerVariableDataHolderBase<ContainerType>&, const Matrix&, const ContainerVariableDataHolderBase<ContainerType>&);

INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::NodesContainerType)
INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ConditionsContainerType)
INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE(ModelPart::ElementsContainerType)

#undef INSTANTIATE_UTILITY_METHOD_FOR_CONTAINER_TYPE

}